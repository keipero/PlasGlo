r"""保存済みのZDPlasKin窒素データを原本を変更せずに比較する。

確認済みのローカルソースを出力先へコピーし、比較用ネイティブ計算をビルドする。物理設定の変更は電子間衝突の無効化だけとする。別に、元データの96
反応・25断面積をBOLOS/SciPyで計算する。通常のKomuro機構は書き換えない。
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from scipy.integrate import solve_ivp

from kinetics import KB, KB_EV, simulate
from bolos_coupling import BolosSwarm, CROSS_SECTIONS, normalize_eedf
from bolos import grid, parser as collision_parser, solver


def read_csv(path):
    r"""ネイティブ密度履歴を読み、cm^-3からm^-3へ変換する。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。

    Returns
    -------
    tuple
        粒子名、時刻配列 [s]、密度配列 [m^-3]。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     path = Path(directory) / "densities.csv"
    ...     text = chr(10).join(["time_s,E", "0,1", "1,2", ""])
    ...     _ = path.write_text(text, encoding="utf-8")
    ...     names, time, density = read_csv(path)
    >>> names, density.tolist()
    (['E'], [[1000000.0, 2000000.0]])
    """
    with Path(path).open(encoding="utf-8-sig", newline="") as fp:
        rows = csv.reader(fp)
        names = next(rows)
        data = np.array([[float(v) for v in row] for row in rows if row])
    if (
        names[0] != "time_s"
        or not np.isfinite(data).all()
        or not np.all(np.diff(data[:, 0]) > 0)
    ):
        raise ValueError("Expected finite, strictly increasing time_s data.")
    return (
        names[1:],
        data[:, 0],
        data[:, 1:].T * 1e6,
    )  # Archive explicitly writes cm^-3.


def number(text):
    r"""Fortran倍精度指数をPythonの指数表記に変換する。

    Parameters
    ----------
    text : str
        解析対象の文字列。

    Returns
    -------
    str
        d/Dによる指数をeによる指数へ変換した数値式。

    Examples
    --------
    >>> number("1.d-10")
    '1.e-10'
    """
    return re.sub(r"(?i)(\d+(?:\.\d*)?|\.\d+)[d]([+-]?\d+)", r"\1e\2", text)


def driver_conditions(path):
    r"""元のネイティブ実行ドライバーから計算条件を読み取る。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。

    Returns
    -------
    dict
        圧力、温度、換算電界、初期電子密度、終了時刻。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> from py3zdplaschem.paths import SAMPLES_DIR
    >>> values = driver_conditions(
    ...     SAMPLES_DIR / "nitrogen/user_code.F90"
    ... )
    >>> values["temperature_k"]
    300.0
    """
    source = Path(path).read_text(encoding="utf-8")
    found = {}
    for name in (
        "pressure_pa",
        "temperature_k",
        "field_td",
        "seed_density",
        "end_time_s",
    ):
        match = re.search(
            r"(?im)^\s*double precision, parameter :: "
            + name
            + r"\s*=\s*([^!\n]+)",
            source,
        )
        if not match:
            raise ValueError(f"Cannot verify driver parameter {name}")
        found[name] = float(number(match[1]).strip())
    return found


class Expression:
    r"""確認済みFortran速度式のための制限付き算術評価器。

    eval/exec、属性参照、ファイル操作、任意の関数呼び出しは許可しない。
    """

    functions = {
        "exp": math.exp,
        "log": math.log,
        "sqrt": math.sqrt,
        "min": min,
        "max": max,
    }

    def __init__(self, text):
        r"""制限付き算術式を構文木へ変換する。

        Parameters
        ----------
        text : str
            解析対象の文字列。

        Examples
        --------
        >>> expr = Expression("2.d0 * Tg")
        >>> expr({"tg": 300.0})
        600.0
        """
        self.node = ast.parse(number(text).lower().strip(), mode="eval").body

    def __call__(self, values):
        r"""指定された変数の値で構文解析済みの式を評価する。

        Parameters
        ----------
        values : mapping
            変数名をキー、数値を値とする式評価用の辞書。

        Returns
        -------
        float
            制限付き算術式の計算結果。

        Examples
        --------
        >>> Expression("sqrt(tg)")({"tg": 900.0})
        30.0
        """

        def visit(node):
            r"""許可された算術構文木のノードを1つ評価する。

            Parameters
            ----------
            node : ast.AST
                制限付き式評価器で処理する算術構文木のノード。

            Returns
            -------
            float
                構文木ノードが表す数値。

            Raises
            ------
            ValueError
                入力が要求する形状、単位、モデルの制約を満たさない場合。
            """
            if isinstance(node, ast.Constant) and isinstance(
                node.value, (int, float)
            ):
                return node.value
            if isinstance(node, ast.Name):
                return values[node.id]
            if isinstance(node, ast.UnaryOp) and isinstance(
                node.op, (ast.UAdd, ast.USub)
            ):
                return visit(node.operand) * (
                    -1 if isinstance(node.op, ast.USub) else 1
                )
            if isinstance(node, ast.BinOp):
                a, b = visit(node.left), visit(node.right)
                if isinstance(node.op, ast.Add):
                    return a + b
                if isinstance(node.op, ast.Sub):
                    return a - b
                if isinstance(node.op, ast.Mult):
                    return a * b
                if isinstance(node.op, ast.Div):
                    return a / b
                if isinstance(node.op, ast.Pow):
                    return a**b
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id in self.functions
                and not node.keywords
            ):
                return self.functions[node.func.id](
                    *(visit(a) for a in node.args)
                )
            raise ValueError(f"Unsupported rate expression: {ast.dump(node)}")

        return float(visit(self.node))


def solve_archive_swarm(path, field, temperature, cells=400, *, target="N2"):
    r"""保存コードの純気体EEDFと電子衝突係数を計算する。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。
    field : float
        換算電界E/N [Td]。
    temperature : float
        一定の気体温度 [K]。
    cells : int, optional
        BOLOSで使用するエネルギー格子のセル数。
    target : str, optional
        断面積データで使用される純気体の基底状態の標的名。

    Returns
    -------
    tuple
        SIチャンネル係数、EEDF診断、解析した衝突定義。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    """
    with Path(path).open(encoding="utf-8") as fp:
        collisions = collision_parser.parse(fp)
    b = solver.BoltzmannSolver(grid.LinearGrid(0, 60, cells))
    processes = b.load_collisions(collisions)
    b.set_density(target, 1.0)
    b.EN = field * solver.TOWNSEND
    b.kT = temperature * solver.KB / solver.ELECTRONVOLT
    with np.errstate(all="ignore"):
        b.init()
        coarse = b.grid
        f = b.converge(b.maxwell(2), maxn=1000, rtol=1e-6)
        b.grid = grid.QuadraticGrid(0, 15 * b.mean_energy(f), cells)
        b.init()
        f, it, error = b.converge(
            b.grid.interpolate(f, coarse), maxn=1000, rtol=1e-6, full=True
        )
        f, minimum, tail_repaired = normalize_eedf(b, f)
        rates = {
            f"{p.target_name}->{p.product}".upper(): float(b.rate(f, p))
            for p in processes
            if p.product
        }
        if not np.isfinite(list(rates.values())).all() or any(
            value < 0 for value in rates.values()
        ):
            raise RuntimeError("Invalid archived-dataset rate coefficients")
        reduced_mobility = float(b.mobility(f))
        if not np.isfinite(reduced_mobility) or reduced_mobility <= 0:
            raise RuntimeError("Invalid BOLOS electron reduced mobility")
    return (
        rates,
        dict(
            mean_energy_eV=float(b.mean_energy(f)),
            iterations=int(it),
            error=float(error),
            collision_count=len(processes),
            uncorrected_eedf_min=minimum,
            eedf_tail_repaired=tail_repaired,
            electron_reduced_mobility_1_V_m_s=reduced_mobility,
        ),
        collisions,
    )


def side(text):
    r"""旧形式の反応式の片側にある同一粒子種の出現数を数える。

    Parameters
    ----------
    text : str
        解析対象の文字列。

    Returns
    -------
    collections.Counter
        大文字表記の粒子名をキー、出現数を値とする辞書。

    Examples
    --------
    >>> side("N2 + N2")["N2"]
    2
    """
    return Counter(t.strip().upper() for t in text.strip().split(" + "))


class ArchiveMechanism:
    def __init__(
        self,
        path,
        names,
        conditions,
        swarm_rates,
        energy,
        *,
        expected_counts=(21, 96),
    ):
        r"""反応機構と速度係数のメタデータを構築する。

        Parameters
        ----------
        path : str or pathlib.Path
            読み込む入力ファイルのパス。
        names : sequence of str
            密度配列の行と同じ順序の粒子種名。
        conditions : dict or object
            気体条件、電界設定、旧モデルの設定情報。
        swarm_rates : dict of str to float
            SI単位系の電子衝突チャンネル係数。
        energy : float or None
            平均電子エネルギー [eV]。Teに依存しない式ではNoneを許容する。
        expected_counts : tuple of int or None, optional
            期待する粒子種数と反応数。Noneの場合は件数の検査を行わない。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。
        """
        self.names, self.index = names, {
            s.upper(): i for i, s in enumerate(names)
        }
        self.n0 = conditions["pressure_pa"] / (
            KB * conditions["temperature_k"]
        )
        variables = {
            "tgas": conditions["temperature_k"],
            "en": conditions["field_td"],
            "te": 2 * energy / (3 * KB_EV),
        }
        reactions = []
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            text = line.strip()
            if text.startswith("$") and "=" in text:
                name, expression = text[1:].split("=", 1)
                variables[name.strip().lower()] = Expression(expression)(
                    variables
                )
            if "=>" in text and not text.startswith("#"):
                equation, rate = text.split("!", 1)
                reactions.append((equation.strip(), rate.strip()))
        ns, nr = len(names), len(reactions)
        self.stoich = np.zeros((ns, nr), dtype=int)
        self.powers = np.zeros((ns, nr), dtype=int)
        self.third = np.zeros(nr, dtype=int)
        self.k = np.empty(nr)
        self.orders = np.empty(nr, dtype=int)
        self.equations = []
        for j, (equation, expression) in enumerate(reactions):
            left, right = [side(x) for x in equation.split("=>")]
            if left["ANY_NEUTRAL"] != right["ANY_NEUTRAL"]:
                raise ValueError("Unbalanced ANY_NEUTRAL")
            self.third[j] = left["ANY_NEUTRAL"]
            self.orders[j] = sum(left.values())
            for name in left.keys() | right.keys():
                if name == "ANY_NEUTRAL":
                    continue
                i = self.index[name]
                self.powers[i, j] = left[name]
                self.stoich[i, j] = right[name] - left[name]
            if expression.upper().startswith("BOLSIG "):
                key = expression[7:].replace(" ", "").upper()
                self.k[j] = swarm_rates[key]
            else:
                self.k[j] = Expression(expression)(variables) * 1e-6 ** (
                    self.orders[j] - 1
                )
            self.equations.append(equation)
        if expected_counts is not None and (ns, nr) != expected_counts:
            raise ValueError(
                f"Unexpected archived mechanism: {ns} species, {nr} reactions"
            )
        self.scaled = self.k * self.n0 ** (self.orders - 1)
        normalized_names = [s.upper() for s in names]
        self.charge = np.array(
            [
                -1 if s == "E" else int(s.endswith("^+"))
                for s in normalized_names
            ]
        )
        self.atoms = np.array(
            [
                0 if s == "E" else int(re.match(r"(?:AR|N)(\d*)", s)[1] or 1)
                for s in normalized_names
            ]
        )
        self.neutral = self.charge == 0
        np.testing.assert_array_equal(self.atoms @ self.stoich, 0)
        np.testing.assert_array_equal(self.charge @ self.stoich, 0)
        if not np.isfinite(self.k).all() or (self.k < 0).any():
            raise ValueError("Invalid replay coefficients")

    def terms(self, y):
        r"""反応物密度因子と旧形式の中性第三体密度和を構築する。

        Parameters
        ----------
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        tuple
            非負の規格化密度、反応物の積因子、中性種密度和。
        """
        q = np.maximum(y, 0)
        products = q[:, None] ** self.powers
        return q, products, q[self.neutral].sum()

    def rhs(self, t, y):
        r"""規格化した粒子密度の時間微分を評価する。

        Parameters
        ----------
        t : float
            時刻 [s]。
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        numpy.ndarray
            規格化した粒子密度の時間微分 [s^-1]。
        """
        _, products, m = self.terms(y)
        return self.stoich @ (
            self.scaled * products.prod(axis=0) * m**self.third
        )

    def jac(self, t, y):
        r"""規格化された反応源項の解析的Jacobianを評価する。

        Parameters
        ----------
        t : float
            時刻 [s]。
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        numpy.ndarray
            形状が (粒子種数, 粒子種数) の解析的Jacobian [s^-1]。
        """
        q, products, m = self.terms(y)
        dr = np.zeros((len(self.k), len(y)))
        for i in range(len(y)):
            if y[i] < 0:
                continue
            active = self.powers[i] > 0
            if active.any():
                p = products[:, active].copy()
                p[i] = self.powers[i, active] * q[i] ** (
                    self.powers[i, active] - 1
                )
                dr[active, i] = (
                    self.scaled[active]
                    * p.prod(axis=0)
                    * m ** self.third[active]
                )
            if self.neutral[i]:
                active = self.third > 0
                dr[active, i] += (
                    self.scaled[active]
                    * products[:, active].prod(axis=0)
                    * self.third[active]
                    * m ** (self.third[active] - 1)
                )
        return self.stoich @ dr


def replay_archive(path, names, time, initial, conditions, out):
    r"""元のサンプル反応機構をBOLOS係数で積分する。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。
    names : sequence of str
        密度配列の行と同じ順序の粒子種名。
    time : numpy.ndarray
        厳密に単調増加する出力時刻 [s]。
    initial : numpy.ndarray
        SI単位系の各粒子種の初期密度 [m^-3]。
    conditions : dict or object
        気体条件、電界設定、旧モデルの設定情報。
    out : pathlib.Path
        比較結果とレポートを保存するディレクトリ。

    Returns
    -------
    tuple
        密度履歴、反応機構、電子群診断、係数、衝突定義。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    """
    rates, swarm, collisions = solve_archive_swarm(
        path / "bolsigdb.dat",
        conditions["field_td"],
        conditions["temperature_k"],
    )
    model = ArchiveMechanism(
        path / "kinet.inp", names, conditions, rates, swarm["mean_energy_eV"]
    )
    solution = solve_ivp(
        model.rhs,
        (time[0], time[-1]),
        initial / model.n0,
        t_eval=time,
        method="BDF",
        jac=model.jac,
        rtol=1e-7,
        atol=1e-4 / model.n0,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    density = solution.y * model.n0
    atom = model.atoms @ density
    conservation = float(np.max(abs(atom / atom[0] - 1)))
    np.savez_compressed(
        out / "same_mechanism_bolos.npz",
        time_s=time,
        densities_m3=density,
        species=np.array(names),
    )
    (out / "same_mechanism_bolos.json").write_text(
        json.dumps(
            dict(
                swarm=swarm,
                max_atom_relative_error=conservation,
                rhs_evaluations=solution.nfev,
                electron_electron_collisions=False,
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    with (out / "same_mechanism_coefficients.csv").open(
        "w", newline="", encoding="utf-8"
    ) as fp:
        writer = csv.writer(fp)
        writer.writerow(["reaction", "order", "coefficient_SI"])
        writer.writerows(zip(model.equations, model.orders, model.k))
    return density, model, swarm, rates, collisions


def native_control(reference, out, disable_ee=True):
    r"""比較用にコピーしたZDPlasKinをビルドして実行する。

    Parameters
    ----------
    reference : pathlib.Path
        元のネイティブ実行フォルダ。変更するのは比較用コピーだけ。
    out : pathlib.Path
        比較結果とレポートを保存するディレクトリ。
    disable_ee : bool, optional
        コピーしたネイティブ計算で電子間衝突を無効にするかどうか。

    Returns
    -------
    tuple
        粒子名、ネイティブ時刻 [s]、SI密度履歴。

    Raises
    ------
    FileNotFoundError
        必要な元ソース、コンパイラー、入力ファイルが存在しない場合。
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    """
    root = reference.parents[2]
    dvode = (
        root / "ZDPlasKin_2.0a_Windows/ZDPlasKin_2.0a_Windows/dvode_f90_m.F90"
    )
    compiler = Path(
        "C:/rtools44/x86_64-w64-mingw32.static.posix/bin/gfortran.exe"
    )
    if not dvode.is_file() or not compiler.is_file():
        raise FileNotFoundError(
            "Local original DVODE source / compiler unavailable for controlled"
            " comparison"
        )
    target = out / (
        "zdplaskin_ee_off" if disable_ee else "zdplaskin_original_rebuild"
    )
    target.mkdir(exist_ok=True)
    for name in (
        "zdplaskin_m.F90",
        "kinet.inp",
        "bolsigdb.dat",
        "bolsig_x86_64_g.dll",
    ):
        shutil.copy2(reference / name, target / name)
    shutil.copy2(dvode, target / "dvode_f90_m.F90")
    source = (reference / "user_code.F90").read_text(encoding="utf-8")
    needle = "  call ZDPlasKin_init()"
    if source.count(needle) != 1:
        raise ValueError("Cannot safely locate native control initialization")
    if disable_ee:
        source = source.replace(
            needle,
            needle + "\n  call ZDPlasKin_set_config(BOLSIG_EE_FRAC=1.d100)",
        )
    (target / "user_code.F90").write_text(source, encoding="utf-8")
    env = os.environ.copy()
    env["PATH"] = str(compiler.parent) + os.pathsep + env.get("PATH", "")
    temp = target / "tmp"
    temp.mkdir(exist_ok=True)
    env["TMP"] = env["TEMP"] = str(temp.resolve())
    args = [
        str(compiler),
        "-std=legacy",
        "-O2",
        "-ffree-line-length-none",
        "dvode_f90_m.F90",
        "zdplaskin_m.F90",
        "user_code.F90",
        "bolsig_x86_64_g.dll",
        "-o",
        "run_control.exe",
    ]
    compiled = subprocess.run(
        args, cwd=target, env=env, capture_output=True, text=True, timeout=120
    )
    (target / "compile.log").write_text(
        compiled.stdout + compiled.stderr, encoding="utf-8"
    )
    if compiled.returncode:
        raise RuntimeError(
            f"Native control compile failed: {compiled.stderr[-2000:]}"
        )
    result = subprocess.run(
        [str((target / "run_control.exe").resolve())],
        cwd=target,
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )
    (target / "stdout.log").write_text(
        result.stdout + result.stderr, encoding="utf-8"
    )
    if result.returncode:
        raise RuntimeError(
            f"Native control failed: {result.stdout[-2000:]}"
            f" {result.stderr[-2000:]}"
        )
    if disable_ee and "electron-electron collisions ON" in result.stdout:
        raise RuntimeError(
            "Electron-electron collisions remained enabled in native control"
        )
    return read_csv(target / "densities.csv")


def main():
    r"""コマンドライン比較を実行し、検証用ファイルを保存する。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    FileNotFoundError
        必要な元ソース、コンパイラー、入力ファイルが存在しない場合。
    """
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("reference_csv", type=Path)
    cli.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent
        / "validation_zdplaskin_nitrogen",
    )
    args = cli.parse_args()
    ref = args.reference_csv.resolve()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    names, time, density = read_csv(ref)
    conditions = driver_conditions(ref.parent / "user_code.F90")
    np.testing.assert_allclose(
        time[-1], conditions["end_time_s"], rtol=1e-12, atol=0
    )
    print(
        "Rebuilding the unmodified original to verify the supplied CSV...",
        flush=True,
    )
    check_names, check_time, check_density = native_control(
        ref.parent, out, disable_ee=False
    )
    if check_names != names:
        raise ValueError("Rebuilt original species mismatch")
    np.testing.assert_allclose(time, check_time, rtol=1e-14, atol=1e-25)
    reproduction_error = float(
        np.max(
            np.max(abs(check_density - density), axis=1) / density.max(axis=1)
        )
    )
    if reproduction_error > 1e-5:
        raise RuntimeError(
            "Unmodified original failed to reproduce source CSV:"
            f" {reproduction_error}"
        )
    print(
        "Building ZDPlasKin native control with electron-electron collisions"
        " disabled...",
        flush=True,
    )
    off_names, off_time, off_density = native_control(ref.parent, out)
    if off_names != names:
        raise ValueError("Native species mismatch")
    np.testing.assert_allclose(time, off_time, rtol=1e-14, atol=1e-25)
    print(
        "Replaying the same 96 reactions and 25 cross sections using"
        " BOLOS/SciPy...",
        flush=True,
    )
    replay, model, swarm, rates, collisions = replay_archive(
        ref.parent, names, time, density[:, 0], conditions, out
    )
    current_path = out / "current_komuro/densities.npz"
    if not current_path.is_file():
        raise FileNotFoundError(
            "Run the documented current_komuro case first."
        )
    with np.load(current_path) as data:
        current_names = data["species"].tolist()
        current = data["densities_m3"]
        current_time = data["time_s"]
    np.testing.assert_allclose(time, current_time, rtol=1e-12, atol=1e-25)
    build_report(
        out,
        ref,
        names,
        time,
        density,
        off_density,
        replay,
        model,
        swarm,
        rates,
        collisions,
        current_names,
        current,
        conditions,
        reproduction_error,
    )
    print(f"Comparison saved: {out}", flush=True)


def build_report(
    out,
    ref,
    names,
    time,
    reference,
    off,
    replay,
    model,
    swarm,
    rates,
    collisions,
    current_names,
    current,
    conditions,
    reproduction_error,
):
    r"""比較指標、図、入力の出典情報を保存する。

    Parameters
    ----------
    out : pathlib.Path
        比較結果とレポートを保存するディレクトリ。
    ref : pathlib.Path
        出典照合に使用する元のZDPlasKin密度CSV。
    names : sequence of str
        密度配列の行と同じ順序の粒子種名。
    time : numpy.ndarray
        厳密に単調増加する出力時刻 [s]。
    reference : numpy.ndarray
        元のネイティブ計算の密度履歴 [m^-3]。
    off : numpy.ndarray
        電子間衝突を無効にしたネイティブ計算の密度履歴。
    replay : numpy.ndarray
        条件をそろえたPython計算の密度履歴 [m^-3]。
    model : MassActionKernel
        粒子種の順序と反応源項を定義するモデル。
    swarm : dict
        条件をそろえたBOLOS計算の診断情報。
    rates : dict or numpy.ndarray
        SI単位系の電子衝突チャンネル係数。
    collisions : sequence of dict
        比較のメタデータに残す、解析済み衝突定義。
    current_names : sequence of str
        旧空気モデルとの比較で使用する粒子種の順序。
    current : numpy.ndarray
        旧空気モデルの密度履歴 [m^-3]。
    conditions : dict or object
        気体条件、電界設定、旧モデルの設定情報。
    reproduction_error : float
        ネイティブ再ビルド結果のピーク規格化最大誤差。
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    index = {s: i for i, s in enumerate(names)}
    ci = {s: i for i, s in enumerate(current_names)}
    pairs = [
        ("E", ["e"], "electron"),
        ("N2", ["N2"], "ground N2"),
        ("N2^+", ["N2+"], "N2+"),
        ("N2(C3)", ["N2*(C)"], "C state: grouping differs"),
        ("N2(A3)", ["N2*(A1)", "N2*(A2)"], "A states: grouping differs"),
        ("N2(B3)", ["N2*(B)"], "B states: grouping differs"),
        ("N", ["N(^4S)"], "ground N"),
        ("N2(V1)", ["N2*(v1)"], "N2(v1)"),
        ("N4^+", ["N4+"], "N4+"),
    ]
    ns = time * 1e9
    colors = {
        "reference": "#2369bd",
        "current": "#d46913",
        "off": "#7a7a7a",
        "replay": "#299451",
    }
    fig, axes = plt.subplots(
        3, 3, figsize=(13, 10), sharex=True, layout="constrained"
    )
    metrics = []
    for ax, (name, ours, title) in zip(axes.flat, pairs):
        i = index[name]
        project = current[[ci[s] for s in ours]].sum(axis=0)
        for values, label, color, style in (
            (reference[i], "ZDPlasKin supplied", colors["reference"], "-"),
            (project, "Current Komuro model", colors["current"], "-"),
            (
                off[i],
                "ZDPlasKin same 96 reactions, e-e off",
                colors["off"],
                "--",
            ),
            (
                replay[i],
                "BOLOS same 96 reactions, e-e off",
                colors["replay"],
                ":",
            ),
        ):
            ax.plot(ns, values, label=label, color=color, ls=style, lw=1.5)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel(r"Density [m$^{-3}$]")
        ax.grid(alpha=0.25)
        for case, values, target in (
            ("current_vs_supplied", project, reference[i]),
            ("matched_bolos_vs_ee_off", replay[i], off[i]),
            ("ee_off_vs_supplied", off[i], reference[i]),
        ):
            peak = float(target.max())
            ip = int(values.argmax())
            rp = int(target.argmax())
            metrics.append(
                dict(
                    case=case,
                    species=name,
                    reference_peak_m3=peak,
                    compared_peak_m3=float(values[ip]),
                    reference_peak_ns=float(ns[rp]),
                    compared_peak_ns=float(ns[ip]),
                    peak_relative_difference=float(values[ip] / peak - 1),
                    reference_final_m3=float(target[-1]),
                    compared_final_m3=float(values[-1]),
                    final_ratio=(
                        float(values[-1] / target[-1])
                        if target[-1] > 0
                        else None
                    ),
                    nrmse_by_reference_peak=float(
                        np.sqrt(np.mean((values - target) ** 2)) / peak
                    ),
                )
            )
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, fontsize=9)
    fig.suptitle(
        "Pure N2, 1000 Td continuous, 1 atm, 300 K, initial electrons = 1e14"
        " m^-3\nSupplied record ends at 1 ns; current and archived reaction"
        " mechanisms differ",
        fontsize=12,
    )
    fig.savefig(out / "comparison.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), layout="constrained")
    for ax, name in zip(axes, ("E", "N2(C3)")):
        i = index[name]
        for values, label, color in (
            (reference[i], "ZDPlasKin supplied", colors["reference"]),
            (off[i], "ZDPlasKin e-e off", colors["off"]),
            (replay[i], "BOLOS same mechanism, e-e off", colors["replay"]),
        ):
            ax.semilogy(
                ns,
                np.where(values > 0, values, np.nan),
                label=label,
                color=color,
            )
        ax.set(
            xlim=(0, 0.15),
            xlabel="Time [ns]",
            ylabel=r"Density [m$^{-3}$]",
            title=name,
        )
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle("Early avalanche: controlled same-mechanism comparison")
    fig.savefig(out / "early_avalanche.png", dpi=160)
    plt.close(fig)
    with (out / "comparison_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as fp:
        writer = csv.DictWriter(fp, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    e = index["E"]
    n2 = index["N2"]
    c = index["N2(C3)"]
    gamma = float(
        np.log(reference[e, 1] / reference[e, 0]) / (time[1] - time[0])
    )
    interval = float(time[1] - time[0])
    kion = gamma / reference[n2, 0]
    kc = float(
        reference[c, 1]
        * gamma
        / (reference[e, 0] * reference[n2, 0] * np.expm1(gamma * interval))
    )
    linked_c = sum(
        value
        for equation, value in zip(model.equations, model.k)
        if equation.startswith("e + N2 => e + N2(C3)")
    )
    ee_fraction = reference[e] / np.maximum(reference[n2], 1)
    onset = np.flatnonzero(ee_fraction >= 1e-5)
    current_density = float(
        conditions["pressure_pa"] / (KB * conditions["temperature_k"])
    )
    report = dict(
        conditions=conditions,
        reference_csv=str(ref),
        reference_rows=len(time),
        original_rebuild_max_error_normalized_by_species_peak=(
            reproduction_error
        ),
        reference_density_unit=(
            "cm^-3, converted to m^-3 by multiplying by 1e6"
        ),
        actual_end_ns=float(time[-1] * 1e9),
        electron_electron_onset_estimated_ns=(
            float(ns[onset[0]]) if len(onset) else None
        ),
        initial_transport_comparison=dict(
            reference_ionization_estimate_m3_s=kion,
            bolos_siglo_ionization_m3_s=rates["N2->N2^+"],
            reference_C_lumped_estimate_m3_s=kc,
            bolos_siglo_C_lumped_m3_s=linked_c,
        ),
        C_loss_initial_lifetime_ns=dict(
            current=1e9 / (2.8e7 + 2e-17 * current_density),
            archive=1e9 / (2.45e7 + 1e-17 * current_density),
        ),
        archived_swarm_bolos=swarm,
        metrics=metrics,
    )
    add_all_species_verification(report, out, names, off, replay)
    provenance = {}
    for name in (
        "densities.csv",
        "user_code.F90",
        "kinet.inp",
        "bolsigdb.dat",
        "zdplaskin_m.F90",
        "stdout.log",
    ):
        p = ref.parent / name
        provenance[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    report["source_sha256"] = provenance
    (out / "report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_markdown(out, report)


def add_all_species_verification(report, out, names, off, replay):
    r"""全粒子種のピーク密度で規格化した誤差を比較レポートへ追加する。

    Parameters
    ----------
    report : dict
        保存または更新する比較指標と出典情報。
    out : pathlib.Path
        比較結果とレポートを保存するディレクトリ。
    names : sequence of str
        密度配列の行と同じ順序の粒子種名。
    off : numpy.ndarray
        電子間衝突を無効にしたネイティブ計算の密度履歴。
    replay : numpy.ndarray
        条件をそろえたPython計算の密度履歴 [m^-3]。
    """
    peaks = off.max(axis=1)
    entries = []
    for i, name in enumerate(names):
        entries.append(
            dict(
                species=name,
                max_absolute_difference_over_reference_peak=float(
                    abs(replay[i] - off[i]).max() / peaks[i]
                ),
                nrmse_over_reference_peak=float(
                    np.sqrt(np.mean((replay[i] - off[i]) ** 2)) / peaks[i]
                ),
                peak_relative_difference=float(replay[i].max() / peaks[i] - 1),
            )
        )
    report["matched_all_species_verification"] = entries
    report["matched_all_species_max_pointwise_error_over_reference_peak"] = (
        max(e["max_absolute_difference_over_reference_peak"] for e in entries)
    )
    with (out / "matched_all_species_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as fp:
        writer = csv.DictWriter(fp, fieldnames=list(entries[0]))
        writer.writeheader()
        writer.writerows(entries)
    coarse = out / "current_interval_check/densities.npz"
    fine = out / "current_komuro/densities.npz"
    if coarse.is_file() and fine.is_file():
        with np.load(coarse) as a, np.load(fine) as b:
            refinement = []
            for name in ("e", "N2", "N2*(C)", "N4+"):
                i = b["species"].tolist().index(name)
                refinement.append(
                    dict(
                        species=name,
                        max_absolute_difference_over_fine_peak=float(
                            abs(
                                a["densities_m3"][i] - b["densities_m3"][i]
                            ).max()
                            / b["densities_m3"][i].max()
                        ),
                    )
                )
        report["current_interval_check_0p005_vs_0p0025_ns"] = refinement


def write_markdown(out, report):
    r"""比較概要をMarkdownのレポートに書き込む。

    Parameters
    ----------
    out : pathlib.Path
        比較結果とレポートを保存するディレクトリ。
    report : dict
        保存または更新する比較指標と出典情報。
    """
    rebuild_error = report[
        "original_rebuild_max_error_normalized_by_species_peak"
    ]
    matched_error = report[
        "matched_all_species_max_pointwise_error_over_reference_peak"
    ]
    lines = [
        "# 純N2：ZDPlasKinとBOLOS＋0Dの比較",
        "",
        "元の実行フォルダは変更せず、比較用の計算だけを本フォルダで実行しました。",
        "元CSVの密度はcm^-3です。比較では全密度を1e6倍してm^-3へ統一しています。",
        "",
        "## 確認した条件",
        "",
        (
            "- 純N2、101325 Pa、300 K。初期電子・N2+密度は各1e8 cm^-3 = 1e14"
            " m^-3。"
        ),
        "- E/N=1000 Tdを連続印加。アフターグローなし。",
        (
            "- フォルダ名や古いコメントは10 nsですが、"
            "CSV・実行コードの設定値・ログは1 ns。10001時刻。"
        ),
        (
            "- ZDPlasKinは21種・96反応、SIGLO N2の25断面積。"
            "BOLSIGの標的は基底状態N2のみ。"
        ),
        (
            "- 現行Komuroモデルは44種・134反応、同梱Phelps断面積。"
            "今回の比較実行では純N2を指定。"
        ),
        "",
        "## 比較の構成",
        "",
        "1. supplied：利用者が指定した元CSV。電子―電子衝突が自動的に有効になる設定。",
        "2. current：現行Komuroモデル。外部条件だけを元CSVに合わせた計算。",
        (
            "3. ZD e-e off：元の96反応・断面積・初期条件をそのまま使い、"
            "電子―電子衝突だけを無効にした対照計算。"
        ),
        (
            "4. matched BOLOS：同じ96反応・25断面積をBOLOS＋SciPy BDFで計算。"
            "電子―電子衝突なし。"
        ),
        "matched BOLOSは比較専用です。現行の反応機構は書き換えていません。",
        (
            "native対照計算は元ソースのコピーをコンパイルしました。"
            "変更点は初期化後の `BOLSIG_EE_FRAC=1.d100` の指定だけです。"
        ),
        "",
        (
            "変更前の元ソースも別に再ビルドして元CSVと照合しました。"
            "全21種のピーク密度で規格化した最大誤差は"
            f" {rebuild_error:.6g} です。"
        ),
        "",
        (
            "同じ96反応でのBOLOSとZD e-e offの全21種・全10001時刻の比較では、"
            "各種の参照ピーク密度で規格化した最大の絶対差は"
            f" {100 * matched_error:.4g}%でした。"
            "微量密度での相対差ではありません。"
        ),
        "",
        "## 数値比較",
        "",
        "| 比較 | 粒子 | 最大密度の差 [%] | 終端密度の比 | 全軌道NRMSE [%] |",
        "|---|---|---:|---:|---:|",
    ]
    for row in report["metrics"]:
        if row["species"] not in ("E", "N2", "N2(C3)", "N4^+"):
            continue
        ratio = (
            "—" if row["final_ratio"] is None else f"{row['final_ratio']:.5g}"
        )
        lines.append(
            f"| {row['case']} | {row['species']} |"
            f" {100 * row['peak_relative_difference']:.4g} | {ratio} |"
            f" {100 * row['nrmse_by_reference_peak']:.4g} |"
        )
    lines += [
        "",
        (
            "NRMSEは `sqrt(mean((比較値-参照値)^2))/参照値の最大値`。"
            "比較前の時刻シフトや密度の係数合わせは行っていません。"
        ),
        "",
        "## N2(C3)の定義と寿命の違い",
        "",
        (
            "- ZD機構はC3への励起にE3・a''1への励起を直接合算します。"
            "現行機構はCとEを別々に追跡します。"
        ),
        (
            "- A・Bの状態への振り分けも同一ではありません。"
            "A状態の比較は現行A1+A2を参考表示したもので、"
            "完全に等価な状態ではありません。"
        ),
        "- ZD機構のC3の損失：放射2.45e7 s^-1、N2失活1e-17 m^3/s。",
        "- 現行Cの損失：放射2.8e7 s^-1、N2失活が2経路で合計2e-17 m^3/s。",
        (
            "- 初期純N2条件の実効寿命：ZD機構"
            f" {report['C_loss_initial_lifetime_ns']['archive']:.5g} ns、現行"
            f" {report['C_loss_initial_lifetime_ns']['current']:.5g} ns。"
            "空気80:20の0.609 nsとは条件が異なります。"
        ),
        "",
        "## その他の機構・データの違い",
        "",
        (
            "- ZD機構はN+・N3+、原子とのイオン変換、振動緩和VT、"
            "原子再結合などを含みます。"
        ),
        (
            "- 現行は回転励起を独立したN2*(J)密度として追跡しますが、"
            "ZD機構はROTをEEDFの衝突に含めるだけです。"
        ),
        (
            "- 現行PhelpsデータにはN2+(B2SIGMA)の電離チャンネルが追加されています。"
            "元SIGLOデータは全電離の1チャンネルです。"
        ),
        (
            "- 初期電離率が低い区間ではBOLSIG側の電子―電子衝突は無効ですが、"
            "その後は有効になります。BOLOSにはこの項がありません。"
        ),
        "",
        "## 解釈の範囲",
        "",
        (
            "元ZDPlasKinのログも、"
            "BOLSIGに標的として設定していない種の割合が大きくなると警告しています。"
            "現行モデル側にも同じ制約があります。"
        ),
        (
            "両方とも1000 Tdを連続印加すると基底状態N2が強く消費されます。"
            "数値的な相互比較と、"
            "実放電に対する物理的妥当性の確認は区別する必要があります。"
        ),
        "",
        (
            "外部条件を揃えただけのcurrentとsuppliedの差を、"
            "そのままBOLOSの誤差と判断することはできません。"
            "同じ機構・断面積で比較したmatched BOLOSとZD e-e"
            " offの差を併記しています。"
        ),
        "",
        "## 保存ファイル",
        "",
        "- comparison.png：主要9系列の線形軸比較。",
        "- early_avalanche.png：同じ96反応での初期成長の比較。",
        (
            "- comparison_metrics.csv / report.json：数値指標、条件、"
            "初期実効寿命、元ファイルのSHA-256。"
        ),
        "- same_mechanism_bolos.npz：比較用96反応のBOLOS結果。",
        (
            "- zdplaskin_ee_off/：電子―電子衝突を無効化した対照計算のソース・ログ・出力。"
        ),
        "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
