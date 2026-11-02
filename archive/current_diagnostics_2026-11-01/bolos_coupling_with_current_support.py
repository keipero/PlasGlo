r"""定常BOLOSのローカル接続。BOLSIG出力は読み込まない。

同梱Phelps断面積は基底状態N2/O2の標的だけを含む。他の粒子種は化学反応として追跡するが、未提供の電子散乱を仮定したり、その密度を基底状
態の割合へ混ぜたりしない。
"""

from __future__ import annotations
from importlib import import_module

from pathlib import Path
import sys

import numpy as np

# Import the local clone, independently of the current working
# directory.
BOLOS_DIR = Path(__file__).resolve().parent / "Py3bolos"
if not (BOLOS_DIR / "bolos" / "__init__.py").is_file():
    raise ModuleNotFoundError(
        f"Py3bolosのクローンが見つかりません: {BOLOS_DIR}。"
        "このフォルダにkeipero/Py3bolosをクローンしてください。",
        name="bolos",
    )
_loaded_bolos = sys.modules.get("bolos")
if (
    _loaded_bolos is not None
    and Path(_loaded_bolos.__file__).resolve().parent != BOLOS_DIR / "bolos"
):
    raise ImportError(
        "別のbolosが読み込み済みです。Pythonカーネルを再起動してください。"
    )
if str(BOLOS_DIR) not in sys.path:
    sys.path.insert(0, str(BOLOS_DIR))
grid = import_module("bolos").grid
parser = import_module("bolos").parser
solver = import_module("bolos").solver

CROSS_SECTIONS = (
    Path(__file__).resolve().parent / "Phelps/Phelps_N2-O2_20261026.txt"
)
THREE_BODY_TARGET = "__O2_three_body__"

# The C identifiers in the chemistry refer to this particular Phelps
# dataset.
# Verify products as well as thresholds: several channels share a
# threshold.
CHANNELS = [
    ("N2", "EFFECTIVE", None, 0),
    *[
        ("N2", "EXCITATION", p, t)
        for p, t in [
            ("N2(rot)", 0.020),
            ("N2(v1)", 0.29),
            ("N2(v1res)", 0.291),
            ("N2(v2)", 0.59),
            ("N2(v3)", 0.88),
            ("N2(v4)", 1.17),
            ("N2(v5)", 1.47),
            ("N2(v6)", 1.76),
            ("N2(v7)", 2.06),
            ("N2(v8)", 2.35),
            ("N2(A3,v0-4)", 6.17),
            ("N2(A3,v5-9)", 7.00),
            ("N2(B3)", 7.35),
            ("N2(W3)", 7.36),
            ("N2(A3,v10-)", 7.80),
            ("N2(B'3)", 8.16),
            ("N2(a'1)", 8.40),
            ("N2(a1)", 8.55),
            ("N2(w1)", 8.89),
            ("N2(C3)", 11.03),
            ("N2(E3)", 11.87),
            ("N2(a''1)", 12.25),
            ("N2(SUM)", 13.00),
        ]
    ],
    ("N2", "IONIZATION", "N2^+", 15.60),
    ("N2", "IONIZATION", "N2^+(B2SIGMA)", 18.80),
    ("O2", "ATTACHMENT", None, 0),
    ("O2", "ATTACHMENT", "O^-+O", 0),
    ("O2", "EFFECTIVE", None, 0),
    *[
        ("O2", "EXCITATION", p, t)
        for p, t in [
            ("O2(rot)", 0.020),
            ("O2(v1)", 0.19),
            ("O2(v1res)", 0.19),
            ("O2(v2)", 0.38),
            ("O2(v2res)", 0.38),
            ("O2(v3)", 0.57),
            ("O2(v4)", 0.75),
            ("O2(a1)", 0.977),
            ("O2(b1)", 1.627),
            ("O2(4.5eV)", 4.50),
            ("O2(6.0eV)", 6.00),
            ("O2(8.4eV)", 8.40),
            ("O2(9.97eV)", 9.97),
        ]
    ],
    ("O2", "IONIZATION", "O2^+", 12.06),
]


def normalize_eedf(b, distribution):
    r"""有意な負値を拒否し、丸め誤差の範囲でEEDF末尾を延長する。

    Parameters
    ----------
    b : bolos.solver.BoltzmannSolver
        エネルギー格子と規格化規約を提供するBOLOSソルバー。
    distribution : numpy.ndarray
        ソルバーのエネルギー格子上のEEDF値。

    Returns
    -------
    tuple
        規格化したEEDF、補正前の最小値、末尾補正の有無。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    """
    f = distribution.copy()
    minimum = float(np.min(f))
    if not np.isfinite(f).all() or minimum < -1e-12 * np.max(f):
        raise RuntimeError("BOLOSのEEDFに非有限値または有意な負値があります。")
    # A hard floor in an unresolved tail produces a cliff in log(F).
    # Continue the last decreasing slope resolved above relative 1e-12.
    tail_repaired = False
    if minimum <= 0:
        reliable = np.flatnonzero(f > np.max(f) * 1e-12)
        if len(reliable) >= 2:
            previous, last = reliable[-2:]
            slope = (np.log(f[last]) - np.log(f[previous])) / (
                b.cenergy[last] - b.cenergy[previous]
            )
            if slope < 0 and last < len(f) - 1:
                f[last + 1:] = f[last] * np.exp(
                    np.maximum(
                        slope * (b.cenergy[last + 1:] - b.cenergy[last]), -690
                    )
                )
                tail_repaired = True
    f = b._normalized(np.maximum(f, np.max(f) * 1e-300))
    return f, minimum, tail_repaired


class BolosSwarm:
    r"""固定した換算電界と現在のN2/O2密度でEEDFを再計算する。

    E/NのNはn(N2)+n(O2)。C27の断面積は密度規格化されており、EEDFの損失重みはx(O2)*n(O2)[cm^-3]*si
    gma27となる。補助標的にこの重みを与え、重みを掛けない積分値を三体係数k3
    [m^6/s]へ換算する。O2密度がゼロでも除算を生じさせず、キャッシュも変更しない。
    """

    def __init__(
        self, path=CROSS_SECTIONS, cells=400, rtol=1e-6, max_iterations=1000
    ):
        r"""衝突データを読み込んで検証し、BOLOSを初期化する。

        Parameters
        ----------
        path : str or pathlib.Path, optional
            読み込む入力ファイルのパス。
        cells : int, optional
            BOLOSで使用するエネルギー格子のセル数。
        rtol : float, optional
            収束判定または時間積分の相対許容誤差。
        max_iterations : int, optional
            1回の収束計算で許容するBOLOSの最大反復回数。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。
        """
        if not isinstance(cells, (int, np.integer)) or cells < 50:
            raise ValueError("BOLOSの格子数は50以上の整数にしてください。")
        if not np.isfinite(rtol) or not 0 < rtol < 1:
            raise ValueError("BOLOSのrtolは0より大きく1未満にしてください。")
        if (
            not isinstance(max_iterations, (int, np.integer))
            or max_iterations < 1
        ):
            raise ValueError("BOLOSの最大反復回数は正の整数にしてください。")
        self.path, self.cells = Path(path), cells
        self.rtol, self.max_iterations = rtol, max_iterations
        with self.path.open(encoding="utf-8") as fp:
            collisions = parser.parse(fp)
        if len(collisions) != len(CHANNELS):
            raise ValueError("C1〜C43に対応するPhelps断面積43過程が必要です。")
        self.labels = {}
        for cid, (data, expected) in enumerate(zip(collisions, CHANNELS), 1):
            target, kind, product, threshold = expected
            if (data["target"], data["kind"], data.get("product")) != (
                target,
                kind,
                product,
            ) or not np.isclose(
                data.get("threshold", 0), threshold, rtol=0, atol=5e-4
            ):
                raise ValueError(
                    f"C{cid}の断面積の種・過程・生成物・閾値が一致しません。"
                )
            category = {
                "EFFECTIVE": "Effective (momentum)",
                "EXCITATION": "Excitation",
                "IONIZATION": "Ionization",
                "ATTACHMENT": "Attachment",
            }[kind]
            self.labels[cid] = f"{target} {category}"
            if kind in ("EXCITATION", "IONIZATION"):
                self.labels[cid] += (
                    f" {threshold:.3f} eV"
                    if threshold == 0.020
                    else f" {threshold:.2f} eV"
                )
        self.bsolver = solver.BoltzmannSolver(grid.LinearGrid(0, 60, cells))
        # C27 is a three-body process, so it must not be subtracted as
        # an
        # ordinary two-body inelastic channel from O2 EFFECTIVE data.
        ordinary = self.bsolver.load_collisions(
            collisions[:26] + collisions[27:]
        )
        attachment = dict(collisions[26], target=THREE_BODY_TARGET)
        self.attachment = self.bsolver.add_process(**attachment)
        self.processes = ordinary[:26] + [self.attachment] + ordinary[26:]
        self.eedf = None
        self.last_field_td = None
        self.solve_count = 0

    def solve(self, field_td, densities, gas_temperature=300.0):
        r"""宣言した標的の密度に対する定常EEDFを求める。

        Parameters
        ----------
        field_td : float
            換算電界E/N [Td]。
        densities : dict of str to float
            基底状態N2とO2の数密度 [m^-3]。
        gas_temperature : float, optional
            一定の気体温度 [K]。

        Returns
        -------
        dict
            SIチャンネル係数、平均エネルギー [eV]、組成比、診断情報。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。
        RuntimeError
            ソルバーが収束しないか、計算した物理量が不正な場合。
        """
        if not np.isfinite(field_td) or field_td < 0:
            raise ValueError("E/Nは0以上の有限値にしてください。")
        if not np.isfinite(gas_temperature) or gas_temperature <= 0:
            raise ValueError("気体温度は正の有限値にしてください。")
        values = np.array(
            [densities.get(name, 0.0) for name in ("N2", "O2")], dtype=float
        )
        if (
            not np.isfinite(values).all()
            or np.any(values < 0)
            or values.sum() <= 0
        ):
            raise ValueError(
                "BOLOSに渡すN2/O2密度が不正、または両方ともゼロです。"
            )
        total = float(values.sum())
        fractions = values / total
        b = self.bsolver
        b.set_density("N2", fractions[0])
        b.set_density("O2", fractions[1])
        b.set_density(THREE_BODY_TARGET, fractions[1] * values[1] * 1e-6)
        b.kT = gas_temperature * solver.KB / solver.ELECTRONVOLT
        b.EN = field_td * solver.TOWNSEND
        try:
            # BOLOS evaluates unused zero-energy boundaries and tiny
            # tails.
            # Validate the resulting distribution/rates below instead of
            # allowing hundreds of identical floating-point warnings.
            with np.errstate(
                divide="ignore",
                invalid="ignore",
                over="ignore",
                under="ignore",
            ):
                if field_td == 0:
                    # Rebuild the energy grid when switching off. A hot
                    # grid
                    # and hot initial guess converge very slowly at zero
                    # field.
                    if self.last_field_td != 0:
                        b.grid = grid.QuadraticGrid(0, 15 * b.kT, self.cells)
                        guess = b.maxwell(b.kT)
                    else:
                        guess = self.eedf
                    b.init()
                elif self.eedf is None or self.last_field_td == 0:
                    b.grid = grid.LinearGrid(0, 60, self.cells)
                    b.init()
                    coarse = b.grid
                    f = b.converge(
                        b.maxwell(2.0),
                        maxn=self.max_iterations,
                        rtol=self.rtol,
                    )
                    b.grid = grid.QuadraticGrid(
                        0, 15 * b.mean_energy(f), self.cells
                    )
                    b.init()
                    guess = b.grid.interpolate(f, coarse)
                else:
                    b.init()
                    guess = self.eedf
                f, iterations, error = b.converge(
                    guess,
                    maxn=self.max_iterations,
                    rtol=self.rtol,
                    full=True,
                    # At zero field, the default pseudo-time increment
                    # can
                    # falsely accept an almost unchanged thermal initial
                    # guess.
                    delta0=1e20 if field_td == 0 else 1e14,
                )
                f, minimum, tail_repaired = normalize_eedf(b, f)
                raw_rates = np.array([b.rate(f, p) for p in self.processes])
                energy = float(b.mean_energy(f))
                reduced_mobility = float(b.mobility(f))
        except solver.ConvergenceError as exc:
            raise RuntimeError(
                f"BOLOSが収束しません：E/N={field_td:g} Td, "
                f"N2/O2={fractions.tolist()}, N={total:.6g} m^-3"
            ) from exc
        if (
            not np.isfinite(raw_rates).all()
            or np.any(raw_rates < -1e-12 * np.max(raw_rates))
            or not np.isfinite(energy)
            or energy <= 0
            or not np.isfinite(reduced_mobility)
            or reduced_mobility <= 0
        ):
            raise RuntimeError(
                "BOLOSのEEDF・反応係数・平均エネルギーが不正です。"
            )
        raw_rates = np.maximum(raw_rates, 0.0)
        self.eedf = f.copy()
        self.last_field_td = float(field_td)
        self.solve_count += 1
        rates = np.r_[0.0, raw_rates]
        rates[
            27
        ] *= 1e-6  # m^3 cm^3/s -> m^6/s; density appears in chemistry twice.
        return dict(
            rates=rates,
            energy_eV=energy,
            iterations=int(iterations),
            electron_reduced_mobility_1_V_m_s=reduced_mobility,
            error=float(error),
            uncorrected_eedf_min=minimum,
            target_density_m3=total,
            eedf_tail_repaired=tail_repaired,
            fractions={"N2": float(fractions[0]), "O2": float(fractions[1])},
        )
