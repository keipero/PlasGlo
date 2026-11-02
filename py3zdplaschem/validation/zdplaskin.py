r"""通常の計算から分離した、条件をそろえたネイティブ比較。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from ..constants import KB, KB_EV
from ..paths import OUTPUTS_DIR
from ..solvers.simulation import simulate
from ..chemistry.sample_mechanism import SampleMechanism
from ..electrons.sample_swarm import solve_sample_swarm
from ..io.reference_data import read_csv, driver_conditions


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
    rates, swarm, collisions = solve_sample_swarm(
        path / "bolsigdb.dat",
        conditions["field_td"],
        conditions["temperature_k"],
    )
    model = SampleMechanism(
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

    Examples
    --------
    プロジェクトのルートから実行します::

        python verification/zdplaskin.py
    """
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("reference_csv", type=Path)
    cli.add_argument(
        "--output",
        type=Path,
        default=OUTPUTS_DIR / "validation/zdplaskin_nitrogen",
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
