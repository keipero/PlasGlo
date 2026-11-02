r"""共通経路の窒素計算を再実行したZDPlasKinで検証する。
"""

from importlib import import_module

from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np

MODEL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODEL))
load_case = import_module("py3zdplaschem").load_case
simulate = import_module("py3zdplaschem").simulate
save_result = import_module("py3zdplaschem").save_result
read_csv = import_module("py3zdplaschem.io.reference_data").read_csv
native_control = import_module(
    "py3zdplaschem.validation.zdplaskin"
).native_control


def metrics(names, predicted, reference):
    r"""粒子種ごとのピーク密度で規格化した時間履歴の差を計算する。

    Parameters
    ----------
    names : sequence of str
        密度配列の行と同じ順序の粒子種名。
    predicted : numpy.ndarray
        形状が (粒子種数, 時刻数) の予測密度 [m^-3]。
    reference : numpy.ndarray
        形状が (粒子種数, 時刻数) の参照密度 [m^-3]。

    Returns
    -------
    list of dict
        粒子種ごとの絶対誤差とピーク規格化誤差。
    """
    peaks = np.max(np.abs(reference), axis=1)
    scale = np.where(peaks > 0, peaks, 1.0)
    error = predicted - reference
    return [
        dict(
            species=name,
            reference_peak_m3=float(peaks[i]),
            max_absolute_difference_m3=float(np.max(np.abs(error[i]))),
            max_error_over_reference_peak=float(
                np.max(np.abs(error[i])) / scale[i]
            ),
            nrmse_over_reference_peak=float(
                np.sqrt(np.mean(error[i] ** 2)) / scale[i]
            ),
        )
        for i, name in enumerate(names)
    ]


def main():
    r"""共通経路の窒素計算を再実行したZDPlasKinの対照と比較する。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。
    AssertionError
        必要な数値比較または出典照合に失敗した場合。

    Examples
    --------
    プロジェクトのルートから実行します::

        python verification/validate_nitrogen.py
    """
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--reference",
        type=Path,
        default=MODEL.parent.parent
        / "ZDplasKin_analysis/runs/nitrogen_1000Td_10ns/20261021_234859_884",
    )
    cli.add_argument(
        "--output",
        type=Path,
        default=MODEL / "outputs/verification/gas_generalization/nitrogen",
    )
    cli.add_argument("--cells", type=int, default=1600)
    args = cli.parse_args()
    out, reference = args.output.resolve(), args.reference.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sources = [
        reference / name
        for name in (
            "kinet.inp",
            "bolsigdb.dat",
            "user_code.F90",
            "zdplaskin_m.F90",
            "bolsig_x86_64_g.dll",
            "densities.csv",
        )
    ]

    def hashes():
        r"""比較の前後にネイティブ計算の元ファイルをハッシュ化する。

        Returns
        -------
        dict
            元ファイルの絶対パスをキー、SHA-256を値とする辞書。
        """
        return {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        }

    before = hashes()
    print("Running the common nitrogen case...", flush=True)
    case = load_case("nitrogen").with_updates(electrons={"cells": args.cells})
    result = simulate(case)
    names = result.mechanism.species
    assert len(names) == 21 and len(result.mechanism.reactions) == 96
    save_result(
        result,
        out / "python",
        species=["E", "N2", "N2^+", "N2(C3)", "N4^+", "N"],
        balance_species=["N2(C3)"],
    )
    supplied_names, supplied_time, supplied = read_csv(
        reference / "densities.csv"
    )
    print(
        "Building/running ZDPlasKin with the original electron-electron"
        " settings...",
        flush=True,
    )
    rebuilt_names, rebuilt_time, rebuilt = native_control(
        reference, out, disable_ee=False
    )
    print(
        "Building/running ZDPlasKin with electron-electron collisions"
        " disabled...",
        flush=True,
    )
    matched_names, matched_time, matched = native_control(
        reference, out, disable_ee=True
    )
    for other in (supplied_names, rebuilt_names, matched_names):
        if [name.upper() for name in other] != names:
            raise ValueError("Native/Python species order mismatch.")
    for time in (supplied_time, rebuilt_time, matched_time):
        np.testing.assert_allclose(result.time, time, rtol=1e-12, atol=1e-25)
    comparisons = dict(
        python_vs_ee_off=metrics(names, result.densities, matched),
        python_vs_original=metrics(names, result.densities, supplied),
        rebuilt_original_vs_supplied=metrics(names, rebuilt, supplied),
    )
    maximum = {
        key: max(row["max_error_over_reference_peak"] for row in rows)
        for key, rows in comparisons.items()
    }
    assert hashes() == before, "Original reference inputs or outputs changed."
    report = dict(
        cells=args.cells,
        species_count=21,
        reaction_count=96,
        time_points=len(result.time),
        conditions=dict(
            field_td=1000.0,
            temperature_K=300.0,
            pressure_Pa=101325.0,
            end_ns=1.0,
            initial_electron_density_m3=1e14,
            initial_ion_density_m3=1e14,
        ),
        source_sha256=before,
        original_files_unchanged=True,
        max_error_over_reference_peak=maximum,
        comparisons=comparisons,
        atom_error=result.diagnostics["max_atom_relative_error"],
        net_charge_error_m3=result.diagnostics["max_net_charge_error_m3"],
        comparison_note=(
            "The e-e-off native run is the matched comparison. The original"
            " run enables e-e collisions automatically; its difference is"
            " reported separately. No time shifts, coefficient fitting or"
            " density rescaling."
        ),
        limitations=result.diagnostics["limitations"],
    )
    # Use the before-refactor baseline captured for this task: identical
    # source
    # conditions and 400 cells. Older notebooks may use a different seed
    # density.
    baseline = (
        MODEL / "outputs/verification/gas_generalization/baseline/nitrogen.npz"
    )
    if baseline.is_file():
        regression = (
            result
            if args.cells == 400
            else simulate(case.with_updates(electrons={"cells": 400}))
        )
        with np.load(baseline) as old:
            np.testing.assert_array_equal(old["species"], np.array(names))
            np.testing.assert_array_equal(old["time_s"], regression.time)
            scale = np.maximum(
                np.max(np.abs(old["densities_m3"]), axis=1), 1.0
            )
            delta = float(
                np.max(
                    np.abs(regression.densities - old["densities_m3"])
                    / scale[:, None]
                )
            )
        report["refactor_regression_cells"] = 400
        report["refactor_error_over_previous_species_peak"] = delta
        if delta > 1e-10:
            raise AssertionError(
                f"Nitrogen changed during generalization: {delta}"
            )
    for key, rows in comparisons.items():
        with (out / f"{key}.csv").open(
            "w", encoding="utf-8", newline=""
        ) as fp:
            writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), layout="constrained")
    for name in ("E", "N2^+", "N2(C3)", "N4^+"):
        i = names.index(name)
        (line,) = axes[0].semilogy(
            result.time * 1e9,
            np.maximum(result.densities[i], 1.0),
            label=f"{name}: Python",
        )
        axes[0].semilogy(
            result.time * 1e9,
            np.maximum(matched[i], 1.0),
            "--",
            color=line.get_color(),
            label=f"{name}: ZD e-e off",
        )
        scale = max(np.max(np.abs(matched[i])), 1.0)
        axes[1].plot(
            result.time * 1e9,
            100 * (result.densities[i] - matched[i]) / scale,
            label=name,
        )
    axes[0].set(
        ylabel="Density [m$^{-3}$]",
        title=(
            f"Nitrogen: common runner vs native ZDPlasKin ({args.cells} cells)"
        ),
    )
    axes[1].set(
        xlabel="Time [ns]", ylabel="Difference / reference species peak [%]"
    )
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
    fig.savefig(out / "comparison.png", dpi=160)
    plt.close(fig)
    (out / "report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    lines = [
        "# 窒素：共通計算経路とZDPlasKinの比較",
        "",
        (
            "21種・96反応、1000 Td、300 K、101325 Pa、1"
            f" ns、{len(result.time)}時刻、BOLOS {args.cells}セル。"
        ),
        "",
        (
            "密度はm⁻³に統一し、粒子種・時刻を照合。"
            "最大差は各粒子種の参照ピーク密度で規格化しています。"
        ),
        "",
        "| 比較 | 全粒子種・全時刻の最大差 [%] |",
        "|---|---:|",
    ]
    labels = dict(
        python_vs_ee_off="Python対ZDPlasKin（電子間衝突なし）",
        python_vs_original="Python対元のZDPlasKin結果",
        rebuilt_original_vs_supplied="元設定の再実行対保存済み結果",
    )
    lines += [
        f"| {labels[key]} | {100 * value:.8g} |"
        for key, value in maximum.items()
    ]
    lines += [
        "",
        (
            "電子間衝突なしの対照がモデル条件を合わせた比較です。"
            "元の計算は電子間衝突が自動的に有効になるため、別に示しています。"
        ),
        (
            "元入力・元出力のハッシュが比較前後で一致しました。"
            "係数・密度の調整や時刻シフトは行っていません。"
        ),
        (
            "この比較は指定条件での数値検証であり、"
            "他のガス種・条件の物理的妥当性を保証しません。"
        ),
        "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(
        json.dumps(
            {
                key: value
                for key, value in report.items()
                if key not in ("comparisons", "source_sha256")
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
