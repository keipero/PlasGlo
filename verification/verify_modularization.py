r"""再編後の計算を同じ環境で保存した再編前の密度履歴と比較する。
"""

from importlib import import_module

from pathlib import Path
import hashlib
import json
import sys

import numpy as np

MODEL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODEL))
simulate = import_module("plasglo").simulate
VERIFICATION_OUTPUT_DIR = import_module(
    "plasglo.paths"
).VERIFICATION_OUTPUT_DIR
run_case = import_module("samples.gases.run_samples").run_case

OUT = VERIFICATION_OUTPUT_DIR / "modularization"
BASELINE = OUT / "baseline"


def compare(name, names, time, density):
    r"""計算した時間履歴を保存済みの数値基準と比較する。

    Parameters
    ----------
    name : str
        この処理で使用する粒子種・ケース・出力などの識別名。
    names : sequence of str
        密度配列の行と同じ順序の粒子種名。
    time : numpy.ndarray
        厳密に単調増加する出力時刻 [s]。
    density : numpy.ndarray
        基準と比較する、形状が (粒子種数, 時刻数) の実密度 [m^-3]。

    Returns
    -------
    dict
        基準データのパス、ピーク規格化最大誤差、合否を含む辞書。

    Raises
    ------
    AssertionError
        必要な数値比較または出典照合に失敗した場合。
    """
    baseline = BASELINE / f"{name}.npz"
    current = (
        VERIFICATION_OUTPUT_DIR / "gas_generalization/baseline" / f"{name}.npz"
    )
    if name == "nitrogen" and current.is_file():
        baseline = current
    with np.load(baseline) as old:
        np.testing.assert_array_equal(
            np.char.upper(old["species"]), np.char.upper(np.array(names))
        )
        np.testing.assert_allclose(time, old["time_s"], rtol=1e-14, atol=1e-25)
        scale = np.maximum(np.max(np.abs(old["densities_m3"]), axis=1), 1.0)
        peak_error = float(
            np.max(
                np.max(np.abs(density - old["densities_m3"]), axis=1) / scale
            )
        )
        if peak_error > 1e-10:
            raise AssertionError(
                f"{name} changed after module split: {peak_error}"
            )
    return dict(
        baseline=str(baseline),
        max_difference_over_species_peak=peak_error,
        unchanged_within_tolerance=True,
    )


manifest = json.loads((BASELINE / "manifest.json").read_text(encoding="utf-8"))
for relative, expected in manifest["input_hashes"].items():
    assert (
        hashlib.sha256((MODEL / relative).read_bytes()).hexdigest() == expected
    ), relative
archive = MODEL / "archive/current_diagnostics_2026-11-01"
for item in json.loads(
    (archive / "source_manifest.json").read_text(encoding="utf-8")
):
    assert (
        hashlib.sha256((archive / item["archived"]).read_bytes()).hexdigest()
        == item["sha256"]
    )
assert not (MODEL / "current_diagnostics.py").exists()
assert not (MODEL / "test_current_diagnostics.py").exists()

report = dict(
    inputs_unchanged=True,
    archived_source_hashes_verified=True,
    current_feature_detached=True,
    cases={},
)
print("Air pulse...", flush=True)
air = simulate(**manifest["air_settings"])
report["cases"]["air"] = compare(
    "air", air.mechanism.species, air.time, air.densities
)
print("Argon sample...", flush=True)
names, time, density, _, _ = run_case("argon")
report["cases"]["argon"] = compare("argon", names, time, density)
print("Nitrogen common case (400 cells)...", flush=True)
names, time, density, _, _ = run_case("nitrogen")
report["cases"]["nitrogen"] = compare("nitrogen", names, time, density)
# Exercise previous imports to ensure they re-export exactly the
# canonical types.
legacy_mechanism = import_module("kinetics").Mechanism
Mechanism = import_module("plasglo.chemistry.mechanism").Mechanism
legacy_swarm = import_module("bolos_coupling").BolosSwarm
BolosSwarm = import_module("plasglo.electrons.bolos_backend").BolosSwarm

assert legacy_mechanism is Mechanism and legacy_swarm is BolosSwarm
report["legacy_imports_compatible"] = True
(OUT / "report.json").write_text(
    json.dumps(report, indent=2), encoding="utf-8"
)
print(json.dumps(report, indent=2))
