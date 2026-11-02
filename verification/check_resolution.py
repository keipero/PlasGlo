r"""窒素・アルゴンの電離係数のエネルギー格子依存性を調べる。
"""

from importlib import import_module

from pathlib import Path
import json
import sys

MODEL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODEL))
VERIFICATION_OUTPUT_DIR = import_module(
    "plasglo.paths"
).VERIFICATION_OUTPUT_DIR
SAMPLE_OUTPUT_DIR = import_module("plasglo.paths").SAMPLE_OUTPUT_DIR
solve_sample_swarm = import_module(
    "plasglo.electrons.sample_swarm"
).solve_sample_swarm


def main():
    r"""エネルギー格子を細かくしたときのサンプル係数の変化を調べる。

    Examples
    --------
    プロジェクトのルートから実行します::

        python verification/check_resolution.py
    """
    results = {}
    for sample, target, field, temperature, channel in [
        ("argon", "Ar", 51.0, 301.0, "AR->AR^+"),
        ("nitrogen", "N2", 1000.0, 300.0, "N2->N2^+"),
    ]:
        rows = []
        for cells in [400, 800, 1600]:
            rates, state, _ = solve_sample_swarm(
                MODEL / "samples/gases" / sample / "bolsigdb.dat",
                field,
                temperature,
                cells=cells,
                target=target,
            )
            rows.append(
                dict(cells=cells, ionization_m3_s=rates[channel], **state)
            )
        results[sample] = dict(
            rows=rows,
            ionization_change_400_to_1600=rows[-1]["ionization_m3_s"]
            / rows[0]["ionization_m3_s"]
            - 1,
            mean_energy_change_400_to_1600=rows[-1]["mean_energy_eV"]
            / rows[0]["mean_energy_eV"]
            - 1,
        )
    output = VERIFICATION_OUTPUT_DIR
    output.mkdir(parents=True, exist_ok=True)
    (output / "grid_resolution.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
