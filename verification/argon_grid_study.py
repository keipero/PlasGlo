r"""Arサンプルのセル数・上限エネルギー・反復条件の影響を分離する。
"""

from importlib import import_module

from pathlib import Path
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp

MODEL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODEL))
VERIFICATION_OUTPUT_DIR = import_module(
    "py3zdplaschem.paths"
).VERIFICATION_OUTPUT_DIR
SAMPLE_OUTPUT_DIR = import_module("py3zdplaschem.paths").SAMPLE_OUTPUT_DIR
grid = import_module("py3zdplaschem.electrons.bolos_backend").grid
parser = import_module("py3zdplaschem.electrons.bolos_backend").parser
solver = import_module("py3zdplaschem.electrons.bolos_backend").solver
normalize_eedf = import_module(
    "py3zdplaschem.electrons.bolos_backend"
).normalize_eedf
SampleMechanism = import_module(
    "py3zdplaschem.chemistry.sample_mechanism"
).SampleMechanism
read_csv = import_module("py3zdplaschem.io.reference_data").read_csv

SAMPLE = MODEL / "samples/gases/argon"
OUT = VERIFICATION_OUTPUT_DIR


def swarm(cells, maximum_eV=None, rtol=1e-6):
    r"""アルゴンのエネルギー格子依存性を1条件分計算する。

    Parameters
    ----------
    cells : int, optional
        BOLOSで使用するエネルギー格子のセル数。
    maximum_eV : float or None, optional
        固定した上限エネルギー [eV]。Noneの場合は適応的に決める。
    rtol : float, optional
        収束判定または時間積分の相対許容誤差。

    Returns
    -------
    dict
        格子条件、平均エネルギー、電離係数、収束診断。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    """
    with (SAMPLE / "bolsigdb.dat").open(encoding="utf-8") as handle:
        collisions = parser.parse(handle)
    b = solver.BoltzmannSolver(grid.LinearGrid(0, 60, cells))
    processes = b.load_collisions(collisions)
    b.set_density("Ar", 1.0)
    b.EN = 51 * solver.TOWNSEND
    b.kT = 301 * solver.KB / solver.ELECTRONVOLT
    with np.errstate(all="ignore"):
        b.init()
        f = b.converge(b.maxwell(2), maxn=1000, rtol=rtol)
        old_grid = b.grid
        maximum = 15 * b.mean_energy(f) if maximum_eV is None else maximum_eV
        b.grid = grid.QuadraticGrid(0, maximum, cells)
        b.init()
        f, iterations, error = b.converge(
            b.grid.interpolate(f, old_grid), maxn=1000, rtol=rtol, full=True
        )
        f, minimum, repaired = normalize_eedf(b, f)
        rates = {
            f"{p.target_name}->{p.product}".upper(): float(b.rate(f, p))
            for p in processes
            if p.product
        }
    if not np.isfinite(list(rates.values())).all() or any(
        value < 0 for value in rates.values()
    ):
        raise RuntimeError("Invalid coefficients")
    return dict(
        cells=cells,
        maximum_eV=float(maximum),
        rtol=rtol,
        ionization_m3_s=rates["AR->AR^+"],
        mean_energy_eV=float(b.mean_energy(f)),
        error=float(error),
        iterations=int(iterations),
        uncorrected_eedf_min=minimum,
        tail_repaired=repaired,
    )


def density(row, settings, names, time):
    r"""指定した電子群係数でアルゴンの化学反応を積分する。

    Parameters
    ----------
    row : dict
        平均エネルギーと電離係数を含む電子群計算結果。
    settings : dict
        サンプル条件と初期密度のメタデータ。
    names : sequence of str
        密度配列の行と同じ順序の粒子種名。
    time : numpy.ndarray
        厳密に単調増加する出力時刻 [s]。

    Returns
    -------
    numpy.ndarray
        形状が (粒子種数, 時刻数) の粒子数密度 [m^-3]。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    """
    mechanism = SampleMechanism(
        SAMPLE / "kinet.inp",
        names,
        settings["conditions"],
        {"AR->AR^+": row["ionization_m3_s"]},
        row["mean_energy_eV"],
        expected_counts=(3, 2),
    )
    initial = np.array(
        [settings["initial_densities_cm3"].get(s, 0) * 1e6 for s in names]
    )
    solution = solve_ivp(
        mechanism.rhs,
        (0, time[-1]),
        initial / mechanism.n0,
        method="BDF",
        jac=mechanism.jac,
        t_eval=time,
        rtol=1e-10,
        atol=1e-5 / mechanism.n0,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y * mechanism.n0


def main():
    r"""アルゴンの格子数・上限エネルギー・収束条件の感度を分けて調べる。

    Examples
    --------
    プロジェクトのルートから実行します::

        python verification/argon_grid_study.py
    """
    settings = json.loads((SAMPLE / "case.json").read_text(encoding="utf-8"))
    raw_names, time, native = read_csv(OUT / "argon_native/densities.csv")
    names = [s.upper() for s in raw_names]
    electron = names.index("E")
    old_rows = json.loads(
        (OUT / "grid_resolution.json").read_text(encoding="utf-8")
    )["argon"]["rows"]
    results = {
        "adaptive_cutoff_cells": old_rows.copy(),
        "fixed_cutoff": [],
        "cutoff_checks": [],
        "iteration_checks": [],
    }
    for cells in [3200, 6400]:
        print(f"Adaptive cutoff, {cells} cells", flush=True)
        results["adaptive_cutoff_cells"].append(swarm(cells))
    for cells in [1600, 3200, 6400]:
        print(f"Fixed 100 eV, {cells} cells, rtol=1e-8", flush=True)
        results["fixed_cutoff"].append(swarm(cells, 100.0, 1e-8))
    # For a quadratic grid, N proportional to sqrt(Emax) approximately
    # preserves cell widths near a fixed threshold while changing the
    # cutoff.
    for maximum in [80.0, 140.0]:
        cells = round(6400 * np.sqrt(maximum / 100.0))
        print(f"Cutoff {maximum:g} eV, {cells} cells", flush=True)
        results["cutoff_checks"].append(swarm(cells, maximum, 1e-8))
    results["iteration_checks"].append(swarm(6400, 100.0, 1e-6))
    native_peak = float(np.max(native[electron]))
    for group, rows in results.items():
        previous = None
        for row in rows:
            predicted = density(row, settings, names, time)
            row["electron_peak_m3"] = float(predicted[electron].max())
            row["electron_error_over_native_peak"] = float(
                np.max(np.abs(predicted[electron] - native[electron]))
                / native_peak
            )
            if previous is not None:
                row["ionization_relative_change_from_previous"] = (
                    row["ionization_m3_s"] / previous["ionization_m3_s"] - 1
                )
                row[
                    "electron_trajectory_change_from_previous_over_fine_peak"
                ] = float(
                    np.max(
                        np.abs(
                            predicted[electron] - previous_density[electron]
                        )
                    )
                    / predicted[electron].max()
                )
            previous, previous_density = row, predicted
    (OUT / "argon_grid_study.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    for group, rows in results.items():
        print(group, flush=True)
        for row in rows:
            print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
