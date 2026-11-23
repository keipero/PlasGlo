r"""振動準位の有無で、純 CO2 に最初の数サイクルを与えたときの効率を、E/N ごとに比べる。

各 E/N で、純 CO2 から数サイクル（パルス＋パルス間）を計算し、次を出す。

- 分解した CO2 の割合 X（振動準位の CO2 は、まだ CO2 として数える）。
- 振動に蓄えられたエネルギー（準位エネルギーの和）。
- 投入エネルギー（電子密度の時間積分と、EEDF から求めた入力電力）。
- 実効コスト = 投入エネルギー / 分解した CO2 の数 [eV/分子]。

使い方::

    python docs/co2_draft/nrp_vib_scan.py \
        outputs/co2_mechanism_vib/case_co2_vib.json \
        outputs/co2_mechanism/case_co2.json --td 30 50 100 200
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import add_vibration as vib  # noqa: E402
import nrp_efficiency as eff  # noqa: E402
import nrp_pulse_train as npt  # noqa: E402
from plasglo import load_case  # noqa: E402
from plasglo.electrons import mixture  # noqa: E402


def vib_energy(state):
    r"""振動準位に蓄えられたエネルギーの密度 [eV/m^3]。

    Parameters
    ----------
    state : dict of str to float
        化学種の密度 [m^-3]。

    Returns
    -------
    float
        準位エネルギーと密度の積の和。

    Examples
    --------
    >>> round(vib_energy({"CO2v1": 2.0, "CO2": 5.0}), 3)
    0.582
    """
    return sum(vib.energy(n) * state.get(vib.level(n), 0.0)
               for n in range(1, vib.N_LEVELS + 1))


def total_co2(state):
    r"""基底と振動準位を合わせた CO2 の密度 [m^-3]。

    Parameters
    ----------
    state : dict of str to float
        化学種の密度 [m^-3]。

    Returns
    -------
    float
        CO2 と ``CO2v*`` の合計。

    Examples
    --------
    >>> total_co2({"CO2": 3.0, "CO2v2": 1.0, "CO": 5.0})
    4.0
    """
    return state.get("CO2", 0.0) + sum(
        state.get(vib.level(n), 0.0) for n in range(1, vib.N_LEVELS + 1))


def run_cycle(case, state, args, swarm, charges):
    r"""1サイクルを計算し、投入エネルギーと終状態を返す。

    Parameters
    ----------
    case : plasglo.Case
        ケース。
    state : dict of str to float
        サイクル開始時の密度 [m^-3]。
    args : argparse.Namespace
        ``td``（この E/N）、``pulse_ns``、``period_us``、``seed_ne``。
    swarm : plasglo.electrons.mixture.ElectronSwarm
        入力電力の計算用。
    charges : dict of str to int
        化学種の電荷。

    Returns
    -------
    tuple
        (終状態, 投入エネルギー [eV/m^3])。
    """
    state = npt.ensure_seed(state, args.seed_ne)
    dens = {k: state.get(k, 0.0) for k in ("CO2", "CO", "O2")}
    dens["CO2"] += sum(state.get(vib.level(n), 0.0)
                       for n in range(1, vib.N_LEVELS + 1))
    power = (eff.efficiency_at(swarm, args.td, dens)["power"]
             if args.td > 0 else 0.0)
    pulse = npt.run_phase(case, state, args.td, args.pulse_ns, 41, 5.0,
                          rtol=args.rtol)
    names = [getattr(s, "name", s) for s in pulse.mechanism.species]
    ne = np.asarray(pulse.densities)[names.index("E")]
    t = np.asarray(pulse.time)
    n_target = sum(dens.values())
    energy = float(np.trapezoid(ne, t)) * n_target * power
    mid = npt.neutralize(npt.final_state(pulse), charges)
    inter_ns = args.period_us * 1e3 - args.pulse_ns
    edges = [e for e in npt.GAP_EDGES_NS if e < inter_ns] + [inter_ns]
    state = mid
    for lo, hi in zip(edges[:-1], edges[1:]):
        result = npt.run_phase(case, state, 0.0, hi - lo, 21, (hi - lo) / 4,
                               rtol=args.rtol)
        state = npt.neutralize(npt.final_state(result), charges)
    return state, energy


def main():
    r"""E/N ごとに、振動あり・なしの効率を計算して表示し、JSON に書く。"""
    logging.disable(logging.WARNING)
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("case_vib")
    ap.add_argument("case_novib")
    ap.add_argument("--td", type=float, nargs="+", default=[50, 100, 200])
    ap.add_argument("--cycles", type=int, default=5)
    ap.add_argument("--pulse-ns", type=float, default=20.0)
    ap.add_argument("--period-us", type=float, default=100.0)
    ap.add_argument("--seed-ne", type=float, default=1e20)
    ap.add_argument("--out", default="vib_scan.json")
    ap.add_argument(
        "--atol", type=float, default=1e10,
        help="密度の絶対許容誤差 [m^-3]（振動準位は密度が小さく、負値が出る）")
    ap.add_argument("--rtol", type=float, default=1e-9)
    args = ap.parse_args()
    npt.ATOL_DENSITY = args.atol
    results = []
    for label, path in (("振動あり", args.case_vib), ("振動なし", args.case_novib)):
        case = load_case(path)
        s = case.electron_settings()
        p = s["profile_data"]
        swarm = mixture.ElectronSwarm(
            p["cross_sections"], p["targets"], p["channels"], cells=s["cells"])
        species = case.data["mechanism"]["species"]
        names = [x["name"] for x in species]
        charges = {x["name"]: x["charge"] for x in species}
        n0 = case.data["initial_densities_m3"]["CO2"]
        # 電子の種（CO2 -> CO2^+ + E）の再結合は、電界がなくても CO を作る。
        # 電子 1 個あたり CO が 1 個できる人為的な分解なので、電界 0 の
        # 1 サイクルで基準を求め、あとで引く。
        zero = argparse.Namespace(**{**vars(args), "td": 0.0})
        base = {n: 0.0 for n in names}
        base["CO2"] = n0
        base, _ = run_cycle(case, base, zero, swarm, charges)
        x_seed = 1 - total_co2(base) / n0
        print(f"{label}: 電子の種による基準の X = {100 * x_seed:.5f}%/サイクル",
              flush=True)
        for td in args.td:
            args.td_current = td
            run = argparse.Namespace(**{**vars(args), "td": td})
            state = {n: 0.0 for n in names}
            state["CO2"] = n0
            total, rows, t0 = 0.0, [], time.time()
            for k in range(1, args.cycles + 1):
                state, e = run_cycle(case, state, run, swarm, charges)
                total += e
                x = 1 - total_co2(state) / n0 - k * x_seed
                co_made = state["CO"] + state.get("C", 0.0)
                v_e = vib_energy(state)
                cost = total / (x * n0) if x > 0 else float("inf")
                rows.append(dict(cycle=k, x=x, energy=total, vib=v_e,
                                 cost=cost))
                print(f"{label} E/N={td:5.0f} Td 周期{k:2d}: X={100 * x:7.4f}% "
                      f"E_vib={v_e:.2e} 投入={total:.2e} eV/m3 "
                      f"コスト={cost:7.2f} eV/分子 CO={co_made:.2e} "
                      f"({time.time() - t0:.0f}s)", flush=True)
            results.append(dict(model=label, td=td, rows=rows))
    Path(args.out).write_text(json.dumps(results, ensure_ascii=False),
                              encoding="utf-8")


if __name__ == "__main__":
    main()
