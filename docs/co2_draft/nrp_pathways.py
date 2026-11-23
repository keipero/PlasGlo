r"""飽和したパルス列の状態から1サイクル計算し、化学種ごとの反応経路を出す。

``nrp_pulse_train.py`` の記録の最後の状態から、直接計算を1サイクル行い、
指定した化学種ごとに、全反応の積算寄与（パルス中＋パルス間）を出す。
記録（上位6反応のみ）と違い、全反応を対象にする。

使い方::

    python docs/co2_draft/nrp_pathways.py outputs/co2_nrp/loki_1e20 \
        --case outputs/co2_mechanism/case_co2.json --seed-ne 1e20
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nrp_pulse_train as npt  # noqa: E402
from plasglo import load_case  # noqa: E402

SPECIES = ("CO2", "CO", "O2", "O3", "O", "O(1D)", "CO(a3Pi)", "CO2(e1)")


def cycle_pathways(case, state, args, species):
    r"""1サイクル（パルス＋パルス間）を計算し、化学種ごとの寄与を足し合わせる。

    Parameters
    ----------
    case : plasglo.Case
        ケース。
    state : dict of str to float
        サイクル開始時の密度 [m^-3]。
    args : argparse.Namespace
        ``td``、``pulse_ns``、``period_us``、``seed_ne``、``pulse_update_ns``。
    species : tuple of str
        解析する化学種。

    Returns
    -------
    dict
        化学種から、``npt.summarize`` と同じ形（ただし上位は全反応）の辞書。
    """
    charges = {s["name"]: s["charge"]
               for s in case.data["mechanism"]["species"]}
    state = npt.ensure_seed(state, args.seed_ne)
    pulse = npt.run_phase(case, state, args.td, args.pulse_ns, 81,
                          args.pulse_update_ns)
    parts = {sp: [npt.contributions(pulse, sp)] for sp in species}
    state = npt.neutralize(npt.final_state(pulse), charges)
    inter_ns = args.period_us * 1e3 - args.pulse_ns
    edges = [e for e in npt.GAP_EDGES_NS if e < inter_ns] + [inter_ns]
    for lo, hi in zip(edges[:-1], edges[1:]):
        result = npt.run_phase(case, state, 0.0, hi - lo, 21, (hi - lo) / 4)
        for sp in species:
            parts[sp].append(npt.contributions(result, sp))
        state = npt.neutralize(npt.final_state(result), charges)
    return {sp: npt.summarize(p, top=10_000) for sp, p in parts.items()}


def main():
    r"""飽和状態の1サイクルの、化学種ごとの反応経路を表示・保存する。"""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", help="nrp_pulse_train.py の出力フォルダ")
    ap.add_argument("--case", required=True)
    ap.add_argument("--td", type=float, default=100.0)
    ap.add_argument("--pulse-ns", type=float, default=20.0)
    ap.add_argument("--period-us", type=float, default=100.0)
    ap.add_argument("--seed-ne", type=float, default=1e20)
    ap.add_argument("--pulse-update-ns", type=float, default=5.0)
    ap.add_argument("--bolos-cache-tol", type=float, default=1e-3)
    ap.add_argument("--top", type=int, default=8)
    args = ap.parse_args()
    if args.bolos_cache_tol > 0:
        import bolos_cache

        bolos_cache.install(args.bolos_cache_tol)
    path = Path(args.run) / "cycles.jsonl"
    with path.open(encoding="utf-8-sig") as fp:
        records = [json.loads(line) for line in fp if line.strip()]
    last = [r for r in records if "pulse" in r][-1]
    case = npt.apply_variant(load_case(args.case), "loki", Path(args.run))
    names = [s["name"] for s in case.data["mechanism"]["species"]]
    state = {n: last["state_end_gap"].get(n, 0.0) for n in names}
    out = cycle_pathways(case, state, args, SPECIES)
    (Path(args.run) / "pathways_saturated.json").write_text(
        json.dumps(dict(cycle=last["cycle"], species=out),
                   ensure_ascii=False, indent=1),
        encoding="utf-8",
    )
    print(f"記録の最後（サイクル {last['cycle']}）の状態から1サイクル計算\n")
    for sp, s in out.items():
        prod = sum(i["m3"] for i in s["production"])
        loss = sum(i["m3"] for i in s["loss"])
        print(f"== {sp}: 生成 {prod:.3e}, 消滅 {loss:.3e}, "
              f"正味 {s['change']:+.3e} m^-3")
        for label, key, total in (("生成", "production", prod),
                                  ("消滅", "loss", loss)):
            for item in s[key][:args.top]:
                print(f"  {label} {item['id']} {item['eq']:34s} "
                      f"{item['m3']:.3e} ({100 * item['m3'] / total:5.1f}%)")


if __name__ == "__main__":
    main()
