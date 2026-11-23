r"""飽和状態での、反応速度に対する感度（スクリーニング）を求める。

``nrp_pulse_train.py`` の最後の状態（飽和付近）から、1サイクルを直接計算する。
反応ごとに速度を ``factor`` 倍（と 1/``factor`` 倍）して、1サイクルの
分解率の増分 ``g`` の変化 ``dg`` を求める。飽和は ``g = 0`` の点なので、
飽和時の分解率の変化の目安は ``dg / |dg/dX|`` である。``dg/dX`` は、記録の
増分 ``dX`` と ``X`` の傾き（飽和への接近）から求める。

これは、遅い種（CO、O2、O3）の再調整を含まない、1次の近似である。上位の
反応は、飽和まで走らせ直して確かめる。電子衝突（断面積）の反応は、速度を
倍率で変えられないので、対象外である。

使い方::

    python docs/co2_draft/nrp_sensitivity.py outputs/co2_nrp/loki_1e20 \
        --case outputs/co2_mechanism/case_co2.json
"""

import argparse
import csv
import json
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bolos_cache  # noqa: E402
import nrp_pulse_train as npt  # noqa: E402
from plasglo import load_case  # noqa: E402


def slope_dg_dx(records, x_lo=0.50):
    r"""増分 dX と X の傾き dg/dX を、記録から直線回帰で求める。

    Parameters
    ----------
    records : list of dict
        ``cycles.jsonl`` の記録（``x_end`` と ``dx`` を持つもの）。
    x_lo : float, optional
        回帰に使う X の下限。

    Returns
    -------
    float
        ``dg/dX``（飽和へ近づくときは負）。

    Examples
    --------
    >>> r = [{"x_end": 0.5 + 0.01 * i, "dx": 1e-3 - 1e-4 * i}
    ...      for i in range(6)]
    >>> round(slope_dg_dx(r), 6)
    -0.01
    """
    pts = [(r["x_end"], r["dx"]) for r in records
           if "dx" in r and "x_end" in r and r["x_end"] >= x_lo]
    x, g = np.array(pts).T
    return float(np.polyfit(x, g, 1)[0])


def one_cycle(case, state, args, n0):
    r"""1サイクルを計算し、分解率の増分と終状態を返す。

    Parameters
    ----------
    case : plasglo.Case
        ケース。
    state : dict of str to float
        サイクル開始時の密度 [m^-3]。
    args : argparse.Namespace
        ``td``、``pulse_ns``、``period_us``、``seed_ne``。
    n0 : float
        初期の CO2 密度 [m^-3]。

    Returns
    -------
    tuple
        (分解率の増分 g、終状態)。
    """
    species = case.data["mechanism"]["species"]
    charges = {s["name"]: s["charge"] for s in species}
    start = state["CO2"]
    state = npt.ensure_seed(state, args.seed_ne)
    pulse = npt.run_phase(case, state, args.td, args.pulse_ns, 41, 5.0)
    state = npt.neutralize(npt.final_state(pulse), charges)
    inter_ns = args.period_us * 1e3 - args.pulse_ns
    edges = [e for e in npt.GAP_EDGES_NS if e < inter_ns] + [inter_ns]
    for lo, hi in zip(edges[:-1], edges[1:]):
        result = npt.run_phase(case, state, 0.0, hi - lo, 21, (hi - lo) / 4)
        state = npt.neutralize(npt.final_state(result), charges)
    return (start - state["CO2"]) / n0, state


def scaled_case(src, rows, header, rid, factor, folder):
    r"""1つの反応の速度を ``factor`` 倍したケースを書き、読み込む。

    Parameters
    ----------
    src : dict
        元のケース（JSON を読んだもの）。
    rows : list of dict
        機構の CSV の行。
    header : list of str
        CSV の列名。
    rid : str
        倍率をかける反応ID。
    factor : float
        倍率。
    folder : pathlib.Path
        一時ファイルの置き場所。

    Returns
    -------
    plasglo.Case
        倍率をかけたケース。
    """
    new = [dict(r) for r in rows]
    for r in new:
        if r["id"] == rid:
            r["value"] = f"({factor})*({r['value']})"
    csv_path = folder / "mechanism_scaled.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=header)
        writer.writeheader()
        writer.writerows(new)
    case = json.loads(json.dumps(src))
    case["mechanism"]["reactions"] = csv_path.name
    path = folder / "case_scaled.json"
    path.write_text(json.dumps(case, ensure_ascii=False), encoding="utf-8")
    return load_case(path)


def main():
    r"""反応ごとの感度を計算し、表と JSON に書く。"""
    logging.disable(logging.WARNING)
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run")
    ap.add_argument("--case", required=True)
    ap.add_argument("--factor", type=float, default=2.0)
    ap.add_argument("--td", type=float, default=100.0)
    ap.add_argument("--pulse-ns", type=float, default=20.0)
    ap.add_argument("--period-us", type=float, default=100.0)
    ap.add_argument("--seed-ne", type=float, default=1e20)
    ap.add_argument("--only", nargs="*", help="対象の反応ID（省略で全部）")
    args = ap.parse_args()
    bolos_cache.install(1e-3)
    run = Path(args.run)
    folder = run / "sens_tmp"
    folder.mkdir(exist_ok=True)
    src_path = Path(args.case)
    src = json.loads(src_path.read_text("utf-8-sig"))
    csv_path = src_path.parent / src["mechanism"]["reactions"]
    with csv_path.open(encoding="utf-8-sig") as fp:
        rows = list(csv.DictReader(fp))
        header = list(rows[0])
    with (run / "cycles.jsonl").open(encoding="utf-8-sig") as fp:
        records = [json.loads(line) for line in fp if line.strip()]
    last = [r for r in records if "pulse" in r][-1]
    n0 = src["initial_densities_m3"]["CO2"]
    names = [s["name"] for s in src["mechanism"]["species"]]
    state0 = {n: last["state_end_gap"].get(n, 0.0) for n in names}
    slope = slope_dg_dx([r for r in records if "dx" in r])
    print(f"dg/dX = {slope:.3e}（サイクル {last['cycle']} 付近）", flush=True)
    base_case = scaled_case(src, rows, header, "", 1.0, folder)
    g0, _ = one_cycle(base_case, dict(state0), args, n0)
    print(f"基準 g0 = {g0:.4e}/サイクル", flush=True)
    results = []
    targets = [r for r in rows if r["kind"] == "expression"
               and (not args.only or r["id"] in args.only)]
    for r in targets:
        entry = dict(id=r["id"], equation=r["equation"], g0=g0)
        for label, f in (("up", args.factor), ("down", 1 / args.factor)):
            case = scaled_case(src, rows, header, r["id"], f, folder)
            g, _ = one_cycle(case, dict(state0), args, n0)
            entry[label] = g
            entry[f"dx_{label}"] = (g - g0) / abs(slope)
        results.append(entry)
        print(f"{r['id']} {r['equation']:34s} ×{args.factor:g}: "
              f"ΔX_sat={100 * entry['dx_up']:+.3f}%pt  "
              f"×1/{args.factor:g}: {100 * entry['dx_down']:+.3f}%pt",
              flush=True)
    (run / "sensitivity.json").write_text(
        json.dumps(dict(slope=slope, g0=g0, factor=args.factor,
                        results=results), ensure_ascii=False),
        encoding="utf-8")


if __name__ == "__main__":
    main()
