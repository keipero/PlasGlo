r"""``nrp_pulse_train.py`` の記録（cycles.jsonl）から、飽和の兆候と反応別の収支を出す。

- サイクルごとの分解率 X と、1サイクルあたりの増分 dX（飽和すると 0 に近づく）。
- 選んだサイクルの、CO2 の正味の収支（パルス中＋パルス間）の反応別の内訳。
  記録は、各位相の上位6反応だけなので、それ以外は「その他」にまとめる。

使い方::

    python docs/co2_draft/nrp_report.py outputs/co2_nrp/loki_1e20 \
        --cycles 1 10 20 40
"""

import argparse
import json
from pathlib import Path


def load(path):
    r"""cycles.jsonl を読み込む。

    Parameters
    ----------
    path : pathlib.Path
        実行結果のフォルダ。

    Returns
    -------
    list of dict
        サイクルごとの記録。
    """
    with (path / "cycles.jsonl").open(encoding="utf-8-sig") as fp:
        return [json.loads(line) for line in fp if line.strip()]


def net_by_reaction(record):
    r"""1サイクルの CO2 の正味の収支を、反応ごとに足し合わせる。

    Parameters
    ----------
    record : dict
        1サイクルの記録。

    Returns
    -------
    tuple
        反応IDから (反応式, 正味の寄与 [m^-3]) への辞書と、記録に含まれない
        反応（上位6以外）の正味の寄与の合計。
    """
    net, shown = {}, 0.0
    other = 0.0
    for phase in ("pulse", "gap"):
        p = record[phase]
        listed = 0.0
        for item in p["production"]:
            eq, value = net.get(item["id"], (item["eq"], 0.0))
            net[item["id"]] = (eq, value + item["m3"])
            listed += item["m3"]
        for item in p["loss"]:
            eq, value = net.get(item["id"], (item["eq"], 0.0))
            net[item["id"]] = (eq, value - item["m3"])
            listed -= item["m3"]
        other += (p["total_production"] - p["total_loss"]) - listed
    return net, other


def main():
    r"""飽和の兆候の表と、選んだサイクルの反応別の収支を表示する。"""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", help="nrp_pulse_train.py の出力フォルダ")
    ap.add_argument("--cycles", type=int, nargs="*", default=[])
    ap.add_argument(
        "--every",
        type=int,
        default=1,
        help="表に出す行の間隔（直接計算の記録だけ数える）",
    )
    args = ap.parse_args()
    records = load(Path(args.run))
    direct = [r for r in records if "pulse" in r]
    first = direct[0]
    # 初期の CO2 は、サイクル1の終了時の値に、そのサイクルの純損失を足す。
    n_start = first["state_end_gap"]["CO2"] - (
        first["pulse"]["change"] + first["gap"]["change"]
    )
    print(f"{'cycle':>6s} {'X[%]':>9s} {'dX[%]':>9s} {'O':>9s} {'CO':>9s} "
          f"{'O2':>9s} {'O3':>9s} {'ne_end':>9s}")
    previous, count = 0.0, 0
    for r in records:
        s = r["state_end_gap"]
        x = r.get("x_end", 1 - s["CO2"] / n_start)
        dx = r.get("dx", x - previous)
        previous = x
        if "jump" in r:
            print(f"{r['cycle']:6d} {100 * x:9.4f} {100 * dx:9.5f}   "
                  f"<- ジャンプ（{r['jump']} サイクル分）")
            continue
        count += 1
        if count % args.every:
            continue
        print(f"{r['cycle']:6d} {100 * x:9.4f} {100 * dx:9.5f} "
              f"{s['O']:9.2e} {s['CO']:9.2e} {s['O2']:9.2e} "
              f"{s.get('O3', 0.0):9.2e} {r['state_end_pulse']['E']:9.2e}")
    for k in args.cycles:
        # 指定したサイクル以下で、最も近い直接計算の記録を使う。
        record = [r for r in direct if r["cycle"] <= k][-1]
        net, other = net_by_reaction(record)
        change = record["pulse"]["change"] + record["gap"]["change"]
        print(f"\n== サイクル {record['cycle']}：CO2 の正味の変化 "
              f"{change:.3e} m^-3（{100 * change / n_start:+.5f}%）")
        for rid, (eq, value) in sorted(net.items(), key=lambda kv: kv[1][1]):
            print(f"  {rid} {eq:34s} {value:+.3e}  "
                  f"({100 * value / n_start:+.5f}%)")
        print(f"  その他（上位6以外）{'':17s} {other:+.3e}")


if __name__ == "__main__":
    main()
