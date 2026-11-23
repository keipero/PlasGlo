r"""``nrp_sensitivity.py`` の結果と、飽和まで走らせた検証から、感度の図を描く。

トルネード図：反応ごとに、速度を 2 倍・1/2 倍にしたときの、飽和時の分解率の
変化の目安（1次近似）。検証の計算（``verify_<ID>_x2`` フォルダ）があれば、
その結果（飽和まで走らせた分解率 - 基準の飽和値）を点で重ねる。

使い方::

    python docs/co2_draft/nrp_sens_plot.py outputs/co2_nrp
"""

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = ["Yu Gothic", "Meiryo", "sans-serif"]


def final_x(folder):
    r"""出力フォルダの最後の記録の分解率 X を返す。

    Parameters
    ----------
    folder : pathlib.Path
        ``nrp_pulse_train.py`` の出力フォルダ。

    Returns
    -------
    tuple
        (X、サイクル番号、最後の dX)。記録がなければ None。
    """
    path = folder / "cycles.jsonl"
    if not path.exists():
        return None
    with path.open(encoding="utf-8-sig") as fp:
        records = [json.loads(line) for line in fp if line.strip()]
    records = [r for r in records if "x_end" in r]
    if not records:
        return None
    last = records[-1]
    return last["x_end"], last["cycle"], last["dx"]


def main():
    r"""トルネード図を ``sensitivity_tornado.png`` に書く。"""
    base = Path(sys.argv[1])
    run = base / "loki_1e20"
    doc = json.loads((run / "sensitivity.json").read_text("utf-8"))
    rows = sorted(doc["results"], reverse=True,
                  key=lambda r: max(abs(r["dx_up"]), abs(r["dx_down"])))[:14]
    rows = rows[::-1]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(10, 7))
    up = [100 * r["dx_up"] for r in rows]
    down = [100 * r["dx_down"] for r in rows]
    ax.barh(y + 0.2, up, 0.38, color="#c0392b", label="速度 ×2（1次近似）")
    ax.barh(y - 0.2, down, 0.38, color="#2471a3", label="速度 ×1/2（1次近似）")
    base_x = final_x(run)
    plotted = False
    for i, r in enumerate(rows):
        res = final_x(base / f"verify_{r['id']}_x2")
        if res and base_x:
            dx = 100 * (res[0] - base_x[0])
            ax.plot(dx, i + 0.2, "k*", ms=12,
                    label=None if plotted else "検証（×2、実際に走らせた）")
            ax.annotate(f"{dx:+.1f}", (dx, i + 0.2), xytext=(6, -3),
                        textcoords="offset points", fontsize=8)
            plotted = True
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['id']}  {r['equation']}" for r in rows],
                       fontsize=9)
    ax.axvline(0, color="k", lw=0.8)
    ax.set(xlabel="飽和時の分解率の変化 ΔX_sat [%pt]",
           title="飽和した分解率の、反応速度への感度（E/N=100 Td、10 kHz）")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    out = base / "sensitivity_tornado.png"
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == "__main__":
    main()
