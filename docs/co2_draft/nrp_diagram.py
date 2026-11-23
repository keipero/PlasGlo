r"""``nrp_pathways.py`` の結果（pathways_saturated.json）から、反応経路図を描く。

飽和した1サイクルの、化学種間の流れ [m^-3 / サイクル] を、矢印の太さ
（流れの対数）で表す。赤は CO2 の分解、青は CO2 の再生成、灰色は CO2 を
直接は変えない往復や O 系の流れ。

使い方::

    python docs/co2_draft/nrp_diagram.py outputs/co2_nrp/loki_1e20
"""

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

plt.rcParams["font.family"] = ["Yu Gothic", "Meiryo", "sans-serif"]

NODES = {
    "CO2": (1.0, 4.6), "CO2(e1)": (1.0, 8.0),
    "CO(a3Pi)": (5.5, 8.0), "CO": (5.5, 4.6),
    "O(1D)": (9.8, 8.0), "O": (9.8, 5.0),
    "O2": (9.8, 1.2), "O3": (5.5, 1.2),
}
LABEL = {"CO(a3Pi)": "CO(a³Π)", "O(1D)": "O(¹D)", "CO2(e1)": "CO₂(e1)",
         "CO2": "CO₂", "O2": "O₂", "O3": "O₃"}
RED, BLUE, GRAY = "#c0392b", "#2471a3", "#7f8c8d"


def flux(data, species, ids):
    r"""指定した反応の、ある化学種への積算寄与の大きさ [m^-3] を足し合わせる。

    Parameters
    ----------
    data : dict
        ``pathways_saturated.json`` の ``species``。
    species : str
        化学種。
    ids : list of str
        反応ID。

    Returns
    -------
    float
        積算量の合計 [m^-3]。

    Examples
    --------
    >>> d = {"A": {"production": [{"id": "R1", "m3": 2.0}], "loss": []}}
    >>> flux(d, "A", ["R1", "R9"])
    2.0
    """
    s = data[species]
    return sum(i["m3"] for k in ("production", "loss") for i in s[k]
               if i["id"] in ids)


def edges(data):
    r"""化学種間の流れを、辞書のリストにして返す。

    各流れは、進行方向の左側に ``off`` だけずらして描く（逆向きの流れと
    重ならない）。``tp`` はラベルの位置（始点から終点への割合）。

    Parameters
    ----------
    data : dict
        ``pathways_saturated.json`` の ``species``。

    Returns
    -------
    list of dict
        流れ（始点 a、終点 b、量 v、ラベル t、色 c、ずらし off、位置 tp）。
    """
    f = flux
    return [
        dict(a="CO2", b="CO2(e1)", v=f(data, "CO2(e1)", ["R027"]),
             t="R027 振動励起", c=GRAY, off=0.2, tp=0.5),
        dict(a="CO2(e1)", b="CO2", v=f(data, "CO2(e1)", ["R028"]),
             t="R028 失活", c=GRAY, off=0.2, tp=0.5),
        dict(a="CO2", b="CO(a3Pi)", v=f(data, "CO2", ["R003"]),
             t="R003 解離\n→CO(a)+O", c=RED, off=0.2, tp=0.5),
        dict(a="CO(a3Pi)", b="CO2", v=f(data, "CO2", ["R008", "R009"]),
             t="R008,R009\n再生成", c=BLUE, off=0.2, tp=0.85),
        dict(a="CO2", b="CO", v=f(data, "CO2", ["R002"]),
             t="R002 解離\n→CO+O(¹D)", c=RED, off=0.2, tp=0.62),
        dict(a="CO", b="CO2",
             v=f(data, "CO2", ["R014", "R015", "R016", "R025"]),
             t="R014,R015,R016,R025\n（O と反応して再生成）", c=BLUE,
             off=0.2, tp=0.5),
        dict(a="CO", b="CO(a3Pi)", v=f(data, "CO(a3Pi)", ["R030"]),
             t="R030 電子励起", c=GRAY, off=0.25, tp=0.5),
        dict(a="CO(a3Pi)", b="CO",
             v=f(data, "CO(a3Pi)", ["R010", "R006", "R007", "R011", "R013"]),
             t="消光 R010,R006,\nR007,R011,R013", c=GRAY, off=0.25, tp=0.78),
        dict(a="CO(a3Pi)", b="CO", v=f(data, "CO2", ["R012"]),
             t="R012 +CO₂→2CO+O\n（CO₂の分解）", c=RED, off=1.6, tp=0.3),
        dict(a="O(1D)", b="O", v=f(data, "O(1D)", ["R018", "R046", "R044"]),
             t="R018 等 失活", c=GRAY, off=0.25, tp=0.5),
        dict(a="O2", b="O", v=f(data, "O", ["R033", "R034"]),
             t="R033,R034\n電子衝突解離", c=GRAY, off=0.25, tp=0.5),
        dict(a="O", b="O3", v=f(data, "O3", ["R053", "R054", "R055"]),
             t="R053-055\nO₃生成", c=GRAY, off=0.25, tp=0.5),
        dict(a="O3", b="O2", v=f(data, "O3", ["R049", "R050", "R051"]),
             t="R049-051\nO+O₃→2O₂", c=GRAY, off=-0.25, tp=0.5, side=-1),
    ]


def draw(ax, e):
    r"""1本の流れを、矢印とラベルで描く。

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        描画先。
    e : dict
        ``edges`` の要素。
    """
    (x1, y1), (x2, y2) = NODES[e["a"]], NODES[e["b"]]
    length = math.hypot(x2 - x1, y2 - y1)
    ux, uy = (x2 - x1) / length, (y2 - y1) / length
    nx, ny = -uy, ux  # 進行方向の左
    trim = 0.8
    p1 = (x1 + ux * trim + nx * e["off"], y1 + uy * trim + ny * e["off"])
    p2 = (x2 - ux * trim + nx * e["off"], y2 - uy * trim + ny * e["off"])
    width = 1.6 + 2.0 * math.log10(max(e["v"] / 1e20, 1.0))
    ax.add_patch(FancyArrowPatch(
        p1, p2, arrowstyle="-|>", mutation_scale=12 + width * 2,
        lw=width, color=e["c"], alpha=0.85, shrinkA=0, shrinkB=0,
    ))
    side = e.get("side", 1)
    px = p1[0] + (p2[0] - p1[0]) * e["tp"] + nx * 0.5 * side
    py = p1[1] + (p2[1] - p1[1]) * e["tp"] + ny * 0.5 * side
    ax.text(px, py, f"{e['t']}\n{e['v'] / 1e21:.2f}e21", fontsize=8,
            color=e["c"], ha="center", va="center",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))


def main():
    r"""pathways_saturated.json を読み、pathway_diagram.png を書く。"""
    run = Path(sys.argv[1])
    doc = json.loads((run / "pathways_saturated.json").read_text("utf-8"))
    data = doc["species"]
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(-0.8, 11.0)
    ax.set_ylim(-0.4, 9.0)
    ax.axis("off")
    for e in edges(data):
        draw(ax, e)
    for name, (x, y) in NODES.items():
        ax.add_patch(FancyBboxPatch(
            (x - 0.6, y - 0.32), 1.2, 0.64, boxstyle="round,pad=0.05",
            fc="#fdebd0" if name == "CO2" else "#eef3f7", ec="#34495e",
            lw=1.6, zorder=3))
        ax.text(x, y, LABEL.get(name, name), ha="center", va="center",
                fontsize=12, weight="bold", zorder=4)
    ax.set_title(
        f"飽和状態（サイクル {doc['cycle']}）の1サイクルの反応経路 "
        "[m⁻³/サイクル]\n赤：CO₂分解　青：CO₂再生成　灰：往復・O系　"
        "太さ：流れの対数", fontsize=12)
    out = run / "pathway_diagram.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
