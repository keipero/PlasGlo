r"""``nrp_vib_scan.py`` の結果（vib_scan.json）から、振動あり・なしの比較図を描く。

左：実効エネルギーコスト（投入エネルギー / 分解した CO2 の数）と、EEDF だけで
求めたコスト（``nrp_efficiency.py``、純 CO2、化学なし）。
右：振動あり / なし のコストの比 - 1 [%]。サイクル終了時に残る振動エネルギー
の、投入エネルギーに対する割合も重ねる。

使い方::

    python docs/co2_draft/nrp_vib_plot.py outputs/co2_nrp
"""

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = ["Yu Gothic", "Meiryo", "sans-serif"]


def last_cycle(results, model):
    r"""あるモデルの、最後のサイクルの値を E/N ごとに取り出す。

    Parameters
    ----------
    results : list of dict
        ``vib_scan.json`` の中身。
    model : str
        ``振動あり`` または ``振動なし``。

    Returns
    -------
    dict
        E/N から、最後のサイクルの記録への辞書。

    Examples
    --------
    >>> rows = [{"cycle": 1}, {"cycle": 2}]
    >>> r = [{"model": "a", "td": 1, "rows": rows}]
    >>> last_cycle(r, "a")[1]["cycle"]
    2
    """
    return {r["td"]: r["rows"][-1] for r in results if r["model"] == model}


def main():
    r"""比較図を ``vib_comparison.png`` に書く。"""
    folder = Path(sys.argv[1])
    results = json.loads((folder / "vib_scan.json").read_text("utf-8"))
    with_vib = last_cycle(results, "振動あり")
    without = last_cycle(results, "振動なし")
    tds = sorted(with_vib)
    cost_v = np.array([with_vib[t]["cost"] for t in tds], dtype=float)
    cost_n = np.array([without[t]["cost"] for t in tds], dtype=float)
    ok = np.isfinite(cost_v) & (cost_v < 1e4)
    td = np.array(tds, dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    ax.plot(td[ok], cost_n[ok], "o-", label="振動なし（化学入り）")
    ax.plot(td[ok], cost_v[ok], "s--", label="振動あり（化学入り）")
    eff_path = folder / "loki_1e20" / "efficiency_map.json"
    if eff_path.exists():
        doc = json.loads(eff_path.read_text("utf-8"))
        e_td = np.array(doc["td"])
        e_cost = [r["cost_eV"] for r in doc["table"][0]["rows"]]
        ax.plot(e_td, e_cost, ":", color="gray",
                label="EEDF のみ（純CO2、再生成なし）")
    ax.axhline(2.93, color="k", lw=0.6)
    ax.text(td.min(), 3.2, "解離のエンタルピー 2.93 eV", fontsize=8)
    ax.set(xscale="log", yscale="log", xlabel="E/N [Td]",
           ylabel="解離1回あたりのエネルギーコスト [eV/分子]",
           title="振動準位の有無によるコスト（純CO2から3サイクル）",
           xlim=(15, 200), ylim=(2, 2e3))
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    ax = axes[1]
    diff = 100 * (cost_v / cost_n - 1)
    ax.plot(td[ok], diff[ok], "o-", color="C3", label="コストの差 [%]")
    ax.axhline(0, color="k", lw=0.6)
    ax.set(xscale="log", xlabel="E/N [Td]", ylabel="コストの差 [%]\n"
           "（振動あり / なし − 1）",
           title="振動の効果と、残る振動エネルギー", xlim=(15, 200))
    twin = ax.twinx()
    stored = [with_vib[t]["vib"] / with_vib[t]["energy"] * 1e3 for t in tds]
    twin.plot(td, stored, "^:", color="C2",
              label="残る振動エネルギー / 投入 [‰]")
    twin.set_ylabel("残る振動エネルギー / 投入 [‰]", color="C2")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    twin.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    out = folder / "vib_comparison.png"
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == "__main__":
    main()
