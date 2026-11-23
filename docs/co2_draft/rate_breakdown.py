r"""Biagi の純CO2の電子衝突係数をカテゴリ別に出し、解離の内訳を見積もる。

Biagi は、CO2の全過程（弾性を除く91本）をチャンネルとして宣言し、分類ごとの反応係数を
E/Nごとに表示する。E/Nごとに ``ElectronSwarm`` を作り直す（格子の履歴依存を
避けるため。README 参照）。
"""

import argparse
import collections
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_co2_channel_draft as b  # noqa: E402
from plasglo.electrons.mixture import ElectronSwarm  # noqa: E402

ORDER = [
    "vibration",
    "diss_1Delta",
    "diss_1Pi",
    "diss_dipole",
    "diss_excitation",
    "excitation_other",
    "attachment",
    "ion_CO2+",
    "ion_CO2+_A",
    "ion_CO2+_B",
    "ion_O+",
    "ion_CO+",
    "ion_other",
]


def main():
    r"""カテゴリ別の係数表と、解離の内訳・感度の目安を表示する。"""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--source-dir",
        help="断面積ファイルのフォルダ（既定：環境変数 CO2_XS_DIR）",
    )
    b.CO2_SET = "biagi"
    b.set_source_dir(ap.parse_args().source_dir)
    path = b.SOURCE_DIR / b.CO2_FILES["biagi"]
    with path.open(encoding="utf-8") as fp:
        collisions = b.parser.parse(fp)
    channels, category = [], {}
    for k, c in enumerate(collisions):
        if c["kind"] == "ELASTIC":
            continue
        cid = f"X{k:03d}"
        channels.append(
            dict(
                id=cid,
                target="CO2",
                kind=c["kind"],
                product=c.get("product"),
                threshold_eV=float(c.get("threshold", 0.0)),
            )
        )
        category[cid] = b.classify_co2(c)[0]
    groups = collections.defaultdict(list)
    for cid, name in category.items():
        groups[name].append(cid)
    n_gas = 101325 / (1.380649e-23 * 300)
    print("純CO2、300 K、101325 Pa。反応係数 [m^3/s]")
    names = " ".join(f"{o[:9]:>9s}" for o in ORDER)
    print(f"{'E/N':>5s} {'<e>':>5s} " + names)
    for td in (50.0, 100.0, 200.0, 300.0, 500.0):
        swarm = ElectronSwarm([path], {"CO2": "CO2"}, channels, cells=400)
        state = swarm.solve(td, {"CO2": n_gas}, 300.0)
        sums = {
            g: sum(state["rates"][i] for i in ids) for g, ids in groups.items()
        }
        print(
            f"{td:5.0f} {state['energy_eV']:5.2f} "
            + " ".join(f"{sums.get(o, 0.0):9.2e}" for o in ORDER)
        )
        diss = sum(sums[g] for g in ORDER if g.startswith("diss"))
        parts = {
            g: 100 * sums[g] / diss for g in ORDER if g.startswith("diss")
        }
        print(
            f"      解離の合計 {diss:.2e}  内訳[%] "
            + "  ".join(f"{g[5:]}={v:4.1f}" for g, v in parts.items())
            + f"  | 他の電子励起/解離 = {sums['excitation_other'] / diss:5.1f}"
        )


if __name__ == "__main__":
    main()
