r"""CO2の解離の反応係数を、IST-Lisbon（更新版）と Biagi で比べる。

純CO2（300 K、101325 Pa）で、E/Nごとに ``ElectronSwarm`` を作り直して計算する
（格子の履歴依存を避けるため）。

- IST-Lisbon：解離2チャンネル（7 eV の15%、10.5 eV の全て）の合計。
- Biagi：ラベル付きの解離（¹Δ、¹Π、中性双極子、解離励起）の合計 ``k_lab`` と、
  その他の電子励起の合計 ``k_other``。IST-Lisbon に一致させるのに要る、
  ``k_other`` の解離分岐比 ``f = (k_IST - k_lab) / k_other`` も出す。

測定値は、Liu et al. 2025（Plasma Sources Sci. Technol. 34, 035003）の図2(b)
（Morillo-Candas ら 2020）を目で読んだ値で、±0.3×10^-16 m^3/s 程度の誤差を持つ。
"""

import argparse
import logging
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_co2_channel_draft as b  # noqa: E402
from plasglo.electrons.mixture import ElectronSwarm  # noqa: E402

UNIT = 1e-16
MEASURED = {60: 0.75, 70: 0.70, 80: 1.5, 90: 2.0, 100: 3.5, 105: 4.7}
N_GAS = 101325 / (1.380649e-23 * 300)


def channels_of(path, classify):
    r"""CO2の全過程（弾性を除く）をチャンネルとして宣言し、分類も返す。

    Parameters
    ----------
    path : pathlib.Path
        断面積ファイル。
    classify : callable
        1過程から分類名を返す関数（``b.classify_co2`` など）。

    Returns
    -------
    tuple
        チャンネルのリストと、チャンネルIDから分類名への辞書。
    """
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
        category[cid] = classify(c)[0]
    return channels, category


def sums(path, channels, category, td):
    r"""E/Nを指定して、分類ごとの反応係数の合計を求める。

    Parameters
    ----------
    path : pathlib.Path
        断面積ファイル。
    channels : list of dict
        チャンネルの宣言。
    category : dict
        チャンネルIDから分類名への辞書。
    td : float
        換算電界 [Td]。

    Returns
    -------
    tuple
        平均エネルギー [eV] と、分類ごとの反応係数の辞書 [m^3/s]。
    """
    swarm = ElectronSwarm([path], {"CO2": "CO2"}, channels, cells=400)
    state = swarm.solve(float(td), {"CO2": N_GAS}, 300.0)
    out = {}
    for cid, name in category.items():
        out[name] = out.get(name, 0.0) + state["rates"][cid]
    return state["energy_eV"], out


def main():
    r"""IST-Lisbon と Biagi の解離の反応係数を、E/Nごとに並べて表示する。"""
    logging.disable(logging.CRITICAL)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--source-dir",
        help="断面積ファイルのフォルダ（既定：環境変数 CO2_XS_DIR）",
    )
    folder = b.set_source_dir(ap.parse_args().source_dir)
    ist = folder / b.CO2_FILES["ist"]
    biagi = folder / b.CO2_FILES["biagi"]
    if not biagi.is_file():
        raise SystemExit(f"{biagi} がありません。")
    ist_ch, ist_cat = channels_of(ist, b.classify_co2_ist)
    bi_ch, bi_cat = channels_of(biagi, b.classify_co2)
    print("純CO2、300 K、101325 Pa。単位 1e-16 m^3/s")
    print(
        f"{'E/N':>4s} {'測定':>5s} | {'IST <e>':>7s} {'k_diss':>8s} | "
        f"{'Bi <e>':>6s} {'k_lab':>7s} {'k_other':>8s} {'f(等価)':>8s}"
    )
    for td in (50, 60, 70, 80, 90, 100, 105, 150, 200, 300, 500):
        e1, s1 = sums(ist, ist_ch, ist_cat, td)
        e2, s2 = sums(biagi, bi_ch, bi_cat, td)
        k_ist = (s1["diss_7eV"] + s1["diss_10p5eV"]) / UNIT
        labeled = ("diss_1Delta", "diss_1Pi", "diss_dipole", "diss_excitation")
        lab = sum(s2[n] for n in labeled) / UNIT
        other = s2["excitation_other"] / UNIT
        meas = MEASURED.get(td)
        print(
            f"{td:4d} {'' if meas is None else f'{meas:5.2f}':>5s} | "
            f"{e1:7.3f} {k_ist:8.3f} | {e2:6.3f} {lab:7.3f} {other:8.3f} "
            f"{(k_ist - lab) / other:8.3f}"
        )


if __name__ == "__main__":
    main()
