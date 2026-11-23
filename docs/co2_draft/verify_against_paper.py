r"""LoKI-GM の反応機構を、原典の論文の表と数値で照合する。

照合する論文：

- T. C. Dias, C. Fromentin, L. L. Alves, "A reaction mechanism for
  oxygen plasmas", Plasma Sources Sci. Technol. 32, 084003 (2023)
  （arXiv:2305.03674）の Table 1（電子の再結合）と
  Table 2（重い粒子の反応 R23〜R67）。``oxygen_novib.chem`` と比べる。
- P. Koelman et al., "A Comprehensive Chemical Model for the Splitting
  of CO2 in Non-Equilibrium Plasmas", Plasma Process. Polym. 14,
  1600155 (2017) の Table A3・A4。``CO2_novib.chem`` の「%Koelman」と
  注記された反応と比べる。
- B. Loenders, R. Michiels, A. Bogaerts の Supporting Information
  （プラズマ触媒による CO2 分解の微視的反応速度の研究）の Table S1。
  CO(a3P)・O(1D) の反応と比べる（``LOENDERS`` のコメント参照）。

``PAPER``・``KOELMAN`` の式は、論文の表を目で読んで書き写したもの
（Tg は K、速度係数は m^3/s と m^6/s）。LoKI 側は、該当行を
``loki_to_plasglo`` と同じ変換（``rate_expression``）で評価する。変換で
取り込まなかった反応（励起 O2 が反応物など）も、式そのものを比べる。

使い方::

    python docs/co2_draft/verify_against_paper.py \
        --loki-dir <.chem のフォルダ>
"""

import argparse
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import loki_to_plasglo as conv  # noqa: E402
from plasglo.chemistry.sample_mechanism import Expression  # noqa: E402

exp = math.exp
A, B, HZ = "O2(a1Dg)", "O2(b1Sg+)", "O2(A3Su+_C3Du_c1Su-)"
X, O3X, O3E = "O2(X)", "O3(X)", "O3(exc)"
P, D = "O(3P)", "O(1D)"
ION, ANION = "O(+,gnd)", "O(-,gnd)"
O2P = "O2(+,X)"

# 論文の番号 -> (LoKI の反応式, 論文の式 f(Tg))
PAPER = {
    "R23": (f"{A} + {P} -> {X} + {P}", lambda t: 7e-22),
    "R24": (
        f"{A} + {P} + {X} -> {X} + {P} + {X}",
        lambda t: 3e-44,
    ),
    "R25": (
        f"2{A} -> {B} + {X}",
        lambda t: 1.81e-24 * (t / 300) ** 3.8 * exp(700 / t),
    ),
    "R26": (f"{A} + {X} -> {X} + {X}", lambda t: 2.2e-24 * (t / 300) ** 0.8),
    "R27": (f"{B} + {P} -> {A} + {P}", lambda t: 4e-20),
    "R28": (
        f"{B} + {P} -> {X} + {P}",
        lambda t: 4e-20 + 1e-16 * exp(-3700 / t),
    ),
    "R29": (
        f"{B} + {P} -> {X} + {D}",
        lambda t: 3.39e-17 * (300 / t) ** 0.1 * exp(-4201 / t),
    ),
    "R30": (f"{HZ} + {P} -> {P} + {X}", lambda t: 4.95e-18),
    "R31": (f"{HZ} + {P} -> {D} + {A}", lambda t: 2.7e-18),
    "R32": (f"{HZ} + {P} -> {D} + {B}", lambda t: 1.35e-18),
    "R33": (f"{HZ} + {X} -> {X} + {X}", lambda t: 0.23e-19),
    "R34": (f"{HZ} + {X} -> {A} + {X}", lambda t: 1.86e-19),
    "R35": (f"{HZ} + {X} -> {B} + {X}", lambda t: 0.21e-19),
    "R36": (f"{D} + {P} -> {P} + {P}", lambda t: 8e-18),
    "R37": (
        f"2{P} + {X} -> {X} + {X}",
        lambda t: 0.5 * 3.81e-42 / t * exp(-170 / t),
    ),
    "R38": (
        f"2{P} + {X} -> {A} + {X}",
        lambda t: 0.33 * 3.81e-42 / t * exp(-170 / t),
    ),
    "R39": (
        f"2{P} + {X} -> {B} + {X}",
        lambda t: 0.17 * 3.81e-42 / t * exp(-170 / t),
    ),
    "R40": (f"3{P} -> {X} + {P}", lambda t: 2.5e-43 * t**-0.63),
    "R41": (f"{D} + {X} -> {P} + {X}", lambda t: 7e-18 * exp(67 / t)),
    "R42": (f"{D} + {X} -> {P} + {A}", lambda t: 8e-19 * exp(67 / t)),
    "R43": (f"{D} + {X} -> {P} + {B}", lambda t: 2.56e-17 * exp(67 / t)),
    "R44": (f"{B} + {O3X} -> 2{X} + {P}", lambda t: 1.5e-17),
    "R45": (f"{D} + {O3X} -> 2{X}", lambda t: 1.2e-16),
    "R46": (f"{D} + {O3X} -> {X} + 2{P}", lambda t: 1.2e-16),
    "R47": (
        f"{A} + {O3X} -> 2{X} + {P}",
        lambda t: 5.2e-17 * exp(-2840 / t),
    ),
    "R48": (
        f"{P} + {O3X} -> 2{X}",
        lambda t: 0.5 * 1.8e-17 * exp(-2300 / t),
    ),
    "R49": (
        f"{P} + {O3X} -> {A} + {X}",
        lambda t: 0.33 * 1.8e-17 * exp(-2300 / t),
    ),
    "R50": (
        f"{P} + {O3X} -> {B} + {X}",
        lambda t: 0.17 * 1.8e-17 * exp(-2300 / t),
    ),
    "R51": (f"{O3E} + {P} -> {O3X} + {P}", lambda t: 2e-19),
    "R52": (f"{O3E} + {X} -> {O3X} + {X}", lambda t: 3e-21),
    "R53": (
        f"{A} + {O3E} -> 2{X} + {P}",
        lambda t: 5 * 5.2e-17 * exp(-(2840 - 1553) / t),
    ),
    "R54": (
        f"{P} + {O3E} -> 2{X}",
        lambda t: 8e-18 * exp(-(2060 - 1553) / t),
    ),
    "R55": (
        f"2{P} + {X} -> {O3X} + {P}",
        lambda t: 2.1e-46 * exp(345 / t),
    ),
    "R56": (
        f"{P} + {X} + {X} -> {O3X} + {X}",
        lambda t: 0.33 * 6.4e-47 * exp(663 / t),
    ),
    "R57": (
        f"{P} + {X} + {X} -> {O3E} + {X}",
        lambda t: 0.67 * 6.4e-47 * exp(663 / t),
    ),
    "R58": (
        f"{P} + {X} + {O3X} -> 2{O3X}",
        lambda t: 1.66e-46 * exp(t / 300),
    ),
    "R59": (f"{ANION} + {X} -> {O3X} + e", lambda t: 1e-18),
    "R60": (f"{ANION} + {A} -> {O3X} + e", lambda t: 0.75 * 1.9e-16),
    "R61": (f"{ANION} + {P} -> {X} + e", lambda t: 1.3e-15),
    "R62": (f"{ANION} + {B} -> {P} + {X} + e", lambda t: 6.9e-16),
    "R63": (f"{ION} + {ANION} -> 2{P}", lambda t: 2.8e-13),
    "R64": (f"{ION} + {O3X} -> {O2P} + {X}", lambda t: 1e-16),
    "R65": (
        f"{ION} + {X} -> {O2P} + {P}",
        lambda t: 2e-17 * (300 / t) ** 0.5,
    ),
    "R66": (
        f"{ION} + {A} -> {O2P} + {P}",
        lambda t: 2e-17 * (300 / t) ** 0.5,
    ),
    "R67": (
        f"{O2P} + {ANION} -> {X} + {P}",
        lambda t: 9.6e-14 * (300 / t) ** 0.5,
    ),
}
# 電子温度に依存する再結合（Table 1）。引数は Te [eV]。
PAPER_TE = {
    "R21": (f"e + {O2P} -> 2{P}", lambda te: 0.36 * 2.2e-14 * te**-0.5),
    "R22": (f"e + {O2P} -> {P} + {D}", lambda te: 0.64 * 2.2e-14 * te**-0.5),
}


# Koelman et al., Plasma Process. Polym. 14, 1600155 (2017) の Table A3
# （中性）と Table A4（イオン）。ページを画像にして、目で読み取った値。
# LoKI の CO2_novib.chem で「%Koelman」と注記された反応に対応する。
CO2X, COX, CX, O2X = "CO2(X)", "CO(X)", "C(X)", "O2(X)"
CO2P, COP = "CO2(+,X)", "CO(+,X)"
KOELMAN = {
    "N4 (M=CO2)": (
        f"{P} + {COX} + {CO2X} -> {CO2X} + {CO2X}",
        lambda t: 2 * 8.2e-46 * exp(-1510 / t),
    ),
    "N4 (M=CO)": (
        f"{P} + {COX} + {COX} -> {CO2X} + {COX}",
        lambda t: 8.2e-46 * exp(-1510 / t),
    ),
    "N4 (M=O2)": (
        f"{P} + {COX} + {O2X} -> {CO2X} + {O2X}",
        lambda t: 8.2e-46 * exp(-1510 / t),
    ),
    "N8": (f"{O2X} + {CX} -> {COX} + {P}", lambda t: 3.0e-17),
    "I2": (f"{ION} + {CO2X} -> {O2P} + {COX}", lambda t: 8.1e-16),
    "I5": (f"{COP} + {CO2X} -> {CO2P} + {COX}", lambda t: 1.0e-15),
    "I12": (f"{ANION} + {COX} -> {CO2X} + e", lambda t: 5.5e-16),
    "I23": (f"{CO2P} + {P} -> {O2P} + {COX}", lambda t: 1.64e-16),
    "I24": (f"{CO2P} + {P} -> {ION} + {CO2X}", lambda t: 9.62e-17),
    "I25": (f"{CO2P} + {X} -> {O2P} + {CO2X}", lambda t: 5.3e-17),
    "I30": (f"{COP} + {X} -> {O2P} + {COX}", lambda t: 1.2e-16),
}


# B. Loenders, R. Michiels, A. Bogaerts, "How plasma reaction conditions
# affect the optimal catalyst: A microkinetic study of plasma-catalytic
# CO2 splitting" の Supporting Information, Table S1（単位 cm^3/s。
# ここでは m^3/s に換算）。CO(a3P) の反応は Cenian et al. 1994
# （Contrib. Plasma Phys. 34, 25）、O(1D) + CO2 は Dunlea &
# Ravishankara 2004 の値。LoKI は CO(a3P) の反応に、別の出典
# （Wysong 2000、Skrzypkowski 1999）を使っている。
CM3 = 1e-6
LOENDERS = {
    "S1 O(1D)+CO2": (
        f"{D} + {CO2X} -> {P} + {CO2X}",
        lambda t: 0.74e-10 * CM3 * exp(133 / t),
    ),
    "S1 CO+ + CO2": (
        f"{COP} + {CO2X} -> {CO2P} + {COX}",
        lambda t: 1e-9 * CM3,
    ),
    "S1 CO(a)+O2 -> CO+O2": (
        f"CO(a3P) + {O2X} -> {COX} + {O2X}",
        lambda t: 5.0e-11 * CM3,
    ),
    "S1 CO(a)+O2 -> CO+2O": (
        f"CO(a3P) + {O2X} -> {COX} + 2{P}",
        lambda t: 5.0e-11 * CM3,
    ),
    "S1 CO(a)+O2 -> CO2+O": (
        f"CO(a3P) + {O2X} -> {CO2X} + {P}",
        lambda t: 3.0e-11 * CM3,
    ),
    "S1 CO(a)+CO -> 2CO": (
        f"CO(a3P) + {COX} -> 2{COX}",
        lambda t: 1.4e-10 * CM3,
    ),
    "S1 CO(a)+CO -> CO2+C": (
        f"CO(a3P) + {COX} -> {CO2X} + {CX}",
        lambda t: 1.4e-12 * CM3,
    ),
    "S1 CO(a)+CO2 -> CO+CO2": (
        f"CO(a3P) + {CO2X} -> {COX} + {CO2X}",
        lambda t: 1.5e-11 * CM3,
    ),
    "S1 CO(a)+CO2 -> 2CO+O": (
        f"CO(a3P) + {CO2X} -> 2{COX} + {P}",
        lambda t: 1.5e-11 * CM3,
    ),
}


def loki_value(item, tg, te_ev):
    r"""LoKI の1行の速度係数を、変換と同じ式で評価する。

    Parameters
    ----------
    item : dict
        ``loki_to_plasglo.parse_chem`` の1行。
    tg : float
        気体温度 [K]。
    te_ev : float
        電子温度 [eV]。

    Returns
    -------
    float
        速度係数 [m^3/s または m^6/s]。
    """
    expression = Expression(conv.rate_expression(item["func"], item["params"]))
    return expression(
        {"tgas": tg, "tg": tg, "te": te_ev / conv.KB_EV}
    )


def compare(loki_file, paper, paper_te=None):
    r"""論文の各反応を、LoKI の行と照合する。

    Parameters
    ----------
    loki_file : pathlib.Path
        LoKI の ``.chem`` ファイル。
    paper : dict
        論文の番号から、(LoKI の反応式, 式 f(Tg)) への辞書。
    paper_te : dict or None, optional
        電子温度（eV）に依存する反応の辞書。

    Returns
    -------
    list of dict
        反応ごとの照合結果（見つからない場合は ``found`` が False）。
    """
    index = {}
    for item in conv.parse_chem(loki_file):
        left, right = conv.split_reaction(item["reaction"])
        index.setdefault(conv.canonical(left, right), []).append(item)
    out = []
    table = [(k, v, False) for k, v in paper.items()]
    table += [(k, v, True) for k, v in (paper_te or {}).items()]
    for number, (text, formula), electron in table:
        left, right = conv.split_reaction(text)
        hits = index.get(conv.canonical(left, right), [])
        row = dict(number=number, text=text, found=len(hits) == 1)
        if len(hits) == 1:
            worst = 0.0
            points = (
                [(300.0, te) for te in (0.5, 1.0, 3.0)]
                if electron
                else [(tg, 1.0) for tg in (300.0, 600.0, 1500.0)]
            )
            for tg, te in points:
                want = formula(te if electron else tg)
                got = loki_value(hits[0], tg, te)
                worst = max(worst, abs(got - want) / abs(want))
            row["worst"] = worst
            row["k300"] = (formula(1.0 if electron else 300.0),
                           loki_value(hits[0], 300.0, 1.0))
            row["line"] = hits[0]["line"]
        else:
            row["hits"] = len(hits)
        out.append(row)
    return out


def main():
    r"""照合を実行し、結果の表と、一致しない反応を表示する。"""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--loki-dir", required=True)
    args = ap.parse_args()
    loki = Path(args.loki_dir)
    papers = [
        ("Dias 2023（酸素。Table 1・2）", "oxygen_novib.chem", PAPER, PAPER_TE),
        ("Koelman 2017（Table A3・A4）", "CO2_novib.chem", KOELMAN, None),
        (
            "Loenders らの SI（Table S1。CO(a)・O(1D) の反応）",
            "CO2_novib.chem",
            LOENDERS,
            None,
        ),
    ]
    for title, name, paper, paper_te in papers:
        rows = compare(loki / name, paper, paper_te)
        print(f"\n## {title} と {name}")
        print(f"{'番号':11s} {'k(300K) 論文':>13s} {'k(300K) LoKI':>13s} "
              f"{'最大相対差':>10s}  LoKI行  反応（LoKI の表記）")
        for r in rows:
            if r["found"]:
                print(f"{r['number']:11s} {r['k300'][0]:13.3e} "
                      f"{r['k300'][1]:13.3e} {r['worst']:10.1e}  "
                      f"{r['line']:5d}   {r['text']}")
            else:
                print(f"{r['number']:11s} {'LoKI に該当行なし':>26s}  "
                      f"({r.get('hits', 0)}件)  {r['text']}")
        found = [r for r in rows if r["found"]]
        bad = [r for r in found if r["worst"] > 1e-6]
        print(f"照合：{len(found)}/{len(rows)} 反応を LoKI の行に対応づけ、"
              f"最大相対差 {max(r['worst'] for r in found):.1e}。"
              f"1e-6 超：{len(bad)} 反応")
        for r in bad:
            print(f"  差あり：{r['number']} {r['text']}  "
                  f"300 K で 論文 {r['k300'][0]:.3e}、LoKI {r['k300'][1]:.3e}"
                  f"（{(r['k300'][1] / r['k300'][0] - 1) * 100:+.1f}%）")


if __name__ == "__main__":
    main()
