r"""LoKI-GM の CO2・CO・O2 の反応機構を、PlasGlo の入力形式へ変換する。

LoKI-GM（https://github.com/LoKI-Suite/LoKI-GM、GPL-3.0）の次の3ファイルを読み、
PlasGlo の反応CSVとケースJSONを書き出す。

- ``Code/Input/CO2/CO2_novib.chem``
  （Silva et al., Plasma Sources Sci. Technol. 29, 125020 (2020) に基づく。
  O(1D)、CO(a3Pi)、CO2(e1) を含む）
- ``Code/Input/CO/CO_simple_X.chem``
- ``Code/Input/Oxygen/oxygen_novib.chem``
  （Dias et al., Plasma Sources Sci. Technol. 32, 084003 (2023) に基づく）

**LoKI-GM のデータ（GPL-3.0）も、変換後の反応表も、このリポジトリには含めない。**
このスクリプトは、ユーザーがダウンロードした元ファイルを読み、ユーザーが指定した
フォルダ（既定は、Git が追跡しない ``outputs/``）へ出力する。

変換の約束（LoKI-GM のソース ``Code/RateCoeffFunctions`` で確認）：

- 速度係数は SI 単位 [m^3/s]、三体反応は [m^6/s]（``length_unit`` は ``m``）。
- ``arrheniusGasTemp``：a*exp(b/Tg)
- ``powerGasTemp``：a*Tg^b（300 K で割らない）
- ``modifiedArrheniusGasTemp``：a*Tg^b*exp(c/Tg)
- ``powerGasElecTemp``：a*Tg^b*Te^c（Te は eV）
- ``powerElectronTemp``：a*Te^b（Te は eV）
- ``offsetArrheniusGasTemp``：a0 + a*exp(b/Tg)
- ``expGasTemp``：a*exp(Tg/b)
- ``eedf``：電子衝突。断面積のチャンネルに対応付ける。

取り込まない反応（理由は ``conversion_report.md`` に出す）：壁・輸送、追跡しない種を
反応物にもつもの、断面積データのない電子衝突、反応の前後で同じになるもの。
追跡しない励起種（O2(a1Dg)、O2(b1Sg+)、O2(A3Su+_C3Du_c1Su-)、O3(exc)）が
**生成物**のときは、基底状態（O2(X)、O3(X)）に集約する。これは仮定で、
励起種の反応性（例：O2(a1Dg) + O3）は失われる。

使い方::

    python docs/co2_draft/loki_to_plasglo.py --loki-dir <3ファイルのフォルダ> \
        --source-dir <断面積のフォルダ> --out outputs/co2_mechanism
"""

import argparse
import collections
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_co2_channel_draft as build  # noqa: E402

KB_EV = 8.617333262e-5  # Boltzmann定数 [eV/K]。te [K] から Te [eV] への換算
CHEM_FILES = ("CO2_novib.chem", "CO_simple_X.chem", "oxygen_novib.chem")

# LoKI の種名 -> PlasGlo の種名
RENAME = {
    "e": "E",
    "CO2(X)": "CO2",
    "CO(X)": "CO",
    "O2(X)": "O2",
    "O3(X)": "O3",
    "O(3P)": "O",
    "O(1D)": "O(1D)",
    "C(X)": "C",
    "CO(a3P)": "CO(a3Pi)",
    "CO2(e1)": "CO2(e1)",
    "CO2(+,X)": "CO2^+",
    "CO(+,X)": "CO^+",
    "O(+,gnd)": "O^+",
    "O2(+,X)": "O2^+",
    "O(-,gnd)": "O^-",
}
# 追跡しない励起種を、生成物のときだけ基底状態に集約する。
ALIAS_PRODUCTS = {
    "O2(a1Dg)": "O2(X)",
    "O2(b1Sg+)": "O2(X)",
    "O2(A3Su+_C3Du_c1Su-)": "O2(X)",
    "O3(exc)": "O3(X)",
}
# PlasGlo の種名 -> (元素組成, 電荷, 励起種か)
SPECIES = {
    "E": ({}, -1, False),
    "CO2": ({"C": 1, "O": 2}, 0, False),
    "CO": ({"C": 1, "O": 1}, 0, False),
    "O2": ({"O": 2}, 0, False),
    "O3": ({"O": 3}, 0, False),
    "O": ({"O": 1}, 0, False),
    "O(1D)": ({"O": 1}, 0, True),
    "C": ({"C": 1}, 0, False),
    "CO(a3Pi)": ({"C": 1, "O": 1}, 0, True),
    "CO2(e1)": ({"C": 1, "O": 2}, 0, True),
    "CO2^+": ({"C": 1, "O": 2}, 1, False),
    "CO^+": ({"C": 1, "O": 1}, 1, False),
    "O^+": ({"O": 1}, 1, False),
    "O2^+": ({"O": 2}, 1, False),
    "O^-": ({"O": 1}, -1, False),
}
# 'gas'（全気体の密度）に使う、第三体の重み
THIRD_BODY = ("CO2", "CO", "O2", "O3", "O", "C")

# 電子衝突：(標的, 種類, 生成物に含まれる文字列 または None, しきい値 または None)
# 反応は、反応物・生成物の集合で照合する（``canonical``）。
ELECTRON_SELECTORS = {
    "e + CO2(X) -> e + e + CO2(+,X)": ("CO2", "IONIZATION", None, 13.8),
    "e + CO2(X) -> e + CO(X) + O(1D)": ("CO2", "EXCITATION", "O(1D)", 7.0),
    "e + CO2(X) -> e + CO(a3P) + O(3P)": (
        "CO2",
        "EXCITATION",
        "CO(a3P)",
        10.5,
    ),
    "e + CO2(X) -> CO(X) + O(-,gnd)": ("CO2", "ATTACHMENT", None, None),
    "e + CO2(X) -> e + CO2(e1)": ("CO2", "EXCITATION", "CO2(e1)", 7.0),
    "e + CO(X) -> C(X) + O(-,gnd)": ("CO", "ATTACHMENT", None, None),
    "e + CO(X) -> e + CO(a3P)": ("CO", "EXCITATION", "CO(a3P)", 6.006),
    "e + CO(X) -> e + C(X) + O(3P)": ("CO", "EXCITATION", "C(X)", 11.1),
    "e + CO(X) -> e + e + CO(+,X)": ("CO", "IONIZATION", None, 14.01),
    "e + O2(X) -> e + 2O(3P)": ("O2", "EXCITATION", "6.0eV", 6.0),
    "e + O2(X) -> e + O(3P) + O(1D)": ("O2", "EXCITATION", "8.4eV", 8.4),
    "e + O2(X) -> O(-,gnd) + O(3P)": ("O2", "ATTACHMENT", "O^-", None),
    "e + O2(X) -> e + e + O2(+,X)": ("O2", "IONIZATION", None, 12.06),
}
# 運動量移行・三体付着など、化学反応には使わないが、EEDFには必要なチャンネル。
O2_THREE_BODY = dict(
    species="O2",
    power=1,
    density_unit="cm^-3",
    coefficient_scale=1e-06,
    virtual_target="__O2_three_body__",
)


def parse_chem(path):
    r"""LoKI の ``.chem`` ファイルを、反応の辞書のリストに読み込む。

    Parameters
    ----------
    path : pathlib.Path
        ``.chem`` ファイル。

    Returns
    -------
    list of dict
        ``line``、``file``、``reaction``、``func``、``params``、``comment``
        を持つ辞書。コメント行と空行は除く。

    Examples
    --------
    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     p = Path(tmp) / "x.chem"
    ...     line = "A + B -> C | constantRateCoeff | 1e-16 | % ref"
    ...     _ = p.write_text(line)
    ...     r = parse_chem(p)[0]
    >>> r["reaction"], r["func"], r["params"], r["comment"]
    ('A + B -> C', 'constantRateCoeff', ['1e-16'], 'ref')
    """
    out = []
    text = path.read_text(encoding="utf-8")
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("%"):
            continue
        body, _, comment = line.partition("%")
        fields = [f.strip() for f in body.split("|")]
        if len(fields) < 2 or not fields[0]:
            continue
        params = (
            [p.strip() for p in fields[2].split(",")]
            if len(fields) > 2 and fields[2]
            else []
        )
        out.append(
            dict(
                line=number,
                file=path.name,
                reaction=fields[0],
                func=fields[1],
                params=params,
                comment=comment.strip(),
            )
        )
    return out


def parse_side(text):
    r"""反応式の片側を、(係数, 種名) のリストにする。

    Parameters
    ----------
    text : str
        ``A + 2B`` のような、反応式の片側。

    Returns
    -------
    list of tuple
        (係数 float, LoKI の種名) のリスト。

    Examples
    --------
    >>> parse_side("2O2(a1Dg) + e")
    [(2.0, 'O2(a1Dg)'), (1.0, 'e')]
    """
    terms = []
    for token in text.split(" + "):
        match = re.fullmatch(r"(\d+(?:\.\d+)?)?(.+)", token.strip())
        coef, name = match.groups()
        terms.append((float(coef) if coef else 1.0, name))
    return terms


def split_reaction(text):
    r"""反応式を、反応物と生成物のリストに分ける（``->`` と ``<->`` を同じ扱い）。

    Parameters
    ----------
    text : str
        LoKI の反応式。

    Returns
    -------
    tuple
        反応物のリスト、生成物のリスト（``parse_side`` の形式）。

    Examples
    --------
    >>> split_reaction("e + CO(X) <-> e + CO(a3P)")[1]
    [(1.0, 'e'), (1.0, 'CO(a3P)')]
    """
    left, right = re.split(r"\s*<?->\s*", text.strip(), maxsplit=1)
    return parse_side(left), parse_side(right)


def canonical(reactants, products):
    r"""反応物・生成物の集合を、順序によらない文字列にする。

    Parameters
    ----------
    reactants, products : list of tuple
        ``parse_side`` の形式。

    Returns
    -------
    str
        種名を並べ替え、係数を展開した ``A + B -> C`` の文字列。

    Examples
    --------
    >>> canonical([(1.0, "e"), (2.0, "A")], [(1.0, "B")])
    'A + A + e -> B'
    """

    def side(terms):
        r"""片側の種名を、係数を展開して並べ替え、`` + `` でつなぐ。

        Parameters
        ----------
        terms : list of tuple
            ``parse_side`` の形式。

        Returns
        -------
        str
            種名を並べ替えてつないだ文字列。
        """
        names = []
        for coef, name in terms:
            names += [name] * int(coef)
        return " + ".join(sorted(names))

    return f"{side(reactants)} -> {side(products)}"


def format_side(terms):
    r"""反応式の片側を、PlasGlo の種名で、同じ種を係数にまとめて書く。

    Parameters
    ----------
    terms : list of tuple
        ``parse_side`` の形式。LoKI の種名（``gas`` は第三体 ``M``）。

    Returns
    -------
    str
        ``2E + CO2^+`` のような文字列。

    Examples
    --------
    >>> format_side([(1.0, "e"), (1.0, "e"), (1.0, "CO2(+,X)")])
    '2E + CO2^+'
    """
    counts = collections.OrderedDict()
    for coef, name in terms:
        species = "M" if name == "gas" else RENAME[name]
        counts[species] = counts.get(species, 0) + int(coef)
    return " + ".join(
        (str(n) if n != 1 else "") + species for species, n in counts.items()
    )


def rate_expression(func, params):
    r"""LoKI の速度係数の関数を、PlasGlo の式（SI、te は K）に変換する。

    Parameters
    ----------
    func : str
        LoKI の関数名。
    params : list of str
        係数の文字列。``^`` は ``**`` に置き換える。

    Returns
    -------
    str
        PlasGlo の ``expression`` の式。``tgas`` [K]、``te`` [K] を使う。

    Raises
    ------
    ValueError
        未対応の関数、または係数の数が合わない場合。

    Examples
    --------
    >>> rate_expression("arrheniusGasTemp", ["7.4e-17", "133"])
    '(7.4e-17)*exp((133)/tgas)'
    """
    p = [f"({x.replace('^', '**')})" for x in params]
    te_ev = f"(te*{KB_EV})"
    forms = {
        "constantRateCoeff": (1, "{0}"),
        "arrheniusGasTemp": (2, "{0}*exp({1}/tgas)"),
        "powerGasTemp": (2, "{0}*tgas**{1}"),
        "modifiedArrheniusGasTemp": (3, "{0}*tgas**{1}*exp({2}/tgas)"),
        "powerGasElecTemp": (3, "{0}*tgas**{1}*" + te_ev + "**{2}"),
        "powerElectronTemp": (2, "{0}*" + te_ev + "**{1}"),
        "offsetArrheniusGasTemp": (3, "{0}+{1}*exp({2}/tgas)"),
        "expGasTemp": (2, "{0}*exp(tgas/{1})"),
    }
    if func not in forms:
        raise ValueError(f"未対応の関数：{func}")
    count, template = forms[func]
    if len(p) != count:
        raise ValueError(f"{func}: 係数が{count}個必要（{len(p)}個）")
    return template.format(*p)


def find_channel(collisions, selector):
    r"""断面積から、電子衝突のチャンネルを1つ見つける。

    Parameters
    ----------
    collisions : list of dict
        Py3bolos のパーサが返す過程のリスト。
    selector : tuple
        (標的, 種類, 生成物に含まれる文字列 または None, しきい値 または None)。

    Returns
    -------
    dict or None
        ``ElectronSwarm`` のチャンネル宣言。一意に見つからなければ None。
    """
    target, kind, product, threshold = selector
    found = [
        c
        for c in collisions
        if c["target"] == target
        and c["kind"] == kind
        and (product is None or product in (c.get("product") or ""))
        and (
            threshold is None
            or abs(float(c.get("threshold", 0.0)) - threshold) < 5e-4
        )
    ]
    if len(found) != 1:
        return None
    c = found[0]
    return dict(
        target=target,
        kind=kind,
        product=c.get("product"),
        threshold_eV=float(c.get("threshold", 0.0)),
    )


def convert(loki_dir, source_dir):
    r"""3つの ``.chem`` を読み、反応表・種・チャンネルと、変換の記録を作る。

    Parameters
    ----------
    loki_dir : pathlib.Path
        ``CHEM_FILES`` の3ファイルがあるフォルダ。
    source_dir : pathlib.Path
        断面積ファイル（``build_co2_channel_draft`` が使う4ファイル）のフォルダ。

    Returns
    -------
    dict
        ``rows``（反応CSVの行）、``species``、``channels``、``used``（使った
        断面積ファイル）、``log``（取り込まなかった反応と理由）、``sources``
        （LoKI の入力ファイルのSHA-256）。
    """
    collisions, used = [], []
    for name in build.sources().values():
        path = source_dir / name
        used.append(path)
        with path.open(encoding="utf-8") as fp:
            collisions += build.parser.parse(fp)
    selectors = {
        canonical(*split_reaction(k)): (k, v)
        for k, v in ELECTRON_SELECTORS.items()
    }
    rows, log, channels, sources = [], collections.defaultdict(list), {}, {}
    for name in CHEM_FILES:
        path = loki_dir / name
        sources[name] = build.sha256(path)
        for item in parse_chem(path):
            text = item["reaction"]
            where = f"{item['file']}:{item['line']}"
            reactants, products = split_reaction(text)
            names = [n for _, n in reactants + products]
            if "wall" in names:
                log["壁・輸送（PlasGlo には無い）"].append((where, text))
                continue
            if any(c != int(c) for c, _ in reactants + products):
                log["非整数の係数"].append((where, text))
                continue
            gas = "gas" in names
            untracked = [
                n for _, n in reactants if n not in RENAME and n != "gas"
            ]
            if untracked:
                log["追跡しない種が反応物（励起種など）"].append((where, text))
                continue
            mapped_products, bad = [], False
            for coef, n in products:
                n = ALIAS_PRODUCTS.get(n, n)
                if n not in RENAME and n != "gas":
                    bad = True
                    break
                mapped_products.append((coef, n))
            if bad:
                log["追跡しない種が生成物"].append((where, text))
                continue
            left, right = canonical(reactants, mapped_products).split(" -> ")
            if left == right:
                log["反応の前後で同じ（EEDFのみ）"].append((where, text))
                continue
            if item["func"] == "eedf":
                hit = selectors.get(canonical(reactants, products))
                if hit is None:
                    log["電子衝突：断面積データなし"].append((where, text))
                    continue
                channel = find_channel(collisions, hit[1])
                if channel is None:
                    log["電子衝突：断面積の過程を特定できない"].append(
                        (where, text)
                    )
                    continue
                cid = (
                    f"{channel['target']}->{channel['kind'][:3]}_"
                    f"{channel['threshold_eV']:.3f}_"
                    f"{(channel['product'] or '').replace(' ', '')}"
                )
                channel = dict(id=cid, **channel, label=hit[0])
                channels[cid] = channel
                kind, value = "electron", cid
            else:
                try:
                    kind, value = "expression", rate_expression(
                        item["func"], item["params"]
                    )
                except ValueError as exc:
                    log[f"未対応：{exc}"].append((where, text))
                    continue

            rows.append(
                dict(
                    id=f"R{len(rows) + 1:03d}",
                    equation=(
                        f"{format_side(reactants)} -> "
                        f"{format_side(mapped_products)}"
                    ),
                    kind=kind,
                    value=value,
                    parameters="",
                    length_unit="m",
                    source_file=f"LoKI-GM/{item['file']}",
                    line=item["line"],
                    reference=item["comment"]
                    + ("（M：全中性種）" if gas else ""),
                    loki=(item["func"], item["params"], text),
                )
            )
    used_names = {
        n for r in rows for n in re.split(r" \+ | -> ", r["equation"])
    }
    used_names = {re.sub(r"^\d+", "", n) for n in used_names} - {"M"}
    species = [n for n in SPECIES if n in used_names]
    return dict(
        rows=rows,
        species=species,
        channels=list(channels.values()),
        used=used,
        log=log,
        sources=sources,
        collisions=collisions,
    )


def write_outputs(result, out):
    r"""反応CSV、ケースJSON、変換の記録を、``out`` に書き出す。

    Parameters
    ----------
    result : dict
        ``convert`` の戻り値。
    out : pathlib.Path
        出力フォルダ（なければ作る）。

    Returns
    -------
    pathlib.Path
        ケースJSONのパス。
    """
    import csv

    out.mkdir(parents=True, exist_ok=True)
    columns = [
        "id",
        "equation",
        "kind",
        "value",
        "parameters",
        "length_unit",
        "source_file",
        "line",
        "reference",
    ]
    table = out / "mechanism_co2.csv"
    with table.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result["rows"])
    # EEDFに必要な、化学反応に使わないチャンネル（三体付着）を加える。
    channels = list(result["channels"])
    three = find_channel(result["collisions"], ("O2", "ATTACHMENT", None, 0.0))
    if three is None:
        for c in result["collisions"]:
            if c["target"] == "O2" and c["kind"] == "ATTACHMENT" and not c.get(
                "product"
            ):
                three = dict(
                    target="O2",
                    kind="ATTACHMENT",
                    product=None,
                    threshold_eV=float(c.get("threshold", 0.0)),
                )
    if three is not None:
        channels.append(
            dict(
                id="O2->ATT_3BODY",
                **three,
                label="O2 三体付着（EEDFのみ）",
                density_normalization=O2_THREE_BODY,
            )
        )
    n0 = 101325 / (1.380649e-23 * 300)
    case = {
        "schema_version": 1,
        "name": "co2_loki_draft",
        "mechanism": {
            "schema_version": 1,
            "species": [
                dict(
                    name=n,
                    composition=SPECIES[n][0],
                    charge=SPECIES[n][1],
                    **(
                        {"electron": True}
                        if n == "E"
                        else {"excited": SPECIES[n][2]}
                    ),
                )
                for n in result["species"]
            ],
            "reactions": "mechanism_co2.csv",
            "third_bodies": {
                "M": {n: 1.0 for n in THIRD_BODY if n in result["species"]}
            },
            "variables": [],
            "notes": (
                "LoKI-GM（GPL-3.0）の CO2_novib.chem、CO_simple_X.chem、"
                "oxygen_novib.chem から変換した叩き台。再配布しない。"
            ),
        },
        "gas": {
            "temperature_K": 300.0,
            "pressure_Pa": 101325.0,
            "mole_fractions": {"CO2": 1.0},
        },
        "initial_densities_m3": {"E": 1e12, "CO2^+": 1e12, "CO2": n0},
        "field": {"td": 100.0, "duration_ns": 20.0},
        "electrons": {
            "profile": {
                # ケースJSONは、パスを自分の場所から解決するため、絶対パスで書く。
                "cross_sections": [str(p.resolve()) for p in result["used"]],
                "targets": {"CO2": "CO2", "CO": "CO", "O2": "O2"},
                "channels": channels,
                "source_sha256": {
                    str(p.resolve()): build.sha256(p) for p in result["used"]
                },
            },
            "update_mode": "coupled",
            "update_ns": 0.5,
            "cells": 400,
            "rtol": 1e-06,
        },
        "solver": {
            "end_ns": 50.0,
            "points": 51,
            "rtol": 1e-08,
            "atol_density_m3": 0.0001,
        },
        "diagnostics_species": [],
        "notes": "叩き台。条件は例であり、実験条件ではない。",
    }
    path = out / "case_co2.json"
    path.write_text(
        json.dumps(case, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# 変換の記録",
        "",
        f"- 取り込んだ反応：{len(result['rows'])}",
        f"- 種：{len(result['species'])}（{', '.join(result['species'])}）",
        f"- 電子衝突のチャンネル：{len(result['channels'])}",
        "- LoKI の入力ファイルの SHA-256：",
    ]
    lines += [f"  - `{k}`: `{v}`" for k, v in result["sources"].items()]
    lines += ["", "## 取り込まなかった反応", ""]
    for reason, items in sorted(result["log"].items()):
        lines.append(f"### {reason}（{len(items)}）")
        lines += [f"- `{w}` {t}" for w, t in items]
        lines.append("")
    (out / "conversion_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return path


def verify_expressions(rows):
    r"""変換した式が、LoKI の定義を直接計算した値と一致することを確かめる。

    Parameters
    ----------
    rows : list of dict
        ``convert`` が作った反応の行。

    Returns
    -------
    float
        全反応・全温度での、相対差の最大値。
    """
    import math

    from plasglo.chemistry.sample_mechanism import Expression

    def direct(func, params, tg, te_ev):
        r"""LoKI の定義式を、変換を介さず直接計算する。

        Parameters
        ----------
        func : str
            LoKI の関数名。
        params : list of str
            係数の文字列。
        tg : float
            気体温度 [K]。
        te_ev : float
            電子温度 [eV]。

        Returns
        -------
        float
            速度係数。
        """
        # 係数の文字列は、ダウンロードしたデータなので eval せず、
        # PlasGlo の制限付き式評価器で数値にする。
        a = [Expression(x.replace("^", "**"))({}) for x in params]
        if func == "constantRateCoeff":
            return a[0]
        if func == "arrheniusGasTemp":
            return a[0] * math.exp(a[1] / tg)
        if func == "powerGasTemp":
            return a[0] * tg ** a[1]
        if func == "modifiedArrheniusGasTemp":
            return a[0] * tg ** a[1] * math.exp(a[2] / tg)
        if func == "powerGasElecTemp":
            return a[0] * tg ** a[1] * te_ev ** a[2]
        if func == "powerElectronTemp":
            return a[0] * te_ev ** a[1]
        if func == "offsetArrheniusGasTemp":
            return a[0] + a[1] * math.exp(a[2] / tg)
        if func == "expGasTemp":
            return a[0] * math.exp(tg / a[1])
        raise ValueError(func)

    worst = 0.0
    for row in rows:
        if row["kind"] != "expression":
            continue
        func, params, _ = row["loki"]
        expr = Expression(row["value"])
        for tg in (300.0, 600.0, 1500.0):
            for te_ev in (0.5, 1.0, 3.0):
                got = expr({"tgas": tg, "tg": tg, "te": te_ev / KB_EV})
                want = direct(func, params, tg, te_ev)
                worst = max(worst, abs(got - want) / abs(want))
    return worst


def main():
    r"""変換を実行し、``--verify`` なら式の検証も行う。"""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--loki-dir", required=True, help="3つの .chem のフォルダ")
    ap.add_argument("--source-dir", help="断面積のフォルダ（CO2_XS_DIR）")
    ap.add_argument("--out", default="outputs/co2_mechanism")
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()
    build.CO2_SET = "ist"
    folder = build.set_source_dir(args.source_dir)
    result = convert(Path(args.loki_dir), folder)
    path = write_outputs(result, Path(args.out))
    print(f"反応 {len(result['rows'])}、種 {len(result['species'])}、"
          f"チャンネル {len(result['channels'])}")
    for reason, items in sorted(result["log"].items()):
        print(f"  取り込まず：{reason}：{len(items)}")
    print("ケース:", path)
    if args.verify:
        print(f"式の検証：最大相対差 {verify_expressions(result['rows']):.2e}")


if __name__ == "__main__":
    main()
