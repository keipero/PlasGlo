r"""反応入力、粒子種の元素組成、化学量論の解析。
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re
from ..paths import CHEMISTRY


def parse_side(text):
    r"""反応式の一方の側を解析し、化学種とその係数をカウントする。

    Parameters
    ----------
    text : str
        反応式の片側の文字列。

    Returns
    -------
    terms : collections.Counter
        化学種をキー、係数を値とするCounterオブジェクト。

    Examples
    --------
    >>> parse_side("N2 + O2")
    Counter({'N2': 1, 'O2': 1})
    >>> parse_side("2N2 + 3O2")
    Counter({'O2': 3, 'N2': 2})
    """
    terms = Counter()
    # 前後の空白を削除し、' + 'で分割して各項を処理
    # ex. "N2 + e" -> ["N2", "e"]
    for term in text.strip().split(" + "):
        if term:  # 空白でない場合
            # 係数と化学種を抽出, ex. "2N2" -> (2, "N2"), "e" -> (1, "e")
            match = re.fullmatch(r"(\d+)?(.+)", term.strip())
            if match is None:  # 係数がない場合は1とする
                continue
            terms[match[2]] += int(match[1] or 1)
    return terms


def composition(name):
    r"""化学種の元素組成と電荷を解析する。

    Parameters
    ----------
    name : str
        化学種の名前。

    Returns
    -------
    tuple[int, int, int]
        窒素原子数、酸素原子数、電荷のタプル。

    Examples
    --------
    >>> composition("N2O+")
    (2, 1, 1)
    >>> composition("O2-")
    (0, 2, -1)
    """
    if name in ("e", "M"):
        # (N, O, charge) = (0, 0, -1) for electron, (0, 0, 0) for M
        return (0, 0, -1 if name == "e" else 0)
    # 化学種名から括弧やイオンの符号を除去
    base = re.split(r"[(*]", name)[0].rstrip("+-")

    # 正規表現で元素とその数を抽出し、カウンターに格納
    atoms = Counter()
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", base):
        atoms[element] += int(count or 1)
    # ex. N2O+ -> (2, 1, 1), O2- -> (0, 2, -1)
    return (
        atoms["N"],
        atoms["O"],
        int(name.endswith("+")) - int(name.endswith("-")),
    )


@dataclass
class Reaction:
    r"""反応のデータクラス。

    Attributes
    ----------
    equation : str
        反応式の文字列。 "e + N2 -> e + N2*(v1)"
    category : str
        反応のカテゴリ。 "EXCITATION"
    kind : str
        反応の種類。 "BOLSIG"
    value : str
        反応の値。 "C3 N2 EXCITATION 0.29 eV"
    length_unit : str
        長さの単位。 "m"
    line : int
        反応式が定義されている行番号。
    source_file : str
        反応式が定義されているファイルのパス。
    """

    equation: str
    category: str
    kind: str
    value: str
    length_unit: str
    line: int
    source_file: str = ""


def read_reactions(path=CHEMISTRY):
    r"""反応リストを読み込み、Reactionオブジェクトのリストを返す。

    Parameters
    ----------
    path : str, optional
        読み込む反応リストのパス, by default CHEMISTRY

    Returns
    -------
    list[Reaction]
        反応オブジェクトのリスト

    Raises
    ------
    ValueError
        反応リストがありません。

    Examples
    --------
    >>> reactions = read_reactions()
    >>> len(reactions)
    134
    >>> reactions[0].kind
    'field_table'
    """
    reactions = []
    active, category = False, ""
    for number, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), 1
    ):
        if line == "reaction_list":
            active = True
        elif active and line.startswith("ignored_species"):
            break
        elif active and line.startswith("# "):
            category = line[2:].strip()
        elif active and "->" in line and not line.startswith("#"):
            equation, kind, value, *unit = map(str.strip, line.split(","))
            reactions.append(
                Reaction(
                    equation,
                    category,
                    kind,
                    value,
                    unit[0] if unit else "m",
                    number,
                    str(Path(path).resolve()),
                )
            )
    if not reactions:
        raise ValueError("反応リストがありません。")
    return reactions
