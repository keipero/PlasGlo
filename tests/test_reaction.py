r"""反応式解析と反応リスト読み込みの単体テスト。
"""

from collections import Counter
from pathlib import Path

import pytest
from py3zdplaschem.chemistry.reactions import (
    parse_side,
    composition,
    Reaction,
    read_reactions,
)

# ============================================================
# parse_side()
# ============================================================


@pytest.mark.parametrize(
    ("equation", "expected"),
    [
        ("N2 + O2", Counter({"N2": 1, "O2": 1})),
        ("2N2 + O2", Counter({"N2": 2, "O2": 1})),
        ("N2 + N2", Counter({"N2": 2})),
        ("3O + 2N", Counter({"O": 3, "N": 2})),
        ("N2+ + e", Counter({"N2+": 1, "e": 1})),
        ("O2- + N2+", Counter({"O2-": 1, "N2+": 1})),
        ("", Counter()),
        ("   ", Counter()),
    ],
)
def test_parse_side(equation, expected):
    r"""反応式の片側の解析を検証する。

    Parameters
    ----------
    equation : str
        テスト対象の反応式の片側。
    expected : object
        パラメーター化したテストで期待する解析結果。
    """
    assert parse_side(equation) == expected


def test_parse_side_whitespace():
    r"""反応式の前後の空白を正しく処理することを検証する。

    """
    result = parse_side("  2N2 + O2  ")

    assert result == Counter({"N2": 2, "O2": 1})


# ============================================================
# composition()
# ============================================================


@pytest.mark.parametrize(
    ("species", "expected"),
    [
        ("e", (0, 0, -1)),
        ("M", (0, 0, 0)),
        ("N", (1, 0, 0)),
        ("N2", (2, 0, 0)),
        ("O", (0, 1, 0)),
        ("O2", (0, 2, 0)),
        ("NO", (1, 1, 0)),
        ("N2O", (2, 1, 0)),
        ("NO2", (1, 2, 0)),
        ("N2+", (2, 0, 1)),
        ("O2-", (0, 2, -1)),
        ("NO+", (1, 1, 1)),
        ("N2(A)", (2, 0, 0)),
        ("O2(a)", (0, 2, 0)),
    ],
)
def test_composition(species, expected):
    r"""粒子種の元素組成と電荷の解析を検証する。

    Parameters
    ----------
    species : str
        生成・消失または密度を調べる粒子種名。
    expected : object
        パラメーター化したテストで期待する解析結果。
    """
    assert composition(species) == expected


# ============================================================
# Reaction dataclass
# ============================================================


def test_reaction_creation():
    r"""反応オブジェクトの生成を検証する。

    """
    reaction = Reaction(
        equation="N2 + e -> N2+ + 2e",
        category="Ionization",
        kind="BOLSIG",
        value="test",
        length_unit="m",
        line=10,
    )

    assert reaction.equation == "N2 + e -> N2+ + 2e"
    assert reaction.category == "Ionization"
    assert reaction.kind == "BOLSIG"
    assert reaction.value == "test"
    assert reaction.length_unit == "m"
    assert reaction.line == 10
    assert reaction.source_file == ""


def test_reaction_source_file():
    r"""反応の元ファイル情報を検証する。

    """
    reaction = Reaction(
        equation="N2 -> N + N",
        category="Dissociation",
        kind="constant",
        value="1.0e-15",
        length_unit="m",
        line=5,
        source_file="chemistry.txt",
    )

    assert reaction.source_file == "chemistry.txt"


# ============================================================
# read_reactions()
# ============================================================


def test_read_reactions(tmp_path: Path):
    r"""反応リストの読み込みを検証する。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。
    """
    chemistry = tmp_path / "chemistry.txt"

    chemistry.write_text(
        """
header
reaction_list
# Ionization
N2 + e -> N2+ + 2e, BOLSIG, ionization, m
# Dissociation
O2 + e -> O + O + e, BOLSIG, dissociation, m
ignored_species
""".strip(),
        encoding="utf-8",
    )

    reactions = read_reactions(chemistry)

    assert len(reactions) == 2

    assert reactions[0].equation == "N2 + e -> N2+ + 2e"
    assert reactions[0].category == "Ionization"
    assert reactions[0].kind == "BOLSIG"
    assert reactions[0].value == "ionization"
    assert reactions[0].length_unit == "m"
    assert reactions[0].line == 4
    assert reactions[0].source_file == str(chemistry.resolve())

    assert reactions[1].equation == "O2 + e -> O + O + e"
    assert reactions[1].category == "Dissociation"


def test_read_reactions_default_unit(tmp_path: Path):
    r"""単位省略時に既定の単位を使用することを検証する。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。
    """
    chemistry = tmp_path / "chemistry.txt"

    chemistry.write_text(
        "reaction_list\nN2 -> N + N, constant, 1.0e-15\n",
        encoding="utf-8",
    )

    reactions = read_reactions(chemistry)

    assert len(reactions) == 1
    assert reactions[0].length_unit == "m"


def test_read_reactions_ignored_species(tmp_path: Path):
    r"""ignored_species以降を反応として読み込まないことを検証する。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。
    """
    chemistry = tmp_path / "chemistry.txt"

    chemistry.write_text(
        "reaction_list\n"
        "N2 -> N + N, constant, 1.0e-15\n"
        "ignored_species\n"
        "O2 -> O + O, constant, 1.0e-15\n",
        encoding="utf-8",
    )

    reactions = read_reactions(chemistry)

    assert len(reactions) == 1
    assert reactions[0].equation == "N2 -> N + N"


def test_read_reactions_empty(tmp_path: Path):
    r"""空の反応リストを拒否することを検証する。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。
    """
    chemistry = tmp_path / "chemistry.txt"

    chemistry.write_text(
        "reaction_list\nignored_species\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="反応リストがありません"):
        read_reactions(chemistry)


def test_read_reactions_missing_file(tmp_path: Path):
    r"""存在しない反応ファイルへのエラーを検証する。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。
    """
    with pytest.raises(FileNotFoundError):
        read_reactions(tmp_path / "missing.txt")
