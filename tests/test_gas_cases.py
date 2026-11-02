r"""ガス種に依存しない反応計算、入力検証、混合標的、回帰検証。
"""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from py3zdplaschem import load_case, simulate, save_result, species_balance
from py3zdplaschem.chemistry.network import (
    ReactionNetwork,
    NetworkReaction,
    Species,
)
from py3zdplaschem.electrons.mixture import ElectronSwarm


def write_json(path, data):
    r"""構造化したテストデータをUTF-8のJSONで保存する。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。
    data : dict
        この処理で使用する構造化された入力メタデータ。
    """
    path.write_text(json.dumps(data), encoding="utf-8")


def neutral_case(tmp_path):
    r"""数値検証用に構成したH/He/C/O反応。実物理モデルではない。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。

    Returns
    -------
    Case
        一時的な中性気体の数値テスト用ケース。実物理機構ではない。
    """
    species = [
        dict(name="H2", composition={"H": 2}),
        dict(name="H", composition={"H": 1}),
        dict(name="He", composition={"He": 1}),
        dict(name="CO2", composition={"C": 1, "O": 2}),
    ]
    with (tmp_path / "rates.csv").open(
        "w", newline="", encoding="utf-8"
    ) as fp:
        writer = csv.writer(fp)
        writer.writerow(
            ["id", "equation", "kind", "value", "parameters", "length_unit"]
        )
        writer.writerow(
            [
                "dissociation",
                "H2 + He -> 2H + He",
                "analytic",
                "c1",
                "1e-20",
                "m",
            ]
        )
    mechanism = dict(schema_version=1, species=species, reactions="rates.csv")
    write_json(
        tmp_path / "case.json",
        dict(
            schema_version=1,
            name="manufactured",
            mechanism=mechanism,
            gas=dict(temperature_K=300.0, pressure_Pa=101325.0),
            initial_densities_m3=dict(H2=1e20, He=1e22, CO2=2e20),
            electrons=dict(backend="none"),
            solver=dict(end_ns=1e6, points=11, rtol=1e-9, atol_density_m3=1.0),
        ),
    )
    return load_case(tmp_path / "case.json")


def test_other_elements_have_analytic_solution_and_roundtrip_outputs(
    tmp_path, monkeypatch
):
    r"""他元素で解析解との一致と保存・再読込を検証する。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。
    monkeypatch : pytest.MonkeyPatch
        一時的に作業ディレクトリを変更するpytestのfixture。
    """
    case = neutral_case(tmp_path)
    # Relative file paths do not depend on the current working
    # directory.
    monkeypatch.chdir(tmp_path.parent)
    result = simulate(case)
    model = result.mechanism
    expected = 1e20 * np.exp(-1e-20 * 1e22 * result.time)
    np.testing.assert_allclose(
        result.densities[model.index["H2"]], expected, rtol=2e-8
    )
    np.testing.assert_allclose(
        result.densities[model.index["H"]],
        2 * (1e20 - expected),
        rtol=2e-7,
        atol=1e7,
    )
    assert model.elements == ["C", "H", "He", "O"]
    assert result.diagnostics["max_atom_relative_error"] < 1e-12
    balance = species_balance(result, "H2")
    np.testing.assert_allclose(
        balance["net_m3_s"],
        -1e-20 * result.densities[model.index["H2"]] * 1e22,
    )
    save_result(result, tmp_path / "out", balance_species=["H2"])
    rerun = simulate(load_case(tmp_path / "out/case.json"))
    np.testing.assert_array_equal(rerun.densities, result.densities)
    assert (tmp_path / "out/densities.png").is_file()
    assert json.loads((tmp_path / "out/balances.json").read_text())["H2"]


def test_distinct_third_body_efficiencies_match_finite_difference_jacobian():
    r"""第三体効率が異なる反応のJacobianを有限差分と比較する。
    """
    species = [
        Species("H2", {"H": 2}),
        Species("H", {"H": 1}),
        Species("He", {"He": 1}),
    ]
    reactions = [
        NetworkReaction("forward", "H2 + M -> 2H + M", "expression", "1e-20"),
        NetworkReaction("backward", "2H + X -> H2 + X", "expression", "1e-43"),
    ]
    model = ReactionNetwork(
        species,
        reactions,
        temperature_K=300.0,
        pressure_Pa=101325.0,
        third_bodies={"M": {"H2": 2.0, "He": 0.5}, "X": {"He": 3.0}},
    )
    model.set_swarm({}, None)
    y = np.array([0.2, 0.1, 0.7])
    step = 1e-6
    numerical = np.column_stack(
        [
            (
                model.rhs(0, y + np.eye(3)[i] * step)
                - model.rhs(0, y - np.eye(3)[i] * step)
            )
            / (2 * step)
            for i in range(3)
        ]
    )
    np.testing.assert_allclose(
        model.jacobian(0, y), numerical, rtol=1e-8, atol=1e-4
    )
    np.testing.assert_allclose(
        model.conserved @ model.rhs(0, y), 0.0, atol=1e-8
    )


def test_explicit_isotope_species_is_not_parsed_as_a_stoichiometric_prefix():
    r"""同位体の粒子名を化学量論係数と誤認しないことを検証する。
    """
    model = ReactionNetwork(
        [Species("13C", {"C13": 1}), Species("13C2", {"C13": 2})],
        [
            NetworkReaction(
                "association", "2 13C -> 13C2", "expression", "1e-20"
            )
        ],
        temperature_K=300.0,
        pressure_Pa=101325.0,
    )
    np.testing.assert_array_equal(model.stoich[:, 0], [-2, 1])


def test_electron_name_and_multiple_ion_charge_are_explicit_metadata():
    r"""電子の別名と多価イオンの電荷がメタデータで扱えることを検証する。
    """
    model = ReactionNetwork(
        [
            Species("electron", {}, -1, electron=True),
            Species("He+", {"He": 1}, 1),
            Species("He++", {"He": 1}, 2),
        ],
        [
            NetworkReaction(
                "ionization",
                "electron + He+ -> 2electron + He++",
                "expression",
                "1e-20",
            )
        ],
        temperature_K=300.0,
        pressure_Pa=101325.0,
    )
    model.set_swarm({}, None)
    assert model.electron_species == "electron"
    np.testing.assert_array_equal(model.charge, [-1, 1, 2])
    np.testing.assert_allclose(
        model.conserved @ model.rhs(0.0, np.array([0.1, 0.2, 0.3])), 0.0
    )


@pytest.mark.parametrize(
    "change, message",
    [
        ({"gas": {"mole_fractions": {"H2": 0.6, "He": 0.6}}}, "sum to one"),
        ({"initial_densities_m3": {"missing": 1.0}}, "Unknown initial"),
        ({"initial_densities_m3": {"H2": -1.0}}, "nonnegative"),
        ({"solver": {"points": 2.5}}, "integer"),
        ({"field": {"duration_ns": 0.5}}, "requires coupled"),
    ],
)
def test_invalid_case_is_rejected_before_integration(
    tmp_path, change, message
):
    r"""不正なケースが時間積分の前に拒否されることを検証する。

    Parameters
    ----------
    tmp_path : pathlib.Path
        pytestが用意する独立した一時ディレクトリ。
    change : dict
        入力検証を試すため、意図的に与える不正なケース設定。
    message : str
        不正な入力に対して期待するエラーメッセージの一部。
    """
    case = neutral_case(tmp_path)
    with pytest.raises(ValueError, match=message):
        simulate(case.with_updates(**change))


def test_unbalanced_mechanism_is_rejected():
    r"""保存則を満たさない反応機構が拒否されることを検証する。

    """
    with pytest.raises(ValueError, match="conservation"):
        ReactionNetwork(
            [Species("H2", {"H": 2}), Species("H", {"H": 1})],
            [NetworkReaction("bad", "H2 -> H", "expression", "1.")],
            temperature_K=300.0,
            pressure_Pa=101325.0,
        )


def test_no_electron_channel_is_silently_omitted():
    r"""電子衝突チャンネルの不足を暗黙に無視しないことを検証する。

    """
    case = load_case("nitrogen").with_updates(electrons={"backend": "none"})
    with pytest.raises(ValueError, match="Electron reactions require"):
        simulate(case)


def test_mixed_targets_and_order_independent_channel_mapping():
    r"""混合標的と、定義順に依存しないチャンネル対応を検証する。

    """
    profiles = [
        load_case(name).electron_settings()["profile_data"]
        for name in ("argon", "nitrogen")
    ]
    paths = [p for profile in profiles for p in profile["cross_sections"]]
    channels = [c for profile in profiles for c in profile["channels"]]
    targets = {"Ar": "argon_ground", "N2": "nitrogen_ground"}
    densities = {"argon_ground": 2.2e25, "nitrogen_ground": 2.2e24}
    a = ElectronSwarm(paths, targets, channels, cells=100).solve(
        100.0, densities
    )
    b = ElectronSwarm(
        paths, targets, list(reversed(channels)), cells=100
    ).solve(100.0, densities)
    assert a["energy_eV"] == b["energy_eV"]
    assert a["fractions"]["argon_ground"] == pytest.approx(10 / 11)
    assert a["fractions"]["nitrogen_ground"] == pytest.approx(1 / 11)
    assert a["rates"] == b["rates"]
    assert all(
        np.isfinite(value) and value >= 0 for value in a["rates"].values()
    )
    broken = [dict(channels[0], product="not-in-the-data"), *channels[1:]]
    with pytest.raises(ValueError, match="channel"):
        ElectronSwarm(paths, targets, broken, cells=100)


@pytest.mark.parametrize(
    "name, counts",
    [("nitrogen", (21, 96)), ("argon", (3, 2)), ("air", (44, 134))],
)
def test_bundled_cases_share_the_validated_kernel(name, counts):
    r"""同梱ケースが共通の検証済みカーネルを使うことを検証する。

    Parameters
    ----------
    name : str
        この処理で使用する粒子種・ケース・出力などの識別名。
    counts : tuple of int
        期待する粒子種数と反応数。
    """
    model = load_case(name).build_network()
    assert (len(model.species), len(model.reactions)) == counts
    np.testing.assert_array_equal(model.conserved @ model.stoich, 0)
