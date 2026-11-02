"""反応寄与の単位、係数切替、経路グラフ、保存形式を検証する。"""

import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from plasglo import (
    Result, analyze_pathways, save_pathways, species_balance,
)
from plasglo.chemistry.network import (
    NetworkReaction, ReactionNetwork, Species,
)


def _result():
    """既知の反応速度を持つ三化学種・四反応の結果を作る。

    Returns
    -------
    Result
        速度が順に20、10、12、2 m^-3 s^-1となる密度一定の結果。
    """
    equations = ["A -> B", "A -> C", "B -> C", "C -> A"]
    model = ReactionNetwork(
        [Species(name, {"N": 1}) for name in ("A", "B", "C")],
        [
            NetworkReaction(f"R{j}", equation, "expression", value)
            for j, (equation, value) in enumerate(
                zip(equations, ("2", "1", "4", "1")), 1,
            )
        ],
        temperature_K=300, pressure_Pa=101325,
    )
    model.set_swarm({}, None)
    return Result(
        model, np.array([0.0, 1.0, 2.0]),
        np.array([[10.0] * 3, [3.0] * 3, [2.0] * 3]), {},
    )


class PathwayTests(unittest.TestCase):
    """反応速度の再構成から積算・出力までの診断を検証する。"""

    def test_rates_rankings_and_existing_balance(self):
        """既知の速度と順位が得られ、既存の生成消滅診断と一致する。"""
        result = _result()
        before = result.mechanism.k.copy()
        density = result.densities.copy()
        analysis = analyze_pathways(result, "A")
        np.testing.assert_allclose(
            analysis.rates_m3_s, np.array([[20], [10], [12], [2]]) @
            np.ones((1, 3)),
        )
        np.testing.assert_allclose(
            analysis.integrated_contributions_m3, [-40, -20, 0, 4],
        )
        self.assertEqual(analysis.production[0]["reaction_id"], "R4")
        self.assertEqual(analysis.loss[0]["reaction_id"], "R1")
        self.assertAlmostEqual(analysis.loss[0]["fraction"], 2 / 3)
        balance = species_balance(result, "A")
        np.testing.assert_allclose(
            np.maximum(analysis.contributions_m3_s, 0).sum(axis=0),
            balance["production_m3_s"],
        )
        np.testing.assert_allclose(
            -np.minimum(analysis.contributions_m3_s, 0).sum(axis=0),
            balance["loss_m3_s"],
        )
        np.testing.assert_array_equal(result.mechanism.k, before)
        np.testing.assert_array_equal(result.densities, density)

    def test_switch_between_saved_times_and_at_interval_end(self):
        """未保存の係数切替点で区間を分け、右端の切替を平均しない。"""
        result = _result()
        initial = result.mechanism.k.copy()
        result.swarm_history = [
            dict(time_s=0.0, coefficients_SI=initial),
            dict(time_s=0.75, coefficients_SI=initial * 3),
            dict(time_s=2.0, coefficients_SI=initial * 100),
        ]
        analysis = analyze_pathways(result, "B")
        np.testing.assert_array_equal(analysis.time_s, [0, 0.75, 1, 2])
        np.testing.assert_allclose(analysis.rates_m3_s[0], [20, 60, 60, 2000])
        np.testing.assert_allclose(
            analysis.integrated_rates_m3, np.array([20, 10, 12, 2]) * 4.5,
        )
        partial = analyze_pathways(result, "B", (0.5, 1.5))
        np.testing.assert_allclose(
            partial.integrated_rates_m3, np.array([20, 10, 12, 2]) * 2.5,
        )

    def test_partial_interval_interpolates_densities(self):
        """未保存の区間端点を線形補間して速度を積算する。"""
        result = _result()
        result.densities[0] = [0, 10, 20]
        analysis = analyze_pathways(result, "A", (0.25, 1.75))
        np.testing.assert_allclose(analysis.time_s, [0.25, 1, 1.75])
        self.assertAlmostEqual(analysis.integrated_rates_m3[0], 30)
        self.assertAlmostEqual(analysis.density_change_m3, 15)

    def test_integral_matches_analytic_decay_density_change(self):
        """一次反応の解析解から、積算寄与と密度差の一致を確認する。"""
        result = _result()
        result.time = np.linspace(0, 1, 1001)
        a = 10 * np.exp(-3 * result.time)
        result.densities = np.array([a, 0 * a, 0 * a])
        result.mechanism.k[:] = [2, 1, 0, 0]
        analysis = analyze_pathways(result, "A")
        self.assertAlmostEqual(
            analysis.integrated_contributions_m3.sum()
            / analysis.density_change_m3, 1, places=5,
        )

    def test_graph_preserves_catalyst_and_third_body_rate(self):
        """電子の両側参加と第三体効率を保ち、架空の移動辺を作らない。"""
        model = ReactionNetwork(
            [Species("A", {"N": 1}), Species("B", {"N": 1}),
             Species("E", {}, charge=-1, electron=True)],
            [
                NetworkReaction("exc", "A + E -> B + E", "expression", "2"),
                NetworkReaction("col", "A + M -> B + M", "expression", "1"),
            ],
            temperature_K=300, pressure_Pa=101325,
            third_bodies={"M": {"A": 2, "B": 0.5}},
        )
        model.set_swarm({}, None)
        result = Result(
            model, np.array([0.0, 1.0]),
            np.array([[10, 10], [4, 4], [3, 3]]), {},
        )
        analysis = analyze_pathways(result, "E")
        np.testing.assert_allclose(analysis.integrated_rates_m3, [60, 220])
        self.assertEqual(analysis.production, [])
        self.assertEqual(analysis.loss, [])
        edges = analysis.graph["edges"]
        self.assertIn(
            dict(source="species:E", target="reaction:exc",
                 stoichiometry=1, integrated_m3=60.0), edges,
        )
        self.assertIn(
            dict(source="reaction:exc", target="species:E",
                 stoichiometry=1, integrated_m3=60.0), edges,
        )
        self.assertNotIn(
            "species:M", [n["id"] for n in analysis.graph["nodes"]],
        )

    def test_zero_rates_have_empty_rankings_and_graph(self):
        """全反応が停止している場合もゼロ除算せず空の順位を返す。"""
        result = _result()
        result.mechanism.k[:] = 0
        analysis = analyze_pathways(result, "B")
        self.assertEqual(analysis.production, [])
        self.assertEqual(analysis.loss, [])
        self.assertEqual(analysis.graph, dict(nodes=[], edges=[]))

    def test_graph_and_contribution_use_stoichiometric_multiplicity(self):
        """二量体化の粒子数係数を、正味寄与とグラフの辺に反映する。"""
        model = ReactionNetwork(
            [Species("A", {"N": 1}), Species("B", {"N": 2})],
            [NetworkReaction("dim", "2A -> B", "expression", "1")],
            temperature_K=300, pressure_Pa=101325,
        )
        model.set_swarm({}, None)
        result = Result(
            model, np.array([0.0, 1.0]), np.array([[3, 3], [0, 0]]), {},
        )
        analysis = analyze_pathways(result, "A")
        np.testing.assert_allclose(analysis.integrated_rates_m3, [9])
        np.testing.assert_allclose(analysis.integrated_contributions_m3, [-18])
        edge = analysis.graph["edges"][0]
        self.assertEqual(edge["stoichiometry"], 2)
        self.assertAlmostEqual(edge["integrated_m3"], 18)

    def test_invalid_density_and_time_grid_are_rejected(self):
        """密度の形状不一致、非有限値、不正な時刻を拒否する。"""
        for time, density in (
            ([0, 0, 1], np.ones((3, 3))),
            ([0, 1, np.inf], np.ones((3, 3))),
            ([0, 1, 2], np.ones((3, 2))),
            ([0, 1, 2], np.full((3, 3), np.nan)),
        ):
            result = _result()
            result.time = np.array(time)
            result.densities = density
            with self.assertRaises(ValueError):
                analyze_pathways(result, "A")

    def test_invalid_ranges_and_history_are_rejected(self):
        """範囲外の積算、未知の化学種、不正な係数履歴を拒否する。"""
        result = _result()
        for bounds in ((-1, 1), (1, 3), (1, 1), (1, 0), (0, np.nan)):
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                analyze_pathways(result, "A", bounds)
        with self.assertRaises(ValueError):
            analyze_pathways(result, "unknown")
        for history in (
            [dict(time_s=1, coefficients_SI=result.mechanism.k)],
            [dict(time_s=0, coefficients_SI=np.zeros(2))],
            [dict(time_s=0, coefficients_SI=np.full(4, np.nan))],
            [dict(time_s=0, coefficients_SI=result.mechanism.k)] * 2,
        ):
            result.swarm_history = history
            with self.assertRaises(ValueError):
                analyze_pathways(result, "A")

    def test_saved_arrays_rankings_and_graph(self):
        """保存した配列と順位表、グラフの参照・重みを確認する。"""
        analysis = analyze_pathways(_result(), "A")
        with tempfile.TemporaryDirectory() as directory:
            save_pathways(analysis, directory)
            path = Path(directory)
            with np.load(path / "reaction_rates.npz") as arrays:
                np.testing.assert_allclose(
                    arrays["integrated_rates_m3"], [40, 20, 24, 4],
                )
                self.assertEqual(arrays["species"].item(), "A")
            with (path / "reaction_ranking.csv").open(
                encoding="utf-8", newline="",
            ) as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["direction"], "production")
            self.assertEqual(rows[0]["reaction_id"], "R4")
            graph = json.loads(
                (path / "pathway_graph.json").read_text(encoding="utf-8"),
            )
            ids = {n["id"] for n in graph["nodes"]}
            self.assertTrue(all(
                e["source"] in ids and e["target"] in ids
                for e in graph["edges"]
            ))
            self.assertAlmostEqual(graph["integrated_net_m3"], -56)


if __name__ == "__main__":
    unittest.main()
