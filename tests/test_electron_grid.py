r"""E/Nを順に変えたときの、BOLOSエネルギー格子の履歴依存を検査する。

窒素は平均エネルギーがE/Nで大きく変わる（10 Tdで約1 eV、1000 Tdで約16 eV）
ため、格子の上限が以前のE/Nに固定される不具合を検出しやすい。
"""

import unittest

import numpy as np

from plasglo import load_case
from plasglo.electrons.mixture import ElectronSwarm

DENSITY = {"N2": 2.4e25}


def make_swarm():
    r"""窒素の断面積で、100セルの `ElectronSwarm` を作る。

    Returns
    -------
    ElectronSwarm
        新規（履歴のない）電子群計算。

    Examples
    --------
    >>> make_swarm().cells
    100
    """
    profile = load_case("nitrogen").electron_settings()["profile_data"]
    return ElectronSwarm(
        profile["cross_sections"],
        profile["targets"],
        profile["channels"],
        cells=100,
    )


class ElectronGridHistoryTest(unittest.TestCase):
    r"""格子の上限が、以前のE/Nではなく現在の平均エネルギーに従うこと。"""

    def assert_same_state(self, state, reference):
        r"""平均エネルギーと全チャンネルの係数が一致することを確かめる。

        Parameters
        ----------
        state : dict
            検査対象の計算結果。
        reference : dict
            履歴のない新規インスタンスの計算結果。

        Examples
        --------
        >>> True
        True
        """
        self.assertAlmostEqual(
            state["energy_eV"], reference["energy_eV"],
            delta=1e-6 * reference["energy_eV"],
        )
        for name, value in reference["rates"].items():
            self.assertTrue(
                np.isclose(state["rates"][name], value, rtol=1e-5, atol=0.0),
                f"{name}: {state['rates'][name]} != {value}",
            )

    def test_increasing_field_matches_fresh_solution(self):
        r"""低E/Nの後に高E/Nを解いても、新規の解と一致する。"""
        reference = make_swarm().solve(300.0, DENSITY, 300.0)
        swarm = make_swarm()
        swarm.solve(10.0, DENSITY, 300.0)
        state = swarm.solve(300.0, DENSITY, 300.0)
        self.assert_same_state(state, reference)
        self.assertEqual(swarm.regrid_count, 2)

    def test_decreasing_field_matches_fresh_solution(self):
        r"""高E/Nの後に低E/Nを解いても、新規の解と一致する。"""
        reference = make_swarm().solve(10.0, DENSITY, 300.0)
        swarm = make_swarm()
        swarm.solve(1000.0, DENSITY, 300.0)
        state = swarm.solve(10.0, DENSITY, 300.0)
        self.assert_same_state(state, reference)

    def test_constant_field_does_not_regrid(self):
        r"""一定のE/Nを繰り返しても、格子を作り直さない。"""
        swarm = make_swarm()
        swarm.solve(100.0, DENSITY, 300.0)
        count, upper = swarm.regrid_count, swarm.bsolver.grid.x1
        swarm.solve(100.0, DENSITY, 300.0)
        self.assertEqual(swarm.regrid_count, count)
        self.assertEqual(swarm.bsolver.grid.x1, upper)

    def test_zero_field_does_not_use_regrid(self):
        r"""電界0の計算では、格子の作り直しを数えない。"""
        swarm = make_swarm()
        swarm.solve(0.0, DENSITY, 300.0)
        self.assertEqual(swarm.regrid_count, 0)

    def test_grid_upper_limit_follows_mean_energy(self):
        r"""格子の上限が、平均エネルギーの目安の範囲に収まる。"""
        swarm = make_swarm()
        for td in (10.0, 100.0, 400.0, 20.0):
            state = swarm.solve(td, DENSITY, 300.0)
            target = swarm.adaptive_factor * state["energy_eV"]
            self.assertGreaterEqual(
                swarm.bsolver.grid.x1, target / swarm.regrid_ratio
            )
            self.assertLessEqual(
                swarm.bsolver.grid.x1, target * swarm.regrid_ratio
            )


if __name__ == "__main__":
    unittest.main()
