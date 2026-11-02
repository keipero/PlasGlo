r"""依存クローンの接続、密度依存付着、電界切替を検証する。
"""

from pathlib import Path
import unittest

import numpy as np

from py3zdplaschem.electrons.bolos_backend import BOLOS_DIR, BolosSwarm
from bolos import solver
from py3zdplaschem import Mechanism, read_reactions, simulate


class BolosConnectionTests(unittest.TestCase):
    def test_solver_is_loaded_from_the_clone(self):
        r"""ソルバーが指定クローンから読み込まれることを検証する。

        """
        self.assertEqual(
            Path(solver.__file__).resolve(), BOLOS_DIR / "bolos/solver.py"
        )

    def test_three_body_attachment_matches_eedf_loss_frequency(self):
        r"""三体付着のEEDF損失周波数と化学反応速度の一致を検証する。

        """
        swarm = BolosSwarm()
        density = {"N2": 1.96e25, "O2": 4.9e24}
        state = swarm.solve(100, density)
        raw = swarm.bsolver.rate(swarm.eedf, swarm.attachment)
        # Physical loss frequency must be k3*n(O2)^2, with no extra
        # third body.
        chemical_loss = state["rates"][27] * density["O2"] ** 2
        weight = density["O2"] / sum(density.values()) * density["O2"] * 1e-6
        collision_loss = raw * weight * sum(density.values())
        self.assertAlmostEqual(chemical_loss / collision_loss, 1.0, places=12)
        model = Mechanism(read_reactions(), None, 100)
        model.set_swarm(
            state["rates"],
            state["energy_eV"],
            swarm.labels,
            three_body_si=True,
        )
        index = next(
            i
            for i, reaction in enumerate(model.reactions)
            if reaction.value.startswith("C27 ")
        )
        self.assertAlmostEqual(
            model.k[index] / state["rates"][27], 1.0, places=12
        )

    def test_field_scan_and_pure_gas_endpoints(self):
        r"""換算電界の走査と純気体の両端条件を検証する。

        """
        for field, fraction in [
            (1, 0.8),
            (100, 0.8),
            (2000, 0.8),
            (100, 1.0),
            (100, 0.0),
        ]:
            with self.subTest(field=field, fraction=fraction):
                swarm = BolosSwarm()
                state = swarm.solve(
                    field,
                    {"N2": fraction * 2.45e25, "O2": (1 - fraction) * 2.45e25},
                )
                self.assertTrue(np.isfinite(state["rates"]).all())
                self.assertTrue((state["rates"] >= 0).all())
                self.assertGreater(state["energy_eV"], 0)
                self.assertLessEqual(state["error"], swarm.rtol)

    def test_off_grid_pulse_switch_and_conservation(self):
        r"""出力格子と一致しない電界切替と保存則を検証する。

        """
        result = simulate(
            field_td=1000,
            field_duration_ns=0.05,
            end_ns=0.1,
            seed_density=1e13,
            points=17,
            bolos_update_ns=0.03,
        )
        switch = [
            h
            for h in result.swarm_history
            if np.isclose(h["time_s"], 0.05e-9, rtol=0, atol=1e-24)
        ]
        self.assertEqual(len(switch), 1)
        self.assertEqual(switch[0]["field_td"], 0)
        self.assertEqual(result.swarm_history[0]["field_td"], 1000)
        self.assertEqual(result.swarm_history[-1]["field_td"], 0)
        self.assertEqual(result.diagnostics["reaction_count"], 134)
        self.assertLess(result.diagnostics["max_atom_relative_error"], 1e-9)
        self.assertTrue(np.isfinite(result.densities).all())
        self.assertLess(
            switch[0]["energy_eV"], result.swarm_history[0]["energy_eV"]
        )


if __name__ == "__main__":
    unittest.main()
