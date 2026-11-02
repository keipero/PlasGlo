r"""モデルを変更せず、係数履歴から反応源項を再構成する。
"""

import unittest

import numpy as np

from py3zdplaschem import Mechanism, Reaction, Result, species_balance


class BalanceTests(unittest.TestCase):
    def test_balance_uses_correct_history_at_switch_without_mutating_model(
        self,
    ):
        r"""電界切替時の係数履歴が正しく選ばれ、モデルを変更しないことを検証する。
        """
        reaction = Reaction(
            "N2*(C) -> N2", "Radiative decay", "c1", "3e7", "m", 1
        )
        model = Mechanism([reaction], None, 1000)
        model.set_swarm(np.zeros(44), 1.0, {})
        densities = np.zeros((len(model.species), 3))
        densities[model.index["N2*(C)"]] = 1e15
        result = Result(
            model,
            np.array([0, 0.5e-9, 1e-9]),
            densities,
            {},
            [
                dict(time_s=0, coefficients_SI=np.array([3e7])),
                dict(time_s=1e-9, coefficients_SI=np.array([6e7])),
            ],
        )
        before = model.k.copy()
        balance = species_balance(result, "N2*(C)")
        np.testing.assert_allclose(balance["loss_m3_s"], [3e22, 3e22, 6e22])
        np.testing.assert_allclose(balance["production_m3_s"], 0)
        np.testing.assert_allclose(
            balance["loss_lifetime_s"], [1 / 3e7, 1 / 3e7, 1 / 6e7]
        )
        np.testing.assert_array_equal(model.k, before)


if __name__ == "__main__":
    unittest.main()
