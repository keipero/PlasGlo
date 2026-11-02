r"""電流の符号・単位・変位成分・不足入力を検証するテスト。
"""

import tempfile
from pathlib import Path
import unittest

import numpy as np
from scipy.constants import elementary_charge as e, epsilon_0

from current_diagnostics import (
    calculate_current,
    reduced_field_to_electric_field,
    save_current_result,
)


class CurrentTests(unittest.TestCase):
    def calculate(self, **options):
        r"""明示した上書き条件で保存用の電流テストを計算する。

        Parameters
        ----------
        **options : dict
            保存用の電流テスト条件を上書きするキーワード引数。

        Returns
        -------
        CurrentResult
            テスト条件に対する電流履歴と合計値の完全性を示す情報。
        """
        settings = dict(
            time_s=np.array([0.0, 1e-9, 2e-9]),
            densities_m3=np.array([[1e16] * 3, [1e16] * 3]),
            species=["E", "N2^+"],
            charges=[-1, 1],
            electric_field_V_m=1e4,
            electron_reduced_mobility_1_V_m_s=1e23,
            reference_density_m3=1e24,
            ion_mobilities_m2_V_s={"N2^+": 1e-4},
        )
        settings.update(options)
        return calculate_current(**settings)

    def test_neutral_plasma_has_additive_conduction_and_area_scaling(self):
        r"""中性プラズマの伝導成分の加算と電極面積に対する比例を検証する。
        """
        result = self.calculate(area_m2=2e-6)
        expected = e * 1e16 * 1e4 * (0.1 + 1e-4)
        np.testing.assert_allclose(
            result.traces["total_current_A"], expected * 2e-6, rtol=1e-14
        )
        np.testing.assert_array_equal(
            result.traces["displacement_current_density_A_m2"], 0.0
        )

    def test_field_polarity_reverses_conventional_current(self):
        r"""電界の極性を反転すると慣用電流も反転することを検証する。

        """
        positive = self.calculate()
        negative = self.calculate(electric_field_V_m=-1e4)
        np.testing.assert_allclose(
            negative.traces["total_current_density_A_m2"],
            -positive.traces["total_current_density_A_m2"],
        )

    def test_reduced_mobility_and_field_density_cancel_in_electron_drift(self):
        r"""電子ドリフトで換算移動度と基準密度の因子が相殺することを検証する。
        """
        n = np.array([1e24, 5e23, 2e23])
        result = self.calculate(
            reference_density_m3=n,
            electric_field_V_m=reduced_field_to_electric_field(100.0, n),
        )
        np.testing.assert_allclose(
            result.traces["electron_current_density_A_m2"],
            e * 1e16 * 1e23 * 100e-21,
        )

    def test_linear_ramp_displacement_on_nonuniform_time_grid(self):
        r"""不均一時刻格子における線形電界の変位電流を検証する。

        """
        time = np.array([0, 0.2e-9, 2e-9])
        result = self.calculate(
            time_s=time, electric_field_V_m=1e4 + 3e12 * time
        )
        np.testing.assert_allclose(
            result.traces["displacement_current_density_A_m2"],
            epsilon_0 * 3e12,
            rtol=1e-12,
        )

    def test_explicit_derivative_and_gap_voltage(self):
        r"""明示した電界微分と電極間電圧の関係を検証する。

        """
        result = self.calculate(
            electric_field_derivative_V_m_s=-1e12, electrode_gap_m=0.001
        )
        np.testing.assert_allclose(
            result.traces["displacement_current_density_A_m2"],
            -epsilon_0 * 1e12,
        )
        np.testing.assert_allclose(result.traces["voltage_V"], 10.0)

    def test_missing_ion_mobility_does_not_claim_total(self):
        r"""イオン移動度が不足すると完全な合計値を返さないことを検証する。

        """
        result = self.calculate(ion_mobilities_m2_V_s={}, area_m2=1e-6)
        self.assertFalse(result.total_available)
        self.assertEqual(result.metadata["missing_ion_mobilities"], ["N2^+"])
        self.assertNotIn("total_current_A", result.traces)
        self.assertNotIn("total_current_density_A_m2", result.traces)
        self.assertIn("known_components_current_A", result.traces)

    def test_unpopulated_ion_does_not_prevent_total(self):
        r"""密度ゼロのイオンの移動度不足が完全性を妨げないことを検証する。

        """
        result = self.calculate(
            densities_m3=[[1e16] * 3, [0.0] * 3], ion_mobilities_m2_V_s={}
        )
        self.assertTrue(result.total_available)

    def test_negative_ion_and_multiple_charge_magnitudes(self):
        r"""負イオンと複数の電荷数の扱いを検証する。

        """
        result = self.calculate(
            species=["E", "X-"],
            charges=[-1, -2],
            ion_mobilities_m2_V_s={"X-": 1e-4},
        )
        np.testing.assert_allclose(
            result.traces["ion_current_density_A_m2"],
            e * 2 * 1e16 * 1e-4 * 1e4,
        )

    def test_invalid_physical_inputs_are_rejected(self):
        r"""不正な物理入力が拒否されることを検証する。

        """
        for options in [
            dict(area_m2=0),
            dict(reference_density_m3=0),
            dict(time_s=[0, 0, 1]),
            dict(electron_reduced_mobility_1_V_m_s=-1),
            dict(ion_mobilities_m2_V_s={"N2^+": -1}),
            dict(densities_m3=[[-1] * 3, [1e16] * 3]),
        ]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.calculate(**options)

    def test_saved_output_preserves_incomplete_status(self):
        r"""保存結果に未完の電流評価という状態が残ることを検証する。

        """
        result = self.calculate(ion_mobilities_m2_V_s={})
        with tempfile.TemporaryDirectory() as temp:
            save_current_result(result, temp)
            with np.load(Path(temp) / "current.npz") as data:
                self.assertNotIn("total_current_density_A_m2", data.files)
            self.assertTrue((Path(temp) / "current_metadata.json").is_file())


if __name__ == "__main__":
    unittest.main()
