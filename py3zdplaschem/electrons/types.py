"""電子群計算の状態と履歴に共通する辞書の型。"""

from typing import TypedDict

import numpy as np


class SwarmState(TypedDict):
    """電子群計算が返す状態。

    Attributes
    ----------
    rates : dict of str to float
        チャンネル識別子ごとの反応速度係数。
    energy_eV : float
        平均電子エネルギー [eV]。
    iterations : int
        EEDFの収束に要した反復数。
    error : float
        EEDFの収束誤差。
    uncorrected_eedf_min : float
        末尾補正前のEEDFの最小値。
    eedf_tail_repaired : bool
        EEDFの末尾を補正したかどうか。
    target_density_m3 : float
        衝突標的の総数密度 [m^-3]。
    fractions : dict of str to float
        各衝突標的のモル分率。
    """

    rates: dict[str, float]
    energy_eV: float
    iterations: int
    error: float
    uncorrected_eedf_min: float
    eedf_tail_repaired: bool
    target_density_m3: float
    fractions: dict[str, float]


class SwarmHistory(SwarmState):
    """電子群の状態に積分時点の条件を加えた履歴レコード。

    Attributes
    ----------
    time_s : float
        積分開始からの時刻 [s]。
    field_td : float
        換算電界 [Td]。
    electron_temperature_K : float
        電子温度 [K]。
    target_densities_m3 : dict of str to float
        各衝突標的の数密度 [m^-3]。
    unsupported_neutral_fraction : float
        衝突データに含まれない中性粒子の割合。
    coefficients_SI : numpy.ndarray
        化学反応に用いたSI単位系の速度係数。
    """

    time_s: float
    field_td: float
    electron_temperature_K: float
    target_densities_m3: dict[str, float]
    unsupported_neutral_fraction: float
    coefficients_SI: np.ndarray
