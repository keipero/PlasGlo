r"""固定反応速度係数を使うサンプル反応機構の時間積分。
"""

import numpy as np
from scipy.integrate import solve_ivp


def integrate_prescribed_rates(
    mechanism, time_s, initial_density_m3, *, rtol=1e-8, atol_density_m3=1e-4
):
    r"""固定された反応速度係数で化学反応の時間発展を積分する。

    サンプル用の反応機構に対して、SI単位系の初期密度を機構内部の
    規格化密度へ変換し、BDF法で常微分方程式を積分する。反応速度係数は
    積分中に更新しない。絶対許容誤差も、密度の規格化に合わせて変換する。

    Parameters
    ----------
    mechanism : object
        `rhs` と `jac` を持つ反応機構。`n0` は規格化に使う基準密度を
        [m^-3] で表す。
    time_s : array-like
        出力時刻 [s]。1次元で、有限値かつ厳密に単調増加する必要がある。
    initial_density_m3 : array-like
        各化学種の初期数密度 [m^-3]。`mechanism` が想定する化学種の順序
        で指定する。
    rtol : float, optional
        BDF積分の相対許容誤差。デフォルトは1e-8。
    atol_density_m3 : float, optional
        SI単位系での絶対許容誤差 [m^-3]。規格化して積分器へ渡す。
        デフォルトは1e-4。

    Returns
    -------
    scipy.integrate.OdeResult
        `scipy.integrate.solve_ivp` が返す積分結果。`y` は基準密度で
        規格化された値なので、SI単位系へ戻すには `result.y * mechanism.n0`
        とする。

    Raises
    ------
    ValueError
        `time_s` が1次元でない、要素が2未満、有限値でない、または
        厳密に単調増加していない場合。

    Notes
    -----
    この関数は、サンプル機構のように電界・気体温度・反応速度係数を
    あらかじめ計算して固定するケースに適している。時間積分中に電子群
    分布や反応速度係数を更新する必要がある場合は、通常の
    `simulate` を使用する。

    Examples
    --------
    >>> from plasglo import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> result = simulate(case)
    >>> model = result.mechanism
    >>> solved = integrate_prescribed_rates(
    ...     model, result.time, result.densities[:, 0],
    ... )
    >>> solved.success
    True
    """
    time = np.asarray(time_s, dtype=float)
    if (
        time.ndim != 1
        or len(time) < 2
        or not np.isfinite(time).all()
        or np.any(np.diff(time) <= 0)
    ):
        raise ValueError("Time must be finite and strictly increasing.")
    return solve_ivp(
        mechanism.rhs,
        (time[0], time[-1]),
        np.asarray(initial_density_m3) / mechanism.n0,
        method="BDF",
        jac=mechanism.jac,
        t_eval=time,
        rtol=rtol,
        atol=atol_density_m3 / mechanism.n0,
    )
