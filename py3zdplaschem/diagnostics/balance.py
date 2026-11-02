r"""反応速度係数の履歴を使った生成・消滅速度の診断。
"""

from __future__ import annotations

import numpy as np


def species_balance(result, species):
    r"""指定した化学種の生成速度と消滅速度を再構成する。

    電界印加中の軌跡には最終時刻の `model.k` だけを適用せず、各時刻で
    有効だった反応速度係数を使う。BOLOSとの結合境界では、境界後の係数が
    その区間の右辺を表す。この診断処理は反応機構や積分済みの解を変更しない。

    Parameters
    ----------
    result : Result
        `simulate` が返す計算結果。`mechanism`、`time`、`densities`、
        `swarm_history` を参照する。
    species : str
        生成・消滅速度を求める化学種名。`result.mechanism.index` に
        登録されている必要がある。

    Returns
    -------
    dict[str, numpy.ndarray]
        時刻ごとの診断結果を格納する辞書。キーは次の通り。

        ``production_m3_s``
            生成速度 [m^-3 s^-1]。
        ``loss_m3_s``
            消滅速度 [m^-3 s^-1]。
        ``net_m3_s``
            正味の生成速度 [m^-3 s^-1]。生成速度から消滅速度を引いた値。
        ``loss_lifetime_s``
            消滅速度に基づく寿命 [s]。密度または消滅速度が正でない時刻は
            `numpy.nan`。

    Raises
    ------
    ValueError
        指定した化学種が登録されていない場合、または出力時刻より前の
        反応速度係数履歴がない場合。

    Notes
    -----
    `swarm_history` がある場合は各出力時刻以前で最後に更新された係数を
    使う。履歴がない場合は `result.mechanism.k` を全時刻で使う。
    反応速度は内部の規格化密度をSI単位系へ戻して計算する。

    Examples
    --------
    >>> from py3zdplaschem import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> result = simulate(case)
    >>> model = result.mechanism
    >>> balance = species_balance(result, "AR^+")
    >>> balance["loss_lifetime_s"].shape == result.time.shape
    True
    """
    model = result.mechanism
    # 対象化学種の行から、その化学種が増減する反応だけを選ぶ。
    if species not in model.index:
        raise ValueError(f"未知の粒子種：{species}")
    row = model.stoich[model.index[species]]
    active = row != 0
    production = np.zeros(len(result.time))
    loss = np.zeros(len(result.time))

    # 各出力時刻に対応する、直前の反応速度係数更新位置を求める。
    if result.swarm_history:
        history_times = np.array([h["time_s"] for h in result.swarm_history])
        positions = (
            np.searchsorted(history_times, result.time, side="right") - 1
        )
        if np.any(positions < 0):
            raise ValueError("出力時刻より前の反応係数履歴がありません。")

    for j in range(len(result.time)):
        # BOLOS履歴がある場合は時刻ごとの係数、ない場合は固定係数を使う。
        k = (
            result.swarm_history[positions[j]]["coefficients_SI"]
            if result.swarm_history
            else model.k
        )

        # SI密度を規格化し、質量作用則から各反応の反応速度を計算する。
        q = np.maximum(result.densities[:, j] / model.density, 0.0)
        rates = model.physical_reaction_rates(q, k)[active]

        # 化学量論係数を反応速度に掛け、生成分と消滅分に分ける。
        contributions = row[active] * rates
        production[j] = contributions[contributions > 0].sum()
        loss[j] = -contributions[contributions < 0].sum()

    # 消滅速度が正の時刻だけ、密度 / 消滅速度から寿命を求める。
    density = result.densities[model.index[species]]
    lifetime = np.full(len(result.time), np.nan)
    np.divide(density, loss, out=lifetime, where=(density > 0) & (loss > 0))
    return dict(
        production_m3_s=production,
        loss_m3_s=loss,
        net_m3_s=production - loss,
        loss_lifetime_s=lifetime,
    )
