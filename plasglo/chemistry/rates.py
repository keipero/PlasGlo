r"""重粒子反応と電子温度依存反応の解析的速度式。
"""

from __future__ import annotations

import numpy as np


def analytic_rate(kind, values, tg, te, td):
    r"""解析的な反応速度係数を計算する。

    `kind` で指定された速度式に、空白区切りの係数文字列 `values` と
    気体温度・電子温度・換算電界を代入する。係数の単位は反応機構の
    入力データで定義されている単位を使用する。

    対応している速度式は、定数、気体温度 `Tg` または電子温度 `Te` に
    依存する式、および換算電界 `Td` に依存する式である。`c1`、`c2`
    などは `values` に指定した係数を表す。

    Parameters
    ----------
    kind : str
        使用する速度式の名前。例えば ``'c1'``、
        ``'c1*exp(-c2/Tg)'``、``'c1*(300/Te)**c2'`` などを指定する。
    values : str
        `kind` の式に使う係数を空白で区切った文字列。例えば
        ``'3.0e-10'`` や ``'1.0e-12 2.5e3'`` を指定する。
    tg : float
        気体温度 `Tg` [K]。
    te : float
        電子温度 `Te` [K]。
    td : float
        換算電界 `Td` [Td]。

    Returns
    -------
    float
        指定した条件における反応速度係数。

    Raises
    ------
    ValueError
        `kind` が未対応の場合、または `values` の係数を浮動小数点数へ
        変換できない場合。

    Examples
    --------
    定数の速度係数を計算する。

    >>> analytic_rate('c1', '3.0e-10', 300.0, 2.0, 100.0)
    3e-10

    気体温度に依存する速度係数を計算する。

    >>> round(
    ...     analytic_rate(
    ...         "c1*(300/Tg)**c2", "1e-10 2", 600.0, 2.0, 100.0
    ...     ),
    ...     20,
    ... )
    2.5e-11
    """
    c = [*map(float, values.split()), 0.0, 0.0, 0.0]
    a, b, d, f = c[:4]
    if kind == "c1":
        return a
    if kind == "c1*exp(-c2/Tg)":
        return a * np.exp(-b / tg)
    if kind == "c1*(Tg/300)**c2*exp(-c3/Tg)":
        return a * (tg / 300) ** b * np.exp(-d / tg)
    if kind == "c1*(300/Tg)**c2":
        return a * (300 / tg) ** b
    if kind == "c1*(300/Te)**c2":
        return a * (300 / te) ** b
    if kind == "c1*(Tg/c2)**c3":
        return a * (tg / b) ** d
    if kind == "c1*exp(-(c2/(c3+Td))**c4)":
        return a * np.exp(-((b / (d + td)) ** f))
    if kind == "c1*exp(-(Td/c2)**c3)":
        return a * np.exp(-((td / b) ** d))
    if kind == "10**(c1+c2*(Tg-300))":
        return 10 ** (a + b * (tg - 300))
    if kind == "c1*(300/Tg)**c2*exp(-c3/Tg)":
        return a * (300 / tg) ** b * np.exp(-d / tg)
    if kind == "c1*Tg**c2*exp(-c3/Tg)":
        return a * tg**b * np.exp(-d / tg)
    raise ValueError(f"未対応の反応係数形式：{kind}")
