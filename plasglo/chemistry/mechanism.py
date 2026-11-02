r"""空気機構の質量作用則による反応計算と解析的Jacobian。
"""

from __future__ import annotations

import re
import numpy as np
from ..constants import KB, KB_EV
from ..electrons.swarm_table import SwarmData
from .reactions import parse_side, composition
from .rates import analytic_rate
from .kernel import MassActionKernel


class Mechanism(MassActionKernel):
    r"""dn/dt = S R。密度を初期気体密度で規格化して積分する。

    """

    def __init__(
        self,
        reactions,
        swarm,
        field_td,
        pressure_pa=101325.0,
        *,
        gas_temperature_k=300.0,
        n2_fraction=0.8,
    ):
        r"""反応機構を初期化する。

        反応式から化学種、化学量論係数、反応次数を構築し、指定された
        圧力・気体温度・N2モル分率から初期気体密度を計算する。
        `swarm` が指定されている場合は、`field_td` における電子衝突
        速度係数と平均電子エネルギーも初期化する。

        Parameters
        ----------
        reactions : sequence
            `Reaction` オブジェクトのシーケンス。各反応式は `->` で
            反応物と生成物を区切る。
        swarm : SwarmData or object or None
            電子群速度データ。`SwarmData` を渡す場合は、N2:O2=80:20、
            気体温度300 Kの固定テーブルだけが使用できる。`None` の場合は
            電子群データによる初期化を行わない。
        field_td : float
            換算電界 E/N [Td]。0以上の有限値を指定する。
        pressure_pa : float, optional
            気体圧力 [Pa]。正の有限値を指定する。デフォルトは101325.0。
        gas_temperature_k : float, optional
            気体温度 [K]。正の有限値を指定する。デフォルトは300.0。
        n2_fraction : float, optional
            N2のモル分率。0以上1以下を指定する。O2のモル分率は
            `1 - n2_fraction` で計算される。デフォルトは0.8。

        Raises
        ------
        ValueError
            `pressure_pa` が正の有限値でない場合。
        ValueError
            `field_td` が0以上の有限値でない場合。
        ValueError
            `gas_temperature_k` が正の有限値でない場合。
        ValueError
            `n2_fraction` が0以上1以下でない場合。
        ValueError
            固定 `SwarmData` に対して、気体温度またはN2モル分率が
            対応する条件と異なる場合。
        ValueError
            反応式で第三体 `M` が保存されていない場合。
        ValueError
            反応式で元素または電荷が保存されていない場合。

        Examples
        --------
        反応リストを読み込み、電子群データを使わずに初期化する。

        >>> from plasglo import Mechanism, read_reactions
        >>> model = Mechanism(read_reactions(), None, field_td=100.0)
        >>> model.initial_fractions
        {'N2': 0.8, 'O2': 0.19999999999999996}

        圧力と気体温度、N2モル分率を指定して初期化する。

        >>> model = Mechanism(
        ...     read_reactions(),
        ...     None,
        ...     field_td=100.0,
        ...     pressure_pa=2.0e5,
        ...     gas_temperature_k=400.0,
        ...     n2_fraction=0.7,
        ... )
        >>> model.initial_fractions
        {'N2': 0.7, 'O2': 0.30000000000000004}
        """
        if not np.isfinite(pressure_pa) or pressure_pa <= 0:
            raise ValueError("圧力は正の有限値にしてください。")
        if not np.isfinite(field_td) or field_td < 0:
            raise ValueError("E/Nは0以上の有限値にしてください。")
        if not np.isfinite(gas_temperature_k) or gas_temperature_k <= 0:
            raise ValueError("気体温度は正の有限値にしてください。")
        if not np.isfinite(n2_fraction) or not 0 <= n2_fraction <= 1:
            raise ValueError(
                "N2のモル分率は0〜1にしてください。O2は残りの割合です。"
            )
        if isinstance(swarm, SwarmData) and (
            gas_temperature_k != 300.0 or n2_fraction != 0.8
        ):
            raise ValueError(
                "固定BOLSIGテーブルはN2:O2=80:20、300 K専用です。"
                "BOLOSを使用してください。"
            )

        self.reactions, self.field_td = reactions, field_td
        self.active_field_td = field_td
        self.tg, self.pressure_pa = float(gas_temperature_k), pressure_pa
        self.initial_fractions = {
            "N2": float(n2_fraction),
            "O2": float(1 - n2_fraction),
        }
        self.density = pressure_pa / (KB * self.tg)
        self.te, self.energy = np.nan, np.nan
        sides = [
            (
                parse_side(r.equation.split("->")[0]),
                parse_side(r.equation.split("->")[1]),
            )
            for r in reactions
        ]
        names = {
            s for left, right in sides for s in left.keys() | right.keys()
        } - {"M"}
        self.species = ["e", "N2", "O2"] + sorted(names - {"e", "N2", "O2"})
        self.index = {name: i for i, name in enumerate(self.species)}
        ns, nr = len(self.species), len(reactions)
        self.stoich = np.zeros((ns, nr), dtype=int)
        self.powers = np.zeros((ns, nr), dtype=int)
        self.third_body = np.zeros(nr, dtype=int)
        self.k = np.zeros(nr)
        self.orders = np.zeros(nr, dtype=int)

        # 反応式の各側を解析し、化学種の係数と反応次数を計算する。
        for j, (reaction, (left, right)) in enumerate(zip(reactions, sides)):
            self.orders[j] = sum(left.values())
            if left["M"] != right["M"]:
                raise ValueError("第三体Mが保存されていません。")
            self.third_body[j] = left["M"]
            for name in self.species:
                i = self.index[name]
                self.powers[i, j] = left[name]
                self.stoich[i, j] = right[name] - left[name]
        self.conserved = np.array(
            [composition(name) for name in self.species]
        ).T
        if np.any(self.conserved @ self.stoich):
            raise ValueError("元素または電荷を保存しない反応があります。")
        self.collider_weights = np.array(
            [float(s in ("N2", "O2")) for s in self.species]
        )
        self.charge = self.conserved[-1]
        self.scaled_k = np.zeros(nr)
        if swarm is not None:
            rates, energy = swarm.interpolate(field_td)
            self.set_swarm(rates, energy, swarm.labels)

    def set_swarm(
        self, rates, energy, labels, three_body_si=False, field_td=None
    ):
        r"""電子群データから反応速度係数を設定する。

        電子平均エネルギーから電子温度を計算し、反応の `kind` に応じて
        速度係数を求める。`field_table` 反応では `rates` と `labels` の
        衝突過程番号・ラベルを照合し、それ以外の反応では解析式を評価する。
        設定後は、気体密度で規格化した反応速度に使う `scaled_k` も更新する。

        Parameters
        ----------
        rates : array-like
            衝突過程ごとの反応速度係数。`field_table` 反応で指定された
            `C` 番号を添字として参照できる配列を指定する。
        energy : float
            平均電子エネルギー [eV]。正の有限値を指定する。この値から
            電子温度を計算する。
        labels : mapping
            衝突過程番号をラベルに対応付ける辞書。`field_table` 反応の
            `value` と一致していなければならない。
        three_body_si : bool, optional
            3体付着反応の係数がSI単位系に変換済みかどうか。`True` の
            場合はそのまま使用し、`False` の場合は該当する係数を変換する。
            デフォルトは `False`。
        field_td : float, optional
            解析式で使う換算電界 E/N [Td]。`None` の場合は、コンストラクタ
            で指定した値 `self.field_td` を使う。

        Raises
        ------
        ValueError
            `energy` が正の有限値でない場合。
        ValueError
            `field_table` 反応の衝突過程指定が `C番号 ラベル` の形式でない場合。
        ValueError
            `labels` に登録された衝突過程ラベルが反応の指定と一致しない場合。
        ValueError
            解析式で未知の長さ単位が指定されている場合。
        ValueError
            計算された反応係数に非有限値または負の値が含まれる場合。

        Examples
        --------
        反応機構を作成した後、電子群テーブルからデータを取得して設定する。

        >>> from plasglo import (
        ...     Mechanism,
        ...     read_reactions,
        ...     read_swarm,
        ... )
        >>> model = Mechanism(read_reactions(), None, field_td=100.0)
        >>> swarm = read_swarm()
        >>> rates, energy = swarm.interpolate(100.0)
        >>> model.set_swarm(rates, energy, swarm.labels)

        `field_td` を指定すると、解析式で使う換算電界を上書きできる。

        >>> model.set_swarm(rates, energy, swarm.labels, field_td=200.0)
        """
        if not np.isfinite(energy) or energy <= 0:
            raise ValueError("平均電子エネルギーは正の有限値にしてください。")
        te = 2 * energy / (3 * KB_EV)
        td = self.field_td if field_td is None else field_td
        coefficients = np.zeros(len(self.reactions))
        for j, reaction in enumerate(self.reactions):
            if reaction.kind == "field_table":
                match = re.fullmatch(r"C(\d+) (.+)", reaction.value)
                if match is None:
                    raise ValueError(
                        f"電子衝突過程の指定が不正：{reaction.value}"
                    )
                cid, label = int(match[1]), " ".join(match[2].split())
                if labels.get(cid) != label:
                    raise ValueError(f"衝突過程の対応不一致：{reaction.value}")
                coefficients[j] = rates[cid] * (
                    1e-6 if cid == 27 and not three_body_si else 1.0
                )
            else:
                coefficients[j] = analytic_rate(
                    reaction.kind, reaction.value, self.tg, te, td
                )
                if reaction.length_unit == "cm":
                    coefficients[j] *= 1e-6 ** (self.orders[j] - 1)
                elif reaction.length_unit != "m":
                    raise ValueError("未知の長さ単位です。")
        if not np.all(np.isfinite(coefficients)) or np.any(coefficients < 0):
            raise ValueError("不正な反応係数です。")
        self.k = coefficients
        self.energy, self.te = float(energy), float(te)
        self.active_field_td = float(td)
        self.scaled_k = self.k * self.density ** (self.orders - 1)

    def initial_state(self, seed_density):
        r"""指定した電子密度を持つ初期状態を作成する。

        状態ベクトルは気体密度で規格化され、`self.species` の順序で
        格納される。N2とO2は指定されたモル分率に応じて中性種と1価の
        正イオンに分配され、電子密度は `seed_density` になる。

        Parameters
        ----------
        seed_density : float
            初期電子密度 [m^-3]。0より大きく、初期気体密度より小さい
            値を指定する。

        Returns
        -------
        numpy.ndarray
            `self.species` と同じ長さの、規格化された初期密度ベクトル。

        Raises
        ------
        ValueError
            `seed_density` が有限値でない、0以下、または初期気体密度以上
            の場合。

        Examples
        --------
        >>> from plasglo import Mechanism, read_reactions
        >>> model = Mechanism(read_reactions(), None, field_td=100.0)
        >>> state = model.initial_state(1.0e13)
        >>> bool(state[model.index['e']] > 0)
        True
        """
        if (
            not np.isfinite(seed_density)
            or not 0 < seed_density < self.density
        ):
            raise ValueError(
                "初期電子密度は0より大きく気体密度より小さくしてください。"
            )
        y = np.zeros(len(self.species))
        seed = seed_density / self.density
        for name, fraction in self.initial_fractions.items():
            y[self.index[name]] = fraction * (1 - seed)
            y[self.index[name + "+"]] = fraction * seed
        y[self.index["e"]] = seed
        return y

    def _products(self, y):
        r"""反応速度の計算に使う非負の状態量を作成する。

        陰解法の試行値にわずかな負値が含まれていても、反応速度が負に
        ならないように各成分を0以上へ置き換える。元の `y` 自体は変更しない。

        Parameters
        ----------
        y : numpy.ndarray
            `self.species` の順序で並んだ規格化密度ベクトル。

        Returns
        -------
        tuple
            `(q, products, m)` のタプル。`q` は非負化した状態ベクトル、
            `products` は各反応物のべき乗項、`m` はN2とO2の合計密度。

        Examples
        --------
        >>> from plasglo import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> clipped, products, collider = model._products(
        ...     model.initial_state()
        ... )
        >>> clipped.shape, products.shape
        ((3,), (3, 2))
        """
        # 陰解法の試行値が僅かに負でも、負の反応速度を作らない。
        # 保存性を保つため、解や微分値そのものの事後クリップはしない。
        return super()._products(y)

    def reaction_rates(self, y):
        r"""各反応の規格化反応速度を計算する。

        反応係数、反応物濃度のべき乗、第三体濃度を掛け合わせて、
        各反応の反応速度を求める。`y` に負値が含まれる場合は、計算時に
        0として扱う。

        Parameters
        ----------
        y : numpy.ndarray
            `self.species` の順序で並んだ規格化密度ベクトル。

        Returns
        -------
        numpy.ndarray
            反応ごとの規格化反応速度。長さは `len(self.reactions)`。

        Examples
        --------
        >>> from plasglo import Mechanism, read_reactions
        >>> model = Mechanism(read_reactions(), None, field_td=100.0)
        >>> rates = model.reaction_rates(model.initial_state(1.0e13))
        >>> rates.shape == (len(model.reactions),)
        True
        """
        return super().reaction_rates(y)

    def rhs(self, time, y):
        r"""化学反応による状態ベクトルの時間微分を計算する。

        各反応速度に化学量論係数を掛けて合計し、常微分方程式
        ``dy/dt = S @ reaction_rates(y)`` の右辺を返す。`time` は
        `scipy.integrate.solve_ivp` との互換性のために受け取るが、
        この計算では使用しない。

        Parameters
        ----------
        time : float
            時刻。反応速度自体は時刻に直接依存しない。
        y : numpy.ndarray
            `self.species` の順序で並んだ規格化密度ベクトル。

        Returns
        -------
        numpy.ndarray
            `y` と同じ長さの状態ベクトルの時間微分。

        Examples
        --------
        >>> from plasglo import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> model.rhs(0.0, model.initial_state()).shape
        (3,)
        """
        return super().rhs(time, y)

    def jacobian(self, time, y):
        r"""状態方程式の解析的なヤコビアンを計算する。

        `rhs(time, y)` を各状態変数で偏微分した行列を返す。第三体反応
        では、N2とO2の濃度が第三体濃度として反応速度に与える微分も
        含める。`time` は `solve_ivp` との互換性のために受け取るが、
        この計算では使用しない。

        Parameters
        ----------
        time : float
            時刻。ヤコビアン自体は時刻に直接依存しない。
        y : numpy.ndarray
            `self.species` の順序で並んだ規格化密度ベクトル。

        Returns
        -------
        numpy.ndarray
            行と列がそれぞれ状態成分に対応するヤコビアン。
            形状は `(len(self.species), len(self.species))`。

        Examples
        --------
        >>> from plasglo import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> model.jacobian(0.0, model.initial_state()).shape
        (3, 3)
        """
        return super().jacobian(time, y)
