r"""ガス種に依存しない規格化質量作用則と解析的Jacobian。
"""

from __future__ import annotations

import numpy as np


class MassActionKernel:
    r"""粒子種情報と第三体の重みを入力とする共通数値カーネル。

    collider_weightsの形状は (粒子種数,) または (反応数, 粒子種数)。重みは第三体効率であり、追加の追跡粒子種ではな
    い。reaction_ratesはR/N0、physical_reaction_ratesはSI単位系のRを返す。
    """

    def _products(self, y):
        r"""非負の反応物密度因子と第三体の有効密度を構築する。

        Parameters
        ----------
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        tuple
            非負の規格化密度、反応物の積因子、第三体密度和。

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
        q = np.maximum(y, 0.0)
        products = q[:, None] ** self.powers
        weights = self.collider_weights
        if weights.ndim == 1:
            selected = weights != 0
            m = np.sum(q[selected] * weights[selected])
        else:
            if not hasattr(self, "_collider_groups"):
                unique, inverse = np.unique(
                    weights, axis=0, return_inverse=True
                )
                self._collider_groups = (
                    [(w != 0, w[w != 0]) for w in unique],
                    inverse,
                )
            groups, inverse = self._collider_groups
            m = np.array([np.sum(q[selected] * w) for selected, w in groups])[
                inverse
            ]
        return q, products, m

    def reaction_rates(self, y):
        r"""基準気体密度で割った反応速度を計算する。

        Parameters
        ----------
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        numpy.ndarray
            基準密度で割った反応源項 [s^-1]。

        Examples
        --------
        >>> from plasglo import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> model.reaction_rates(model.initial_state()).shape
        (2,)
        """
        _, products, m = self._products(y)
        return self.scaled_k * np.prod(products, axis=0) * m**self.third_body

    def physical_reaction_rates(self, y, coefficients=None):
        r"""SI単位系の実反応速度を計算する。

        Parameters
        ----------
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。
        coefficients : array_like or None, optional
            SI反応係数。Noneの場合はモデルの現在の係数を使用する。

        Returns
        -------
        numpy.ndarray
            反応順に並べた実反応速度 [m^-3 s^-1]。

        Examples
        --------
        >>> from plasglo import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> model.physical_reaction_rates(model.initial_state()).shape
        (2,)
        """
        _, products, m = self._products(y)
        k = self.k if coefficients is None else np.asarray(coefficients)
        return (
            k
            * self.density**self.orders
            * np.prod(products, axis=0)
            * m**self.third_body
        )

    def rhs(self, time, y):
        r"""規格化した粒子密度の時間微分を評価する。

        Parameters
        ----------
        time : float
            時刻 [s]。自律系の反応速度ではこの引数を使用しない。
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        numpy.ndarray
            規格化した粒子密度の時間微分 [s^-1]。

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
        return self.stoich @ self.reaction_rates(y)

    def jacobian(self, time, y):
        r"""規格化された反応源項の解析的Jacobianを評価する。

        Parameters
        ----------
        time : float
            時刻 [s]。自律系の反応速度ではこの引数を使用しない。
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        numpy.ndarray
            形状が (粒子種数, 粒子種数) の解析的Jacobian [s^-1]。

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
        q, products, m = self._products(y)
        dr = np.zeros((len(self.k), len(y)))
        for i in range(len(y)):
            if y[i] < 0:
                continue
            active = self.powers[i] > 0
            if np.any(active):
                p = products[:, active].copy()
                p[i] = self.powers[i, active] * q[i] ** (
                    self.powers[i, active] - 1
                )
                collider = m if np.ndim(m) == 0 else m[active]
                dr[active, i] = (
                    self.scaled_k[active]
                    * np.prod(p, axis=0)
                    * collider ** self.third_body[active]
                )
            weights = (
                np.full(len(self.k), self.collider_weights[i])
                if self.collider_weights.ndim == 1
                else self.collider_weights[:, i]
            )
            active = (self.third_body > 0) & (weights != 0)
            if np.any(active):
                collider = m if np.ndim(m) == 0 else m[active]
                dr[active, i] += (
                    self.scaled_k[active]
                    * np.prod(products[:, active], axis=0)
                    * self.third_body[active]
                    * weights[active]
                    * collider ** (self.third_body[active] - 1)
                )
        return self.stoich @ dr

    def jac(self, time, y):
        r"""規格化された反応源項の解析的Jacobianを評価する。

        Parameters
        ----------
        time : float
            時刻 [s]。自律系の反応速度ではこの引数を使用しない。
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        numpy.ndarray
            形状が (粒子種数, 粒子種数) の解析的Jacobian [s^-1]。

        Examples
        --------
        >>> from plasglo import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> model.jac(0.0, model.initial_state()).shape
        (3, 3)
        """
        return self.jacobian(time, y)
