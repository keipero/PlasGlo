r"""窒素・アルゴンサンプルの反応式と質量作用則の計算。

入力構文はZDPlasKinに由来するが、このモジュールはPythonだけで計算する。
"""

from __future__ import annotations

import ast
from collections import Counter
import math
from pathlib import Path
import re
import numpy as np
from ..constants import KB, KB_EV
from .kernel import MassActionKernel


def number(text):
    r"""Fortran倍精度指数をPythonの指数表記に変換する。

    Parameters
    ----------
    text : str
        解析対象の文字列。

    Returns
    -------
    str
        d/Dによる指数をeによる指数へ変換した数値式。

    Examples
    --------
    >>> number("1.d-10")
    '1.e-10'
    """
    return re.sub(r"(?i)(\d+(?:\.\d*)?|\.\d+)[d]([+-]?\d+)", r"\1e\2", text)


class Expression:
    r"""確認済みFortran速度式のための制限付き算術評価器。

    eval/exec、属性参照、ファイル操作、任意の関数呼び出しは許可しない。
    """

    functions = {
        "exp": math.exp,
        "log": math.log,
        "sqrt": math.sqrt,
        "min": min,
        "max": max,
    }

    def __init__(self, text):
        r"""制限付き算術式を構文木へ変換する。

        Parameters
        ----------
        text : str
            解析対象の文字列。

        Examples
        --------
        >>> expr = Expression("2.d0 * Tg")
        >>> expr({"tg": 300.0})
        600.0
        """
        self.node = ast.parse(number(text).lower().strip(), mode="eval").body

    def __call__(self, values):
        r"""指定された変数の値で構文解析済みの式を評価する。

        Parameters
        ----------
        values : mapping
            変数名をキー、数値を値とする式評価用の辞書。

        Returns
        -------
        float
            制限付き算術式の計算結果。

        Examples
        --------
        >>> Expression("sqrt(tg)")({"tg": 900.0})
        30.0
        """

        def visit(node):
            r"""許可された算術構文木のノードを1つ評価する。

            Parameters
            ----------
            node : ast.AST
                制限付き式評価器で処理する算術構文木のノード。

            Returns
            -------
            float
                構文木ノードが表す数値。

            Raises
            ------
            ValueError
                入力が要求する形状、単位、モデルの制約を満たさない場合。
            """
            if isinstance(node, ast.Constant) and isinstance(
                node.value, (int, float)
            ):
                return node.value
            if isinstance(node, ast.Name):
                return values[node.id]
            if isinstance(node, ast.UnaryOp) and isinstance(
                node.op, (ast.UAdd, ast.USub)
            ):
                return visit(node.operand) * (
                    -1 if isinstance(node.op, ast.USub) else 1
                )
            if isinstance(node, ast.BinOp):
                a, b = visit(node.left), visit(node.right)
                if isinstance(node.op, ast.Add):
                    return a + b
                if isinstance(node.op, ast.Sub):
                    return a - b
                if isinstance(node.op, ast.Mult):
                    return a * b
                if isinstance(node.op, ast.Div):
                    return a / b
                if isinstance(node.op, ast.Pow):
                    return a**b
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id in self.functions
                and not node.keywords
            ):
                return self.functions[node.func.id](
                    *(visit(a) for a in node.args)
                )
            raise ValueError(f"Unsupported rate expression: {ast.dump(node)}")

        return float(visit(self.node))


def side(text):
    r"""旧形式の反応式の片側にある同一粒子種の出現数を数える。

    Parameters
    ----------
    text : str
        解析対象の文字列。

    Returns
    -------
    collections.Counter
        大文字表記の粒子名をキー、出現数を値とする辞書。

    Examples
    --------
    >>> side("N2 + N2")["N2"]
    2
    """
    return Counter(t.strip().upper() for t in text.strip().split(" + "))


class SampleMechanism(MassActionKernel):
    r"""サンプルの反応入力から係数固定の密度方程式を構築する。
    """

    def __init__(
        self,
        path,
        names,
        conditions,
        swarm_rates,
        energy,
        *,
        expected_counts=(21, 96),
    ):
        r"""反応機構と速度係数のメタデータを構築する。

        Parameters
        ----------
        path : str or pathlib.Path
            読み込む入力ファイルのパス。
        names : sequence of str
            密度配列の行と同じ順序の粒子種名。
        conditions : dict or object
            気体条件、電界設定、旧モデルの設定情報。
        swarm_rates : dict of str to float
            SI単位系の電子衝突チャンネル係数。
        energy : float or None
            平均電子エネルギー [eV]。Teに依存しない式ではNoneを許容する。
        expected_counts : tuple of int or None, optional
            期待する粒子種数と反応数。Noneの場合は件数の検査を行わない。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。
        """
        self.names, self.index = names, {
            s.upper(): i for i, s in enumerate(names)
        }
        self.n0 = conditions["pressure_pa"] / (
            KB * conditions["temperature_k"]
        )
        variables = {
            "tgas": conditions["temperature_k"],
            "en": conditions["field_td"],
            "te": 2 * energy / (3 * KB_EV),
        }
        reactions = []
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            text = line.strip()
            if text.startswith("$") and "=" in text:
                name, expression = text[1:].split("=", 1)
                variables[name.strip().lower()] = Expression(expression)(
                    variables
                )
            if "=>" in text and not text.startswith("#"):
                equation, rate = text.split("!", 1)
                reactions.append((equation.strip(), rate.strip()))
        ns, nr = len(names), len(reactions)
        self.stoich = np.zeros((ns, nr), dtype=int)
        self.powers = np.zeros((ns, nr), dtype=int)
        self.third = np.zeros(nr, dtype=int)
        self.k = np.empty(nr)
        self.orders = np.empty(nr, dtype=int)
        self.equations = []
        for j, (equation, expression) in enumerate(reactions):
            left, right = [side(x) for x in equation.split("=>")]
            if left["ANY_NEUTRAL"] != right["ANY_NEUTRAL"]:
                raise ValueError("Unbalanced ANY_NEUTRAL")
            self.third[j] = left["ANY_NEUTRAL"]
            self.orders[j] = sum(left.values())
            for name in left.keys() | right.keys():
                if name == "ANY_NEUTRAL":
                    continue
                i = self.index[name]
                self.powers[i, j] = left[name]
                self.stoich[i, j] = right[name] - left[name]
            if expression.upper().startswith("BOLSIG "):
                key = expression[7:].replace(" ", "").upper()
                self.k[j] = swarm_rates[key]
            else:
                self.k[j] = Expression(expression)(variables) * 1e-6 ** (
                    self.orders[j] - 1
                )
            self.equations.append(equation)
        if expected_counts is not None and (ns, nr) != expected_counts:
            raise ValueError(
                f"Unexpected sample mechanism: {ns} species, {nr} reactions"
            )
        self.scaled = self.k * self.n0 ** (self.orders - 1)
        normalized_names = [s.upper() for s in names]
        self.charge = np.array(
            [
                -1 if s == "E" else int(s.endswith("^+"))
                for s in normalized_names
            ]
        )
        self.atoms = np.array(
            [
                0 if s == "E" else int(re.match(r"(?:AR|N)(\d*)", s)[1] or 1)
                for s in normalized_names
            ]
        )
        self.neutral = self.charge == 0
        self.collider_weights = self.neutral.astype(float)
        self.third_body = self.third
        self.density = self.n0
        np.testing.assert_array_equal(self.atoms @ self.stoich, 0)
        np.testing.assert_array_equal(self.charge @ self.stoich, 0)
        if not np.isfinite(self.k).all() or (self.k < 0).any():
            raise ValueError("Invalid sample rate coefficients")

    def terms(self, y):
        r"""反応物密度因子と旧形式の中性第三体密度和を構築する。

        Parameters
        ----------
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Returns
        -------
        tuple
            非負の規格化密度、反応物の積因子、中性種密度和。
        """
        q = np.maximum(y, 0)
        products = q[:, None] ** self.powers
        return q, products, q[self.neutral].sum()

    def rhs(self, t, y):
        r"""規格化した粒子密度の時間微分を評価する。

        Parameters
        ----------
        t : float
            時刻 [s]。
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
        return super().rhs(t, y)

    def jac(self, t, y):
        r"""規格化された反応源項の解析的Jacobianを評価する。

        Parameters
        ----------
        t : float
            時刻 [s]。
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
        return super().jacobian(t, y)

    @property
    def scaled_k(self):
        r"""規格化した質量作用則の係数を取得する。

        Returns
        -------
        numpy.ndarray
            反応次数qに対してN0**(q-1)を掛けた係数 [s^-1]。
        """
        return self.scaled

    @scaled_k.setter
    def scaled_k(self, value):
        r"""規格化した質量作用則の係数を設定する。

        Parameters
        ----------
        value : numpy.ndarray
            反応次数qに対してk*N0**(q-1)で規格化した係数 [s^-1]。
        """
        self.scaled = value
