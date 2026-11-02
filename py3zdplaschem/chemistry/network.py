r"""検証済みの粒子種・反応データと共通数値カーネルの接続。
"""

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
import re
import numpy as np

from .kernel import MassActionKernel
from .rates import analytic_rate
from ..constants import KB, KB_EV


@dataclass(frozen=True)
class Species:
    name: str
    composition: dict[str, int]
    charge: int = 0
    electron: bool = False
    excited: bool = False

    @classmethod
    def from_dict(cls, data):
        r"""明示した粒子種のメタデータを検証して構築する。

        Parameters
        ----------
        data : dict
            粒子名、元素組成、電荷、電子・励起種の区別。

        Returns
        -------
        Species
            元素組成、電荷、種の区別を検証した粒子種定義。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> item = Species.from_dict(
        ...     {"name": "He+", "composition": {"He": 1}, "charge": 1}
        ... )
        >>> item.charge
        1
        """
        allowed = {"name", "composition", "charge", "electron", "excited"}
        if set(data) - allowed:
            raise ValueError(f"Unknown species fields: {set(data) - allowed}")
        item = cls(**data)
        if (
            not item.name
            or not isinstance(item.name, str)
            or " + " in item.name
        ):
            raise ValueError("Species needs a nonempty, unambiguous name.")
        if not isinstance(item.composition, dict):
            raise ValueError(
                f"{item.name}: composition must map elements to atom counts."
            )
        for element, count in item.composition.items():
            if (
                not isinstance(element, str)
                or not element
                or isinstance(count, bool)
                or not isinstance(count, int)
                or count <= 0
            ):
                raise ValueError(
                    f"{item.name}: invalid element/count {element}: {count}"
                )
        if isinstance(item.charge, bool) or not isinstance(item.charge, int):
            raise ValueError(f"{item.name}: charge must be an integer.")
        if not isinstance(item.electron, bool) or not isinstance(
            item.excited, bool
        ):
            raise ValueError("Electron/excited flags must be boolean.")
        if item.electron and (item.composition or item.charge != -1):
            raise ValueError("An electron has no atoms and charge -1.")
        if not item.electron and not item.composition:
            raise ValueError(
                f"{item.name}: declare its elemental composition."
            )
        return item


@dataclass(frozen=True)
class NetworkReaction:
    id: str
    equation: str
    kind: str
    value: str
    parameters: str = ""
    length_unit: str = "m"
    source_file: str = ""
    line: int = 0
    reference: str = ""
    temperature_min_K: float | None = None
    temperature_max_K: float | None = None


class ReactionNetwork(MassActionKernel):
    r"""宣言した任意の元素を扱う共通の密度方程式。

    電子衝突係数はチャンネルIDで指定する。解析的係数はrates.pyで評価し、任意式は旧入力でも使う制限付き算術パーサーで評価する。読み
    込み時にすべての係数をSIへ変換する。
    """

    def __init__(
        self,
        species,
        reactions,
        *,
        temperature_K,
        pressure_Pa,
        field_td=0.0,
        third_bodies=None,
        variables=(),
        initial_densities=None,
        initial_fractions=None,
    ):
        r"""反応機構と速度係数のメタデータを構築する。

        Parameters
        ----------
        species : sequence of Species
            状態ベクトル順に並べた粒子種の明示的メタデータ。
        reactions : sequence of Reaction
            反応式と速度係数の指定を持つ反応定義。
        temperature_K : float
            一定の気体温度 [K]。
        pressure_Pa : float
            初期状態を規格化するための基準圧力 [Pa]。
        field_td : float, optional
            換算電界E/N [Td]。
        third_bodies : dict, optional
            第三体グループ名をキー、粒子種ごとの効率を値とする辞書。
        variables : sequence of pairs, optional
            依存関係の順に並べた、名前付き反応速度式。
        initial_densities : dict, optional
            各粒子種の初期数密度 [m^-3]。
        initial_fractions : dict, optional
            診断で参照する初期背景気体のモル分率。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> from py3zdplaschem import load_case
        >>> model = load_case("argon").build_network()
        >>> model.elements
        ['Ar']
        """
        species = [Species.from_dict(vars(item)) for item in species]
        self.species_data = species
        self.species = [s.name for s in species]
        if not self.species or len(set(self.species)) != len(self.species):
            raise ValueError("Species names must be nonempty and unique.")
        if sum(s.electron for s in species) > 1:
            raise ValueError("Declare at most one electron species.")
        self.index = {s: i for i, s in enumerate(self.species)}
        self.names = self.species  # compatibility with sample readers
        self.elements = sorted({e for s in species for e in s.composition})
        self.charge = np.array([s.charge for s in species])
        self.neutral = self.charge == 0
        self.electron_species = next(
            (s.name for s in species if s.electron), None
        )
        self.conserved = np.array(
            [[s.composition.get(e, 0) for s in species] for e in self.elements]
            + [self.charge.tolist()],
            dtype=int,
        )
        self.atoms = self.conserved[:-1].sum(axis=0)
        if (
            not np.isfinite(temperature_K)
            or temperature_K <= 0
            or not np.isfinite(pressure_Pa)
            or pressure_Pa <= 0
            or not np.isfinite(field_td)
            or field_td < 0
        ):
            raise ValueError(
                "Temperature/pressure must be positive and E/N nonnegative."
            )
        self.tg, self.pressure_pa = float(temperature_K), float(pressure_Pa)
        self.field_td = self.active_field_td = float(field_td)
        self.density = self.n0 = pressure_Pa / (KB * temperature_K)
        self.initial_fractions = dict(initial_fractions or {})
        self._initial_densities = dict(initial_densities or {})
        self.reactions = list(reactions)
        if not reactions or len({r.id for r in reactions}) != len(reactions):
            raise ValueError("Reaction IDs must be nonempty and unique.")
        ns, nr = len(species), len(reactions)
        self.stoich = np.zeros((ns, nr), dtype=int)
        self.powers = np.zeros((ns, nr), dtype=int)
        self.third_body = np.zeros(nr, dtype=int)
        self.orders = np.zeros(nr, dtype=int)
        self.collider_weights = np.zeros((nr, ns))
        groups = third_bodies or {}
        if set(groups) & set(self.species):
            raise ValueError(
                "Third-body group names must differ from species names."
            )
        for name, weights in groups.items():
            if not weights or set(weights) - set(self.species):
                raise ValueError(f"{name}: unknown or empty collider group.")
            if any(not np.isfinite(w) or w < 0 for w in weights.values()):
                raise ValueError(
                    f"{name}: collider efficiencies must be finite and"
                    " nonnegative."
                )
        self._expressions = {}
        # Import lazily to keep the expression adapter independent of
        # network I/O.
        from .sample_mechanism import Expression

        self.variables = [
            (name.lower(), Expression(text)) for name, text in variables
        ]
        for j, reaction in enumerate(reactions):
            if not reaction.id or reaction.kind not in (
                "electron",
                "analytic",
                "expression",
            ):
                raise ValueError(
                    f"{reaction.id}: unsupported rate type {reaction.kind}."
                )
            if reaction.length_unit not in ("m", "cm"):
                raise ValueError(
                    f"{reaction.id}: length_unit must be m or cm."
                )
            if reaction.kind == "electron" and reaction.length_unit != "m":
                raise ValueError(
                    "Electron channel coefficients must already be SI."
                )
            if reaction.kind == "expression":
                self._expressions[j] = Expression(reaction.value)
            if reaction.equation.count(" -> ") != 1:
                raise ValueError(
                    f"{reaction.id}: use a single spaced -> arrow."
                )

            def parse_side(text):
                r"""宣言済み粒子種と整数の化学量論係数を解析する。

                Parameters
                ----------
                text : str
                    解析対象の文字列。

                Returns
                -------
                collections.Counter
                    粒子名をキー、正の化学量論係数を値とする辞書。

                Raises
                ------
                ValueError
                    入力が要求する形状、単位、モデルの制約を満たさない場合。

                Examples
                --------
                >>> parse_side("2N2 + O2")["N2"]
                2
                """
                terms = Counter()
                known = set(self.species) | set(groups)
                for term in text.split(" + "):
                    term = term.strip()
                    if term in known:
                        terms[term] += 1
                    else:
                        match = re.fullmatch(r"(\d+)\s*(.+)", term)
                        if (
                            match is None
                            or match[2] not in known
                            or int(match[1]) < 1
                        ):
                            raise ValueError(
                                f"{reaction.id}: unknown or invalid species"
                                f" term {term!r}."
                            )
                        terms[match[2]] += int(match[1])
                return terms

            left, right = [
                parse_side(side) for side in reaction.equation.split(" -> ")
            ]
            unknown = (
                (set(left) | set(right)) - set(self.species) - set(groups)
            )
            if unknown:
                raise ValueError(f"{reaction.id}: unknown species {unknown}.")
            self.orders[j] = sum(left.values())
            used = [g for g in groups if left[g] or right[g]]
            if len(used) > 1:
                raise ValueError(
                    "Use at most one collider group per reaction."
                )
            for group in used:
                if left[group] != right[group]:
                    raise ValueError(
                        f"{reaction.id}: collider group is not conserved."
                    )
                self.third_body[j] = left[group]
                for name, efficiency in groups[group].items():
                    self.collider_weights[j, self.index[name]] = efficiency
            for name in self.species:
                i = self.index[name]
                self.powers[i, j] = left[name]
                self.stoich[i, j] = right[name] - left[name]
            if self.orders[j] < 1 or any(
                value <= 0 for value in left.values()
            ):
                raise ValueError(
                    f"{reaction.id}: positive reactant stoichiometry is"
                    " required."
                )
            lo, hi = reaction.temperature_min_K, reaction.temperature_max_K
            if (
                (
                    lo is not None
                    and (not np.isfinite(lo) or lo <= 0 or temperature_K < lo)
                )
                or (
                    hi is not None
                    and (not np.isfinite(hi) or hi <= 0 or temperature_K > hi)
                )
                or (lo is not None and hi is not None and lo > hi)
            ):
                raise ValueError(
                    f"{reaction.id}: gas temperature is outside the declared"
                    " validity range."
                )
        violations = self.conserved @ self.stoich
        if np.any(violations):
            bad = [
                reactions[j].id
                for j in np.flatnonzero(np.any(violations != 0, axis=0))
            ]
            raise ValueError(f"Element/charge conservation violated: {bad}")
        # Identical collider definitions can use the scalar path (legacy
        # arithmetic).
        if np.all(self.collider_weights == self.collider_weights[0]):
            self.collider_weights = self.collider_weights[0].copy()
        self.k = self.scaled_k = np.zeros(nr)
        self.energy = self.te = np.nan
        self.equations = [r.equation for r in reactions]

    @property
    def scaled(self):
        r"""規格化した質量作用則の係数を取得する。

        Returns
        -------
        numpy.ndarray
            反応次数qに対してN0**(q-1)を掛けた係数 [s^-1]。

        Examples
        --------
        >>> from py3zdplaschem import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> model.scaled.shape
        (2,)
        """
        return self.scaled_k

    @scaled.setter
    def scaled(self, value):
        r"""規格化した質量作用則の係数を設定する。

        Parameters
        ----------
        value : numpy.ndarray
            反応次数qに対してk*N0**(q-1)で規格化した係数 [s^-1]。
        """
        self.scaled_k = value

    def set_swarm(
        self, rates, energy, labels=None, three_body_si=True, field_td=None
    ):
        r"""電子衝突係数と解析的反応係数を整合する形で更新する。

        Parameters
        ----------
        rates : dict or numpy.ndarray
            SI単位系の電子衝突チャンネル係数。
        energy : float or None
            平均電子エネルギー [eV]。Teに依存しない式ではNoneを許容する。
        labels : mapping, optional
            旧形式の衝突過程との対応確認に使うラベル。
        three_body_si : bool, optional
            与えた三体反応係数が既にSI単位系かどうか。
        field_td : float, optional
            換算電界E/N [Td]。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> from py3zdplaschem import load_case, simulate
        >>> case = load_case("argon").with_updates(
        ...     solver={"end_ns": 0.01, "points": 3},
        ...     electrons={"cells": 100},
        ... )
        >>> result = simulate(case)
        >>> model = result.mechanism
        >>> state = result.swarm_history[0]
        >>> model.set_swarm(state["rates"], state["energy_eV"])
        >>> bool((model.k >= 0).all())
        True
        """
        if energy is not None and (not np.isfinite(energy) or energy <= 0):
            raise ValueError(
                "Mean electron energy must be positive and finite."
            )
        td = self.field_td if field_td is None else field_td
        self.energy = float(energy) if energy is not None else np.nan
        self.te = 2 * self.energy / (3 * KB_EV)
        variables = {
            "tgas": self.tg,
            "tg": self.tg,
            "te": self.te,
            "en": td,
            "td": td,
        }
        for name, expression in self.variables:
            if name in variables:
                raise ValueError(
                    f"Rate variable shadows a built-in variable: {name}"
                )
            variables[name] = expression(variables)
        coefficients = []
        for j, reaction in enumerate(self.reactions):
            if reaction.kind == "electron":
                try:
                    k = rates[reaction.value]
                except KeyError as exc:
                    raise ValueError(
                        f"{reaction.id}: missing electron channel"
                        f" {reaction.value}"
                    ) from exc
            elif reaction.kind == "analytic":
                k = analytic_rate(
                    reaction.value, reaction.parameters, self.tg, self.te, td
                )
            else:
                k = self._expressions[j](variables)
            if reaction.length_unit == "cm":
                k *= 1e-6 ** (self.orders[j] - 1)
            coefficients.append(k)
        self.k = np.asarray(coefficients)
        if not np.isfinite(self.k).all() or np.any(self.k < 0):
            raise ValueError(
                "Reaction coefficients must be finite and nonnegative."
            )
        self.scaled_k = self.k * self.density ** (self.orders - 1)
        if not np.isfinite(self.scaled_k).all():
            raise ValueError(
                "Normalized coefficients overflow; check units and reaction"
                " order."
            )
        self.active_field_td = float(td)

    def initial_state(self, seed_density=None):
        r"""明示的な初期密度から規格化された初期状態を構築する。

        Parameters
        ----------
        seed_density : None, optional
            旧APIとの対応用引数。必ずNoneとし、初期密度はケースで指定する。

        Returns
        -------
        numpy.ndarray
            基準気体密度で割った各粒子種の初期密度。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> from py3zdplaschem import load_case
        >>> model = load_case("argon").build_network()
        >>> model.initial_state().shape
        (3,)
        """
        if seed_density is not None:
            raise ValueError(
                "Set initial_densities_m3 explicitly in the case."
            )
        unknown = set(self._initial_densities) - set(self.species)
        if unknown:
            raise ValueError(f"Unknown initial species: {unknown}")
        initial = np.array(
            [self._initial_densities.get(name, 0.0) for name in self.species]
        )
        if (
            not np.isfinite(initial).all()
            or np.any(initial < 0)
            or initial.sum() <= 0
        ):
            raise ValueError(
                "Initial densities must be finite, nonnegative and nonempty."
            )
        self.initial_density_m3 = initial
        return initial / self.density
