r"""明示した任意の衝突標的とチャンネルに対する定常BOLOS計算。
"""

from __future__ import annotations

from pathlib import Path
import sys
from .types import SwarmState
from ..paths import BOLOS_DIR
import numpy as np

if not (BOLOS_DIR / "bolos" / "__init__.py").is_file():
    raise ModuleNotFoundError(
        f"Py3bolosのクローンが見つかりません: {BOLOS_DIR}。"
        "このフォルダにkeipero/Py3bolosをクローンしてください。",
        name="bolos",
    )
if str(BOLOS_DIR) not in sys.path:
    sys.path.insert(0, str(BOLOS_DIR))
from bolos import grid, parser, solver
import bolos as _bolos

loaded_file = getattr(_bolos, "__file__", None)
if (
    loaded_file is None
    or Path(loaded_file).resolve().parent != BOLOS_DIR / "bolos"
):
    raise ImportError("bolosインポートエラー")


def normalize_eedf(b, distribution):
    r"""有意な負値を拒否し、丸め誤差の範囲でEEDF末尾を延長する。

    Parameters
    ----------
    b : bolos.solver.BoltzmannSolver
        エネルギー格子と規格化規約を提供するBOLOSソルバー。
    distribution : numpy.ndarray
        ソルバーのエネルギー格子上のEEDF値。

    Returns
    -------
    tuple
        規格化したEEDF、補正前の最小値、末尾補正の有無。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。

    Examples
    --------
    >>> from plasglo import load_case
    >>> from plasglo.electrons.mixture import ElectronSwarm
    >>> profile = load_case("argon").electron_settings()["profile_data"]
    >>> backend = ElectronSwarm(
    ...     profile["cross_sections"], profile["targets"],
    ...     profile["channels"], cells=100,
    ... )
    >>> _ = backend.solve(51.0, {"AR": 2.5e25}, 301.0)
    >>> f, minimum, repaired = normalize_eedf(
    ...     backend.bsolver, backend.eedf
    ... )
    >>> bool((f >= 0).all())
    True
    """
    f = distribution.copy()
    minimum = float(np.min(f))
    if not np.isfinite(f).all() or minimum < -1e-12 * np.max(f):
        raise RuntimeError("BOLOSのEEDFに非有限値または有意な負値があります。")
    # A hard floor in an unresolved tail produces a cliff in log(F).
    # Continue the last decreasing slope resolved above relative 1e-12.
    tail_repaired = False
    if minimum <= 0:
        reliable = np.flatnonzero(f > np.max(f) * 1e-12)
        if len(reliable) >= 2:
            previous, last = reliable[-2:]
            slope = (np.log(f[last]) - np.log(f[previous])) / (
                b.cenergy[last] - b.cenergy[previous]
            )
            if slope < 0 and last < len(f) - 1:
                f[last + 1:] = f[last] * np.exp(
                    np.maximum(
                        slope * (b.cenergy[last + 1:] - b.cenergy[last]), -690
                    )
                )
                tail_repaired = True
    f = b._normalized(np.maximum(f, np.max(f) * 1e-300))
    return f, minimum, tail_repaired


class ElectronSwarm:
    r"""衝突データ、化学種への対応、特殊な密度重みを入力で指定する。

    E/NのNには宣言した衝突標的の密度和を用いる。チャンネルはファイル内の位置ではなく物理的な対応で照合する。密度規格化された断面積には、
    追加の密度重みとSI化学係数への換算を明示する。
    """

    def __init__(
        self,
        paths,
        targets,
        channels,
        *,
        cells=400,
        rtol=1e-6,
        max_iterations=1000,
        coarse_max_eV=60.0,
        adaptive_factor=15.0,
    ):
        r"""衝突データを読み込んで検証し、BOLOSを初期化する。

        Parameters
        ----------
        paths : sequence of path-like
            まとめて読み込むBOLSIG形式の衝突断面積ファイル。
        targets : dict of str to str
            断面積データの標的名を化学種名へ対応付ける辞書。
        channels : sequence of dict
            チャンネルID、物理過程の照合条件、密度補正の定義。
        cells : int, optional
            BOLOSで使用するエネルギー格子のセル数。
        rtol : float, optional
            収束判定または時間積分の相対許容誤差。
        max_iterations : int, optional
            1回の収束計算で許容するBOLOSの最大反復回数。
        coarse_max_eV : float, optional
            初期の粗い格子の上限エネルギー [eV]。
        adaptive_factor : float, optional
            適応格子の上限を決める、平均エネルギーに対する倍率。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> from plasglo import load_case
        >>> from plasglo.electrons.mixture import ElectronSwarm
        >>> profile = load_case("argon").electron_settings()[
        ...     "profile_data"
        ... ]
        >>> backend = ElectronSwarm(
        ...     profile["cross_sections"], profile["targets"],
        ...     profile["channels"], cells=100,
        ... )
        >>> backend.cells
        100
        """
        if not isinstance(cells, (int, np.integer)) or cells < 50:
            raise ValueError("Use at least 50 integer BOLOS energy cells.")
        if not np.isfinite(rtol) or not 0 < rtol < 1:
            raise ValueError("BOLOS rtol must lie between zero and one.")
        if not isinstance(max_iterations, int) or max_iterations < 1:
            raise ValueError(
                "BOLOS max_iterations must be a positive integer."
            )
        if (
            not np.isfinite(coarse_max_eV)
            or coarse_max_eV <= 0
            or not np.isfinite(adaptive_factor)
            or adaptive_factor <= 0
        ):
            raise ValueError("Energy grid limits must be positive and finite.")
        if not targets or len(set(targets.values())) != len(targets):
            raise ValueError(
                "Declare a unique chemical species for each collision target."
            )
        self.targets = dict(targets)
        self.paths = [Path(p) for p in paths]
        self.path = self.paths[0]
        collisions = []
        for path in self.paths:
            with path.open(encoding="utf-8") as fp:
                collisions.extend(parser.parse(fp))
        if not collisions or {c["target"] for c in collisions} != set(targets):
            raise ValueError(
                "Declared targets must match cross-section targets exactly."
            )
        for target in targets:
            if not any(
                c["target"] == target
                and c["kind"] in ("ELASTIC", "EFFECTIVE", "MOMENTUM")
                for c in collisions
            ):
                raise ValueError(
                    "Missing momentum-transfer data for collision target"
                    f" {target}."
                )
        self.channels = list(channels)
        self.channel_ids = [c["id"] for c in channels]
        if (
            not channels
            or any(not value for value in self.channel_ids)
            or len(set(self.channel_ids)) != len(channels)
        ):
            raise ValueError("Channel IDs must be nonempty and unique.")
        indices, used = [], set()
        normalized = {}
        for channel in channels:
            if set(channel) - {
                "id",
                "target",
                "kind",
                "product",
                "threshold_eV",
                "label",
                "density_normalization",
            }:
                raise ValueError(
                    f'Unknown channel fields for {channel["id"]}.'
                )
            matches = [
                i
                for i, c in enumerate(collisions)
                if c["target"] == channel["target"]
                and c["kind"] == channel["kind"]
                and c.get("product") == channel.get("product")
                and np.isclose(
                    c.get("threshold", 0.0),
                    channel.get("threshold_eV", 0.0),
                    rtol=0,
                    atol=5e-4,
                )
            ]
            if len(matches) != 1 or matches[0] in used:
                raise ValueError(
                    f'Missing, ambiguous or duplicate channel: {channel["id"]}'
                )
            i = matches[0]
            indices.append(i)
            used.add(i)
            if channel.get("density_normalization"):
                spec = dict(channel["density_normalization"])
                if set(spec) - {
                    "species",
                    "power",
                    "density_unit",
                    "coefficient_scale",
                    "virtual_target",
                }:
                    raise ValueError("Unknown density-normalization fields.")
                if (
                    spec["density_unit"] not in ("m^-3", "cm^-3")
                    or not isinstance(spec.get("power", 1), int)
                    or spec.get("power", 1) < 1
                    or not np.isfinite(spec["coefficient_scale"])
                    or spec["coefficient_scale"] <= 0
                ):
                    raise ValueError(
                        "Invalid density normalization or SI conversion."
                    )
                spec.setdefault("power", 1)
                spec.setdefault("virtual_target", f"__density_weighted_{i}__")
                if spec["virtual_target"] in targets:
                    raise ValueError(
                        "Virtual targets must differ from real target names."
                    )
                normalized[i] = spec
        virtual_names = [
            item["virtual_target"] for item in normalized.values()
        ]
        if len(set(virtual_names)) != len(virtual_names):
            raise ValueError(
                "Density-weighted channels need distinct virtual target names."
            )
        self.bsolver = solver.BoltzmannSolver(
            grid.LinearGrid(0, coarse_max_eV, cells)
        )
        ordinary_indices = [
            i for i in range(len(collisions)) if i not in normalized
        ]
        ordinary = self.bsolver.load_collisions(
            [collisions[i] for i in ordinary_indices]
        )
        processes = dict(zip(ordinary_indices, ordinary))
        for i, spec in normalized.items():
            processes[i] = self.bsolver.add_process(
                **dict(collisions[i], target=spec["virtual_target"])
            )
        self.processes = [processes[i] for i in indices]
        self.density_channels = [
            (
                self.channel_ids[j],
                channels[j]["target"],
                processes[i],
                normalized[i],
            )
            for j, i in enumerate(indices)
            if i in normalized
        ]
        self.labels = {
            c["id"]: c.get(
                "label",
                f'{c["target"]} {c["kind"]} {c.get("product") or ""}'.strip(),
            )
            for c in channels
        }
        self.collisions = collisions
        self.cells, self.rtol, self.max_iterations = (
            cells,
            rtol,
            max_iterations,
        )
        self.coarse_max_eV, self.adaptive_factor = (
            coarse_max_eV,
            adaptive_factor,
        )
        self.eedf = self.last_field_td = None
        self.solve_count = 0

    def solve(self, field_td, densities, gas_temperature=300.0) -> SwarmState:
        r"""宣言した標的の密度に対する定常EEDFを求める。

        Parameters
        ----------
        field_td : float
            換算電界E/N [Td]。
        densities : dict of str to float
            化学種名をキーとする数密度 [m^-3]。標的と密度補正の相手を含む。
        gas_temperature : float, optional
            一定の気体温度 [K]。

        Returns
        -------
        SwarmState
            SIチャンネル係数、平均エネルギー [eV]、組成比、診断情報。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。
        RuntimeError
            ソルバーが収束しないか、計算した物理量が不正な場合。

        Examples
        --------
        >>> from plasglo import load_case
        >>> from plasglo.electrons.mixture import ElectronSwarm
        >>> profile = load_case("argon").electron_settings()[
        ...     "profile_data"
        ... ]
        >>> backend = ElectronSwarm(
        ...     profile["cross_sections"], profile["targets"],
        ...     profile["channels"], cells=100,
        ... )
        >>> state = backend.solve(51.0, {"AR": 2.5e25}, 301.0)
        >>> state["energy_eV"] > 0
        True
        """
        if (
            not np.isfinite(field_td)
            or field_td < 0
            or not np.isfinite(gas_temperature)
            or gas_temperature <= 0
        ):
            raise ValueError(
                "E/N must be nonnegative and temperature positive (finite)."
            )
        names = list(self.targets)
        values = np.array(
            [densities[self.targets[target]] for target in names]
        )
        if (
            not np.isfinite(values).all()
            or np.any(values < 0)
            or values.sum() <= 0
        ):
            raise ValueError(
                "Collision-target densities must be nonnegative with a"
                " positive total."
            )
        total = float(values.sum())
        fractions = dict(zip(names, (values / total).tolist()))
        b = self.bsolver
        for target in names:
            b.set_density(target, fractions[target])
        for _, target, process, spec in self.density_channels:
            density = densities[spec["species"]]
            if not np.isfinite(density) or density < 0:
                raise ValueError(
                    "Invalid density-normalized collision partner density."
                )
            factor = density * (
                1e-6 if spec["density_unit"] == "cm^-3" else 1.0
            )
            b.set_density(
                spec["virtual_target"],
                fractions[target] * factor ** spec["power"],
            )
        b.kT = gas_temperature * solver.KB / solver.ELECTRONVOLT
        b.EN = field_td * solver.TOWNSEND
        try:
            with np.errstate(
                divide="ignore",
                invalid="ignore",
                over="ignore",
                under="ignore",
            ):
                if field_td == 0:
                    if self.last_field_td != 0:
                        b.grid = grid.QuadraticGrid(
                            0, self.adaptive_factor * b.kT, self.cells
                        )
                        guess = b.maxwell(b.kT)
                    else:
                        guess = self.eedf
                    b.init()
                elif self.eedf is None or self.last_field_td == 0:
                    b.grid = grid.LinearGrid(0, self.coarse_max_eV, self.cells)
                    b.init()
                    coarse = b.grid
                    f = b.converge(
                        b.maxwell(2.0),
                        maxn=self.max_iterations,
                        rtol=self.rtol,
                    )
                    b.grid = grid.QuadraticGrid(
                        0, self.adaptive_factor * b.mean_energy(f), self.cells
                    )
                    b.init()
                    guess = b.grid.interpolate(f, coarse)
                else:
                    b.init()
                    guess = self.eedf
                f, iterations, error = b.converge(
                    guess,
                    maxn=self.max_iterations,
                    rtol=self.rtol,
                    full=True,
                    delta0=1e20 if field_td == 0 else 1e14,
                )
                f, minimum, repaired = normalize_eedf(b, f)
                raw = np.array([b.rate(f, p) for p in self.processes])
                energy = float(b.mean_energy(f))
        except solver.ConvergenceError as exc:
            raise RuntimeError(
                f"BOLOS convergence failed: E/N={field_td},"
                f" fractions={fractions}"
            ) from exc
        if (
            not np.isfinite(raw).all()
            or np.any(raw < -1e-12 * np.max(raw))
            or not np.isfinite(energy)
            or energy <= 0
        ):
            raise RuntimeError(
                "Invalid EEDF, reaction coefficients or mean electron energy."
            )
        rates = dict(zip(self.channel_ids, np.maximum(raw, 0.0).tolist()))
        for channel_id, _, _, spec in self.density_channels:
            rates[channel_id] *= spec["coefficient_scale"]
        self.eedf, self.last_field_td = f.copy(), float(field_td)
        self.solve_count += 1
        return SwarmState(
            rates=rates,
            energy_eV=energy,
            iterations=int(iterations),
            error=float(error),
            uncorrected_eedf_min=minimum,
            eedf_tail_repaired=repaired,
            target_density_m3=total,
            fractions={self.targets[t]: fractions[t] for t in names},
        )
