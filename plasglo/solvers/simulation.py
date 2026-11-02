r"""反応機構のBDF時間積分と電子群輸送計算の結合。
"""

from __future__ import annotations

from dataclasses import dataclass, field
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from pathlib import Path
import sys
from ..chemistry.kernel import MassActionKernel
from ..electrons.types import SwarmHistory
from ..io.cases import Case, load_case
from .frozen import integrate_prescribed_rates


@dataclass
class Result:
    r"""シミュレーション結果を格納するデータクラス。

    Parameters
    ----------
    mechanism : Mechanism
        使用した反応機構。
    time : numpy.ndarray
        出力時刻 [s]。
    densities : numpy.ndarray
        各化学種の数密度 [m^-3]。行は `mechanism.species`、列は
        `time` に対応する。
    diagnostics : dict
        保存則誤差、計算条件、反応数などの診断情報。
    swarm_history : list of SwarmHistory, optional
        BOLOSを使用した場合の電子群計算履歴。デフォルトは空のリスト。
    field_td : numpy.ndarray or None, optional
        各出力時刻における換算電界 [Td]。デフォルトは `None`。
    """

    mechanism: MassActionKernel
    time: np.ndarray
    densities: np.ndarray
    diagnostics: dict
    swarm_history: list[SwarmHistory] = field(default_factory=list)
    field_td: np.ndarray | None = None
    case: Case | None = None


def _integrate_coupled(
    model,
    initial,
    time,
    swarm,
    update_ns,
    rtol,
    atol_density,
    field_duration_ns=None,
):
    r"""BOLOSと化学反応の逐次結合を行う。

    BDFによる各積分区間では反応速度係数を固定し、区間の終点でBOLOSを
    更新する。これにより、解析的な化学反応ヤコビアンと反応速度係数の
    整合性を保ち、Newton反復中や右辺評価中に隠れたBOLOS更新が発生しない。

    Parameters
    ----------
    model : Mechanism
        化学反応機構。
    initial : numpy.ndarray
        規格化された初期状態ベクトル。
    time : numpy.ndarray
        出力時刻 [s]。単調増加する配列を指定する。
    swarm : BolosSwarm
        電子群分布と反応速度係数を計算するオブジェクト。
    update_ns : float
        BOLOSを更新する間隔 [ns]。
    rtol : float
        BDF積分の相対許容誤差。
    atol_density : float
        実密度に対する絶対許容誤差 [m^-3]。
    field_duration_ns : float or None, optional
        電界を印加する時間 [ns]。指定した場合、以後の換算電界を0 Td
        とする。デフォルトは `None`。

    Returns
    -------
    values : numpy.ndarray
        各化学種の規格化密度。形状は
        `(len(initial), len(time))`。
    nfev : int
        BDF右辺関数の評価回数。
    history : list of SwarmHistory
        BOLOS更新ごとの電子エネルギー、反応速度係数、対象密度などの履歴。

    Raises
    ------
    RuntimeError
        状態に非有限値または許容範囲を超える負値がある場合、BOLOSが
        収束しない場合、または積分に失敗した場合。
    ValueError
        BOLOS更新間隔が不正な場合。
    """
    values = np.empty((len(initial), len(time)))
    values[:, 0] = initial
    history: list[SwarmHistory] = []
    nfev = 0
    charge = getattr(model, "charge", model.conserved[-1])
    neutral = np.flatnonzero(charge == 0)
    target_names = list(swarm.targets.values())
    missing = [i for i in neutral if model.species[i] not in target_names]
    off_time = (
        np.inf if field_duration_ns is None else field_duration_ns * 1e-9
    )

    def update(t, y):
        r"""受理した結合区間の終端で電子反応係数を更新する。

        Parameters
        ----------
        t : float
            時刻 [s]。
        y : numpy.ndarray
            基準気体密度で規格化した各粒子種の数密度。

        Raises
        ------
        RuntimeError
            ソルバーが収束しないか、計算した物理量が不正な場合。
        """
        if (
            not np.isfinite(y).all()
            or y.min() * model.density < -100 * atol_density
        ):
            raise RuntimeError(
                "結合更新時に非有限密度または許容誤差を超える負密度があります。"
            )
        q = np.maximum(y, 0.0) * model.density
        densities = {
            name: float(q[model.index[name]]) for name in model.species
        }
        td = model.field_td if t < off_time else 0.0
        try:
            state = swarm.solve(td, densities, model.tg)
        except RuntimeError as exc:
            raise RuntimeError(f"t={t * 1e9:.9g} ns: {exc}") from exc
        model.set_swarm(
            state["rates"],
            state["energy_eV"],
            swarm.labels,
            three_body_si=True,
            field_td=td,
        )
        history.append(
            SwarmHistory(
                time_s=float(t),
                **state,
                field_td=float(td),
                electron_temperature_K=model.te,
                target_densities_m3={
                    name: densities[name] for name in target_names
                },
                unsupported_neutral_fraction=(
                    float(q[missing].sum() / q[neutral].sum())
                    if q[neutral].sum() > 0
                    else 0.0
                ),
                coefficients_SI=model.k.copy(),
            )
        )

    y, start, step = initial.copy(), 0.0, update_ns * 1e-9
    update(start, y)
    boundaries = np.r_[np.arange(step, time[-1], step), time[-1]]
    if 0 < off_time < time[-1]:
        # Align a regular update with the switch exactly, avoiding a
        # duplicate
        # interval of floating-point length when the two times almost
        # coincide.
        near_switch = np.abs(boundaries - off_time) <= 8 * np.finfo(
            float
        ).eps * max(time[-1], step)
        boundaries[near_switch] = off_time
        boundaries = np.r_[boundaries, off_time]
    for stop in np.unique(boundaries):
        if stop <= start:
            raise ValueError("BOLOSの更新間隔が小さすぎます。")
        solution = solve_ivp(
            model.rhs,
            (start, stop),
            y,
            method="BDF",
            jac=model.jacobian,
            dense_output=True,
            rtol=rtol,
            atol=atol_density / model.density,
        )
        if not solution.success:
            raise RuntimeError(
                f"t={start * 1e9:g}〜{stop * 1e9:g} ns: {solution.message}"
            )
        mask = (time > start) & (time <= stop)
        if np.any(mask):
            values[:, mask] = solution.sol(time[mask])
        nfev += solution.nfev
        y, start = solution.y[:, -1], stop
        update(start, y)
    return values, nfev, history


def simulate_case(case, *, _network=None):
    r"""粒子名に依存せず、検証したケースの密度方程式を積分する。

    EEDF係数は初期状態で固定するか、受理した各結合区間の終端で更新する。出力時刻は係数更新の間隔を決めない。

    Parameters
    ----------
    case : Case or path-like
        検証済みケース、同梱ケース名、またはケースJSONのパス。
    _network : ReactionNetwork or None, optional
        準備済みの旧API用ネットワーク。通常のケース実行では指定しない。

    Returns
    -------
    Result
        実密度 [m^-3]と時刻 [s]を含むケースの積分結果。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。

    Examples
    --------
    >>> from plasglo import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> result = simulate(case)
    >>> model = result.mechanism
    >>> result.densities.shape
    (3, 3)
    """
    if not isinstance(case, Case):
        case = load_case(case)
    model = _network if _network is not None else case.build_network()
    field = case.data.get("field", {})
    if set(field) - {"td", "duration_ns"}:
        raise ValueError("Unknown field settings.")
    duration = field.get("duration_ns")
    if duration is not None and (not np.isfinite(duration) or duration < 0):
        raise ValueError("Field duration must be nonnegative and finite.")
    numerics = case.data.get("solver", {})
    if set(numerics) - {"end_ns", "points", "rtol", "atol_density_m3"}:
        raise ValueError("Unknown solver settings.")
    end = numerics.get("end_ns", 1.0)
    points = numerics.get("points", 501)
    rtol, atol = numerics.get("rtol", 1e-6), numerics.get(
        "atol_density_m3", 1e3
    )
    if (
        not np.isfinite(end)
        or end <= 0
        or isinstance(points, bool)
        or not isinstance(points, int)
        or points < 2
    ):
        raise ValueError(
            "End time must be positive; points must be an integer >= 2."
        )
    if (
        not np.isfinite(rtol)
        or not 0 < rtol < 1
        or not np.isfinite(atol)
        or atol <= 0
    ):
        raise ValueError("Invalid BDF tolerance.")
    time = np.linspace(0.0, end * 1e-9, points)
    initial = model.initial_state()
    settings = case.electron_settings()
    backend = settings["backend"]
    mode = settings.get(
        "update_mode", "coupled" if backend == "bolos" else "frozen"
    )
    if mode not in ("frozen", "coupled") or (
        backend != "bolos" and mode != "frozen"
    ):
        raise ValueError("Only BOLOS supports coupled coefficient updates.")
    if duration is not None and mode == "frozen":
        raise ValueError(
            "A field switch requires coupled updates, not frozen coefficients."
        )
    history: list[SwarmHistory] = []
    target_names, swarm = [], None
    if backend == "bolos":
        from ..electrons.mixture import ElectronSwarm

        profile = settings["profile_data"]
        target_names = list(profile["targets"].values())
        needed = set(target_names)
        for channel in profile["channels"]:
            if channel.get("density_normalization"):
                needed.add(channel["density_normalization"]["species"])
        if needed - set(model.species):
            raise ValueError(
                "Unknown electron target/density partner species:"
                f" {needed - set(model.species)}"
            )
        if model.electron_species in target_names:
            raise ValueError(
                "Electron-electron collisions are not supported by this BOLOS"
                " adapter."
            )
        known_channels = {c["id"] for c in profile["channels"]}
        required = {r.value for r in model.reactions if r.kind == "electron"}
        if required - known_channels:
            raise ValueError(
                "Missing electron channel mappings:"
                f" {required - known_channels}"
            )
        swarm = ElectronSwarm(
            profile["cross_sections"],
            profile["targets"],
            profile["channels"],
            cells=settings.get("cells", 400),
            rtol=settings.get("rtol", 1e-6),
            max_iterations=settings.get("max_iterations", 1000),
            coarse_max_eV=profile.get("coarse_max_eV", 60.0),
            adaptive_factor=profile.get("adaptive_factor", 15.0),
        )
        if mode == "coupled":
            update = settings.get("update_ns", 1.0)
            if not np.isfinite(update) or update <= 0:
                raise ValueError(
                    "The BOLOS update interval must be positive and finite."
                )
            normalized, nfev, history = _integrate_coupled(
                model, initial, time, swarm, update, rtol, atol, duration
            )
        else:
            densities = dict(zip(model.species, initial * model.density))
            state = swarm.solve(model.field_td, densities, model.tg)
            model.set_swarm(state["rates"], state["energy_eV"])
            neutral_total = sum(
                densities[s]
                for s in model.species
                if model.charge[model.index[s]] == 0
            )
            missing = sum(
                densities[s]
                for s in model.species
                if model.charge[model.index[s]] == 0 and s not in target_names
            )
            history = [
                SwarmHistory(
                    time_s=0.0,
                    **state,
                    field_td=model.field_td,
                    electron_temperature_K=model.te,
                    target_densities_m3={
                        s: densities[s] for s in target_names
                    },
                    unsupported_neutral_fraction=(
                        missing / neutral_total if neutral_total else 0.0
                    ),
                    coefficients_SI=model.k.copy(),
                )
            ]
    else:
        rates = settings.get("coefficients_SI", {})
        if backend == "none" and any(
            r.kind == "electron" for r in model.reactions
        ):
            raise ValueError(
                "Electron reactions require a BOLOS or prescribed-coefficient"
                " backend."
            )
        model.set_swarm(rates, settings.get("mean_energy_eV"))
    if mode == "frozen":
        solution = integrate_prescribed_rates(
            model,
            time,
            model.initial_density_m3,
            rtol=rtol,
            atol_density_m3=atol,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        normalized, nfev = solution.y, solution.nfev
    density = normalized * model.density
    if not np.isfinite(density).all() or density.min() < -100 * atol:
        raise RuntimeError(
            "Nonfinite or significantly negative densities; check solver"
            " tolerance."
        )
    totals, start = model.conserved @ normalized, model.conserved @ initial
    reference = np.where(start[:-1] > 0, start[:-1], 1.0)
    atom_error = (
        float(
            np.max(np.abs(totals[:-1] - start[:-1, None]) / reference[:, None])
        )
        if len(reference)
        else 0.0
    )
    charge_error = float(
        np.max(np.abs(totals[-1] - start[-1])) * model.density
    )
    neutral_indices = np.flatnonzero(model.charge == 0)
    unsupported_indices = [
        i for i in neutral_indices if model.species[i] not in target_names
    ]
    neutral_total = density[neutral_indices].sum(axis=0)
    unsupported = np.divide(
        density[unsupported_indices].sum(axis=0),
        neutral_total,
        out=np.zeros_like(neutral_total),
        where=neutral_total > 0,
    )
    ion_fraction = float(
        np.max(density[model.charge > 0].sum(axis=0) / model.density)
    )
    background = [model.index[s] for s in model.initial_fractions]
    depletion = (
        float(
            np.max(
                np.abs(normalized[background] - initial[background, None])
                / np.where(initial[background] > 0, initial[background], 1.0)[
                    :, None
                ]
            )
        )
        if background
        else 0.0
    )
    electron = model.index.get(model.electron_species)
    field_values = (
        np.full_like(time, model.field_td)
        if duration is None
        else np.where(time < duration * 1e-9, model.field_td, 0.0)
    )
    coefficients = model.reaction_rates(normalized[:, -1])
    net, turnover = (
        model.stoich @ coefficients,
        np.abs(model.stoich) @ coefficients,
    )
    active = turnover > 0
    limitations = (
        "No transport, flow, gas heating, circuit or independent electron"
        " energy equation. Data coverage and validity ranges limit each"
        " case."
    )
    diagnostics = dict(
        case_name=case.data.get("name", case.source.stem),
        case_file=str(case.source),
        python_executable=sys.executable,
        field_td=model.field_td,
        field_duration_ns=duration,
        final_field_td=float(field_values[-1]),
        end_ns=end,
        gas_temperature_K=model.tg,
        pressure_Pa=model.pressure_pa,
        initial_gas_density_m3=model.density,
        initial_mole_fractions=model.initial_fractions,
        species_count=len(model.species),
        reaction_count=len(model.reactions),
        elements=model.elements,
        rtol=rtol,
        atol_density_m3=atol,
        min_density_m3=float(density.min()),
        max_atom_relative_error=atom_error,
        max_net_charge_error_m3=charge_error,
        max_charge_residual_over_initial_density=charge_error / model.density,
        max_ion_fraction=ion_fraction,
        max_ground_gas_relative_change=depletion,
        max_electron_density_m3=(
            float(density[electron].max()) if electron is not None else 0.0
        ),
        electron_temperature_K=(
            float(model.te) if np.isfinite(model.te) else None
        ),
        max_unsupported_neutral_fraction=(
            float(unsupported.max()) if backend == "bolos" else None
        ),
        unsupported_neutral_species=(
            [model.species[i] for i in unsupported_indices]
            if backend == "bolos"
            else []
        ),
        rhs_evaluations=nfev,
        swarm_backend=backend,
        coefficient_update_mode=mode,
        final_max_relative_source_imbalance=(
            float(np.max(np.abs(net[active]) / turnover[active]))
            if np.any(active)
            else 0.0
        ),
        reaction_source_files=sorted({r.source_file for r in model.reactions}),
        balance_species=case.data.get("diagnostics_species", []),
        assumptions=(
            f"Closed constant volume; fixed gas temperature {model.tg:g} K;"
            " imposed reduced field; explicit initial densities and collider"
            " efficiencies."
        ),
        limitations=limitations,
    )
    if swarm is not None:
        diagnostics.update(
            cross_sections=[str(p) for p in swarm.paths],
            bolos_cells=swarm.cells,
            bolos_rtol=swarm.rtol,
            bolos_solve_count=swarm.solve_count,
            bolos_update_ns=(
                settings.get("update_ns") if mode == "coupled" else None
            ),
            bolos_updates_during_integration=(
                max(len(history) - 2, 0) if mode == "coupled" else 0
            ),
            max_bolos_error=max(h["error"] for h in history),
            initial_mean_energy_eV=history[0]["energy_eV"],
            final_mean_energy_eV=model.energy,
            min_bolos_uncorrected_eedf=min(
                h["uncorrected_eedf_min"] for h in history
            ),
            bolos_eedf_tail_repairs=sum(
                h["eedf_tail_repaired"] for h in history
            ),
            reduced_field_density_species=target_names,
            reduced_field_density_definition="N = "
            + " + ".join(f"n({s})" for s in target_names),
        )
        limitations += (
            " Steady EEDF; no electron-electron or superelastic collisions. "
            "Only declared collision targets contribute to the EEDF. "
            + (
                "First-order coupling: reduce update_ns to check convergence."
                if mode == "coupled"
                else (
                    "Electron coefficients are frozen at the initial"
                    " composition; this approximation must be checked."
                )
            )
        )
        if duration is not None:
            limitations += (
                " The EEDF switches to a steady zero-field solution; transient"
                " cooling is not integrated."
            )
        if unsupported.max() > 0.01 or ion_fraction > 0.01:
            warnings.warn(
                "Collision-data coverage or weak-ionization screen exceeded;"
                " review the case assumptions."
            )
    diagnostics["limitations"] = limitations
    return Result(
        model, time, density, diagnostics, history, field_values, case
    )


def simulate(
    field_td=100.0,
    end_ns=100.0,
    seed_density=1e12,
    pressure_pa=101325.0,
    points=501,
    rtol=1e-6,
    atol_density=1e3,
    mechanism=None,
    *,
    swarm_backend=None,
    bolos_update_ns=1.0,
    bolos_cells=400,
    bolos_rtol=1e-6,
    cross_sections=None,
    field_duration_ns=None,
    gas_temperature_k=300.0,
    n2_fraction=0.8,
    case=None,
):
    r"""明示した任意の気体組成をsimulate(load_case(path))で計算する。

    数値引数を使う旧空気APIも維持する。ケースの条件変更にはcase.with_updatesを使用し、2種類の指定方法を混在させた場合は
    エラーとする。どちらの入口も同じBDF/LFA計算を呼び出す。

    Parameters
    ----------
    field_td : float, optional
        換算電界E/N [Td]。
    end_ns : float, optional
        計算終了時刻 [ns]。
    seed_density : float, optional
        空気用アダプターの初期電子密度 [m^-3]。
    pressure_pa : float, optional
        初期状態を規格化するための基準圧力 [Pa]。
    points : int, optional
        開始・終了時刻を含む出力点数。
    rtol : float, optional
        収束判定または時間積分の相対許容誤差。
    atol_density : float, optional
        実密度に対する積分の絶対許容誤差 [m^-3]。
    mechanism : MassActionKernel, optional
        反応源項と解析的Jacobianを提供するネットワーク。
    swarm_backend : str or None, optional
        旧空気APIの電子係数バックエンド。他の気体はケースで指定する。
    bolos_update_ns : float, optional
        電子計算と化学反応計算の逐次結合間隔 [ns]。
    bolos_cells : int, optional
        BOLOSのエネルギー格子セル数。
    bolos_rtol : float, optional
        BOLOSの相対収束許容誤差。
    cross_sections : path-like or None, optional
        旧空気APIで既定の断面積を置き換えるファイル。
    field_duration_ns : float or None, optional
        電界の印加時間 [ns]。Noneの場合は一定電界を印加し続ける。
    gas_temperature_k : float, optional
        一定の気体温度 [K]。
    n2_fraction : float, optional
        空気用アダプターの初期N2モル分率。残りをO2とする。
    case : Case or path-like, optional
        検証済みケース、同梱ケース名、またはケースJSONのパス。

    Returns
    -------
    Result
        共通形式の密度履歴、反応係数履歴、診断結果。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> from plasglo import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> result = simulate(case)
    >>> model = result.mechanism
    >>> result.densities.shape
    (3, 3)
    """
    options = locals().copy()
    options.pop("case")
    if isinstance(field_td, (Case, str, Path)):
        if case is not None:
            raise ValueError("Specify the case only once.")
        case, options["field_td"] = field_td, 100.0
    if case is not None:
        import inspect

        defaults = inspect.signature(simulate).parameters
        if any(
            value != defaults[name].default for name, value in options.items()
        ):
            raise ValueError(
                "Use case.with_updates(...) instead of numeric air options"
                " with a case."
            )
        return simulate_case(case)

    # These names belong only to the legacy air preset, never the solver
    # core.
    from ..constants import KB
    from ..electrons.swarm_table import read_swarm
    from ..chemistry.network import ReactionNetwork, NetworkReaction, Species

    if mechanism is None:
        if (
            not np.isfinite(pressure_pa)
            or pressure_pa <= 0
            or not np.isfinite(gas_temperature_k)
            or gas_temperature_k <= 0
            or not np.isfinite(n2_fraction)
            or not 0 <= n2_fraction <= 1
        ):
            raise ValueError(
                "Invalid air pressure, temperature or N2 fraction."
            )
        density = pressure_pa / (KB * gas_temperature_k)
        if not np.isfinite(seed_density) or not 0 < seed_density < density:
            raise ValueError(
                "The seed density must lie between zero and the initial gas"
                " density."
            )
        fractions = {"N2": n2_fraction, "O2": 1 - n2_fraction}
        seed = seed_density / density
        initial_m3 = {"e": seed * density}
        for name, fraction in fractions.items():
            initial_m3[name] = fraction * (1 - seed) * density
            initial_m3[name + "+"] = fraction * seed * density
    else:
        model = mechanism
        initial = model.initial_state(seed_density)
        initial_m3 = dict(
            zip(model.species, (initial * model.density).tolist())
        )
        gas_temperature_k, pressure_pa, field_td = (
            model.tg,
            model.pressure_pa,
            model.field_td,
        )
        fractions = model.initial_fractions
    backend = swarm_backend or (
        "prescribed" if mechanism is not None else "bolos"
    )
    if backend not in ("bolos", "table", "prescribed"):
        raise ValueError("swarm_backend must be bolos or table.")
    electrons = dict(
        backend="bolos",
        update_mode="coupled",
        update_ns=bolos_update_ns,
        cells=bolos_cells,
        rtol=bolos_rtol,
    )
    if cross_sections is not None:
        electrons["cross_sections"] = str(Path(cross_sections).resolve())
    network = None
    if mechanism is not None:
        if backend == "bolos":
            raise ValueError(
                "Use an explicit case for a custom mechanism with BOLOS"
                " coupling."
            )
        species = [
            Species(
                name,
                {
                    e: int(count)
                    for e, count in zip(("N", "O"), model.conserved[:-1, i])
                    if count
                },
                int(model.conserved[-1, i]),
                name == "e",
            )
            for i, name in enumerate(model.species)
        ]
        reactions = [
            NetworkReaction(
                f"R{i + 1:03}",
                r.equation,
                "expression",
                repr(float(model.k[i])),
                source_file=r.source_file,
                line=r.line,
            )
            for i, r in enumerate(model.reactions)
        ]
        network = ReactionNetwork(
            species,
            reactions,
            temperature_K=model.tg,
            pressure_Pa=model.pressure_pa,
            field_td=model.field_td,
            third_bodies={"M": {"N2": 1.0, "O2": 1.0}},
            initial_densities=initial_m3,
            initial_fractions=model.initial_fractions,
        )
        electrons = dict(
            backend="none",
            update_mode="frozen",
            mean_energy_eV=(
                float(model.energy) if np.isfinite(model.energy) else None
            ),
        )
    elif backend == "table":
        if gas_temperature_k != 300.0 or n2_fraction != 0.8:
            raise ValueError(
                "The legacy swarm table is restricted to N2:O2=80:20 at 300 K."
            )
        rates, energy = read_swarm().interpolate(field_td)
        coefficients = {
            f"C{i}": float(value) for i, value in enumerate(rates) if i
        }
        coefficients["C27"] *= 1e-6
        electrons = dict(
            backend="prescribed",
            update_mode="frozen",
            coefficients_SI=coefficients,
            mean_energy_eV=float(energy),
        )
    case = load_case("air").with_updates(
        gas=dict(
            temperature_K=gas_temperature_k,
            pressure_Pa=pressure_pa,
            mole_fractions=fractions,
        ),
        initial_densities_m3=initial_m3,
        field=dict(td=field_td, duration_ns=field_duration_ns),
        solver=dict(
            end_ns=end_ns,
            points=points,
            rtol=rtol,
            atol_density_m3=atol_density,
        ),
        electrons=electrons,
    )
    result = simulate_case(case, _network=network)
    if network is not None:
        result.case = (
            None  # An in-memory legacy mechanism is not the air case JSON.
        )
        result.diagnostics.update(case_name="custom_mechanism", case_file=None)
    return result
