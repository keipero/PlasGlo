r"""化学反応計算後に評価する、一様電界の0次元ドリフト・変位電流。

移動度は正の大きさで与え、速度をv_s=sign(q_s)*mu_s*Eとする。シース、空間拡散、外部回路、電極捕集は解かない。密度が正のすべて
のイオンに移動度がある場合だけ合計電流を返す。電流密度Jには形状情報は不要だが、電流Iには電極面積が必要。
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path

import numpy as np
from scipy.constants import elementary_charge, epsilon_0


@dataclass
class CurrentResult:
    traces: dict[str, np.ndarray]
    metadata: dict

    @property
    def total_available(self):
        r"""密度が正のすべてのイオンに移動度データがあるかを確認する。

        Returns
        -------
        bool
            完全な合計電流を明示的に計算できる場合はTrue。
        """
        return self.metadata["all_populated_ions_have_mobility"]


def _series(value, length, name, *, positive=False, nonnegative=False):
    r"""定数または時系列の物理量を配列化し、入力条件を検証する。

    Parameters
    ----------
    value : array_like
        検証・配列化・設定する値。
    length : int
        必要な時間サンプル数。
    name : str
        この処理で使用する粒子種・ケース・出力などの識別名。
    positive : bool, optional
        すべての値に正の制約を課すかどうか。
    nonnegative : bool, optional
        すべての値に非負の制約を課すかどうか。

    Returns
    -------
    numpy.ndarray
        指定したサンプル数を持つ、有限値の1次元配列。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> _series(2.0, 3, "example").tolist()
    [2.0, 2.0, 2.0]
    """
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        array = np.full(length, float(array))
    if array.shape != (length,) or not np.isfinite(array).all():
        raise ValueError(
            f"{name} must be a finite scalar or an array with {length} values."
        )
    if positive and np.any(array <= 0):
        raise ValueError(f"{name} must be positive.")
    if nonnegative and np.any(array < 0):
        raise ValueError(f"{name} must be nonnegative.")
    return array


def reduced_field_to_electric_field(field_td, reference_density_m3):
    r"""明示した基準密度NからE=(E/N)*1e-21*N [V/m]を求める。

    Parameters
    ----------
    field_td : float
        換算電界E/N [Td]。
    reference_density_m3 : array_like
        E/Nとmu_e*Nの両方に使用する正の基準密度 [m^-3]。

    Returns
    -------
    numpy.ndarray
        E=(E/N)*1e-21*Nで求めた電界 [V/m]。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> round(float(reduced_field_to_electric_field(100.0, 1e25)))
    1000000
    """
    density = np.asarray(reference_density_m3, dtype=float)
    field = np.asarray(field_td, dtype=float)
    if (
        not np.isfinite(density).all()
        or np.any(density <= 0)
        or not np.isfinite(field).all()
    ):
        raise ValueError(
            "Reference density must be positive and field must be finite."
        )
    return field * 1e-21 * density


def calculate_current(
    time_s,
    densities_m3,
    species,
    charges,
    electric_field_V_m,
    electron_reduced_mobility_1_V_m_s,
    reference_density_m3,
    *,
    ion_mobilities_m2_V_s=None,
    area_m2=None,
    electrode_gap_m=None,
    electric_field_derivative_V_m_s=None,
    negative_density_tolerance_m3=0.0,
):
    r"""符号付き電界方向に沿う慣用電流を計算する。

    電子の入力はBOLOS.mobility()が返すmu_e*Nとする。イオンは種ごとの実移動度 [m^2/(V s)]を定数または時系列
    で与える。密度が正のイオンで移動度が不足する場合、合計電流を表示せず、既知成分を部分値として残す。電荷中性でも電子と正イオンの慣用電流は
    相殺されない。

    E/Nと電子移動度mu_e=(mu_e*N)/Nには同じNを使用する。変位電流はepsilon_0*dE/dtとする。既知の滑らかな波形
    では時間微分を明示し、それ以外は与えた時刻格子で微分する。理想的な電圧不連続には有限の分解済みピークを定義できない。

    Parameters
    ----------
    time_s : array_like
        厳密に単調増加するサンプル時刻 [s]。
    densities_m3 : array_like
        形状が (粒子種数, 時刻数) の密度 [m^-3]。
    species : sequence of str
        密度配列の行順に並べた、電子E/eを含む粒子名。
    charges : array_like
        粒子種ごとの整数の電荷数。クーロン単位の電荷ではない。
    electric_field_V_m : array_like
        符号付きの一様電界 [V/m]。定数または時系列で指定する。
    electron_reduced_mobility_1_V_m_s : array_like
        正の電子換算移動度mu_e*N [1/(V m s)]。
    reference_density_m3 : array_like
        E/Nとmu_e*Nの両方に使用する正の基準密度 [m^-3]。
    ion_mobilities_m2_V_s : dict or None, optional
        イオン名をキーとする正の移動度の大きさ [m^2/(V s)]。
    area_m2 : float or None, optional
        電極面積 [m^2]。電流密度を電流 [A]に変換するために必要。
    electrode_gap_m : float or None, optional
        印加電圧の推定に使用する電極間距離 [m]。
    electric_field_derivative_V_m_s : array_like or None, optional
        明示する電界時間微分dE/dt [V/(m s)]。Noneでは数値微分を使用する。
    negative_density_tolerance_m3 : float, optional
        診断時にゼロ補正を許容する負密度の範囲 [m^-3]。

    Returns
    -------
    CurrentResult
        電流履歴と、移動度が不足しているイオンなどのメタデータ。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。
    """
    time = np.asarray(time_s, dtype=float)
    if (
        time.ndim != 1
        or len(time) < 2
        or not np.isfinite(time).all()
        or np.any(np.diff(time) <= 0)
    ):
        raise ValueError(
            "time_s must contain at least two finite strictly increasing"
            " times."
        )
    names = [str(name) for name in species]
    if len(set(names)) != len(names):
        raise ValueError("Species names must be unique.")
    charge = np.asarray(charges, dtype=float)
    if (
        charge.shape != (len(names),)
        or not np.isfinite(charge).all()
        or np.any(charge != np.round(charge))
    ):
        raise ValueError(
            "charges must contain one integer charge number per species."
        )
    electrons = [i for i, name in enumerate(names) if name.upper() == "E"]
    if len(electrons) != 1 or charge[electrons[0]] != -1:
        raise ValueError(
            "Exactly one electron species E/e with charge -1 is required."
        )
    density = np.asarray(densities_m3, dtype=float)
    if (
        density.shape != (len(names), len(time))
        or not np.isfinite(density).all()
    ):
        raise ValueError(
            "densities_m3 must have shape (species, time) and finite values."
        )
    if (
        not np.isfinite(negative_density_tolerance_m3)
        or negative_density_tolerance_m3 < 0
    ):
        raise ValueError(
            "negative_density_tolerance_m3 must be finite and nonnegative."
        )
    if np.any(density < -negative_density_tolerance_m3):
        raise ValueError(
            "Negative densities exceed the allowed numerical tolerance."
        )
    clipped = int(np.count_nonzero(density < 0))
    density = np.maximum(
        density, 0.0
    )  # Only this diagnostic; input data is never modified.
    count = len(time)
    field = _series(electric_field_V_m, count, "electric_field_V_m")
    reference = _series(
        reference_density_m3, count, "reference_density_m3", positive=True
    )
    reduced_mobility = _series(
        electron_reduced_mobility_1_V_m_s,
        count,
        "electron_reduced_mobility_1_V_m_s",
        positive=True,
    )
    mobility_e = reduced_mobility / reference
    electron_current = (
        elementary_charge * density[electrons[0]] * mobility_e * field
    )
    if electric_field_derivative_V_m_s is None:
        derivative = (
            np.zeros(count)
            if np.all(field == field[0])
            else np.gradient(field, time, edge_order=2 if count > 2 else 1)
        )
        derivative_method = (
            "sampled field gradient; no unrecorded switching edge is included"
        )
    else:
        derivative = _series(
            electric_field_derivative_V_m_s,
            count,
            "electric_field_derivative_V_m_s",
        )
        derivative_method = "supplied derivative"
    displacement = epsilon_0 * derivative
    traces = dict(
        time_s=time.copy(),
        electric_field_V_m=field.copy(),
        reference_density_m3=reference.copy(),
        electron_mobility_m2_V_s=mobility_e,
        electron_current_density_A_m2=electron_current,
        displacement_current_density_A_m2=displacement,
    )
    provided = (
        {} if ion_mobilities_m2_V_s is None else dict(ion_mobilities_m2_V_s)
    )
    ion_indices = [
        i for i in range(len(names)) if i != electrons[0] and charge[i] != 0
    ]
    ion_names = {names[i] for i in ion_indices}
    if set(provided) - ion_names:
        raise ValueError(
            f"Unknown ion mobility keys: {sorted(set(provided) - ion_names)}"
        )
    ion_current = np.zeros(count)
    missing = []
    for i in ion_indices:
        name = names[i]
        if name not in provided:
            if np.any(density[i] > 0):
                missing.append(name)
            continue
        mobility_i = _series(
            provided[name], count, f"mobility of {name}", nonnegative=True
        )
        traces[f"ion_mobility_m2_V_s:{name}"] = mobility_i.copy()
        contribution = (
            elementary_charge
            * abs(charge[i])
            * density[i]
            * mobility_i
            * field
        )
        traces[f"ion_current_density_A_m2:{name}"] = contribution
        ion_current += contribution
    traces["known_ion_current_density_A_m2"] = ion_current
    traces["known_conduction_current_density_A_m2"] = (
        electron_current + ion_current
    )
    traces["known_components_current_density_A_m2"] = (
        electron_current + ion_current + displacement
    )
    if not missing:
        traces["ion_current_density_A_m2"] = ion_current
        traces["conduction_current_density_A_m2"] = (
            electron_current + ion_current
        )
        traces["total_current_density_A_m2"] = (
            electron_current + ion_current + displacement
        )
    for label, value in [
        ("area_m2", area_m2),
        ("electrode_gap_m", electrode_gap_m),
    ]:
        if value is not None and (
            not np.isscalar(value) or not np.isfinite(value) or value <= 0
        ):
            raise ValueError(
                f"{label} must be a positive finite scalar or None."
            )
    if area_m2 is not None:
        for name, values in list(traces.items()):
            if "current_density_A_m2" in name:
                traces[name.replace("current_density_A_m2", "current_A")] = (
                    values * area_m2
                )
    if electrode_gap_m is not None:
        traces["voltage_V"] = field * electrode_gap_m
    if not all(np.isfinite(values).all() for values in traces.values()):
        raise ValueError(
            "Current calculation overflowed or produced non-finite values."
        )
    metadata = dict(
        area_m2=None if area_m2 is None else float(area_m2),
        electrode_gap_m=(
            None if electrode_gap_m is None else float(electrode_gap_m)
        ),
        all_populated_ions_have_mobility=not bool(missing),
        missing_ion_mobilities=missing,
        supplied_ion_mobilities=sorted(provided),
        zero_ion_mobility_species=[
            name
            for name in sorted(provided)
            if np.all(np.asarray(provided[name]) == 0)
        ],
        clipped_negative_density_values=clipped,
        derivative_method=derivative_method,
        elementary_charge_C=elementary_charge,
        vacuum_permittivity_F_m=epsilon_0,
        model=(
            "Uniform field and cross-section; conventional drift current plus"
            " vacuum displacement current."
        ),
        limitations=(
            "No sheath/electrode collection, spatial diffusion current,"
            " circuit feedback, dielectric barriers or stray capacitance."
            " Mobilities must match the chemistry conditions. This is"
            " post-processing, not an added particle-loss term."
        ),
    )
    return CurrentResult(traces, metadata)


def save_current_result(result, output):
    r"""保存用の電流履歴と診断メタデータを書き出す。

    Parameters
    ----------
    result : Result
        密度履歴、反応機構、有効な反応係数、診断情報。
    output : str or pathlib.Path
        出力ファイルを保存するディレクトリ。
    """
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "current.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(result.traces)
        writer.writerows(zip(*result.traces.values()))
    np.savez_compressed(output / "current.npz", **result.traces)
    (output / "current_metadata.json").write_text(
        json.dumps(
            result.metadata, ensure_ascii=False, indent=2, allow_nan=False
        ),
        encoding="utf-8",
    )


def plot_current_result(result):
    r"""保存用の電流成分と、計算可能な場合の合計電流を描画する。

    Parameters
    ----------
    result : Result
        密度履歴、反応機構、有効な反応係数、診断情報。

    Returns
    -------
    matplotlib.figure.Figure
        面積がある場合は電流、それ以外は電流密度を示す図。
    """
    import matplotlib.pyplot as plt

    integrated = result.metadata["area_m2"] is not None
    suffix = "current_A" if integrated else "current_density_A_m2"
    fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
    curves = [
        (f"electron_{suffix}", "Electron conduction"),
        (f"known_ion_{suffix}", "Known ion conduction"),
        (f"displacement_{suffix}", "Displacement"),
    ]
    curves.append(
        (f"total_{suffix}", "Total")
        if result.total_available
        else (
            f"known_components_{suffix}",
            "Known components (ions incomplete)",
        )
    )
    for key, label in curves:
        ax.plot(result.traces["time_s"] * 1e9, result.traces[key], label=label)
    ax.set(
        xlabel="Time [ns]",
        ylabel="Current [A]" if integrated else "Current density [A/m$^2$]",
        title="Uniform-field 0D current estimate",
    )
    ax.grid(alpha=0.25)
    ax.legend()
    return fig
