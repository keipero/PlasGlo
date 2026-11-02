r"""粒子密度と生成消失のグラフ。
"""

from __future__ import annotations

import numpy as np
from ..chemistry.reactions import composition
from ..diagnostics.balance import species_balance


def plot_species_balance(result, species):
    r"""生成消失と消失寿命を描き、正味の減衰時間とは区別する。

    Parameters
    ----------
    result : Result
        密度履歴、反応機構、有効な反応係数、診断情報。
    species : str
        生成・消失または密度を調べる粒子種名。

    Returns
    -------
    matplotlib.figure.Figure
        生成消失速度と消失寿命を示す図。

    Examples
    --------
    >>> from py3zdplaschem import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> result = simulate(case)
    >>> model = result.mechanism
    >>> import matplotlib.pyplot as plt
    >>> fig = plot_species_balance(result, "AR^+")
    >>> len(fig.axes)
    2
    >>> plt.close(fig)
    """
    import matplotlib.pyplot as plt

    balance = species_balance(result, species)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), layout="constrained")
    time_ns = result.time * 1e9
    for name, label in [
        ("production_m3_s", "Production"),
        ("loss_m3_s", "Loss"),
    ]:
        values = balance[name]
        axes[0].plot(
            time_ns, np.where(values > 0, values, np.nan), label=label
        )
    axes[0].set(
        yscale="log",
        ylabel=r"Rate [m$^{-3}$ s$^{-1}$]",
        title=f"{species}: production and loss",
    )
    axes[0].legend()
    axes[1].plot(time_ns, balance["loss_lifetime_s"] * 1e9)
    axes[1].set(ylabel="Loss lifetime n / L [ns]", xlabel="Time [ns]")
    duration = result.diagnostics.get("field_duration_ns")
    for ax in axes:
        ax.grid(alpha=0.25)
        if duration is not None and 0 <= duration < time_ns[-1]:
            ax.axvline(duration, color="black", ls="--", lw=1)
    return fig


def plot_result(result, species=None):
    r"""粒子種と気体組成のメタデータに基づき密度履歴を描画する。

    Parameters
    ----------
    result : Result
        密度履歴、反応機構、有効な反応係数、診断情報。
    species : sequence of str or None, optional
        描画する粒子種。Noneの場合はメタデータに基づき自動分類する。

    Returns
    -------
    matplotlib.figure.Figure
        電荷と励起状態で分類して表示する密度の図。

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
    >>> import matplotlib.pyplot as plt
    >>> fig = plot_result(result)
    >>> bool(fig.axes)
    True
    >>> plt.close(fig)
    """
    import matplotlib.pyplot as plt

    model = result.mechanism
    backgrounds = set(model.initial_fractions)
    metadata = {item.name: item for item in getattr(model, "species_data", [])}

    def charge(name):
        r"""メタデータまたは旧入力の規約から粒子種の電荷を取得する。

        Parameters
        ----------
        name : str
            この処理で使用する粒子種・ケース・出力などの識別名。

        Returns
        -------
        int
            選択した粒子種の符号付き整数電荷数。
        """
        return (
            metadata[name].charge if name in metadata else composition(name)[2]
        )

    if species is not None:
        selected = [s for s in species if s not in backgrounds]
        background = [s for s in species if s in backgrounds]
        panels = []
        if selected:
            panels.append(("Selected species", selected))
        if background:
            panels.append(("Background gas", background))
        if not panels:
            panels = [("Selected species", [])]
    else:
        charged = [s for s in model.species if charge(s) != 0]
        excited = [
            s
            for s in model.species
            if (metadata[s].excited if s in metadata else "(" in s or "*" in s)
        ]
        neutral = [
            s for s in model.species if charge(s) == 0 and s not in excited
        ]
        panels = [
            ("Charged species", charged),
            ("Excited species / atoms", excited),
            ("Ground-state neutrals", neutral),
        ]
    floor = result.diagnostics.get("atol_density_m3", 1e3)
    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=(10, 3.8 * len(panels)),
        squeeze=False,
        layout="constrained",
    )
    for ax, (title, names) in zip(axes[:, 0], panels):
        for name in names:
            if name not in model.index:
                raise ValueError(f"未知の粒子種：{name}")
            values = result.densities[model.index[name]]
            # 絶対許容誤差以下の微量密度は描かない。CSV/NPZは未加工。
            if np.any(values > floor):
                ax.plot(
                    result.time * 1e9,
                    np.where(values > floor, values, np.nan),
                    label=name,
                )
        ax.set(
            xlabel="Time [ns]",
            ylabel=r"Number density [m$^{-3}$]",
            title=title,
        )
        ax.grid(alpha=0.25)
        duration = result.diagnostics.get("field_duration_ns")
        if duration is not None and 0 <= duration < result.time[-1] * 1e9:
            ax.axvline(duration, color="black", ls="--", lw=1, label="E/N off")
        if ax.lines:
            ax.legend(fontsize=8, ncol=3, loc="best")
    field_label = f"E/N = {model.field_td:g} Td"
    duration = result.diagnostics.get("field_duration_ns")
    if duration is not None:
        field_label += f" until {duration:g} ns, then 0 Td"
    mixture = ", ".join(
        f"{name}={value:g}" for name, value in model.initial_fractions.items()
    )
    fig.suptitle(
        f"0D kinetics: {field_label}, Tg = {model.tg:g} K\n"
        f"{mixture}; "
        f"plotted density > {floor:g} m$^{{-3}}$"
    )
    if (
        result.diagnostics.get("max_unsupported_neutral_fraction") or 0
    ) > 0.01 or result.diagnostics.get("max_ion_fraction", 0) > 0.01:
        axes[0, 0].set_title(
            axes[0, 0].get_title()
            + "\nComposition / ionization screen exceeded: incomplete"
            " electron-collision data",
            fontsize=10,
            color="darkred",
        )
    return fig
