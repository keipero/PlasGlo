r"""密度・係数・電界・電子群状態・生成消失診断の共通出力。
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import re

import numpy as np

from ..diagnostics.balance import species_balance
from ..visualization.plots import plot_result, plot_species_balance


def save_result(result, output, species=None, *, balance_species=None):
    r"""全ガス共通の結果を保存する。生成消失の保存は明示的に選択する。

    Parameters
    ----------
    result : Result
        密度履歴、反応機構、有効な反応係数、診断情報。
    output : str or pathlib.Path
        出力ファイルを保存するディレクトリ。
    species : sequence of str or None, optional
        描画する粒子種。Noneの場合はメタデータに基づき自動分類する。
    balance_species : sequence of str or None, optional
        生成消失を保存する粒子種。Noneではケースの設定を使う。

    Examples
    --------
    >>> from py3zdplaschem import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> result = simulate(case)
    >>> model = result.mechanism
    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     save_result(result, directory)
    ...     written = (Path(directory) / "densities.csv").is_file()
    >>> written
    True
    """
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    model = result.mechanism

    def table(name, header, rows):
        r"""列名を明示した出力CSVを1つ書き込む。

        Parameters
        ----------
        name : str
            この処理で使用する粒子種・ケース・出力などの識別名。
        header : sequence of str
            出力順に並べたCSVの列名。
        rows : iterable
            CSVへ書き込むスカラー値の行。
        """
        with (output / name).open("w", encoding="utf-8", newline="") as fp:
            writer = csv.writer(fp)
            writer.writerow(header)
            writer.writerows(rows)

    def figure(name, fig):
        r"""図を保存し、描画用の資源を解放する。

        Parameters
        ----------
        name : str
            この処理で使用する粒子種・ケース・出力などの識別名。
        fig : matplotlib.figure.Figure
            保存後に閉じるFigureオブジェクト。
        """
        import matplotlib.pyplot as plt

        fig.savefig(output / name, dpi=160)
        plt.close(fig)

    field = (
        result.field_td
        if result.field_td is not None
        else np.full_like(result.time, model.field_td)
    )
    table(
        "densities.csv",
        ["time_s", *model.species],
        np.column_stack((result.time, result.densities.T)),
    )
    np.savez_compressed(
        output / "densities.npz",
        time_s=result.time,
        densities_m3=result.densities,
        species=np.array(model.species),
        field_td=field,
    )
    table("field_history.csv", ["time_s", "field_td"], zip(result.time, field))
    table(
        "coefficients.csv",
        ["reaction", "order", "coefficient_SI"],
        (
            (r.equation, int(model.orders[i]), float(model.k[i]))
            for i, r in enumerate(model.reactions)
        ),
    )
    history = result.swarm_history
    if history:
        first = history[0]["rates"]
        channel_ids = (
            list(first) if isinstance(first, dict) else list(range(len(first)))
        )
        channel_rates = [
            [h["rates"][key] for key in channel_ids] for h in history
        ]
        np.savez_compressed(
            output / "swarm_history.npz",
            time_s=np.array([h["time_s"] for h in history]),
            field_td=np.array([h["field_td"] for h in history]),
            mean_energy_eV=np.array([h["energy_eV"] for h in history]),
            electron_temperature_K=np.array(
                [h["electron_temperature_K"] for h in history]
            ),
            channel_ids=np.array(channel_ids),
            channel_coefficients_SI=np.array(channel_rates),
            reaction_coefficients_SI=np.array(
                [h["coefficients_SI"] for h in history]
            ),
        )
        targets = list(history[0]["target_densities_m3"])
        header = [
            "time_s",
            "mean_energy_eV",
            "electron_temperature_K",
            "N_m3",
            *[f"n_{s}_m3" for s in targets],
            *[f"x_{s}" for s in targets],
            "unsupported_neutral_fraction",
            "iterations",
            "convergence_error",
            "field_td",
        ]
        rows = [
            [
                h["time_s"],
                h["energy_eV"],
                h["electron_temperature_K"],
                h["target_density_m3"],
                *[h["target_densities_m3"][s] for s in targets],
                *[h["fractions"][s] for s in targets],
                h["unsupported_neutral_fraction"],
                h["iterations"],
                h["error"],
                h["field_td"],
            ]
            for h in history
        ]
        table("swarm_history.csv", header, rows)
    (output / "diagnostics.json").write_text(
        json.dumps(
            result.diagnostics, ensure_ascii=False, indent=2, allow_nan=False
        ),
        encoding="utf-8",
    )
    if result.case is not None:
        case_data = result.case.resolved_data()
        (output / "case.json").write_text(
            json.dumps(
                case_data, ensure_ascii=False, indent=2, allow_nan=False
            ),
            encoding="utf-8",
        )
    with (output / "used_reactions.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as fp:
        writer = csv.writer(fp, delimiter="\t")
        writer.writerow(
            [
                "id",
                "reaction",
                "rate_type",
                "k_SI",
                "order",
                "source_line",
                "source_file",
            ]
        )
        writer.writerows(
            (
                getattr(r, "id", f"R{i + 1:03}"),
                r.equation,
                r.kind,
                model.k[i],
                model.orders[i],
                r.line,
                r.source_file,
            )
            for i, r in enumerate(model.reactions)
        )
    figure("densities.png", plot_result(result, species))
    balances = {}
    for name in (
        result.diagnostics.get("balance_species", [])
        if balance_species is None
        else balance_species
    ):
        values = species_balance(result, name)
        safe_name = re.sub(r'[<>:"/\\|?*]', "_", name)
        stem = f"balance_{model.index[name]:03}_{safe_name}"
        table(
            stem + ".csv",
            ["time_s", *values],
            np.column_stack((result.time, *values.values())),
        )
        figure(stem + ".png", plot_species_balance(result, name))
        balances[name] = stem
    if balances:
        (output / "balances.json").write_text(
            json.dumps(balances, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
