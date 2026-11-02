"""係数履歴から反応速度、化学種別の寄与、反応経路を再構成する。"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Literal, TypedDict

import numpy as np


class ReactionContribution(TypedDict):
    """生成側または消滅側で順位付けした反応の積算寄与。"""

    reaction_id: str
    equation: str
    integrated_m3: float
    fraction: float


class PathwayNode(TypedDict):
    """化学種または反応を表すグラフの節点。"""

    id: str
    label: str
    kind: Literal["species", "reaction"]


class PathwayEdge(TypedDict):
    """化学種と反応を結ぶ辺。重みは粒子数密度の積算量とする。"""

    source: str
    target: str
    stoichiometry: int
    integrated_m3: float


class PathwayGraph(TypedDict):
    """原子対応を仮定しない、化学種と反応の二部グラフ。"""

    nodes: list[PathwayNode]
    edges: list[PathwayEdge]


@dataclass
class PathwayAnalysis:
    """指定した化学種と時間区間に対する反応経路解析結果。

    Attributes
    ----------
    species : str
        寄与を解析した化学種。
    time_s : numpy.ndarray
        保存時刻、区間端点、係数更新時刻を合わせた解析時刻 [s]。
    reaction_ids, equations : tuple of str
        各配列の反応軸に対応する識別子と反応式。
    rates_m3_s : numpy.ndarray
        各反応の速度 [m^-3 s^-1]。形状は (反応数, 時刻数)。
    contributions_m3_s : numpy.ndarray
        対象化学種への正味寄与。同じ形状で、生成を正とする。
    integrated_rates_m3 : numpy.ndarray
        区間内の反応速度の積算量 [m^-3]。
    integrated_contributions_m3 : numpy.ndarray
        対象化学種への正味積算寄与 [m^-3]。
    production, loss : list of ReactionContribution
        生成側、消滅側それぞれの積算寄与の降順リスト。
        integrated_m3は両側とも正値、fractionは各側の総量に対する比率。
    density_change_m3 : float
        区間両端の保存密度から得た対象化学種の密度差 [m^-3]。
    graph : PathwayGraph
        積算反応速度が正の全反応のグラフ。辺は反応物側・生成物側の
        係数で重み付けする。両辺に存在する電子等も含む。
        第三体グループは物質移動の節点に含めない。
    """

    species: str
    time_s: np.ndarray
    reaction_ids: tuple[str, ...]
    equations: tuple[str, ...]
    rates_m3_s: np.ndarray
    contributions_m3_s: np.ndarray
    integrated_rates_m3: np.ndarray
    integrated_contributions_m3: np.ndarray
    production: list[ReactionContribution]
    loss: list[ReactionContribution]
    density_change_m3: float
    graph: PathwayGraph


def _rank(ids, equations, values):
    """正の積算量を持つ反応を降順に並べる。

    Parameters
    ----------
    ids, equations : tuple of str
        反応順に並べた識別子と反応式。
    values : numpy.ndarray
        生成量または消滅量を正値で表した配列 [m^-3]。

    Returns
    -------
    list of ReactionContribution
        寄与が正の反応と各側の合計に対する比率。
    """
    total = float(values.sum())
    return [
        ReactionContribution(
            reaction_id=ids[j], equation=equations[j],
            integrated_m3=float(values[j]), fraction=float(values[j] / total),
        )
        for j in np.argsort(-values, kind="stable") if values[j] > 0
    ]


def _graph(model, ids, equations, integrated):
    """反応の左右の係数から二部グラフを作る。

    Parameters
    ----------
    model : MassActionKernel
        化学種、反応物係数、正味の化学量論行列を持つ反応機構。
    ids, equations : tuple of str
        反応順に並べた識別子と反応式。
    integrated : numpy.ndarray
        各反応の積算速度 [m^-3]。

    Returns
    -------
    PathwayGraph
        化学種と反応を節点とし、反応参加量を辺の重みとしたグラフ。
    """
    nodes: list[PathwayNode] = []
    edges: list[PathwayEdge] = []
    participating = set()
    for j in np.flatnonzero(integrated > 0):
        reaction = f"reaction:{ids[j]}"
        nodes.append(PathwayNode(
            id=reaction, label=equations[j], kind="reaction",
        ))
        for side, counts in enumerate((
            model.powers[:, j], model.powers[:, j] + model.stoich[:, j],
        )):
            for i in np.flatnonzero(counts > 0):
                name = model.species[i]
                participating.add(name)
                species = f"species:{name}"
                edges.append(PathwayEdge(
                    source=species if side == 0 else reaction,
                    target=reaction if side == 0 else species,
                    stoichiometry=int(counts[i]),
                    integrated_m3=float(counts[i] * integrated[j]),
                ))
    nodes.extend(
        PathwayNode(id=f"species:{s}", label=s, kind="species")
        for s in model.species if s in participating
    )
    return PathwayGraph(nodes=nodes, edges=edges)


def analyze_pathways(result, species, time_range_s=None) -> PathwayAnalysis:
    """保存した密度と係数履歴から反応別の経路寄与を求める。

    Parameters
    ----------
    result : Result
        simulateが返す密度、反応機構、係数履歴を含む計算結果。
    species : str
        生成・消滅への寄与を順位付けする化学種。
    time_range_s : tuple of float or None, optional
        積算する開始・終了時刻 [s]。Noneでは保存時刻の全区間。
        区間は保存範囲内とし、開始時刻は終了時刻より前とする。

    Returns
    -------
    PathwayAnalysis
        全反応の速度と積算量、対象化学種への寄与、順位、グラフ。
        反応速度は各解析時刻直後に有効な係数で評価する。

    Raises
    ------
    ValueError
        化学種、時間区間、密度配列、係数履歴が不正な場合。

    Notes
    -----
    区間端点と係数更新時刻を解析格子に追加し、未保存時刻の密度は
    線形補間する。係数を各更新区間内で固定して台形積分するため、
    係数の不連続点をまたいで平均しない。積算値の精度は密度の保存
    間隔に依存し、ODE積分器が直接積算した値とは異なる。
    正味寄与はS_ij R_jであり、両辺に同数存在する粒子は順位に含まれない。
    グラフの辺は反応参加量を示し、原子ごとの移動量は意味しない。

    Examples
    --------
    >>> from py3zdplaschem import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> result = simulate(case)
    >>> paths = analyze_pathways(result, "AR^+")
    >>> paths.rates_m3_s.shape == (2, 3)
    True
    >>> paths.production[0]["fraction"] > 0
    True
    """
    model = result.mechanism
    if species not in model.index:
        raise ValueError(f"未知の化学種: {species}")
    saved = np.asarray(result.time, dtype=float)
    densities = np.asarray(result.densities, dtype=float)
    if (
        saved.ndim != 1 or len(saved) < 2
        or not np.isfinite(saved).all() or np.any(np.diff(saved) <= 0)
        or densities.shape != (len(model.species), len(saved))
        or not np.isfinite(densities).all()
    ):
        raise ValueError("有限な密度と、単調増加する2点以上の時刻が必要です。")
    bounds = np.asarray(
        (saved[0], saved[-1]) if time_range_s is None else time_range_s,
        dtype=float,
    )
    if (
        bounds.shape != (2,) or not np.isfinite(bounds).all()
        or bounds[0] < saved[0] or bounds[1] > saved[-1]
        or bounds[0] >= bounds[1]
    ):
        raise ValueError("積算区間は保存範囲内の開始時刻・終了時刻で指定します。")
    start, stop = bounds
    history = result.swarm_history
    updates = np.array([h["time_s"] for h in history], dtype=float)
    if history and (
        not np.isfinite(updates).all() or np.any(np.diff(updates) <= 0)
        or updates[0] > start
    ):
        raise ValueError("係数履歴は区間開始以前から時刻順に並べてください。")
    time = np.unique(np.r_[
        bounds, saved[(saved > start) & (saved < stop)],
        updates[(updates > start) & (updates < stop)],
    ])
    density = np.array([np.interp(time, saved, row) for row in densities])
    normalized = np.maximum(density / model.density, 0.0)
    positions = (
        np.searchsorted(updates, time, side="right") - 1
        if history else np.zeros(len(time), dtype=int)
    )
    coefficients = (
        [np.asarray(h["coefficients_SI"]) for h in history]
        if history else [np.asarray(model.k)]
    )
    nr = len(model.reactions)
    if any(
        k.shape != (nr,) or not np.isfinite(k).all() or np.any(k < 0)
        for k in coefficients
    ):
        raise ValueError("係数は反応数に一致する有限な非負の配列が必要です。")
    rates = np.column_stack([
        model.physical_reaction_rates(normalized[:, j], coefficients[p])
        for j, p in enumerate(positions)
    ])
    integrated = np.zeros(nr)
    for j, dt in enumerate(np.diff(time)):
        # 切替時刻の右端だけ、切替前の係数で左極限を求める。
        right = (
            model.physical_reaction_rates(
                normalized[:, j + 1], coefficients[positions[j]],
            )
            if positions[j] != positions[j + 1] else rates[:, j + 1]
        )
        integrated += 0.5 * dt * (rates[:, j] + right)
    row = model.stoich[model.index[species]]
    contributions = row * integrated
    ids = tuple(
        getattr(r, "id", f"R{j + 1:03}")
        for j, r in enumerate(model.reactions)
    )
    equations = tuple(r.equation for r in model.reactions)
    return PathwayAnalysis(
        species=species, time_s=time, reaction_ids=ids, equations=equations,
        rates_m3_s=rates, contributions_m3_s=row[:, None] * rates,
        integrated_rates_m3=integrated,
        integrated_contributions_m3=contributions,
        production=_rank(ids, equations, np.maximum(contributions, 0)),
        loss=_rank(ids, equations, np.maximum(-contributions, 0)),
        density_change_m3=float(
            density[model.index[species], -1]
            - density[model.index[species], 0]
        ),
        graph=_graph(model, ids, equations, integrated),
    )


def save_pathways(analysis: PathwayAnalysis, output) -> None:
    """反応速度配列、寄与の順位、グラフを専用ディレクトリに保存する。

    Parameters
    ----------
    analysis : PathwayAnalysis
        analyze_pathwaysが返す解析結果。
    output : str or pathlib.Path
        保存先。同名のファイルは上書きする。

    Notes
    -----
    reaction_rates.npzに配列、reaction_ranking.csvに生成側・消滅側の
    順位を保存する。pathway_graph.jsonには節点・辺と解析条件を保存する。
    グラフは全反応を対象とし、順位表は対象化学種の正味寄与を示す。

    Examples
    --------
    >>> import tempfile
    >>> from py3zdplaschem import load_case, simulate
    >>> case = load_case("argon").with_updates(
    ...     solver={"end_ns": 0.01, "points": 3},
    ...     electrons={"cells": 100},
    ... )
    >>> analysis = analyze_pathways(simulate(case), "AR^+")
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     save_pathways(analysis, directory)
    ...     exists = (Path(directory) / "pathway_graph.json").is_file()
    >>> exists
    True
    """
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output / "reaction_rates.npz", time_s=analysis.time_s,
        species=analysis.species, reaction_ids=np.array(analysis.reaction_ids),
        equations=np.array(analysis.equations), rates_m3_s=analysis.rates_m3_s,
        contributions_m3_s=analysis.contributions_m3_s,
        integrated_rates_m3=analysis.integrated_rates_m3,
        integrated_contributions_m3=analysis.integrated_contributions_m3,
        density_change_m3=analysis.density_change_m3,
    )
    with (output / "reaction_ranking.csv").open(
        "w", encoding="utf-8", newline="",
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "direction", "rank", "reaction_id", "equation",
            "integrated_m3", "fraction",
        ])
        for direction, rows in (
            ("production", analysis.production), ("loss", analysis.loss),
        ):
            writer.writerows(
                (direction, rank, row["reaction_id"], row["equation"],
                 row["integrated_m3"], row["fraction"])
                for rank, row in enumerate(rows, 1)
            )
    graph = dict(
        schema_version=1, species=analysis.species,
        time_range_s=[float(analysis.time_s[0]), float(analysis.time_s[-1])],
        weight_unit="m^-3", weight_definition="stoichiometry * integral(R dt)",
        density_change_m3=analysis.density_change_m3,
        integrated_net_m3=float(analysis.integrated_contributions_m3.sum()),
        **analysis.graph,
    )
    (output / "pathway_graph.json").write_text(
        json.dumps(graph, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
