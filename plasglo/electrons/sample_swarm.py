r"""共通BOLOSソルバーを呼ぶ純気体サンプル向けの互換入口。
"""

from pathlib import Path

from .mixture import ElectronSwarm, parser


def solve_sample_swarm(path, field, temperature, cells=400, *, target="N2"):
    r"""共通の電子ソルバーを使用して純気体のEEDFを求める。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。
    field : float
        換算電界E/N [Td]。
    temperature : float
        一定の気体温度 [K]。
    cells : int, optional
        BOLOSで使用するエネルギー格子のセル数。
    target : str, optional
        断面積データで使用される純気体の基底状態の標的名。

    Returns
    -------
    tuple
        SIチャンネル係数、EEDF診断、解析した衝突定義。

    Examples
    --------
    >>> from plasglo.paths import SAMPLES_DIR
    >>> rates, info, collisions = solve_sample_swarm(
    ...     SAMPLES_DIR / "argon/bolsigdb.dat", 51.0, 301.0,
    ...     cells=100, target="Ar",
    ... )
    >>> rates["AR->AR^+"] > 0
    True
    """
    with Path(path).open(encoding="utf-8") as fp:
        collisions = parser.parse(fp)
    channels = [
        dict(
            id=f"channel_{i}",
            target=c["target"],
            kind=c["kind"],
            product=c.get("product"),
            threshold_eV=c.get("threshold", 0.0),
        )
        for i, c in enumerate(collisions)
    ]
    swarm = ElectronSwarm([path], {target: target}, channels, cells=cells)
    state = swarm.solve(field, {target: 1.0}, temperature)
    rates = {
        f'{c["target"]}->{c["product"]}'.replace(" ", "").upper(): state[
            "rates"
        ][channel["id"]]
        for channel, c in zip(channels, collisions)
        if c.get("product")
    }
    report = dict(
        mean_energy_eV=state["energy_eV"],
        collision_count=len(collisions),
        **{
            key: state[key]
            for key in (
                "iterations",
                "error",
                "uncorrected_eedf_min",
                "eedf_tail_repaired",
            )
        },
    )
    return rates, report, collisions
