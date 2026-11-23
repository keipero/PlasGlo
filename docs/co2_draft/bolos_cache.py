r"""``ElectronSwarm.solve`` の結果を、組成が近ければ使い回す（核は変更しない）。

EEDF と電子衝突の反応係数は、E/N、気体温度、標的の組成（割合）、密度規格化
チャンネルの相手の密度だけで決まる（電子間衝突がないため、電子密度には
よらない）。NRP のパルス列では、E/N は 2 通り（パルス中と 0）で、組成は
1サイクルにほとんど動かないので、組成を量子化した鍵で結果を使い回せる。

``install(tol)`` は、``plasglo.electrons.mixture.ElectronSwarm.solve`` を、
キャッシュ付きの関数に差し替える。``tol`` は、標的の割合を量子化する幅
（絶対値）と、密度規格化の相手の密度を量子化する幅（対数）。鍵が同じ範囲の
組成には、最初に計算した値を返す（その範囲の誤差は、割合で最大 ``tol``）。

使い方::

    import bolos_cache
    stats = bolos_cache.install(1e-3)
    ...  # simulate など
    print(stats)  # {'hit': ..., 'miss': ...}
"""

import copy
import math
import sys
from pathlib import Path

MODEL = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MODEL))
from plasglo.electrons import mixture  # noqa: E402


def make_key(swarm, field_td, densities, gas_temperature, tol):
    r"""キャッシュの鍵を作る。

    Parameters
    ----------
    swarm : plasglo.electrons.mixture.ElectronSwarm
        電子群計算（断面積・標的・チャンネルの設定を鍵に含める）。
    field_td : float
        換算電界 [Td]。
    densities : dict of str to float
        化学種の数密度 [m^-3]。
    gas_temperature : float
        気体温度 [K]。
    tol : float
        標的の割合を量子化する幅（絶対値）。密度規格化の相手の密度（対数）は、
        その 10 倍の幅で量子化する。

    Returns
    -------
    tuple
        ハッシュ可能な鍵。
    """
    species = list(swarm.targets.values())
    total = sum(densities[s] for s in species)
    fractions = tuple(round(densities[s] / total / tol) for s in species)
    # 密度規格化の相手（O2 の三体付着）は、EEDF への影響が小さいので、
    # 割合より 10 倍粗く量子化する（O2 は 1 サイクルに 0.2% 増える）。
    partners = tuple(
        round(math.log(max(densities[spec["species"]], 1.0)) / (10 * tol))
        for _, _, _, spec in swarm.density_channels
    )
    setup = (
        tuple(str(p) for p in swarm.paths),
        tuple(sorted(swarm.targets.items())),
        tuple(c["id"] for c in swarm.channels),
        swarm.cells,
    )
    return (setup, round(field_td, 9), round(gas_temperature, 6), fractions,
            partners)


def install(tol=1e-3):
    r"""``ElectronSwarm.solve`` を、キャッシュ付きの関数に差し替える。

    Parameters
    ----------
    tol : float, optional
        量子化の幅。小さいほど正確で、再計算が増える。

    Returns
    -------
    dict
        ``hit``（使い回した回数）と ``miss``（実際に計算した回数）の
        カウンター。差し替えの後も、同じ辞書が更新される。
    """
    original = mixture.ElectronSwarm.solve
    store = {}
    stats = {"hit": 0, "miss": 0}

    def solve(self, field_td, densities, gas_temperature=300.0):
        r"""キャッシュを調べ、なければ元の ``solve`` で計算して、保存する。

        Parameters
        ----------
        field_td : float
            換算電界 [Td]。
        densities : dict of str to float
            化学種の数密度 [m^-3]。
        gas_temperature : float, optional
            気体温度 [K]。

        Returns
        -------
        dict
            元の ``solve`` と同じ形の結果（キャッシュからは、複製を返す）。
        """
        key = make_key(self, field_td, densities, gas_temperature, tol)
        if key in store:
            stats["hit"] += 1
            return copy.deepcopy(store[key])
        stats["miss"] += 1
        state = original(self, field_td, densities, gas_temperature)
        store[key] = copy.deepcopy(state)
        return state

    solve.__doc__ = original.__doc__
    mixture.ElectronSwarm.solve = solve
    return stats
