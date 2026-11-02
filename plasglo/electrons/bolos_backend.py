r"""空気用の旧配列API。数値計算本体はmixture.pyにある。
"""

from importlib import import_module
import json

from pathlib import Path
import numpy as np
from ..paths import BOLOS_DIR, CROSS_SECTIONS
from .mixture import grid, parser, solver, normalize_eedf, ElectronSwarm

THREE_BODY_TARGET = "__O2_three_body__"

# Historical export only; the physical channel definitions live in input
# data.
DATA_DIR = import_module("..paths", __package__).DATA_DIR

CHANNELS = [
    (c["target"], c["kind"], c.get("product"), c.get("threshold_eV", 0.0))
    for c in json.loads(
        (DATA_DIR / "cases/air.json").read_text(encoding="utf-8")
    )["electrons"]["profile"]["channels"]
]


class BolosSwarm:
    r"""同梱空気データとC番号配列を扱う互換アダプター。
    """

    def __init__(
        self, path=CROSS_SECTIONS, cells=400, rtol=1e-6, max_iterations=1000
    ):
        r"""衝突データを読み込んで検証し、BOLOSを初期化する。

        Parameters
        ----------
        path : str or pathlib.Path, optional
            読み込む入力ファイルのパス。
        cells : int, optional
            BOLOSで使用するエネルギー格子のセル数。
        rtol : float, optional
            収束判定または時間積分の相対許容誤差。
        max_iterations : int, optional
            1回の収束計算で許容するBOLOSの最大反復回数。

        Examples
        --------
        >>> backend = BolosSwarm(cells=100)
        >>> backend.cells
        100
        """
        import json
        from .mixture import ElectronSwarm
        from ..paths import DATA_DIR

        profile = json.loads(
            (DATA_DIR / "cases/air.json").read_text(encoding="utf-8")
        )["electrons"]["profile"]
        self._swarm = ElectronSwarm(
            [path],
            profile["targets"],
            profile["channels"],
            cells=cells,
            rtol=rtol,
            max_iterations=max_iterations,
        )
        self.labels = {
            int(key[1:]): value for key, value in self._swarm.labels.items()
        }
        self.attachment = self._swarm.density_channels[0][2]

    def __getattr__(self, name):
        r"""旧APIの属性参照を共通の電子ソルバーへ転送する。

        Parameters
        ----------
        name : str
            内部ソルバーから取得する属性名。

        Returns
        -------
        object
            内部ソルバーから取得した属性。

        Examples
        --------
        >>> backend = BolosSwarm(cells=100)
        >>> backend.cells
        100
        """
        return getattr(self._swarm, name)

    def solve(self, field_td, densities, gas_temperature=300.0):
        r"""宣言した標的の密度に対する定常EEDFを求める。

        Parameters
        ----------
        field_td : float
            換算電界E/N [Td]。
        densities : dict of str to float
            基底状態N2とO2の数密度 [m^-3]。
        gas_temperature : float, optional
            一定の気体温度 [K]。

        Returns
        -------
        dict
            SIチャンネル係数、平均エネルギー [eV]、組成比、診断情報。

        Examples
        --------
        >>> backend = BolosSwarm(cells=100)
        >>> state = backend.solve(100.0, {"N2": 2e25, "O2": 5e24})
        >>> state["rates"].shape
        (44,)
        """
        state = self._swarm.solve(field_td, densities, gas_temperature)
        return dict(
            state,
            rates=np.array(
                [0.0]
                + [state["rates"][key] for key in self._swarm.channel_ids]
            ),
        )
