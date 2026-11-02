r"""旧形式の固定組成電子群テーブルの読み込みと補間。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np
from ..paths import BOLSIG


@dataclass
class SwarmData:
    fields: np.ndarray
    energy: np.ndarray
    rates: np.ndarray
    labels: dict

    def interpolate(self, field):
        r"""換算電界を指定し、固定電子群テーブルを補間する。

        Parameters
        ----------
        field : float
            換算電界E/N [Td]。

        Returns
        -------
        tuple
            衝突チャンネルごとの速度係数と平均電子エネルギー [eV]。

        Raises
        ------
        ValueError
            入力が要求する形状、単位、モデルの制約を満たさない場合。

        Examples
        --------
        >>> table = read_swarm()
        >>> rates, energy = table.interpolate(100.0)
        >>> rates.shape
        (44,)
        """
        if (
            not np.isfinite(field)
            or not self.fields[0] <= field <= self.fields[-1]
        ):
            raise ValueError(
                f"E/Nは{self.fields[0]:g}〜{self.fields[-1]:g} Tdに指定してください。"
            )
        # ゼロ係数を含むため、全テーブルを通常の線形補間とする。
        rates = np.array(
            [np.interp(field, self.fields, row) for row in self.rates.T]
        )
        energy = np.interp(field, self.fields, self.energy)
        return rates, float(energy)


def read_swarm(path=BOLSIG):
    r"""旧形式の固定組成BOLSIG電子群テーブルを読み込む。

    Parameters
    ----------
    path : str or pathlib.Path, optional
        読み込む入力ファイルのパス。

    Returns
    -------
    SwarmData
        表形式の換算電界、電子エネルギー、反応係数。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> table = read_swarm()
    >>> rates, energy = table.interpolate(100.0)
    >>> rates.shape, energy > 0
    ((44,), True)
    """
    fields, energies, rows, labels = [], [], [], {}
    blocks = re.split(
        r"^R\d+\s*$", Path(path).read_text(encoding="utf-8"), flags=re.M
    )[1:]
    for block in blocks:

        def parameter(label):
            r"""BOLSIG表から指定したパラメータを読み取る。

            Parameters
            ----------
            label : str
                BOLSIG出力表のパラメータ見出しと完全一致する文字列。

            Returns
            -------
            float
                指定したパラメータ行から読み取った数値。

            Raises
            ------
            ValueError
                入力が要求する形状、単位、モデルの制約を満たさない場合。
            """
            match = re.search(
                "^" + re.escape(label) + r"\s+(\S+)", block, re.M
            )
            if match is None:
                raise ValueError(f"BOLSIG出力に{label}がありません。")
            return float(match[1])

        if not (
            np.isclose(parameter("Mole fraction N2"), 0.8)
            and np.isclose(parameter("Mole fraction O2"), 0.2)
            and np.isclose(parameter("Gas temperature (K)"), 300.0)
        ):
            raise ValueError(
                "このモデルはN2:O2=80:20、300 Kのテーブル専用です。"
            )
        fields.append(parameter("Electric field / N (Td)"))
        energies.append(parameter("Mean energy (eV)"))
        row = np.full(44, np.nan)
        row[0] = 0.0  # C番号と配列添字を一致させる。
        for match in re.finditer(
            r"^C(\d+)\s+(N2|O2)\s+(.+?)\s+([0-9.Ee+\-]+)\s*$", block, re.M
        ):
            cid = int(match[1])
            label = match[2] + " " + " ".join(match[3].split())
            if cid in labels and labels[cid] != label:
                raise ValueError("実行間で衝突過程名が異なります。")
            labels[cid] = label
            row[cid] = float(match[4])
        rows.append(row)
    x, e, rates = np.asarray(fields), np.asarray(energies), np.asarray(rows)
    if (
        not len(x)
        or not np.all(np.diff(x) > 0)
        or not np.all(np.isfinite(rates))
    ):
        raise ValueError("BOLSIGテーブルの列・順序・数値を確認してください。")
    if np.any(rates < 0) or np.any(e <= 0):
        raise ValueError("負の係数または非正の平均エネルギーがあります。")
    return SwarmData(x, e, rates, labels)
