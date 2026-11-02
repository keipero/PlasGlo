r"""比較用の密度CSVとFortranドライバーの条件を読み込む。

外部ZDPlasKin計算との比較に使用する。
"""

from __future__ import annotations

import csv
from pathlib import Path
import re
import numpy as np
from ..chemistry.sample_mechanism import number


def read_csv(path):
    r"""ネイティブ密度履歴を読み、cm^-3からm^-3へ変換する。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。

    Returns
    -------
    tuple
        粒子名、時刻配列 [s]、密度配列 [m^-3]。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     path = Path(directory) / "densities.csv"
    ...     text = chr(10).join(["time_s,E", "0,1", "1,2", ""])
    ...     _ = path.write_text(text, encoding="utf-8")
    ...     names, time, density = read_csv(path)
    >>> names, density.tolist()
    (['E'], [[1000000.0, 2000000.0]])
    """
    with Path(path).open(encoding="utf-8-sig", newline="") as fp:
        rows = csv.reader(fp)
        names = next(rows)
        data = np.array([[float(v) for v in row] for row in rows if row])
    if (
        names[0] != "time_s"
        or not np.isfinite(data).all()
        or not np.all(np.diff(data[:, 0]) > 0)
    ):
        raise ValueError("Expected finite, strictly increasing time_s data.")
    return (
        names[1:],
        data[:, 0],
        data[:, 1:].T * 1e6,
    )  # Archive explicitly writes cm^-3.


def driver_conditions(path):
    r"""元のネイティブ実行ドライバーから計算条件を読み取る。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。

    Returns
    -------
    dict
        圧力、温度、換算電界、初期電子密度、終了時刻。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> from plasglo.paths import SAMPLES_DIR
    >>> values = driver_conditions(
    ...     SAMPLES_DIR / "nitrogen/user_code.F90"
    ... )
    >>> values["temperature_k"]
    300.0
    """
    source = Path(path).read_text(encoding="utf-8")
    found = {}
    for name in (
        "pressure_pa",
        "temperature_k",
        "field_td",
        "seed_density",
        "end_time_s",
    ):
        match = re.search(
            r"(?im)^\s*double precision, parameter :: "
            + name
            + r"\s*=\s*([^!\n]+)",
            source,
        )
        if not match:
            raise ValueError(f"Cannot verify driver parameter {name}")
        found[name] = float(number(match[1]).strip())
    return found
