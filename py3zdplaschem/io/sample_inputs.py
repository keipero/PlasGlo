r"""サンプルの出典照合とZDPlasKin粒子種リストの読み込み。
"""

import hashlib
import json
from pathlib import Path
import re


def verify_sample_inputs(folder):
    r"""記録したSHA-256でサンプル入力を照合する。

    Parameters
    ----------
    folder : str or pathlib.Path
        出典マニフェストと入力ファイルがあるサンプルディレクトリ。

    Returns
    -------
    list of dict
        元のサンプルファイルの照合済み出典情報。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> from py3zdplaschem.paths import SAMPLES_DIR
    >>> entries = verify_sample_inputs(SAMPLES_DIR / "argon")
    >>> len(entries)
    3
    """
    folder = Path(folder)
    provenance = json.loads(
        (folder / "provenance.json").read_text(encoding="utf-8")
    )
    for item in provenance:
        path = folder / item["local_file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(
                f"Sample input differs from its recorded source: {path}"
            )
    return provenance


def read_species(path):
    r"""元入力のSPECIESブロックから粒子名を読み込む。

    Parameters
    ----------
    path : str or pathlib.Path
        読み込む入力ファイルのパス。

    Returns
    -------
    list of str
        元の順序を保った大文字表記の粒子種名。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> from py3zdplaschem.paths import SAMPLES_DIR
    >>> read_species(SAMPLES_DIR / "argon/kinet.inp")
    ['E', 'AR', 'AR^+']
    """
    text = Path(path).read_text(encoding="utf-8")
    block = re.search(r"(?ms)^SPECIES\s*\n(.*?)^END\s*$", text)
    if block is None:
        raise ValueError("Missing SPECIES block.")
    return (
        " ".join(line.split("#")[0] for line in block[1].splitlines())
        .upper()
        .split()
    )
