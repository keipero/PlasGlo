r"""ローカルPy3bolosで窒素・アルゴンサンプルを実行する。
"""

from __future__ import annotations
from importlib import import_module

import argparse
import json
from pathlib import Path
import sys

MODEL_DIR = Path(__file__).resolve().parents[2]
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))


SAMPLE_OUTPUT_DIR = import_module("plasglo.paths").SAMPLE_OUTPUT_DIR


def run_case(name, output=None, *, cells=400):
    r"""共通ケース実行・出力を使いながらサンプルCLI/APIを維持する。

    Parameters
    ----------
    name : str
        この処理で使用する粒子種・ケース・出力などの識別名。
    output : str or pathlib.Path, optional
        出力ファイルを保存するディレクトリ。
    cells : int, optional
        BOLOSで使用するエネルギー格子のセル数。

    Returns
    -------
    tuple
        粒子名、時刻 [s]、密度 [m^-3]、ネットワーク、診断情報。

    Raises
    ------
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    >>> names, time, density, model, report = run_case(
    ...     "argon", cells=100
    ... )
    >>> len(names), density.shape[0]
    (3, 3)
    """
    from plasglo import load_case, simulate, save_result

    if name not in ("nitrogen", "argon"):
        raise ValueError("Choose nitrogen or argon.")
    result = simulate(load_case(name).with_updates(electrons={"cells": cells}))
    if output is not None:
        save_result(result, output)
    return (
        result.mechanism.species,
        result.time,
        result.densities,
        result.mechanism,
        result.diagnostics,
    )


def main():
    r"""コマンドライン引数で選択した同梱気体サンプルを実行する。

    Examples
    --------
    プロジェクトのルートから実行します::

        python samples/gases/run_samples.py --sample argon
    """
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--sample", choices=["nitrogen", "argon", "all"], default="all"
    )
    cli.add_argument("--output", type=Path, default=SAMPLE_OUTPUT_DIR)
    cli.add_argument(
        "--cells",
        type=int,
        default=400,
        help="BOLOS energy grid cells (at least 50)",
    )
    args = cli.parse_args()
    for name in (
        ["nitrogen", "argon"] if args.sample == "all" else [args.sample]
    ):
        print(f"Running {name}...", flush=True)
        *_, report = run_case(name, args.output / name, cells=args.cells)
        print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
