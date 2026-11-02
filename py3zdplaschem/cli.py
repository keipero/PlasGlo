r"""気体ケースを実行する。旧形式の空気用数値引数にも対応する。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from .paths import OUTPUTS_DIR
from .solvers.simulation import simulate
from .io.results import save_result
from .io.cases import load_case


def main():
    r"""引数を読み取り、気体ケースを実行して結果を保存する。

    Examples
    --------
    プロジェクトのルートから実行します::

        python kinetics.py --case nitrogen
    """
    import matplotlib

    matplotlib.use("Agg")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", help="Case JSON path, or bundled air/nitrogen/argon"
    )
    parser.add_argument(
        "--field-td",
        type=float,
        default=100.0,
        help="一定のE/N [Td]。tableの場合は1〜2000",
    )
    parser.add_argument("--end-ns", type=float, default=100.0)
    parser.add_argument(
        "--field-duration-ns",
        type=float,
        default=None,
        help="E/N印加時間 [ns]。以後は0 Tdで計算継続。省略時は全計算時間に印加",
    )
    parser.add_argument(
        "--seed-density", type=float, default=1e12, help="初期電子密度 [m^-3]"
    )
    parser.add_argument("--pressure-pa", type=float, default=101325.0)
    parser.add_argument(
        "--gas-temperature-k",
        type=float,
        default=300.0,
        help="一定の気体温度 [K]",
    )
    parser.add_argument(
        "--n2-fraction",
        type=float,
        default=0.8,
        help="N2のモル分率 [0〜1]。O2は1から引いた残り",
    )
    parser.add_argument("--points", type=int, default=501)
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol-density", type=float, default=1e3)
    parser.add_argument(
        "--swarm-backend", choices=["bolos", "table"], default="bolos"
    )
    parser.add_argument(
        "--bolos-update-ns",
        type=float,
        default=1.0,
        help="BOLOSと化学計算の結合間隔 [ns]",
    )
    parser.add_argument("--bolos-cells", type=int, default=400)
    parser.add_argument("--bolos-rtol", type=float, default=1e-6)
    parser.add_argument(
        "--cross-sections", type=Path, help="C1〜C43対応のPhelps断面積ファイル"
    )
    parser.add_argument(
        "--species", nargs="+", help="描画する粒子名を空白区切りで指定"
    )
    parser.add_argument(
        "--balance-species", nargs="+", help="生成・消失履歴を保存する粒子名"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.case:
        case = load_case(args.case)
        supplied = {
            value.split("=")[0]
            for value in sys.argv[1:]
            if value.startswith("--")
        }
        if supplied & {"--seed-density", "--n2-fraction", "--swarm-backend"}:
            parser.error(
                "With --case, set initial densities, composition and backend"
                " in the case JSON."
            )
        groups = dict(gas={}, field={}, solver={}, electrons={})
        mapping = {
            "pressure_pa": ("gas", "pressure_Pa"),
            "gas_temperature_k": ("gas", "temperature_K"),
            "field_td": ("field", "td"),
            "field_duration_ns": ("field", "duration_ns"),
            "end_ns": ("solver", "end_ns"),
            "points": ("solver", "points"),
            "rtol": ("solver", "rtol"),
            "atol_density": ("solver", "atol_density_m3"),
            "bolos_cells": ("electrons", "cells"),
            "bolos_rtol": ("electrons", "rtol"),
            "bolos_update_ns": ("electrons", "update_ns"),
            "cross_sections": ("electrons", "cross_sections"),
        }
        for argument, (group, name) in mapping.items():
            if "--" + argument.replace("_", "-") in supplied:
                value = getattr(args, argument)
                groups[group][name] = (
                    str(value.resolve()) if isinstance(value, Path) else value
                )
        result = simulate(case.with_updates(**groups))
        output = args.output or OUTPUTS_DIR / "cases" / case.source.stem
    else:
        result = simulate(
            args.field_td,
            args.end_ns,
            args.seed_density,
            args.pressure_pa,
            args.points,
            args.rtol,
            args.atol_density,
            swarm_backend=args.swarm_backend,
            bolos_update_ns=args.bolos_update_ns,
            bolos_cells=args.bolos_cells,
            bolos_rtol=args.bolos_rtol,
            cross_sections=args.cross_sections,
            field_duration_ns=args.field_duration_ns,
            gas_temperature_k=args.gas_temperature_k,
            n2_fraction=args.n2_fraction,
        )
        output = args.output or OUTPUTS_DIR / "air/default"
    save_result(
        result, output, args.species, balance_species=args.balance_species
    )
    print(json.dumps(result.diagnostics, ensure_ascii=False, indent=2))
    print(f"保存先：{output.resolve()}")
