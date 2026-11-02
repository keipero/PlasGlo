r"""BOLOSと結合した0次元化学反応計算のmarimo入口。
"""

import marimo

__generated_with = "0.25.2"
app = marimo.App()


@app.cell
def _():
    r"""marimoセルを実行し、宣言した値を後続セルへ渡す。

    Returns
    -------
    tuple
        return文の順序で後続セルへ渡す値。
    """
    from pathlib import Path
    import sys
    import marimo as mo
    import numpy as np

    # CLI/Jupyter/marimo share chemistry, coupling and unit conversions.
    model_dir = Path(__file__).resolve().parent.parent
    if str(model_dir) not in sys.path:
        sys.path.insert(0, str(model_dir))
    from plasglo import (
        simulate,
        plot_result,
        save_result,
        read_reactions,
        plot_species_balance,
        species_balance,
    )

    return (
        mo,
        model_dir,
        np,
        plot_result,
        plot_species_balance,
        read_reactions,
        save_result,
        simulate,
        species_balance,
    )


@app.cell
def _(read_reactions):
    r"""marimoセルを実行し、宣言した値を後続セルへ渡す。

    Parameters
    ----------
    read_reactions : callable
        元の空気反応リストを読み込む関数。
    """
    reactions = read_reactions()
    print(f"反応数: {len(reactions)}")
    return


@app.cell
def _(model_dir):
    r"""marimoセルを実行し、宣言した値を後続セルへ渡す。

    Parameters
    ----------
    model_dir : pathlib.Path
        入出力の基準となるプロジェクトのディレクトリ。

    Returns
    -------
    tuple
        return文の順序で後続セルへ渡す値。
    """
    FIELD_TD = 1000.0  # E/N [Td], N = n(N2) + n(O2)
    FIELD_DURATION_NS = 0.05  # E/N印加時間。以後0 Td。Noneなら計算終了まで印加
    END_NS = 3.0
    POINTS = 1201  # 3 nsで0.0025 ns刻み。初期の立上がりも描画する
    SEED_DENSITY = 1e13  # m^-3
    PRESSURE_PA = 101325.0
    GAS_TEMPERATURE_K = 300.0
    N2_FRACTION = 0.8  # N2モル分率。O2モル分率は1-N2_FRACTION
    SWARM_BACKEND = "bolos"  # "table" で旧固定BOLSIGテーブルと比較
    BOLOS_UPDATE_NS = 0.01  # 間隔を短くした計算と比較して収束を確認
    BOLOS_CELLS = 400
    BOLOS_RTOL = 1e-6
    SPECIES = ["e", "N2+", "O2+", "NO", "N2*(C)", "N2", "O2"]  # 全粒子は None
    OUTPUT = model_dir / "outputs/air/default"
    return (
        BOLOS_CELLS,
        BOLOS_RTOL,
        BOLOS_UPDATE_NS,
        END_NS,
        FIELD_DURATION_NS,
        FIELD_TD,
        GAS_TEMPERATURE_K,
        N2_FRACTION,
        OUTPUT,
        POINTS,
        PRESSURE_PA,
        SEED_DENSITY,
        SPECIES,
        SWARM_BACKEND,
    )


@app.cell
def _(
    BOLOS_CELLS,
    BOLOS_RTOL,
    BOLOS_UPDATE_NS,
    END_NS,
    FIELD_DURATION_NS,
    FIELD_TD,
    GAS_TEMPERATURE_K,
    N2_FRACTION,
    POINTS,
    PRESSURE_PA,
    SEED_DENSITY,
    SWARM_BACKEND,
    simulate,
):
    r"""marimoセルを実行し、宣言した値を後続セルへ渡す。

    Parameters
    ----------
    BOLOS_CELLS : int
        BOLOSのエネルギー格子セル数。
    BOLOS_RTOL : float
        BOLOSの相対収束許容誤差。
    BOLOS_UPDATE_NS : float
        電子計算と化学反応計算の結合間隔 [ns]。
    END_NS : float
        計算終了時刻 [ns]。
    FIELD_DURATION_NS : float or None
        電界の印加時間 [ns]。Noneの場合は一定電界とする。
    FIELD_TD : float
        換算電界 [Td]。
    GAS_TEMPERATURE_K : float
        一定の気体温度 [K]。
    N2_FRACTION : float
        空気プリセットにおける窒素の初期モル分率。
    POINTS : int
        要求する出力サンプル数。
    PRESSURE_PA : float
        初期状態を規格化するための基準圧力 [Pa]。
    SEED_DENSITY : float
        初期電子密度 [m^-3]。
    SWARM_BACKEND : str
        設定セルで選択した電子反応係数のバックエンド。
    simulate : callable
        共通の気体ケース実行関数。

    Returns
    -------
    tuple
        return文の順序で後続セルへ渡す値。
    """
    _settings = dict(
        field_td=FIELD_TD,
        end_ns=END_NS,
        seed_density=SEED_DENSITY,
        pressure_pa=PRESSURE_PA,
        swarm_backend=SWARM_BACKEND,
        field_duration_ns=FIELD_DURATION_NS,
        gas_temperature_k=GAS_TEMPERATURE_K,
        n2_fraction=N2_FRACTION,
        bolos_update_ns=BOLOS_UPDATE_NS,
        bolos_cells=BOLOS_CELLS,
        bolos_rtol=BOLOS_RTOL,
        points=POINTS,
    )
    sim_result = simulate(**_settings)
    return (sim_result,)


@app.cell
def _(SPECIES, mo, plot_result, plot_species_balance, sim_result):
    r"""marimoセルを実行し、宣言した値を後続セルへ渡す。

    Parameters
    ----------
    SPECIES : sequence of str or None
        描画する粒子種。Noneの場合はメタデータから自動選択する。
    mo : module
        設定セルから渡されるmarimoの画面表示用モジュール。
    plot_result : callable
        共通の密度描画関数。
    plot_species_balance : callable
        共通の粒子種別生成消失の描画関数。
    sim_result : Result
        後続セルの表示と保存に使用する計算結果。
    """
    fig = plot_result(sim_result, species=SPECIES)
    balance_fig = plot_species_balance(sim_result, "N2*(C)")
    mo.vstack([fig, balance_fig])
    return


@app.cell
def _(FIELD_DURATION_NS, mo, np, sim_result, species_balance):
    r"""marimoセルを実行し、宣言した値を後続セルへ渡す。

    Parameters
    ----------
    FIELD_DURATION_NS : float or None
        電界の印加時間 [ns]。Noneの場合は一定電界とする。
    mo : module
        設定セルから渡されるmarimoの画面表示用モジュール。
    np : module
        設定セルから渡されるNumPyモジュール。
    sim_result : Result
        後続セルの表示と保存に使用する計算結果。
    species_balance : callable
        係数履歴から生成速度と消失速度を再構成する。
    """
    _diagnostics = sim_result.diagnostics
    _balance = species_balance(sim_result, "N2*(C)")
    _outside = (
        _diagnostics.get("max_unsupported_neutral_fraction", 0) > 0.01
        or _diagnostics["max_ion_fraction"] > 0.01
    )
    _message = (
        f'最大電離率：{_diagnostics["max_ion_fraction"]:.3%}。'
        "電子衝突断面積のない中性種の最大割合："
        f'{_diagnostics.get("max_unsupported_neutral_fraction", 0):.3%}。\n\n'
    )
    if (
        FIELD_DURATION_NS is not None
        and 0 <= FIELD_DURATION_NS < sim_result.time[-1] * 1e9
    ):
        _off_index = int(
            np.searchsorted(sim_result.time, FIELD_DURATION_NS * 1e-9)
        )
        _index = sim_result.mechanism.index
        _n2_initial = sim_result.densities[_index["N2"], 0]
        _o2_initial = sim_result.densities[_index["O2"], 0]
        _n2_remaining = (
            sim_result.densities[_index["N2"], _off_index] / _n2_initial
            if _n2_initial > 0
            else float("nan")
        )
        _o2_remaining = (
            sim_result.densities[_index["O2"], _off_index] / _o2_initial
            if _o2_initial > 0
            else float("nan")
        )
        _message += (
            f"遮断直後の最初の出力点（{sim_result.time[_off_index] * 1e9:g}"
            " ns）："
            f"N₂残存率 {_n2_remaining:.3%}、O₂残存率 {_o2_remaining:.3%}、"
            "N₂(C)消失時定数"
            f' {_balance["loss_lifetime_s"][_off_index] * 1e9:.3g} ns。\n\n'
        )
    _message += (
        "BOLOSの電子衝突標的は基底状態N₂/O₂のみです。"
        "電子エネルギー方程式・電子間衝突・超弾性衝突・気体加熱は未実装です。"
        "1%は適用範囲を確認するための目安で、精度保証の閾値ではありません。"
        "消失時定数n/Lには生成による補給を含みません。"
    )
    mo.callout(
        mo.md(_message),
        kind="warn" if _outside else "info",
        title="適用範囲の確認が必要" if _outside else "モデルの仮定と診断",
    )
    return


@app.cell
def _(OUTPUT, SPECIES, save_result, sim_result):
    r"""marimoセルを実行し、宣言した値を後続セルへ渡す。

    Parameters
    ----------
    OUTPUT : pathlib.Path
        生成する計算結果の保存先。
    SPECIES : sequence of str or None
        描画する粒子種。Noneの場合はメタデータから自動選択する。
    save_result : callable
        共通の結果保存関数。
    sim_result : Result
        後続セルの表示と保存に使用する計算結果。
    """
    save_result(sim_result, OUTPUT, species=SPECIES)
    print(f"保存先：{OUTPUT.resolve()}")
    return


if __name__ == "__main__":
    app.run()
