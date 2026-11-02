# Py3ZDPlasChem：0次元プラズマ反応計算

開発用のスタイル・文書検査には `requirements-dev.txt` を使用します。管理対象のPythonコードとノートブックをPEP 8（コード79文字、説明72文字）で検査し、全関数の日本語NumPy形式docstringと引数・戻り値・Examplesを確認できます。

```powershell
env_Py3ZDPlasChem/Scripts/python.exe -m pip install -r Py3ZDPlasChem-master/requirements-dev.txt
env_Py3ZDPlasChem/Scripts/python.exe Py3ZDPlasChem-master/verification/check_style.py
env_Py3ZDPlasChem/Scripts/python.exe Py3ZDPlasChem-master/verification/check_docstrings.py --examples
```

化学量論に基づく粒子密度の時間積分と、Py3bolosによる電子衝突係数の計算を組み合わせるPythonコードです。現在は局所電界近似（LFA）を採用しています。ガス種・元素・電荷・反応・第三体・断面積の対応は入力で定義でき、純N₂・純Ar・N₂/O₂を共通の計算経路で実行します。

## 共通API

```python
from py3zdplaschem import load_case, simulate, save_result

case = load_case("nitrogen")  # またはケースJSONのパス
case = case.with_updates(electrons={"cells": 1600})
result = simulate(case)
save_result(result, "outputs/my_nitrogen")
```

入力設定は「ケースJSON＋反応CSV」の2ファイルです。`data/cases/` のJSONに粒子種・条件・第三体・電子衝突の対応をまとめ、`data/mechanisms/*.csv` に反応一覧を置いています。別途、対象ガスの断面積ファイルが必要です。新しい実行時依存パッケージは追加していません。[入力形式と追加手順](docs/GAS_CASES.md)を参照してください。

## 反応経路解析

計算後の `Result` から、全反応の速度、対象化学種への正味寄与、指定区間での積算寄与と順位、化学種と反応を結ぶグラフ用データを取得できます。

```python
from py3zdplaschem import analyze_pathways, save_pathways

paths = analyze_pathways(result, species="N", time_range_s=(0.0, 1e-9))
print(paths.production[:5])  # 生成側の積算寄与の上位5反応
print(paths.loss[:5])        # 消滅側の積算寄与の上位5反応
save_pathways(paths, "outputs/my_nitrogen/pathways/N")
```

窒素ノートブックにも解析・描画・保存を含めています。対象種は `TARGET_SPECIES`、区間は秒単位の `TIME_RANGE_S` で変更できます。積算は保存密度を用いる後処理なので、短い反応ピークには出力間隔を細かくして収束を確認してください。[解析方法と出力形式](docs/REACTION_PATHWAYS.md)を参照してください。

2026-11-01、機能ごとのモジュール分割と出力先の統一を行いました。電流計算の試作は通常の計算・ノートブックから取り外し、コードを [archive/](archive/current_diagnostics_2026-11-01/README.md) に記録として保存しています。

## 基本実行

窒素の基本ノートブックは [notebooks/run_nitrogen.ipynb](notebooks/run_nitrogen.ipynb) です。Pythonカーネルを選び、上から全セルを実行してください。既定は純N₂、21種・96反応、1000 Td、300 K、101325 Pa、初期電子・N₂⁺各10¹⁴ m⁻³、1 ns、1600エネルギーセルです。条件セルで変更できます。

結果は `outputs/notebooks/nitrogen/` に保存します。密度・反応係数・診断JSON・グラフに加え、`pathways/N/` に反応経路の解析配列・順位CSV・グラフJSON・寄与の図を保存します。対象種を変更した場合は、その名前のディレクトリになります。電流は計算しません。

計算本体は `py3zdplaschem/chemistry/network.py` と `kernel.py`、電子計算は `electrons/mixture.py`、時間積分は `solvers/simulation.py` です。窒素の元入力は `samples/gases/nitrogen/` に保持しています。通常実行にZDPlasKin本体は必要ありません。[コードを読む順序](docs/ARCHITECTURE.md)も参照してください。

必要なパッケージの導入例（ワークスペースのルートで実行）：

```powershell
.\env_Py3ZDPlasChem\Scripts\python.exe -m pip install -r .\Py3ZDPlasChem-master\requirements.txt ipykernel
```

Py3bolosは `Py3bolos/` にあるローカルクローンを読みます。ソース本体には変更を加えていません。tyの探索パスはワークスペース直下の `ty.toml` に設定しています。

## フォルダ構成

```text
Py3ZDPlasChem-master/
├── py3zdplaschem/      # 機能別の計算モジュール
│   ├── chemistry/     # 反応式、化学量論、速度係数、反応カーネル
│   ├── electrons/     # BOLOS接続、定常EEDF、固定swarmテーブル
│   ├── solvers/       # BDF積分と逐次LFA連成
│   ├── diagnostics/   # 生成・消失、履歴からの診断
│   ├── visualization/ # グラフ
│   ├── io/            # 入力照合、読み込み、結果保存
│   └── validation/    # ZDPlasKinとの比較
├── data/              # ケースJSON・反応CSV・粒子種・電子衝突設定・断面積
├── notebooks/         # 基本実行ノートブック
├── examples/          # marimo入口 ZDchem.py
├── samples/gases/     # 窒素・アルゴンの元入力と実行コード
├── tests/             # 有効な自動テスト
├── verification/      # 検証用スクリプトと説明
├── outputs/           # 計算・検証・保存済み試作結果の共通親フォルダ
├── docs/              # モデル方針、構成説明、調査資料
├── archive/           # 実行から切り離した試作コード
└── Py3bolos/          # 利用者がクローンした依存リポジトリ
```

責務と今後の拡張先は [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) を参照してください。`kinetics.py`、`bolos_coupling.py`、`compare_zdplaskin.py` は以前の呼び出しとの互換用の薄い入口として残し、数値計算の本体はパッケージへ移しました。

## CLIと出力先

ワークスペースのルートから実行します。

```powershell
# ZDPlasKinと同じ窒素・アルゴンの入力
python Py3ZDPlasChem-master/samples/gases/run_samples.py --sample all

# 任意のケースJSON、または標準ケース名を指定
python Py3ZDPlasChem-master/kinetics.py --case nitrogen --bolos-cells 1600
python Py3ZDPlasChem-master/kinetics.py --case path/to/my_case.json

# 空気のKomuro機構：44種・134反応
python Py3ZDPlasChem-master/kinetics.py --field-td 100 --end-ns 100

# 空気の基本機構で矩形パルスを計算
python Py3ZDPlasChem-master/kinetics.py --field-td 1000 --field-duration-ns 0.05 --end-ns 3 --bolos-update-ns 0.01 --output Py3ZDPlasChem-master/outputs/air/default
```

| 実行・資料 | 既定の保存先 |
|---|---|
| 窒素ノートブック | `outputs/notebooks/nitrogen/` |
| 窒素・アルゴンのサンプル | `outputs/samples/nitrogen/`、`outputs/samples/argon/` |
| 空気CLI | `outputs/air/default/` |
| marimo `examples/ZDchem.py` | `outputs/air/default/` |
| 接続・格子・数値検証 | `outputs/verification/` |
| 旧窒素比較の保存済みデータ | `outputs/validation/zdplaskin_nitrogen/` |
| 電流試作の保存済み結果 | `outputs/archive/current_diagnostics_2026-11-01/` |

`py3zdplaschem/paths.py` で共通パスを管理しています。CLIの `--output` やノートブックの `OUTPUT_DIR` で別の場所を指定することも可能です。既存出力は同じ内容のまま移動し、移動前後のハッシュを照合しました。

## モデルの範囲

- 現在の密度方程式は体積内の化学反応を解きます。壁面失活・輸送・流入出・気体加熱・外部回路は未導入です。
- 電子エネルギー方程式を独立に積分せず、定常EEDFから反応係数と有効電子温度を求めます。電子間衝突・超弾性衝突は含みません。
- 窒素・アルゴンサンプルと、空気用Komuro機構は別の反応機構です。[サンプルの条件と出典](samples/gases/README.md) を参照してください。
- 保存則の一致は、実放電に対する物理モデルの妥当性を保証しません。特に1000 Tdや強電離・強解離での解釈には注意が必要です。

[モデル方針](docs/modeling_decisions.md)、[再編前の詳細説明](docs/history/README_before_modularization.md) に背景を残しています。

## 検証

```powershell
env_Py3ZDPlasChem/Scripts/python.exe -m unittest discover -s Py3ZDPlasChem-master -p "test_*.py" -v
env_Py3ZDPlasChem/Scripts/python.exe Py3ZDPlasChem-master/verification/verify_modularization.py
```

60件の自動テストと10 subtestsで、反応入力・任意元素の保存則・第三体効率・Jacobian・混合気体の電子計算・電界切替・生成消失診断・反応経路解析を確認しました。

共通経路の窒素21種・96反応を、元入力から再ビルドしたZDPlasKinと比較しました。1600セル・1 ns・10001時刻で、各粒子種の参照ピーク密度で規格化した最大差は、電子間衝突なしの対照に対して約0.32649%です。電子間衝突が自動的に有効になる元設定との差は約11.875%で、別に評価しています。元設定の再実行は保存済み結果と一致しました。[比較結果](outputs/verification/gas_generalization/nitrogen/README.md)・[検証の説明](verification/README.md)を参照してください。
