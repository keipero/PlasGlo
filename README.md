<p align="center">
  <img src="img/plasglo_logo.png" alt="PlasGlo — Plasma Global Modeling Framework" width="640">
</p>

# PlasGlo

**Plasma Global Modeling Framework — Pythonによる0次元プラズマ反応計算**

PlasGloは、化学量論に基づく粒子密度の時間積分と、Py3bolosによる電子衝突反応係数の計算を組み合わせたプラズマ化学シミュレーションコードです。局所電界近似（LFA）を採用し、指定した換算電界 E/N の下で化学種の時間変化を求めます。

Pythonパッケージ名は `plasglo` です。純N₂・純Ar・N₂/O₂のケースを含み、共通APIで計算・可視化・保存・反応経路解析を行えます。

## 主な機能

- **入力による反応機構の定義**：ケースJSONと反応CSVで、化学種・元素組成・電荷・初期条件・第三体効率・電子衝突チャンネルを指定。
- **剛性反応系の時間積分**：共通の質量作用則カーネル、解析的Jacobian、BDF法を使用。
- **電子衝突計算との連成**：定常EEDFから反応係数と有効電子温度を求め、指定間隔で更新。電界の切替時刻で積分区間を分割。
- **診断と結果保存**：元素・電荷の保存、化学種の生成・消滅、密度・反応係数・電子群状態の履歴を確認。
- **反応経路解析**：化学種ごとの反応寄与、区間積算、生成・消滅の順位、グラフ用データを出力。
- **複数の実行入口**：Python API、CLI、Jupyter Notebook、marimo。

## セットアップ

以下のコマンドは、クローンしたリポジトリの直下（このREADMEがあるディレクトリ）で実行してください。Python 3.14.3での動作を確認しています。

### 1. Py3bolosを取得

電子衝突計算には、`Py3bolos/` に配置したサブモジュールを使用します。

```bash
git submodule update --init --recursive
```

ZIPで取得した場合など、Gitサブモジュールを利用できない場合は、[Py3bolos](https://github.com/keipero/Py3bolos)を同じ場所へクローンしてください。

```bash
git clone https://github.com/keipero/Py3bolos Py3bolos
```

### 2. Python環境を準備

```bash
python -m venv .venv
```

利用しているシェルに合わせて仮想環境を有効化します。

```powershell
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
```

```bash
# Windows Git Bash
source .venv/Scripts/activate

# Linux / macOS
source .venv/bin/activate
```

必要なパッケージを導入します。

```bash
python -m pip install -r requirements.txt
```

主な依存パッケージはNumPy、SciPy、Matplotlibです。ノートブックを利用する場合は、追加で実行環境を用意します。

```bash
python -m pip install jupyterlab ipykernel
```

現在はソースツリーから実行する構成です。Python APIの例もリポジトリ直下で実行してください。通常の計算にZDPlasKin本体は必要ありません。

## クイックスタート

### Python API

```python
from plasglo import load_case, simulate, save_result

case = load_case("nitrogen")  # 標準ケース名、またはケースJSONのパス
case = case.with_updates(electrons={"cells": 1600})

result = simulate(case)
save_result(result, "outputs/my_nitrogen")
```

`load_case()` → `simulate()` → `save_result()` が基本の流れです。入力の定義方法は[ケース入力ガイド](docs/GAS_CASES.md)を参照してください。新しいガスを計算するには、対象の反応機構と対応する断面積データを用意する必要があります。

### CLI

```bash
# 窒素ケースを計算
python -m plasglo --case nitrogen --bolos-cells 1600

# 自分のケースJSONを指定
python -m plasglo --case path/to/my_case.json

# 窒素・アルゴンのサンプルを実行
python samples/gases/run_samples.py --sample all

# 空気のKomuro機構：44種・134反応
python -m plasglo --field-td 100 --end-ns 100

# 利用できるオプションを表示
python -m plasglo --help
```

空気の矩形パルスの例です。

```bash
python -m plasglo --field-td 1000 --field-duration-ns 0.05 --end-ns 3 --bolos-update-ns 0.01 --output outputs/air/pulse
```

### Jupyter Notebook

[notebooks/run_nitrogen.ipynb](notebooks/run_nitrogen.ipynb)を開き、準備したPython環境のカーネルを選んで上から実行してください。

既定条件は純N₂、21種・96反応、1000 Td、300 K、101325 Pa、初期電子・N₂⁺密度各10¹⁴ m⁻³、計算時間1 ns、1600エネルギーセルです。条件はノートブック内で変更できます。

密度・反応係数・診断・グラフと、Nの反応経路解析結果を `outputs/notebooks/nitrogen/` に保存します。

## 反応経路解析

計算後の `Result` から、対象化学種の生成・消滅に寄与する反応を調べられます。

```python
from plasglo import analyze_pathways, save_pathways

paths = analyze_pathways(
    result,
    species="N",
    time_range_s=(0.0, 1e-9),
)

print(paths.production[:5])  # 積算生成寄与の上位5反応
print(paths.loss[:5])        # 積算消滅寄与の上位5反応
save_pathways(paths, "outputs/my_nitrogen/pathways/N")
```

反応速度・寄与の配列（NPZ）、順位表（CSV）、化学種と反応を結ぶグラフデータ（JSON）を保存します。ノートブックでは寄与の時間変化と積算寄与の図も表示できます。

積算は保存密度に基づく後処理です。短い反応ピークを解析する場合は、出力間隔を細かくして積算量の収束を確認してください。グラフの辺は反応への参加量を表し、原子ごとの移動経路を特定するものではありません。詳細は[反応経路解析ガイド](docs/REACTION_PATHWAYS.md)を参照してください。

## 入力と出力

| 項目 | 場所・内容 |
|---|---|
| ケース定義 | `data/cases/*.json`：化学種、条件、第三体、電子衝突の対応 |
| 反応機構 | `data/mechanisms/*.csv`：反応式、係数、単位、出典 |
| 電子衝突データ | `data/cross_sections/`：ケースが参照する断面積データ |
| 計算結果 | 密度CSV・NPZ、電界・係数履歴、診断JSON、使用反応一覧、図 |
| 反応経路解析 | 反応速度NPZ、順位CSV、グラフJSON |

既定の出力先は `outputs/` 配下です。CLIの `--output`、APIの `save_result()`、ノートブックの `OUTPUT_DIR` で変更できます。

## モデルの適用範囲

現在のモデルは、**閉鎖・定容・気体温度一定の0次元体積反応モデル**です。

- 電子エネルギー方程式は独立に積分せず、指定した E/N に対する定常EEDFを用います。電界遮断後の電子の過渡冷却は独立に解きません。
- 電子間衝突・超弾性衝突、壁面損失、空間輸送、流入出、気体加熱、外部回路は含みません。
- EEDFに寄与する衝突標的はケース入力と断面積データで定義した範囲に限られ、追跡するすべての生成種が自動的に寄与するわけではありません。
- 窒素・アルゴンサンプルと空気用Komuro機構は、それぞれ異なる反応機構です。

保存則の一致や参照コードとの一致は、実放電に対する物理モデルの妥当性を保証しません。特に高電界・強電離・強解離への適用では、反応機構・断面積・LFAの適用範囲と数値収束を確認してください。

背景は[モデル方針](docs/modeling_decisions.md)、入力と出典は[サンプルの説明](samples/gases/README.md)に記載しています。

## 検証と開発

開発用パッケージを導入し、自動テストを実行できます。

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests -q
```

改名時点で60件の自動テストと10 subtestsが成功しています。反応入力、元素・電荷の保存、第三体効率、Jacobian、混合気体の電子計算、電界切替、生成・消滅診断、反応経路解析を確認しています。

管理対象のPythonコードとノートブックには、PEP 8および日本語のNumPy形式docstringの検査を用意しています。

```bash
python verification/check_style.py
python verification/check_docstrings.py --examples
```

窒素21種・96反応について、元入力から再ビルドしたZDPlasKinと比較しています。1600セル・1 ns・10001時刻で、各化学種の参照ピーク密度で規格化した最大差は、電子間衝突なしの対照に対して約0.32649%でした。電子間衝突が自動的に有効になる元設定との差は約11.875%であり、条件を分けて評価しています。

比較条件・検証履歴・再実行手順は[検証資料](verification/README.md)を参照してください。外部ZDPlasKinとの比較には、別途その実行環境と参照入力が必要です。

## リポジトリ構成

```text
plasglo-master/
├── plasglo/           # 計算パッケージ
│   ├── chemistry/     # 反応式、化学量論、速度係数、反応カーネル
│   ├── electrons/     # 定常EEDF、Py3bolos接続、係数テーブル
│   ├── solvers/       # BDF積分、LFA連成
│   ├── diagnostics/   # 保存則、生成・消滅、反応経路解析
│   ├── io/            # 入力検証、読み込み、結果保存
│   ├── visualization/ # 描画
│   └── validation/    # 参照コードとの比較
├── data/              # ケース、反応機構、断面積
├── notebooks/         # 実行ノートブック
├── examples/          # marimoアプリ
├── samples/gases/     # 窒素・アルゴンの入力と実行例
├── tests/             # 自動テスト
├── verification/      # 検証スクリプトと資料
├── docs/              # 入力仕様、構成、モデル方針
├── img/               # ロゴなどの画像
├── archive/           # 通常実行から切り離した過去の試作
├── outputs/           # 実行時に生成する結果（Git管理対象外）
└── Py3bolos/          # 電子衝突計算の依存サブモジュール
```

実装の責務と読み進め方は[アーキテクチャ](docs/ARCHITECTURE.md)を参照してください。旧版の `py3zdplaschem` を利用していたコードでは、importを `plasglo` に変更してください。
