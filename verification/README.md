# Py3bolos再接続と検証（2026-10-31）

`Py3bolos/bolos/` のクローンへ接続しました。クローンの追跡対象ファイルには変更を加えていません。確認コミットは `cedc937479b17f4d51ce87974b6787d40149d2e3`、実行環境はPython 3.14.3、NumPy 2.4.4、SciPy 1.17.1です。

## 接続・実行テスト

- 当時の自動テスト7件が `env_Py3ZDPlasChem` で成功しました。電流試作の10件はコードとともにarchiveへ退避し、通常テストから取り外しています。
- [窒素ノートブック](../notebooks/run_nitrogen.ipynb) は電流セルを除いた全18セル（コード8セル）を現在の仮想環境の新しいカーネルで実行しました。出力は `../outputs/notebooks/nitrogen/` です。条件は1000 Td・300 K・1 ns・1600セル、元素保存の最大相対誤差は約1.33e-15でした。検証記録は `../outputs/verification/nitrogen_notebook_validation.json` です。
- 読み込むBOLOSがクローン内のものか、三体付着C27のSI変換と消失頻度が一致するかを確認しました。
- 1/100/2000 Td、純N₂・純O₂、出力格子と異なる時刻での電界遮断を確認しました。
- 既存の原子状酸素拡張の保存則・質量作用則・指数減衰・時刻別源項テストも成功しました。
- 空気モデルは1000 Tdを0.05 ns印加、3 nsまで、更新間隔0.01 ns・400セル・136反応でCLI実行しました。301回のBOLOS計算を含め完了し、元素保存の最大相対誤差は 4.44e-15 でした。出力は `../outputs/air/atomic_oxygen/` です。

## ZDPlasKinと同じサンプル入力での比較

各入力は `../samples/gases/` に配置し、コピー元とのSHA-256一致を確認しました。過去の比較出力は再利用せず、外部原本のコピーから対照計算を再ビルドしました。

| 比較 | 全粒子種・全出力時刻の最大差 |
|---|---:|
| 窒素21種・96反応、BOLOSとZDPlasKinの電子間衝突なし対照 | 0.186347% |
| アルゴン3種・2反応、BOLOS係数とZDPlasKin | 3.20542% |
| アルゴンで反応係数をZDPlasKin値に揃えた化学積分のみ | 0.00667141% |

差は粒子種ごとに `max(abs(Python - ZDPlasKin)) / max(abs(ZDPlasKin))` で定義します。微量密度での瞬間相対誤差ではありません。アルゴンの電離係数はBOLOS 7.30450454e-18 m³/s、ZDPlasKin 7.27249949e-18 m³/sです。約0.44%の係数差が密度増幅に影響しています。物理モデルの同一性・全条件での精度を保証する比較ではありません。

400/800/1600セルでも断面積から再計算しました。400→1600セルの電離係数変化は窒素約−0.02195%、アルゴン約−0.23831%です。アルゴンの400セル結果は格子依存が残っています。詳細は `../outputs/verification/grid_resolution.json`。BOLOSはEEDF末尾の傾きや有効断面積の差し引きに関する警告を出すことがあり、警告を理由なく無効化していません。アルゴンの丸め誤差程度の負の末尾には既存の対数傾き補正を共用し、補正前の最小値と補正の有無を診断へ保存します。有意な負値は引き続きエラーです。

## 再実行

追加のアルゴン格子試験では6400セルまで評価しました。現行設定で3200→6400セルの電子密度軌道の変化は約0.060%、上限100 eV固定・rtol=1e-8では約0.026%でした。[格子数・上限・反復条件を分けた結果](ARGON_GRID_STUDY.md)を参照してください。

ワークスペースのルートから実行します。

```powershell
python -m unittest discover -s Py3ZDPlasChem-master -p "test_*.py" -v
python Py3ZDPlasChem-master/samples/gases/run_samples.py --sample all
python Py3ZDPlasChem-master/verification/verify_reconnection.py
python Py3ZDPlasChem-master/verification/check_resolution.py
```

最初の2つはPythonだけで動作します。3つ目のZDPlasKin比較は外部 `ZDplasKin_analysis` 原本とローカルgfortran/BOLSIG DLLを使用します。使用したファイルは変更していません。窒素の元実行は電離率により電子間衝突が有効になるため、そのままの軌道をBOLOSとの一致基準にはしていません。

## 過去の結果の削除

利用者の指示により、既存の空気計算、旧ノートブック出力、アフターグロー解析、旧ZDPlasKin比較、pygmol側の過去のローカル数値出力・ログ等、計134ファイルを削除しました。入力・ソース・文献・検討資料・上流のテスト用参照値とクローン同梱の資料は保持しています。削除記録は `previous_results_cleanup.json`。現在ある出力は再接続後に作成した新しいものです。

## 2026-11-01：電流試作の退避とモジュール化

現在の実行入口は `notebooks/`・`examples/`・`samples/`、計算本体は `py3zdplaschem/`、全ての既定出力は `outputs/` に分離しました。[構成と拡張方針](../docs/ARCHITECTURE.md)を参照してください。

`verify_modularization.py` で窒素・アルゴン・空気パルスを分割前の同一環境の密度履歴と比較し、全ケースの最大差が0であることを確認しました。元入力と退避した試作コードもSHA-256で一致しています。記録は `../outputs/verification/modularization/report.json` と `moves.json` です。

ノートブックの再検証にはnbformat・nbclientが使える検証用Pythonと、ipykernelが入った計算用Pythonを使用します。

```powershell
python Py3ZDPlasChem-master/verification/validate_nitrogen_notebook.py --python env_Py3ZDPlasChem/Scripts/python.exe
env_Py3ZDPlasChem/Scripts/python.exe Py3ZDPlasChem-master/verification/verify_modularization.py
```

電流計算の元コードは `../archive/current_diagnostics_2026-11-01/`、当時の出力は `../outputs/archive/current_diagnostics_2026-11-01/` に保持しています。

## 2026-11-02：通常計算の名称整理

通常計算のモジュールを `chemistry/sample_mechanism.py`・`electrons/sample_swarm.py`、比較入力の読み込みを `io/reference_data.py` に変更しました。入力は `samples/gases/` に移動し、ノートブック・CLI・検証コードの参照を更新しています。外部ZDPlasKinとの比較処理と出典情報にはZDPlasKinの名称を保持しています。[名称の対応と読み順](../docs/ARCHITECTURE.md#窒素アルゴン計算の読み方)を参照してください。

自動テスト7件が成功しました。窒素ノートブックのコード8セルは計算用仮想環境で順次実行し、密度・時刻・粒子名が既存基準と完全一致しました。今回はJupyterカーネル起動が進まなかったため停止し、描画にはAggを使用してPythonとして実行しています。インライン表示の再検証は行っていません。入力・出典記録11ファイルの移動前後のハッシュ一致、両気体の出典ハッシュ、比較モジュールのimport、CLIのヘルプも確認しました。比較CLIの既存の `argparse` import漏れも修正しました。

空気・アルゴンを再計算した `verify_modularization.py` も成功しています。窒素の順次実行の記録は `../outputs/verification/sample_naming_validation.json`、既存基準との比較は `../outputs/verification/modularization/report.json` にあります。

## 2026-11-02：任意の反応追加機能の削除

APIの `chemistry_extension` と `read_reactions()` の `extension`、CLIの `--chemistry-extension`、原子状酸素の追加2反応を削除しました。marimoも基本134反応を一度だけ計算する構成に変更し、比較表示・設定・古いセッション出力を取り除きました。拡張専用テストを削除し、係数履歴から生成消失を再構成するテストは `tests/test_balance.py` へ移しています。

`python -m pytest Py3ZDPlasChem-master/tests -q` は35件と5 subtestsが成功しました。空気の基準密度履歴は、削除前のコードで拡張を使わずに再計算してから更新し、削除後との完全一致を確認しました。変更前の拡張付き基準・manifest・追加反応の入力は `outputs/verification/modularization/baseline/with_atomic_oxygen/` に過去記録として残しています。窒素ノートブックのコード8セルも計算用仮想環境で順次実行し、既存の密度・時刻・粒子名と完全一致しました。

空気CLIは1000 Td・0.05 nsの矩形パルス、終了0.1 ns、結合間隔0.03 nsで実行・保存を確認しました。結果は `../outputs/verification/base_chemistry/`、API引数の削除と窒素の確認記録は同フォルダの `removal_validation.json` です。冒頭の136反応の結果と拡張テストは削除前の記録です。

## 2026-11-02：ガス種に依存しない共通経路

入力をケースJSONと反応CSVに整理しました。元素組成・電荷・初期状態・第三体効率・断面積の標的とチャンネルは入力で定義し、全ケースで同じ反応カーネル、BDF/LFA積分、結果形式、保存・描画を使います。新しい反応データを生成する機能や、新しい物理項は追加していません。[入力形式](../docs/GAS_CASES.md)を参照してください。

### 主検証：ZDPlasKin窒素21種・96反応

利用者指定の窒素データで、1000 Td・300 K・101325 Pa・1 ns、初期電子・N₂⁺各1e14 m⁻³に統一しました。BOLOSは1600セル、比較は全21種・10001時刻です。元入力のコピーからZDPlasKinを再ビルドし、元設定と電子間衝突なしの対照をそれぞれ実行しています。

| 比較 | 各粒子種の参照ピーク密度で規格化した最大差 |
|---|---:|
| 共通Python経路 対 ZDPlasKin電子間衝突なし | 0.326490913% |
| 共通Python経路 対 元のZDPlasKin結果 | 11.875151264% |
| 元設定を再ビルド・再実行 対 保存済み元結果 | 0% |

元設定では電子間衝突が自動的に有効になるため、モデル条件を合わせた対照と区別しています。係数・密度の調整、時刻シフトは行っていません。元入力・元出力のハッシュは比較前後で一致しました。元素保存誤差は約1.33e-15、正味電荷誤差は約1.69e9 m⁻³です。高電離・強解離でのモデル範囲、他ガス種の物理的妥当性は別に評価する必要があります。

改修直前に同じ入力・400セルで保存した基準に対して、共通経路の窒素密度履歴の差は0でした。古いノートブック由来の基準は初期電子密度が1e8 m⁻³で今回の元入力と異なるため、今回の窒素回帰比較には使用していません。

結果は [比較レポート](../outputs/verification/gas_generalization/nitrogen/README.md)、[詳細JSON](../outputs/verification/gas_generalization/nitrogen/report.json)、[図](../outputs/verification/gas_generalization/nitrogen/comparison.png) に保存しています。

```powershell
env_Py3ZDPlasChem/Scripts/python.exe Py3ZDPlasChem-master/verification/validate_nitrogen.py
```

### 実装の確認

`python -m pytest Py3ZDPlasChem-master/tests -q` は50件と5 subtestsが成功しました。H/He/C/Oの人工的な検証用反応で解析解・保存則・結果再読込、反応ごとの第三体効率のJacobianを有限差分で確認しました。同位体の粒子名、多価イオン、電子の別名、入力不整合の拒否、Ar/N₂混合標的のBOLOS計算、チャンネル定義順に依存しない対応付けも確認しています。これらの人工的反応はテスト専用で、物理モデルとして配布していません。

純Arと空気の既存基準との差は、それぞれ約1.01e-15、3.27e-11以下（各粒子のピーク密度で規格化）でした。更新した窒素ノートブックは6コードセルを計算用Pythonで順次実行し、上記1600セルの結果と完全一致しました。Jupyterの新しいカーネルやインライン表示の再検証は今回行っていません。CLI、旧固定テーブル、旧係数固定機構、保存ケースの再読込も確認しています。記録は `outputs/verification/gas_generalization/entrypoints.json` と `outputs/verification/modularization/report.json` です。

## 2026-11-02：PEP 8と日本語NumPy形式docstring

管理対象のPythonファイル55個とノートブック2個（15コードセル）を整形しました。外部のPy3bolos・pygmol、仮想環境、生成済み出力内のコピーは対象外です。コード行は79文字、説明は72文字で検査し、PEP 8が推奨する二項演算子の前での改行と矛盾するW503だけを除外しています。

明示的に定義した全172関数・メソッド・内部関数・テスト関数に、日本語のNumPy形式docstringがあります。標準見出しと型名・識別子は英語のままにし、引数・戻り値・単位・必要な例外を記載しました。83関数にExamplesを設け、パッケージ内の実行可能な使用例226件がすべて成功しました。外部環境が必要な旧MATLABコードでは、その前提と実行対象外の使用例を明記しています。

現在の数値テスト50件と5 subtests、保存用電流コードの既存テスト10件が成功しました。反応カーネル・機構・係数評価の5ファイルではdocstringを除いた構文木が変更前と一致しています。窒素の密度履歴と元のZDPlasKin入力ハッシュ、ノートブック6コードセルの結果も、整形前の検証済み結果と照合しています。

保存用ソースも整形したため、変更前のソースとマニフェストを `outputs/verification/pep8/original_sources.zip` に保存しました。`archive/.../source_manifest.json` の `sha256` は整形後、`pre_documentation_sha256` は変更前の値です。原本の照合を失わないよう、ZIPと対応するハッシュ一覧を残しています。

再検査には `check_style.py` と `check_docstrings.py --examples` を使います。記録は `outputs/verification/pep8/` の `migration.json`・`docstrings.json`・`numerical_verification.json` にあります。
