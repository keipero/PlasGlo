# 共通の0次元反応計算の構成

公開する基本の使い方は `load_case → simulate → save_result` です。ガス種ごとのソルバーは作らず、入力を同じ反応ネットワークに変換します。

| 場所 | 責務 |
|---|---|
| `data/cases/*.json` | 条件・粒子種・元素組成・電荷・第三体・電子衝突の対応 |
| `data/mechanisms/*.csv` | 反応式・係数・単位・出典 |
| `plasglo/io/cases.py` | JSON／CSVの読み込み、パス解決、入力検証 |
| `plasglo/chemistry/network.py` | 元素・電荷の保存確認、反応係数・初期状態 |
| `plasglo/chemistry/kernel.py` | 共通の質量作用則・密度方程式・解析的Jacobian |
| `plasglo/chemistry/rates.py` | 解析的速度係数の式の評価 |
| `plasglo/electrons/mixture.py` | 任意の宣言済み標的の定常BOLOS計算 |
| `plasglo/solvers/simulation.py`、`frozen.py` | 係数固定／逐次LFA結合のBDF積分と共通 `Result` |
| `plasglo/diagnostics/balance.py` | 係数履歴から任意の種の生成・消失を再構成 |
| `plasglo/io/results.py`、`visualization/plots.py` | ガス名に依存しない保存・描画 |

## コードを読む順序

1. `notebooks/run_nitrogen.ipynb`：条件の選択・実行・確認・描画・保存。
2. `data/cases/nitrogen.json` と参照先：入力で何を指定するか。
3. `io/cases.py`：入力から反応ネットワークを作る処理。
4. `chemistry/network.py` → `kernel.py`：係数、化学量論、密度方程式、Jacobian。
5. `electrons/mixture.py` → `solvers/simulation.py`：電子計算と時間積分の接続。

## 数値計算の範囲

閉鎖・定容・気体温度一定、指定E/Nのモデルです。電子係数は初期状態で固定するか、受理したBDF区間の終端ごとに更新します。電界切替では区間を分割し、出力点数と結合間隔は分離します。第三体の対象種と効率は入力で指定し、`M` を特定の気体へ暗黙に結び付けません。

電子エネルギーODE、電子間衝突、超弾性衝突、壁面損失、流れ、気体加熱、外部回路は今回追加していません。新しいガスの反応データや断面積も自動生成しません。

## 既存コードと検証

`Mechanism`、`SampleMechanism` は従来入力向けの互換アダプターで、密度方程式・Jacobianは共通 `MassActionKernel` を使います。`BolosSwarm` と `solve_sample_swarm()` も共通の電子ソルバーを呼びます。既存の数値引数による `simulate(...)` は空気プリセットへの入口として残し、独立した時間積分は持ちません。新しいガス・新しい機構にはケースJSONを使ってください。

比較用の `validation/zdplaskin.py` と元のFortran入力は通常計算から分離しています。`verification/validate_nitrogen.py` が共通経路の窒素を、同じ元入力の新しいZDPlasKin実行と比較します。

入力形式は [GAS_CASES.md](GAS_CASES.md)、現在の検証結果は [verification/README.md](../verification/README.md) にあります。過去の再編・機能削除の記録は同検証資料と `docs/history/`、`archive/` に保持しています。
