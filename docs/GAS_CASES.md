# ガス種を入力で指定する

入力設定は **ケースJSON＋反応CSV** の2ファイルです。電子計算を行う場合は、別途そのガスの断面積を用意します。粒子種・第三体・電子衝突の対応をケースJSONにまとめ、ガス種ごとのPythonソルバーは作りません。JSONとCSVには標準ライブラリを使い、新しい依存パッケージは追加していません。

```python
from py3zdplaschem import load_case, simulate, save_result

case = load_case("nitrogen")  # または自分のケースJSONのパス
result = simulate(case)
save_result(result, "outputs/my_case")
```

標準ケースは `air`（44種・134反応）、`nitrogen`（21種・96反応）、`argon`（3種・2反応）です。新しいガスの反応データや断面積を自動生成する機能はありません。標準ケースのJSONをコピーし、その系の実データに置き換えてください。

## ケースJSON

実例は `data/cases/nitrogen.json` と `argon.json` です。`schema_version` は1です。

| 項目 | 内容 |
|---|---|
| `gas` | `temperature_K`、`pressure_Pa`、任意の `mole_fractions` |
| `initial_densities_m3` | 粒子名から初期密度への辞書。電子・イオンも指定 |
| `field` | E/Nの `td`、任意の `duration_ns`。遮断後は0 Td |
| `mechanism` | `schema_version: 1`、`species`、`reactions`、任意の `third_bodies`・`variables`・出典 |
| `electrons` | `backend`、`update_mode`、数値条件、`profile` |
| `solver` | `end_ns`、`points`、`rtol`、`atol_density_m3` |
| `diagnostics_species` | 生成消失履歴を保存する粒子名の配列。省略可能 |

相対パスはケースJSONから解決します。ただし、CSV内の任意列 `source_file` はそのCSVから解決します。保存された `case.json` は参照先を絶対パスに解決した実行条件です。入力本体はコピーせず、同じ入力ファイルを参照します。

密度はm⁻³、ケース内の時間はns、保存CSVの時間はs、温度はK、圧力はPa、平均電子エネルギーはeVです。圧力は初期規格化密度 `N0=p/(kB*Tg)` を決めます。定容計算中の圧力を一定に拘束するものではありません。

`gas.mole_fractions` を指定すると、中性種の初期値を `x*N0` とします。同じ種を `initial_densities_m3` で指定した場合は明示値が優先します。組成・圧力を変えても明示した初期密度は自動では変わりません。初期電荷は通常ゼロにそろえます。意図的な非中性条件だけは `allow_net_charge: true` で指定できます。

標準窒素の初期電子・N₂⁺密度は元のZDPlasKin入力に合わせて各1e14 m⁻³です。古いノートブックには別の初期密度の保存結果があるため、比較時は条件を照合してください。

```python
case = case.with_updates(
    field={"td": 500},
    solver={"end_ns": 0.5},
    electrons={"cells": 1600},
)
```

`gas`・`field`・`solver`・`electrons` は指定項目だけを上書きします。`initial_densities_m3` など、それ以外の項目は全体を置き換えます。旧空気APIの数値引数とケース指定を同時に混ぜることはできません。

## 粒子種・反応・第三体

`mechanism.species` は次のような配列です。元素組成と電荷を名前から推測しないため、任意元素・多価イオン・同位体の識別名を扱えます。

```json
[
  {"name": "e", "composition": {}, "charge": -1, "electron": true},
  {"name": "He", "composition": {"He": 1}, "charge": 0},
  {"name": "He+", "composition": {"He": 1}, "charge": 1}
]
```

電子は空の組成・電荷-1・`electron: true` とし、励起種には任意で `excited: true` を指定します。Heの例はメタデータの書式例であり、検証済みHe反応機構を同梱するという意味ではありません。

`mechanism.reactions` はCSVのパスです。必須列は `id,equation,kind,value,length_unit`、任意列は `parameters,source_file,line,reference,temperature_min_K,temperature_max_K` です。反応式は ` -> `、化学種の区切りは ` + `、係数は `2H` または `2 H` と書きます。重複ID・未定義の粒子・元素や電荷を保存しない反応は拒否します。

| kind | value | parameters |
|---|---|---|
| `electron` | 電子チャンネルID | 空欄 |
| `analytic` | `rates.py` にある式の種類 | 空白区切りの数値係数 |
| `expression` | 制限された算術式 | 空欄 |

式には `Tgas`／`Tg`、`Te` [K]、`EN`／`Td` [Td]、`exp,log,sqrt,min,max` を使えます。追加変数は `mechanism.variables` に依存順の `[名前, 式]` の配列で定義します。任意のPythonコードや `eval` は実行しません。Teを使う式には、電子計算または明示した平均電子エネルギーが必要です。

`length_unit` は `m` または `cm` です。次数qの係数は `cm` の場合 `1e-6**(q-1)` 倍してSIに変換します。電子係数は変換済みSIとして `m` を指定します。気体温度の適用範囲が分かる場合は温度列と出典を記入します。未記入は精度保証ではありません。

`mechanism.third_bodies` は、例えば `{"M":{"He":1.0,"N2":0.5}}` と定義します。`M` は追跡する粒子ではなく重み付き密度和です。反応式の両辺に同じ数を置きます。異なる効率のグループも別名で指定でき、1反応で使うグループは1つです。旧空気のN₂＋O₂、旧窒素の全中性種和は、それぞれ入力で維持しています。

## 電子計算

`electrons.profile` に `cross_sections`、`targets`、`channels` をまとめます。

- `cross_sections`：BOLSIG形式のファイル、または複数ファイルの配列。
- `targets`：断面積側の標的名から化学種名への対応。各標的に運動量移行データが必要。
- `channels`：任意の `id` と `target,kind,product,threshold_eV`。位置ではなく物理的な対応で一意に照合。

存在しない、曖昧な、重複する対応は拒否します。新しい反応機構ではC番号を使う必要はありません。標準空気だけは元入力との対応のためC番号を維持しています。

E/NのNは、宣言した衝突標的種の現在の密度の合計です。化学計算に含む全中性種が自動的に電子衝突標的になるわけではありません。断面積がない種について、別の気体のデータを代用しません。

密度規格化された特殊な断面積には、対応チャンネルの `density_normalization` で `species,power,density_unit,coefficient_scale` を指定します。補助標的を使ってEEDFの損失頻度と化学係数を整合させます。通常の二体断面積には指定せず、標準空気の三体付着設定を別のデータへ無条件に流用しないでください。

`update_mode: frozen` は初期係数を固定し、`coupled` は受理した区間終端で更新します。電界遮断には `coupled` を使います。出力点数を増やしても結合間隔 `update_ns` は変わりません。`backend: prescribed` は明示したSI係数 `coefficients_SI`、`backend: none` は断面積による電子反応がない機構に使います。Te依存の式には `mean_energy_eV` を明示してください。

## 結果と互換性

全ガスで同じ `Result` とCSV／NPZ／JSON／図を使い、元素数・標的数・チャンネル数に固定値はありません。生成消失履歴は `diagnostics_species` または `save_result(..., balance_species=[...])` で指定します。保存名は `balances.json` に記録します。

旧 `simulate(field_td=..., n2_fraction=...)` と固定swarmテーブルは空気用の互換入口です。旧 `mechanism=` は係数固定の互換用途で、独自機構のBOLOS連成にはケース入力を使います。元入力は参照用に保持し、通常の編集対象はケースJSONとCSVです。原本のハッシュが変われば変換済みデータとの不整合として検出します。

現在の主検証は [窒素とZDPlasKinの比較](../outputs/verification/gas_generalization/nitrogen/README.md) です。枠組みの汎用化と、各ガスの物理モデルの妥当性は別に確認します。壁・輸送・流入出・加熱・回路・電子間衝突・超弾性衝突・独立した電子エネルギー方程式は追加していません。
