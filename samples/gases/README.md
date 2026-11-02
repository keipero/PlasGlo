# 窒素・アルゴンの気体サンプル

Pythonの反応計算とPy3bolosを組み合わせた実行例です。通常の実行にZDPlasKin本体は必要ありません。

通常計算の設定は `data/cases/nitrogen.json`・`argon.json`、反応一覧は `data/mechanisms/nitrogen.csv`・`argon.csv` です（プロジェクトルートからのパス）。このフォルダの `kinet.inp`・断面積・Fortran・`case.json` は元入力の参照記録として保持しています。`run_samples.py` は共通の `load_case → simulate → save_result` を呼ぶ互換入口です。

外部の `ZDplasKin_analysis` 原本から、反応入力 `kinet.inp`、衝突断面積 `bolsigdb.dat`、元の実行設定 `user_code.F90` を改変せずコピーしています。コピー元とSHA-256は各フォルダの `provenance.json` に記録しました。過去の密度出力はここへ持ち込んでいません。`case.json` は元の実行設定をPython用に記録したものです。

| サンプル | 粒子種 | 反応数 | 衝突断面積数 | E/N | Tg | 計算時間 |
|---|---:|---:|---:|---:|---:|---:|
| nitrogen | 21 | 96 | 25 | 1000 Td | 300 K | 1 ns |
| argon | 3 | 2 | 3 | 51 Td | 301 K | 300 ns |

窒素は利用者の保存済みrevision 2機構で、ZDPlasKin公式の窒素機構全体ではありません。原本のフォルダ名とコメントには10 nsが残っていますが、実際の実行設定 `end_time_s=1e-9` と保存済みデータに合わせて1 nsとします。窒素の初期基底N₂密度は101325 Paから決まり、電子・N₂⁺密度は各10⁸ cm⁻³です。アルゴンは初期Ar密度2.5×10¹⁹ cm⁻³、電子・Ar⁺は各1 cm⁻³です。

ワークスペースのルートから実行します。

窒素をノートブックから実行する場合は、notebooks/内の [run_nitrogen.ipynb](../../notebooks/run_nitrogen.ipynb) を開いてください。条件をセル内で変更して、入力照合・計算・診断・描画・保存まで実行できます。

```powershell
python Py3ZDPlasChem-master/samples/gases/run_samples.py --sample all
```

このファイルのあるフォルダであれば `python run_samples.py --sample all` でも実行できます。NumPy・SciPy・Matplotlibと隣接するPy3bolosクローンを使用し、MATLABやFortranコンパイラは通常のPythonサンプル実行には不要です。

`--sample nitrogen` / `--sample argon` で個別に選択できます。既定出力はプロジェクト内の `outputs/samples/` で、密度CSV・NPZ・グラフ、SI反応係数、診断JSONを保存します。密度はm⁻³、時刻は秒で、原本の密度cm⁻³から単位を換算します。反応入力内の二体・三体係数も反応次数に従ってSIへ変換します。

電子エネルギー格子は既定400セルで、`--cells 800` や `--cells 1600` に変更できます。格子比較では、400→1600セルで電離係数が窒素約0.022%、アルゴン約0.238%変化しました。アルゴンの密度は電離による増幅に敏感なので、精度が必要な場合は格子依存も確認してください。

指定E/Nと一定Tgから定常EEDFを求めるLFAです。純ガスの基底状態だけを衝突標的とし、その標的組成は常に1なので、係数を一定として密度を積分します。`ANY_NEUTRAL` は現在の全中性粒子密度を用います。電子間衝突・超弾性衝突・気体加熱・電子エネルギーODEは含みません。

元の窒素実行では電離率に応じて電子間衝突が有効になります。その軌道との差をBOLOSの誤差とみなすことはできないため、検証では元ソースのコピーをビルドし、電子間衝突だけを無効にした対照計算と比較します。反応・断面積の出典と省略項は原本のコメントも参照してください。

アルゴンの `native_csv_driver.F90` は検証用にCSV精度を上げた実行コードです。元の実行コードは `user_code.F90` にそのまま保持しています。
