# 電流計算試作の保管記録

2026-11-01、利用者の指示により電流計算を通常コードと窒素ノートブックから取り外しました。このフォルダはコードの記録用です。通常のimport、ノートブック、自動テストの対象には含めません。

## 保管したもの

- `current_diagnostics.py`：電子・イオン・変位電流を扱った試作コードそのもの。
- `test_current_diagnostics.py`：対応する10件のテスト。
- `CURRENT_MODEL.md`：当時の計算式、入力値、制約の説明。
- `run_nitrogen_with_current.ipynb`：電流セルを含む実行済みノートブック。
- `validate_notebook_with_current.py`：当時の検証コード。
- `bolos_coupling_with_current_support.py`、`compare_zdplaskin_with_current_support.py`：換算移動度を追加した時点の連成側のソース。
- `source_manifest.json`：退避前のパスと保管内容のSHA-256。

各コードは元のバイト列を保っています。現在のディレクトリ構成へ合わせた修正は行っていないため、そのまま実行する入口ではありません。再導入するときの検討資料として使用してください。

当時の出力・検証結果は [outputs/archive/](../../outputs/archive/current_diagnostics_2026-11-01/) に保存しています。現在の窒素ノートブックは [notebooks/run_nitrogen.ipynb](../../notebooks/run_nitrogen.ipynb) で、電流計算を含みません。
