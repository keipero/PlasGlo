# BOLOSの接続状況

2026-10-31、旧 `bolos-python3/` の削除後、利用者がクローンした [keipero/Py3bolos](https://github.com/keipero/Py3bolos) への接続を完了しました。

- 接続先：`Py3ZDPlasChem-master/Py3bolos/bolos/`
- 確認コミット：`cedc937479b17f4d51ce87974b6787d40149d2e3`
- `bolos_coupling.py` はこのクローンを直接読み、別の同名パッケージが読み込み済みならカーネル再起動を求めます。
- N₂/O₂の反応番号対応と三体付着C27の密度規格化・SI変換は連成側に保持しています。
- 空気反応リストは `data/mechanisms/air_chemistry_Komuro_reactions.txt` に原本と同一内容で配置しました。外部afivoフォルダへの依存はありません。
- 窒素・アルゴンのZDPlasKin入力サンプルは `samples/gases/` にあります。
- クローン本体のソース変更や共有Python環境へのBOLOSのインストールは行っていません。

接続・保存則・電界切替などの自動テスト7件と、空気パルス計算、窒素・アルゴンの対照計算を実行しました。[テスト結果と制約](../verification/README.md)を参照してください。

過去のローカル実行結果は利用者の指示で削除しました。今回の削除記録は `verification/previous_results_cleanup.json`、旧BOLOS本体の削除記録は `bolos_removal_manifest.json` です。
