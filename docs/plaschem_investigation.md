# plaschem 採用調査

調査日: 2026-10-31（日本時間）

追加調査（同日）: 同じ著者による独立版 [hanicinecm/pygmol](https://github.com/hanicinecm/pygmol) を確認した。こちらには MIT ライセンス、パッケージ依存定義、方程式文書、テストがある。本メモのライセンス未確認・導入構成についての記述は、最初に指定された martin-hanicinec-ucl/plaschem リポジトリについてのもの。独立版の setup.py は SciPy 1.7.3 固定で、Python 3.14 への導入は未検証。参考候補の比較は reference_repositories.md に記載。

## 結論

現在の E/N 指定・BOLOS 連成・密閉定容系の計算には、chemistry_0d を継続することを推奨する。plaschem は反応機構管理と縮約・感度解析に価値があるが、付属のグローバルモデルへ全面移行しても、現在の目的に対する導入・検証作業が減るとは判断できない。

以下は公開 README・ソース・依存定義と手元のコードの静的調査による評価。plaschem のインストール、テスト実行、フォーク作成、性能比較は行っていない。手元の既存計算コードも変更していない。

## 1. リポジトリの構成と長所

著者の2021年の UCL 博士論文に付随する研究用コード。README は、一般利用向けの独立パッケージとしての完全な文書やサポートを意図していなかったことを明記している。

| 構成 | 役割 |
|---|---|
| plaschem | RP（粒子種）、Reaction（反応）、Chemistry（反応機構）のデータ構造 |
| pygmo_fwork.pygmol | SciPy solve_ivp を使う0次元プラズマモデル |
| pygmo_fwork.global_model / results | 計算のラッパー、ログ、後処理 |
| pygmo_fwork.optichem | グラフによる粒子種の順位付け、骨格機構への縮約、Morris 感度解析 |

反応機構管理には、元素・電荷保存の検査、化学量論行列、粒子・反応の整合的な無効化、XML 入出力、単位変換、LaTeX 表出力がある。反応・粒子オブジェクトには DOI や QDB ID 等の属性もある。データ属性があることは、外部データベースへの現在の接続やデータの正確性を保証するものではない。

中性種だけの系には NeutralChemistry / NeutralModel の別実装もある。

出典: [README](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/README.md)、[反応クラス](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/plaschem/reactions.py)、[Chemistry](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/plaschem/chemistry.py)、[optichem](https://github.com/martin-hanicinec-ucl/plaschem/tree/master/pygmo_fwork/optichem)、[中性種モデル](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/pygmo_fwork/pygmol/neutral_model.py)。

## 2. chemistry_0d と付属モデルの違い

| 項目 | 手元の chemistry_0d | plaschem 付属 pygmol |
|---|---|---|
| 駆動入力 | 指定 E/N、矩形パルスとアフターグロー | 時間依存の吸収電力 P(t) |
| 電子反応係数 | BOLOS の定常 EEDF を逐次更新して積分 | 主に Te に対する Arrhenius 型パラメータ |
| 電子エネルギー | 独立した時間発展方程式なし | 電子エネルギー密度の時間発展を解く |
| 電子密度 | 電子を独立した粒子種として積分 | 重粒子の電荷から準中性を仮定して算出 |
| 気体 | 密閉・定容、Tg 一定 | Tg 一定、流入出・圧力調整源項を含む |
| 壁面 | 現在なし | 拡散損失、付着、壁からの戻り、電子・イオンの壁へのエネルギー損失 |
| 数値積分 | BDF と解析的 Jacobian | 既定 BDF、解析的 Jacobian は未実装 |
| 反応機構 | N2/O2 用の既存反応データを解析 | XML または RP / Reaction から構築 |
| 解析機能 | 保存則、粒子別収支、BOLOS 履歴等 | 上記に加え機構縮約・感度解析の枠組み |

両方とも体積反応の核は質量作用則と化学量論行列による源項計算である。plaschem に変えれば化学量論ベースの方法から別の高精度な方法へ変わる、という関係ではない。差は電子系、境界条件、反応データ、解析・管理機能にある。

付属 pygmol に BOLOS / BOLSIG / LoKI-B を実行して EEDF を更新する仕組みは確認できなかった。Te 依存の係数近似だけで、現在の非 Maxwell 型 EEDF による係数をそのまま再現できるとは限らない。E/N 駆動を維持するには係数更新とエネルギー入力の扱いを設計し直す必要がある。

電子エネルギー方程式は今後の電子冷却モデルの参考になるが、移植だけで、欠けている断面積・緩和経路・気体加熱等が補われるわけではない。

出典: [equations.py](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/pygmo_fwork/pygmol/equations.py)、[model.py](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/pygmo_fwork/pygmol/model.py)、[モデルの数式説明](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/pygmo_fwork/pygmol/documentation/pygmol_doc_content.tex)。手元は README.md、kinetics.py、bolos_coupling.py を参照。

## 3. 採用前に把握すべき具体的な実装

- get_electron_density は n_e = sum(z_i n_i) を使う。空間的な電荷分離を解くモデルではない。
- get_flow_source_rates は feeds が空でも中性種への圧力回復項を加える。回復時間はコード中で 1e-3 s。流入をゼロにするだけでは、一般の密閉系と同じ方程式にはならない。
- get_min_n_correction は 1 m^-3 未満の重粒子密度を、回復時間 1e-10 s の人工源項で押し上げる。これは反応の化学量論から得られる源項ではないため、その部分は元素保存を一般には保証しない。実際の影響量は計算条件に依存し、今回測定していない。
- get_objective_functions(jacobian=True) は NotImplementedError を発生させる。
- Model._solve の solve_ivp 呼び出しには rtol / atol が指定されていない。現在の検証条件に合わせるには、これらを設定できるよう改修する必要がある。
- equations.py 冒頭の Arrhenius 式の説明には指数部が実装と逆の記述がある。実装と Reaction.get_rate_coefficient は exp(-C/T) を使う。コメントだけで数式を移植しないこと。
- 第三体 M は付属モデルでは重粒子総密度。手元の M=N2 基底状態+O2 基底状態とは異なるので、変換時に区別が必要。

これらはコード上の性質であり、著者の対象条件で結果が不正確だと結論するものではない。

出典: [equations.py](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/pygmo_fwork/pygmol/equations.py)、[model.py](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/pygmo_fwork/pygmol/model.py)、[reactions.py](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/plaschem/reactions.py)。

## 4. 導入・保守・ライセンスの確認状況

調査時の GitHub API メタデータでは、plaschem の pushed_at は 2021-09-06。リポジトリは archived ではないが、これを現在も活発に保守されている根拠にはできない。updated_at はコード更新日時と区別した。

README で試験済みなのは Linux / Python 3.9 / pipenv。Pipfile.lock は NumPy 1.21.0、SciPy 1.7.0、pandas 1.2.5、Matplotlib 3.4.2、pyvalem 2.1.0、PyYAML 5.4.1 等に固定されている。Python 3.14 への導入手順として、そのロックをそのまま使えるとは扱えない。

setup.py の version は 0.00 で、install_requires はない。pip install だけで依存一式が揃う構成ではない。README 自体も GitHub からの pip インストールは未試験と述べる。Windows では sh 用テスト起動スクリプトの代替も要る。

単体テストと、保存済み計算結果への回帰テストは存在する。これは長所だが、最新環境での成功や現在の空気・大気圧・サブ ns 条件の検証済みを意味しない。今回、それらのテストは実行していない。

公開ルート一覧に LICENSE ファイルは見当たらず、GitHub API の license も null、確認した setup.py / README に明示的な利用ライセンスは見当たらなかった。フォーク後の改変・配布・他コードへの取り込みを採用方針にする前に、著者に適用ライセンスを確認することが重要。これは調査で確認できた表記の報告であり、利用権限についての法的判断ではない。

出典: [リポジトリ](https://github.com/martin-hanicinec-ucl/plaschem)、[API メタデータ](https://api.github.com/repos/martin-hanicinec-ucl/plaschem)、[Pipfile.lock](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/Pipfile.lock)、[setup.py](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/setup.py)、[計算回帰テスト](https://github.com/martin-hanicinec-ucl/plaschem/blob/master/run_tests/test_pygmol/test_pygmol.py)。

## 5. LoKI-GM との関係

LoKI-GM は MATLAB 製で、電子 Boltzmann 方程式の LoKI-B と化学ソルバ LoKI-C を連成または独立して使える。README は気相・表面の記述、EEDF・輸送係数・反応係数・状態密度の出力、パルス等を説明している。今回確認した README は LoKI-GM_26.07 向けで、推奨する最も古い MATLAB は R2024b。GPL-3.0 の LICENSE がある。API の pushed_at は 2026-10-23 で、plaschem より最近の更新が確認できる。

電子 Boltzmann ソルバを含む総合的なモデルや参考計算系としては有力。ただし MATLAB を必要とし、既存 Python コードへの小さな差し替えにはならない。今回 LoKI-GM の詳細なモデル監査や実行比較は行っていない。

出典: [LoKI-GM README](https://github.com/LoKI-Suite/LoKI-GM/blob/master/README.md)、[API メタデータ](https://api.github.com/repos/LoKI-Suite/LoKI-GM)、[LICENSE.txt](https://github.com/LoKI-Suite/LoKI-GM/blob/master/LICENSE.txt)。

## 6. 採用方法の評価と次の手順

| 方法 | 作業量の見立て | 判断 |
|---|---|---|
| 現行の化学量論+BDF+BOLOS を継続 | 小〜中。足りない物理・反応データの追加が中心 | 現在の用途では推奨 |
| plaschem の反応機構管理だけ接続 | 中。粒子名、状態、係数形式、M、XML/属性の変換が必要 | 機構が増えて管理が難しくなった段階で検討 |
| optichem を現行ソルバにつなぐ | 中〜大。ラッパーが期待する API、結果保存、反応無効化との接続を調査・実装 | 多数の反応の縮約や感度解析が主目的なら価値あり |
| pygmol へ全面移行し現在の条件を再現 | 大。電子系・境界条件・数値制御・検証を改修 | 容易な置き換えとは評価しない |
| LoKI-GM を別の参照計算系として利用 | 中〜大。MATLAB 環境と機構・断面積・条件の対応付けが必要 | 電子系を含むモデル比較に有用 |

作業量はソースからの定性的な推定で、所要日数の見積もりではない。

plaschem をさらに試す場合は、(1) ライセンス確認、(2) 既存環境と分離した環境で依存を現代化、(3) 上流テストの再現、(4) 小さな反応系で反応速度・行列・単位・保存則を比較、(5) 現行44種のデータ変換と同条件比較、の順がよい。反応縮約が目的なら電子モデルまで置き換える必要性を先に検討する。

手元には ZDPlasKin と同一機構・断面積を使った既存の比較記録もあり、今後の変更時の検証資産として残す価値がある。今回その比較結果を再計算・再検証したわけではない。
