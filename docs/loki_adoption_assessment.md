# LoKI-GM の採用判断と Py3ZDPlasChem の開発方針

調査日: 2026-10-31（日本時間）

## 判断

利用者の条件は、MATLAB 環境がなく、Python 環境で計算・連成・改修を続けたいというもの。この条件では、LoKI-GM をそのままフォークするより、既存の Py3ZDPlasChem / chemistry_0d を土台として必要機能を実装・検証する方針が適している。

以下は公開ソース・公式文書と手元のコードを調べた評価。LoKI-GM の実行、移植、数値的な相互検証は実施していない。既存の計算コードは変更していない。

## MATLAB 不要で利用できるか

| 方法 | 調査結果 | 採用上の評価 |
|---|---|---|
| LoKI-GM をフォークしてそのまま実行 | MATLAB 製。公式 README は MATLAB R2024b 以降を推奨 | 現在の環境には合わない |
| Python の MATLAB Engine API から呼ぶ | MATLAB 本体のインストールが必要。Runtime のみでは Engine を実行できない | MATLAB 不要という条件を満たさない |
| MATLAB Runtime を使う | コンパイル済みアプリ・コンポーネントを動かす仕組み。ソースをフォークするだけでは利用できない | 別途ビルド済み成果物と、その作成・更新環境の確保が必要 |
| Octave へ移す | LoKI-GM の Chemistry / Setup は classdef の events / notify 等を使用。Octave 公式資料は events / listeners を未実装としている | そのまま動くとは扱えず、互換化と検証が必要 |
| LoKI-C の Python 移植 | 化学量論による反応源項の核は移せるが、状態階層、反応係数関数、輸送・熱モデル、入力・出力まで関係する | 全体移植は相当の開発。既存コードへ必要機能を加える方が範囲を制御しやすい |

LoKI-B++（LoKI-B-cpp）は MATLAB 不要の C++ 製電子 Boltzmann ソルバで、公式 README に Windows のビルド・CLI 実行手順がある。将来 Python から外部実行する連成案は考えられるが、これは電子系の選択肢であり、LoKI-C の化学・熱・表面モデルを Python で使えるようにするものではない。今回ビルド・Python 連成の実証はしていない。

## ソースから確認した、参考になる点

- `Code/Chemistry.m` の `kinetics` は正逆反応速度を計算し、生成物と反応物の化学量論係数の差を掛けて密度微分を求める。体積反応の核は現在のコードと共通する。
- `.chem` ファイルが反応式、係数関数、パラメータを分離している。反応データの拡張や出典管理の設計を検討する材料になる。
- `detailedBalance` は電子衝突の逆過程の係数、係数関数が用意する逆係数、一部の重粒子二体反応の詳細釣り合いを扱う。任意の反応の逆係数を自動的に生成できる実装ではない。
- 気体加熱、壁との熱交換、拡散・壁反応、表面種を扱う仕組みがある。これらを加えるには、反応以外にも物性値と境界条件が必要になる。
- 定常放電では、圧力・準中性・電子系との整合性を反復して求める。入力した電子密度等を満たす E/N を求める場合があり、E/N を固定して密度の過渡応答を積分する現在の計算と区別が必要。
- パルスの準定常計算では、電子輸送係数や衝突係数等の事前計算テーブルを用いる。`localField` は E/N、`localEnergy` は電子温度を使って補間する。これは現在の組成で BOLOS を逐次再計算する方式と同一ではない。
- アフターグローで時間依存 EEDF を密度と同時に積分する分岐もあるが、`solveEedf` の初期値は false で、ソースに EXPERIMENTAL FEATURE と記されている。通常設定で電子冷却が自動的に解決するとは判断できない。

## Py3ZDPlasChem の完成度を上げる順序案

「完成」の対象は、まず現在扱っている N2/O2・指定 E/N・パルス／アフターグローの0次元計算を基準とする。以下は開発の優先順位案であり、既に機能が完成・検証されたという意味ではない。

1. **再現性の確立**
   反応データの所在、依存環境、実行例、データの版・出典を固定し、別環境で同じ入力から再現できる状態を作る。保存則、解析的 Jacobian、単位変換、BOLOS 結合間隔の既存の検証資産を整理する。
2. **反応機構と適用範囲の監査**
   状態のまとめ方、放射・失活、回転／振動緩和、逆反応、電子衝突標的の不足を整理する。反応を追加するたびに出典、適用範囲、保存則と密度への影響を確認する。
3. **電界オフ後の電子系の改善**
   現状の定常 EEDF への即時切り替えについて、対象時間尺度での誤差を評価する。必要に応じて平均電子エネルギー方程式や時間依存 EEDF を選ぶ。エネルギー式の追加だけで、非 Maxwell 型 EEDF の変化を全て再現できるわけではない。
4. **実験条件に必要な熱・境界モデルの追加**
   気体加熱や壁面損失の重要度を、計算時間、反応速度、拡散・熱緩和時間と比較して判断する。追加するときは密閉・定容モデルとの切り替えを明示する。
5. **同条件の外部比較と利用範囲の明記**
   ZDPlasKin 等との比較は反応機構・断面積・初期条件・近似を合わせる。数値的な一致と実放電への妥当性を区別し、検証した条件範囲を文書化する。

LoKI-GM 全体との機能同等性ではなく、対象条件での再現性、数値収束、保存則、必要な物理の網羅性を完成基準にする。

## 参照元

- [LoKI-GM README](https://github.com/LoKI-Suite/LoKI-GM/blob/master/README.md)
- [Chemistry.m](https://github.com/LoKI-Suite/LoKI-GM/blob/master/Code/Chemistry.m)
- [反応オブジェクト Reaction.m](https://github.com/LoKI-Suite/LoKI-GM/blob/master/Code/Reaction.m)
- [窒素の反応機構例](https://github.com/LoKI-Suite/LoKI-GM/blob/master/Code/Input/Nitrogen/nitrogen_mod.chem)
- [LoKI-C の設定説明](https://github.com/LoKI-Suite/LoKI-GM/blob/master/Documentation/html/Chemistry.html)
- [計算モードの公式説明](https://github.com/LoKI-Suite/LoKI-GM/blob/master/Documentation/html/Getting_started.html)
- [MATLAB Engine の必要条件](https://www.mathworks.com/help/matlab/matlab-engine-for-python.html)
- [MATLAB Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html)
- [Octave classdef の対応状況](https://wiki.octave.org/Classdef)
- [LoKI-B++ README](https://github.com/LoKI-Suite/LoKI-B-cpp/blob/main/README.md)

手元の参照: README.md、kinetics.py、bolos_coupling.py、validation_zdplaskin_nitrogen/README.md。
