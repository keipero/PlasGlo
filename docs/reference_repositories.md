# Py3ZDPlasChem の参考リポジトリ

調査日: 2026-10-31（日本時間）。MATLAB を使わず Python で開発を継続する条件で、公開 README、関連ソース、公式文書、GitHub API メタデータを比較した。インストール・実行・性能比較は未実施。

## 優先する3候補

### 1. ChemPlasKin: 気体・プラズマ反応と熱収支

- [リポジトリ](https://github.com/ShaoX96/ChemPlasKin)、[README](https://github.com/ShaoX96/ChemPlasKin/blob/main/README.md)
- C++ 製。Cantera と CppBOLOS をソースレベルで組み合わせた0次元反応計算。
- 現在の BOLOS + 化学反応という構成に近く、反応熱、プラズマによる加熱、振動エネルギーから並進への緩和を検討する参考として優先する。
- [plasmaReactor.h](https://github.com/ShaoX96/ChemPlasKin/blob/main/src/plasmaReactor.h) に、定圧／定容の切り替え、化学種・温度・エネルギーの状態変数、熱収支、N2 振動エネルギーの緩和等がある。反応ごとのエネルギー配分を移植する場合、明示的な励起状態との二重計上を確認する必要がある。
- [examples](https://github.com/ShaoX96/ChemPlasKin/tree/main/examples) に Air、N2、H2O2He、DBD_H2O2 がある。[反応機構変換ツール](https://github.com/ShaoX96/ChemPlasKin/tree/main/data/PAC_kinetics/ZDPlasKin_kinetics/kineticsParser) も参考候補。
- そのまま使う導入負担は大きめ。README の手順は Cantera 3.0 のソース修正・再ビルドを含み、Windows は未試験で WSL を勧めている。通常の Cantera Python パッケージを入れるだけで ChemPlasKin になるわけではない。
- API の pushed_at は 2026-02-14。[LICENSE](https://github.com/ShaoX96/ChemPlasKin/blob/main/LICENSE) は LGPL 2.1（README には LGPLv2 と記載）。

### 2. Cantera: Python API、反応機構の形式、比較計算

- [リポジトリ](https://github.com/Cantera/cantera)、[公式プラズマパルス例](https://cantera.org/stable/examples/python/reactors/nanosecond_pulse_discharge.html)
- C++ コアと Python API を持つ化学反応・熱力学・輸送ライブラリ。3.2 の公式 Python 例では E/N のガウス波形、EEDF の更新、区間ごとの反応積分が示される。
- 参考にする対象は、入力データと計算コードの分離、YAML 反応機構、反応器 API、EEDF と反応計算の更新、結果保存・診断。
- [v3.2.0 の例のソース](https://github.com/Cantera/cantera/blob/v3.2.0/samples/python/reactors/nanosecond_pulse_discharge.py) はメタン・空気機構、ConstPressureReactor、energy="off" を使用している。現在の N2/O2 密閉定容系とそのまま同条件ではなく、時間依存 EEDF や電子冷却の完全なモデルとして扱わない。
- [電子衝突反応の YAML 定義](https://cantera.org/stable/yaml/reactions.html) には、断面積から係数を求める反応や EEDF 計算用衝突データの記述がある。
- Python から実際に利用・比較できる候補として優先する。ただし、手元の励起種・反応式・係数・熱力学データの対応付けと、Python 3.14 環境の適合性は別途確認が必要。
- API の pushed_at は 2026-10-24。[License.txt](https://github.com/Cantera/cantera/blob/main/License.txt) は BSD 型の3条件を記載。

### 3. 独立版 PyGMol: 電子エネルギー方程式と Python の設計

- [リポジトリ](https://github.com/hanicinecm/pygmol)、[README](https://github.com/hanicinecm/pygmol/blob/master/README.rst)
- 最初に調べた plaschem と同じ著者による、別の独立パッケージ。こちらは [MIT ライセンス](https://github.com/hanicinecm/pygmol/blob/master/LICENSE) を明示している。
- SciPy solve_ivp で重粒子密度と電子エネルギー密度を解く。[方程式の設計説明](https://github.com/hanicinecm/pygmol/blob/master/docs/doc_equations.rst) があり、状態変数は重粒子密度と電子エネルギー密度。電子密度は準中性から算出し、気体温度はパラメータとして一定にする。
- 参考にする対象は、電子エネルギー方程式、Chemistry / PlasmaParameters / Equations / Model の分離、単体・回帰・文書テストの構成。
- 現在の E/N + BOLOS 駆動をそのまま置き換える候補とは区別する。吸収電力、衝突エネルギー損失、粒子数変化との整合性を確認して必要部分を設計する。
- API の pushed_at は 2022-02-25。[setup.py](https://github.com/hanicinecm/pygmol/blob/master/setup.py) は SciPy 1.7.3 固定で、現在の環境への導入・テストは未実施。新しい保守済み環境と同一視しない。

## 後の段階で役立つ候補

| リポジトリ | 参考にする点 | 現在の優先度・制約 |
|---|---|---|
| [CRANE](https://github.com/lcpp-org/crane) | 反応ネットワークの定義、源項の生成、0次元から輸送モデルへの接続 | C++ / MOOSE が必要。Python の小規模0次元コードへ直ちに組み込むには重い。API の pushed_at は 2026-09-24、LGPL-2.1 |
| [RePlaChem](https://github.com/znikolaou/RePlaChem) | DRGEP による重要種・経路の評価と機構縮約 | Fortran 製。ZDPlasKin 形式の機構と専用形式の反応速度履歴が必要で、現在の出力との接続が要る。API の pushed_at は 2026-05-08、GPL-3.0 |
| [CppBOLOS](https://github.com/ShaoX96/CppBOLOS) | 既存 BOLOS と近いアルゴリズムの C++ 実装 | 化学計算全体の代替ではない。性能改善を検討する段階の候補で、速度優位は手元の条件で未測定 |
| [LoKI-B++](https://github.com/LoKI-Suite/LoKI-B-cpp) | 別の電子 Boltzmann ソルバとの比較・連成設計 | C++ / CLI 実行。LoKI-C の代替ではなく、Python 連成には接続実装と係数・断面積の対応付けが要る |

## 開発への使い分け

まず ChemPlasKin の熱・振動緩和モデル、Cantera の Python パルス例、独立版 PyGMol の電子エネルギー式を読む。現行コードの基本構造を維持し、必要な物理を小さな検証可能な単位で追加する。全体の導入に進む場合は、Cantera の通常の Python API が現在の要件を満たすかを先に試す価値がある。

各リポジトリの pushed_at は保守状況の参考指標であり、リリース日、テスト成功、モデルの妥当性を保証するものではない。個々の反応データ・断面積の出典や利用条件はコード本体と別に確認する。
