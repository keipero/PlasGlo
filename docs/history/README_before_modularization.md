# E/Nを指定する0次元プラズマ反応計算

今回の主モデルは局所電界近似（LFA）です。指定E/N・現在の基底N2/O2組成・気体温度からBOLOSで定常EEDFと係数を更新し、電子エネルギー方程式は独立に積分しません。[採用したモデル方針](modeling_decisions.md)を参照してください。

**Py3bolos接続済み（2026-10-31）:** `Py3bolos/bolos/` を直接読み込みます。共有Python環境へのBOLOSのインストールは不要です。反応リストも `mechanisms/` に同梱し、このフォルダ単独で入力を解決します。過去のローカル実行出力は削除し、接続後のテスト結果を作成しました。[接続状況](BOLOS_REPLACEMENT.md)と[テスト結果](verification/README.md)を参照してください。

## ZDPlasKinと同じ窒素・アルゴンサンプル

`samples/zdplaskin/nitrogen/` に21種・96反応、`samples/zdplaskin/argon/` に3種・2反応の入力を置いています。両方に、元の `kinet.inp`、断面積 `bolsigdb.dat`、設定を確認できる `user_code.F90` と出典・ハッシュを保存しています。空気用Komuro機構とは別のサンプルです。

```powershell
python Py3ZDPlasChem-master/samples/zdplaskin/run_samples.py --sample all
```

出力は `output_samples/nitrogen/` と `output_samples/argon/` です。`--sample nitrogen` または `--sample argon` で個別に実行できます。[条件・単位・制約](samples/zdplaskin/README.md)を参照してください。

窒素の基本実行ノートブックは [run_nitrogen.ipynb](run_nitrogen.ipynb) です。NumPy・SciPy・Matplotlib・IPythonが使えるPythonカーネルで開き、上から全セルを実行してください。条件セルでE/N・温度・圧力・初期電子密度・計算時間・格子数を変更できます。既定は1000 Td・300 K・1 ns・1600セルで、全21種の密度とグラフ等を `output_notebook_nitrogen/` に保存します。入力データの照合、保存則の診断、描画も含みます。

電流の後処理セルも追加しています。BOLOSの換算電子移動度から電子伝導電流、入力したイオン移動度からイオン伝導電流、電界の時間変化から変位電流を計算します。断面積が未入力なら電流密度のみを表示し、イオン移動度が不足すると総電流は出力しません。[電流モデル・入力値・制約](CURRENT_MODEL.md)を参照してください。

## ZDchem.py：1000 Td・0.05 nsのアフターグロー比較

`ZDchem.py` の現在の設定は1000 Tdを0.05 ns印加し、3 nsまで計算します。
`CHEMISTRY_EXTENSION = "atomic_oxygen"` により、元の134反応に次の2反応を追加します。
`None` に戻すと元の機構です。`COMPARE_ORIGINAL = True` では同一条件で両機構を計算します。

- N₂(B) + O(³P) → NO + N(⁴S)
- N₂(C) + O(³P) → NO + N(⁴S)

両係数は3×10⁻¹⁶ m³/sです。出典は
[Blanchard, Minesi and Laux (2026), PSST 35 085006, Tables 3–4](https://doi.org/10.1088/1361-6595/ae8935)。
同論文のSections 3.1–3.2は、これらが**直接測定されていない仮定係数**であることを明記しています。
原典に挙げられたShkurenkov et al. (2014), PSST 23 065003のモデルに基づく気体運動論的な推定です。
今回の実装は2026年論文に従い、生成物をNO(X)とN(⁴S)へまとめます。
NOの振動エネルギー・発熱は解かず、元機構で複数状態をまとめたN₂(B)へB状態の係数を適用する近似も含みます。
**1000 Td・高電離率におけるモデル全体の検証済み修正ではなく、欠けていた重粒子反応の影響を調べる拡張です。**
元のafivo反応データは編集せず、追加反応と出典・仮定を `mechanisms/atomic_oxygen.txt` に保存しています。

基底状態N₂/O₂とその他の選択粒子は別のパネルに描画します。
N₂(C)の生成率・消失率と消失時定数 `n/L` も表示・保存します。
これは生成による補給を含む正味の減衰時定数とは異なります。
源項の復元では、各時刻のBOLOS更新履歴を使い、最終時刻の0 Td係数を印加中へ適用しません。
電界を切っても原子の再結合・励起種の緩和は有限の反応速度で進むため、基底気体を初期値へリセットしません。

出力先は `output_atomic_oxygen/`、元機構の比較結果はその `original/` です。
従来の `output/` は上書きしません。
`n2c_balance.csv` / `.png` と、ノートブック実行時の `n2c_mechanism_comparison.png` を追加します。
`used_reactions.tsv` には反応ごとの `source_file` 列を追加しています。

BOLOSの衝突標的は依然として基底状態N₂/O₂のみで、電子間・超弾性衝突、電子エネルギー方程式、気体加熱は未実装です。
断面積のない中性種または電離率が1%を超える場合は図とノートブックへ確認メッセージを表示します。
1%は適用範囲の確認用目安であり、それ以下の精度や物理的妥当性を保証する基準ではありません。
追加反応でN₂(C)の減衰が速くなっても、この制約は残ります。

CLIで同じ拡張を使う場合：

```powershell
python Py3ZDPlasChem-master/kinetics.py --field-td 1000 --field-duration-ns 0.05 --end-ns 3 --seed-density 1e13 --bolos-update-ns 0.01 --points 1201 --chemistry-extension atomic_oxygen --output Py3ZDPlasChem-master/output_atomic_oxygen
```

APIは `simulate(..., chemistry_extension="atomic_oxygen")` です。
API/CLIの既定値は元の機構を維持しています。

Juliaのstreamer計算から独立したPythonプログラムです。空間的に一様な密閉・定容系について、134反応・44粒子種の密度を時間積分します。既定では、同梱のPython 3版BOLOSを現在の基底状態N₂/O₂密度で実行し、EEDF・電子衝突係数・電子温度を逐次更新します。BOLSIGの出力ファイルは読みません。縦軸は密度[m⁻³]、横軸は時間[ns]です。

## 実行

ワークスペースのルートで実行します。Python環境にNumPy・SciPy・Matplotlibが必要です。

```powershell
python -m pip install -r Py3ZDPlasChem-master/requirements.txt
python Py3ZDPlasChem-master/kinetics.py --field-td 100 --end-ns 100 --bolos-update-ns 1 --output Py3ZDPlasChem-master/output_bolos_100Td
```

描画する粒子を絞る場合は、記号を含む名前を引用符で囲みます。CSV・NPZには描画の選択にかかわらず全粒子を保存します。

```powershell
python Py3ZDPlasChem-master/kinetics.py --field-td 200 --end-ns 10 --species e "N2+" "O2+" "O2-" --output Py3ZDPlasChem-master/output_200Td
```

主なオプション：

| オプション | 既定値 | 意味 |
|---|---:|---|
| `--field-td` | 100 | 印加中のE/N[Td]。BOLOSは0以上、tableは1〜2000 |
| `--field-duration-ns` | 省略 | E/N印加時間[ns]。以後は0 Td。省略すると計算終了まで印加 |
| `--swarm-backend` | bolos | BOLOS逐次計算。`table`で旧固定テーブルと比較 |
| `--bolos-update-ns` | 1 | BOLOSと化学計算の結合間隔[ns] |
| `--bolos-cells` | 400 | 電子エネルギー格子数。50以上 |
| `--bolos-rtol` | 1e-6 | EEDFの反復収束許容誤差 |
| `--cross-sections` | 同梱Phelps | C1〜C43に対応する断面積ファイル |
| `--end-ns` | 100 | 計算終了時刻[ns] |
| `--seed-density` | 1e12 | 初期電子密度[m⁻³]。等量の正イオンも配置 |
| `--pressure-pa` | 101325 | 初期気体密度を決める圧力[Pa] |
| `--gas-temperature-k` | 300 | 計算中一定の気体温度[K] |
| `--n2-fraction` | 0.8 | 初期N₂モル分率、0〜1。O₂は残りの割合 |
| `--points` | 501 | 出力点数。積分内部の刻み幅は自動調整 |
| `--rtol` | 1e-6 | 相対許容誤差 |
| `--atol-density` | 1e3 | 密度の絶対許容誤差[m⁻³] |

marimoの `ZDchem.py` からもE/N・結合間隔・描画粒子を変更できます。CLIと同じ `kinetics.py` を呼び出し、BOLOSは隣接する `Py3bolos/` のクローンから読み込みます。marimoを利用する場合は `python -m pip install marimo` が必要です。

E/Nの印加時間と、アフターグローを含む計算時間は別々に設定できます。例えば1000 Tdを0.05 nsだけ印加し、5 nsまで計算する場合：

```powershell
python Py3ZDPlasChem-master/kinetics.py --field-td 1000 --field-duration-ns 0.05 --end-ns 5 --bolos-update-ns 0.05 --gas-temperature-k 300 --n2-fraction 0.8 --output Py3ZDPlasChem-master/output_pulse
```

Python/Jupyter/marimoでは、`FIELD_DURATION_NS`、`END_NS`、`GAS_TEMPERATURE_K`、`N2_FRACTION` を編集します。APIは `simulate(field_td=1000, field_duration_ns=0.05, end_ns=5, gas_temperature_k=300, n2_fraction=0.8)` です。N₂:O₂=70:30なら `n2_fraction=0.7` とします。印加時間を `None` にすると従来の連続印加、0にすると最初から0 Tdです。

電界を切る時刻では、出力点やBOLOS更新間隔に関係なく化学積分を必ず分割します。切った後は0 Tdの定常EEDFを再計算し、電子温度依存・E/N依存の解析式も更新します。放射・衝突失活・その他の化学反応は継続します。電子エネルギー方程式を解いていないので、電界オフ直後の電子冷却過程は再現せず、定常EEDFへ即時に切り替える近似です。0 TdのEEDFを気体温度のMaxwell分布へ強制はしません。超弾性衝突のない現在の断面積モデルの限界に注意してください。

初期気体密度は `p/(kB*Tg)` で決め、初期正イオンの配分も指定したN₂/O₂比に合わせます。温度はBOLOSの気体温度と化学反応の解析式へ反映します。ただし断面積そのものは同梱データを再利用します。特にPhelpsの有効回転励起断面積は300 Kを想定しており、大きく温度を変える場合は断面積・反応係数の適用温度範囲の見直しが必要です。

固定BOLSIGテーブルの `table` モードは80:20・300 K・連続印加専用です。印加時間・組成・温度を変更する計算には `bolos` を使用してください。

## BOLOSと化学計算の結合

1. 現在の `n(N2), n(O2)` から混合比を求め、定常Boltzmann方程式を解きます。
2. EEDFから41個の電子衝突反応係数と平均電子エネルギーを求めます。電子温度を更新し、電子温度に依存する解析式の反応係数も計算し直します。
3. その係数を固定して、指定した結合間隔だけBDFで化学反応を積分します。
4. 確定した区間終端の密度でBOLOSを再実行します。最終時刻も再計算して履歴に保存します。

BOLOSの初回は粗い線形格子で解き、平均エネルギーの15倍までの二次格子で再計算します。以降は同じ格子と前回のEEDFを初期推定として使います。0 Tdへの切り替え時は低エネルギー格子を作り直します。大きな組成変化ではエネルギー格子の収束も確認してください。EEDFの負値は最大値の1e-12を超える大きさならエラーとします。丸め誤差で負値が出る裾は、解像できている末尾の対数勾配で正の指数関数として延長し、再規格化します。単純な正値クリップでlog(EEDF)に崖を作るとBOLOSの積分が異常な係数を返すためです。補正前の最小値と補正回数はdiagnosticsに記録します。0 Tdでは小さすぎる擬似時間刻みによる偽の収束も避けています。

これは一次精度の逐次結合です。化学積分の `rtol` だけでは結合誤差を制御できません。`--bolos-update-ns 0.5`、`0.25` のように短くして、注目粒子の密度が収束することを確認してください。`--points` は出力のサンプリング数で、BOLOSの更新頻度には影響しません。区間中は係数を固定するため、既存の解析的Jacobianはその区間の化学方程式と整合します。

終了時刻でのBOLOS再計算は、その後の積分には使われません。`bolos_updates_during_integration` は初期・最終時刻を除いた更新回数です。計算時間より結合間隔が長ければ積分中の更新は0回となり、組成変化が1%を超えた場合は個別に警告します。サブnsの計算では、既定の1 nsをそのまま使わず、対象の時間尺度に合わせて結合間隔を短くしてください。

このモデルでは `N = n(N2基底状態) + n(O2基底状態)` と定義し、指定したE/Nを維持します。したがって、Nが変わると実電場 `E=(E/N)*N` も変わります。実電場固定の計算ではありません。

**断面積の範囲：** 現在のPhelpsデータは基底状態N₂/O₂に対する43過程のみです。原子・NO・O₃・励起種・イオンの密度は化学計算で追跡しますが、それらを衝突標的にしたEEDFへの寄与は計算できません。励起種を基底状態のN₂/O₂へ合算することもしていません。断面積のない中性種の割合と電離率が1%を超えると警告します。この閾値は精度保証ではありません。全44種の寄与を含めるには、それぞれの弾性・非弾性・超弾性過程の断面積と対応付けの追加が必要です。

旧計算との比較：

```powershell
python Py3ZDPlasChem-master/kinetics.py --swarm-backend table --field-td 100 --end-ns 100 --output Py3ZDPlasChem-master/output_table_100Td
```

Python APIは `simulate(..., swarm_backend="bolos", bolos_update_ns=0.5)` です。既存の `mechanism=` を渡す呼び出しは既定でその固定係数を維持します。カスタムMechanismをBOLOSで更新する場合は `swarm_backend="bolos"` も明示します。

## 出力

- `densities.png`：荷電粒子・励起粒子/原子・基底状態中性粒子の3段グラフ。現在は密度・時間とも線形軸。`--species`指定時は1段。印加時間を指定した場合は電界オフ時刻に破線を表示。
- `densities.csv`：1列目が秒、残りが各粒子の密度[m⁻³]。初期時刻を含みます。
- `densities.npz`：`time_s`、`densities_m3`（粒子×時刻）、`species`。
- `field_history.csv`：出力時刻[s]とE/N[Td]。NPZにも `field_td` を保存します。
- `used_reactions.tsv`：反応式・最終時刻に再計算したSI係数・反応次数・元ファイル行番号。tableでは固定係数。
- `swarm_history.csv`：BOLOS更新時刻[s]、平均電子エネルギー[eV]、電子温度[K]、N₂/O₂密度[m⁻³]と混合比、断面積のない中性種の割合、反復回数・収束誤差、E/N[Td]。BOLOSモードのみ。
- `swarm_history.npz`：更新時刻、平均電子エネルギー、電子温度、`channel_coefficients_SI`（更新時刻×44、C0はダミー）、`reaction_coefficients_SI`（更新時刻×134）。C27はm⁶/s、その他の衝突係数はm³/s。反応係数は反応次数に応じたSI単位。最終時刻の行は終端組成での再計算値で、その後の化学積分には使われません。
- `diagnostics.json`：条件、温度、元素・電荷保存誤差、負密度、組成変化、最終時刻の反応源項の釣り合いなど。

`example_100Td/` は旧固定テーブルによる100 Td・100 nsの実行例です。既存の `output/`、`notebook_output/` も自動では更新しないため、計算モードは各出力のdiagnosticsで確認してください。

`example_bolos_100Td/` は新しいBOLOS結合による100 Td・100 ns、初期電子密度1e12 m⁻³、結合間隔1 nsの実行例です。初期時刻と最終時刻を含めてBOLOSを101回実行しています。

`n2c_review_1000Td/` は、保存済みの1000 Td・0.5 nsの結果についてN₂(C)の生成・消滅と実効寿命を調べた診断です。一定の生成源Pと失活頻度λなら、初期密度0からの応答は `nC=P/λ*(1-exp(-λ*t))` ですが、密閉系でN₂/O₂が枯渇してPとλが急変する場合、この一定条件の応答式は適用できません。急生成の後に長い寿命でゆっくり減る区間も、短い時間軸では平坦に見えます。平坦化だけで生成・消滅の釣り合いやモデルの物理的妥当性を判断しないでください。

## 方程式と「平衡」の意味

各反応rの反応速度と、各粒子sの密度変化を計算します。

\[
R_r=k_r\prod_s n_s^{\alpha_{sr}},\qquad
\frac{dn_s}{dt}=\sum_r\nu_{sr}R_r.
\]

化学量論係数だけでは時間発展は求まりません。反応係数・温度・初期密度が必要です。自然放出は一次反応として、`R=A*n`、寿命は`τ=1/A`です。元データの`N2*(B) -> N2*(A1)`、`N2*(C) -> N2*(B)`、`NO(A) -> NO`には既にA[s⁻¹]が指定されています。リストにない放出経路は補っていません。

衝突失活・クエンチングも元リストにある反応を含みます。例えばB状態の放射とN₂による失活が並行する場合、他の経路や供給を無視した実効寿命は`1/(A+kq*n(N2))`です。寿命を別途掛けると二重計上になります。

この計算は指定したE/Nの印加中と、その後の0 Tdでの反応を時間積分します。印加時間を省略すると全計算時間にわたってエネルギーを与え続けます。熱力学的平衡や詳細釣り合いを解くものではありません。反応源項が釣り合えば定常状態になり得ますが、全反応の逆反応・緩和反応がそろっていないので、必ず定常状態に達するとは限りません。実行完了は定常到達を意味しません。

`final_max_relative_source_imbalance` は各粒子の`|正味源項|/生成消滅の絶対寄与の和`の最大値です。微量生成物も含むため、これだけで定常判定はしません。

## データとモデルの仮定

- 反応式・生成物・解析式：`mechanisms/air_chemistry_Komuro_reactions.txt`。Ono・Komuro 2020を基にしたafivo同梱修正版で、2022年の変更を含みます。
- BOLOS用断面積：`Phelps/Phelps_N2-O2_20261026.txt`。C番号だけでなく標的種・過程種・生成物・閾値を照合します。反応リストの過程名も照合します。
- 比較用の固定BOLSIG出力：`Phelps/N2-O2_80-20_2000Td.dat`。N₂ 80%・O₂ 20%、300 K、1〜2000 Td。`table` を明示した場合のみ使用します。
- 反応リストは元ファイルを改変せず本フォルダへコピーしています。コピー元とSHA-256は `mechanisms/air_chemistry_provenance.json` に記録しています。断面積・反応リスト・Py3bolosの参照先は実行時の作業ディレクトリに依存しません。
- 反応物・生成物を読み取って化学量論行列を構築し、元素・電荷の保存を起動時に検証します。全粒子を追跡し、生成種を無断に省略しません。
- 気体温度は指定値で固定、既定値300 K。電子温度は平均エネルギーから`Te=2*energy_eV/(3*kB_eV)`で求め、再結合係数に使います。電子エネルギー方程式は解きません。
- 初期のN₂/O₂比は指定値、既定値80:20。電子と同数の正イオンを同じ比率で配置し、その分だけ親分子を減らして元素・電荷を保存します。その他の粒子は初期0です。
- N₂・O₂も反応で変化する密閉・定容モデル。密度は初期`N0=p/(kB*Tg)`で規格化して計算します。計算中の圧力を固定するモデルではありません。
- E/Nは印加時間まで指定値、以後0 Td。印加時間省略時は全計算時間で固定します。電子衝突係数と電子温度はBOLOSの各更新時刻で再計算し、区間内だけ固定します。`table` の場合はEEDFを再計算せず、旧来の組成変化・電離率警告を出します。
- Mは基底状態の`n(N2)+n(O2)`と定義します。全励起分子や生成中性粒子を含む第三体モデルではありません。
- 解析式のcm単位を二体では1e-6、三体では1e-12でSIへ変換。Phelps三体付着C27の未加重積分値へ1e-6を掛けて三体係数[m⁶/s]を求め、化学計算では `k3*ne*nO2^2` とします。
- C27はEEDF計算でも第三体である現在のO₂密度[cm⁻³]を掛けます。BOLOS内部では仮想標的の重みを `x(O2)*n(O2)*1e-6` として表し、化学計算側で第三体を二重に掛けないようにします。圧力を変更するとこのEEDF計算にも反映されます。`table` にはこの補正がないため、BOLOSとの差はソルバの差だけではありません。
- 放射は光学的に薄い場合の損失として扱い、放出光の再吸収・光電離・放射輸送は含みません。
- 回転・振動励起種には元リスト内で消費経路がありません。長時間では蓄積し、現実の平衡を再現できません。高電場・長時間では緩和反応、生成種の電子衝突断面積、気体加熱などの追加が必要になります。
- 空間輸送・拡散・壁面損失・電場の自己無撞着計算・外部背景電離は含めません。移動度と拡散係数は空間勾配のない今回の方程式には不要です。

## 数値解法と検証

SciPyのBDF法と解析的Jacobianを使用します。陰解法の反復中に現れる負の試行密度は、反応速度の評価に限り0として扱います。解や密度微分を事後クリップしないため、化学量論上の保存性を維持します。絶対許容誤差の100倍を超える負密度はエラーとし、微小な負値もCSVには残します。グラフでは絶対許容誤差以下の微量密度を非表示にします（既定では1e3 m⁻³以下）。これは密度の物理的な下限ではなく、数値精度以下の値をグラフで強調しないための表示条件です。

```powershell
python -m unittest discover -s Py3ZDPlasChem-master -p "test_*.py" -v
```

現在の7件の自動テストは、拡張反応の保存則・質量作用則・指数減衰、電界切替時の源項、クローンの読込先、三体付着のEEDF損失とSI係数、1/100/2000 Tdと純N₂/純O₂、出力格子と一致しない時刻の電界遮断を確認します。窒素・アルゴンのZDPlasKin比較は `verification/verify_reconnection.py` で別途実行します。

参考：[ZDPlasKin公式説明](https://www.zdplaskin.laplace.univ-tlse.fr/wp-content/uploads/2016/12/readme.pdf)、[SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)。本プログラムはBOLOSと結合した小規模な0次元モデルで、ZDPlasKin全体の機能を再実装するものではありません。
