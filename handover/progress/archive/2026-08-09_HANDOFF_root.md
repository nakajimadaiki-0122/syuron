# Tesla型マイクロチャネル研究 引き継ぎ書

更新日: 2026-08-09

## 0. この引き継ぎ書の目的

VS Code上で研究を継続するときに、過去の議論を再構成せずに次の解析へ入るための状態固定ファイルです。

最重要方針は以下です。

> **熱解析・新構造・生体模倣へ進む前に、既知のTesla valve流動を現在の2Dソルバで再現できるか検証する。**

現在は「先行研究の再現前」です。完了済みなのは解析解による基本ソルバ検証と、2026年7月形状の流動再計算だけです。

---

## 1. 研究の大目的

半導体冷却向けマイクロ流路について、Tesla型構造を出発点として、

1. 既知のTesla流動・伝熱を再現する
2. 現状性能とあるべき姿のギャップを定量化する
3. ボトルネックを特定する
4. そのボトルネックに対してのみ新規構造を提案する

という順序で研究を再構成する。

**禁止事項:** ボトルネック未特定のまま、脈動流・ディンプル・生体模倣・トポロジー最適化・CHTへ飛ばない。

---

## 2. 現在地

### 完了

- 2D直線流路の基本検証
  - `np.roll` により上下境界が意図せず周期化していたバグを発見
  - solid padding追加後に直線流路が収束
  - `u_max/u_avg = 1.4816`
  - NY=12での離散理論値 `1.48276` と一致
- 2026年7月Tesla形状（DXF）の流動再計算
  - Re_Dh = 100
  - forward / reverse両方向を計算
  - 圧力ダイオード性はほぼ1
  - loopは「死水域」ではなく、一定量の流量が存在

### 未完了

- Gamboa et al. (2005) の再現
- Porwal & Thompson (2016/2018) の再現
- Han et al. (2024) の再現
- 7月形状と既知Tesla形状の同条件直接比較
- 熱解析再開
- CHT / R_cond分解
- 新規構造提案

したがって、**先行研究はまだ1つも再現できていない**。

---

## 3. 現在の7月形状 baseline（ユーザー報告値）

条件:

- Re_Dh = 100
- D_h = 2 mm
- 水 300 K
- 12 cells/mm
- lattice inlet velocity U = 0.05
- Ma ≈ 0.087
- inlet: 放物線速度分布
- outlet: 定圧

結果:

| 指標 | Forward | Reverse |
|---|---:|---:|
| Δp [Pa] | 12.466 | 12.525 |
| loop流量比 平均 | 21.3% | 18.8% |
| loop流量比 最大 | 30.9% | 28.7% |
| u_loop / u_main | 0.222 | 0.225 |
| mass imbalance | 5.10% | 5.12% |

圧力ダイオード性:

`Di_p = Δp_reverse / Δp_forward = 1.005`

### 解釈の注意

- 「loopが死んでいる」は撤回。流速が低いことと流量が小さいことは同義ではない。
- 現条件では**方向依存圧損がほぼ発生していない**。
- ただし mass imbalance ≈ 5.1% は大きい。絶対圧損・loop流量比の最終議論前に改善する。
- 12 cells/mmは複雑形状のstaircase誤差確認には粗い可能性がある。24 cells/mmとの比較を行う。
- `Di_p=1.005` だけから「形状が悪い」と断定しない。既知Tesla形状を同じソルバで解くことが必要。

---

## 4. 直近の研究質問

現在の最重要質問は1つです。

> **現在の2Dソルバは、既知の正常なTesla valveの方向依存流動を再現できるか。**

判別:

- Gamboa形状で論文2D CFDに近いDiが出る
  - ソルバはTesla流動を表現可能
  - 7月形状のgeometry起因の可能性が高まる
- Gamboa形状でもDi ≈ 1
  - solver / BC / rasterization / grid / 2D実装を疑う

この判別が終わるまで熱解析に戻らない。

---

## 5. 次の実行順序

### Gate 0: 直線流路validationを完成させる

既に速度分布は合格。追加で以下を確認する。

- `u_max/u_avg`: 理論値から1%以内
- Darcy定義で `f_D * Re ≈ 96`: 数%以内
- `|m_in - m_out| / m_in < 1%`（目標0.5%）
- steady residualが十分収束

**mass flowは必ず `rho * u_n` を積分する。** 単純な速度積分ではなく質量流量で評価する。

### Gate 1: 7月形状のgrid check

Re=100だけでよい。

- 12 cells/mm
- 24 cells/mm

比較:

- Di_p
- loop flow fraction
- Δp

12→24で結論が変わるなら、現在baselineは格子依存。

### Phase 1: Gamboa 2005 2D Tesla benchmark

最優先。

1. `339_1.pdf` Fig.1 + Table 2からoptimized geometryを再構築
2. 熱なし
3. Gamboaと同じ2D定常非圧縮層流条件
4. Re=100 forward/reverse
5. Δp_F, Δp_R, Di_pを取得
6. 問題なければ Re = 50, 100, 200, 500 に拡張
7. Gamboa Fig.7の**2D CFD曲線**と比較

実験値ではなく、まず同じ2D CFDを比較対象にする。

### Phase 2: Porwal 2016/2018

Gamboa後。

目的:

- 低Reの3D GMF Tesla流動を独立確認
- helix/main mass-flow ratioを比較
- multi-stageでのforward/reverse差を理解

### Phase 3: Han 2024 flow-only

まだ熱を入れず、Han geometryについて

- Re=100
- Re=300

で以下を見る。

- streamline / velocity field
- branch flow
- vortex
- Δp_F / Δp_R
- Di_p

### Phase 4: Han vs 7月形状

同一solver・同一Re・同一流体・同一BCに揃えて比較。

主指標:

`phi_loop = m_dot_loop / m_dot_total`

他:

- local velocity
- vorticity
- pressure loss
- branch/rejoin mechanism

### Phase 5: 熱解析再開

流動benchmarkが通った後だけ。

- Nu
- R_conv
- R_cal

### Phase 6: 3D + CHT

必要になった時点でsolidを追加し、R_condを含める。

### Phase 7: 新構造

ボトルネックが確定してからのみ実施。

---

## 6. Gamboa et al. 2005 benchmark仕様

論文: *Improvements in Fixed-Valve Micropump Performance Through Shape Optimization of Valves*, Journal of Fluids Engineering 127(2), 339–346.

### 数値条件

- 2D
- steady
- incompressible
- laminar
- Navier–Stokes
- inlet: uniform velocity
- outlet: zero pressure
- wall: no-slip
- characteristic velocity: inlet mean velocity U
- characteristic length: valve channel width w_v
- `D_h = 2 w_v`
- `Re = rho U D_h / mu`
- original grid: 16 elements across valve channel
- grid doublingで解変化 <4%

### Optimized geometry parameters

全長さは `w_v` で無次元化。

- X2 = 1.60
- n = 0.797
- Y3 = 0.608
- LENOUT = 2.94
- R = 2.35
- alpha = 41.9 deg
- beta = 71.7 deg（独立変数ではなくX2, n, Y3から決まる）
- nX2 = 1.2752

Fig.1の定義:

- X2: forward-flow inlet segment length
- n: scale factor giving coordinate nX2
- Y3: coordinate defining outer tangent location of straight return segment
- R: loop outer radius
- LENOUT: outlet segment length
- alpha: outlet segment angle
- origin: inlet channel centerline上、左端から0.5 channel-widthの位置

重要な設計原理:

> return angle betaを大きくし、loop return flowにmain-channel flowと逆向きの速度成分を与えることが高Diに寄与する。ただしreverse時にloopへ十分な流量が入ることも必要。

### 注意

- Table 2の「Average Di = 1.50」はRe=0–2000での平均であり、Re=100のtarget値ではない。
- Re=100 targetはFig.7の2D CFD曲線から比較する。必要なら後でdigitizeする。
- 論文自身が2D CFDと実験の差を報告しているため、最初のvalidation対象は実験値ではなく2D CFD値。

---

## 7. Porwal 2016/2018 benchmark仕様

- GMF Tesla
- 3D CFD
- square cross-section
- D_H = 1 mm
- Re = 25–200
- N = 10 stages（主解析）
- valve-to-valve distance = 1 mm
- hydrodynamic entrance length `delta = 0.05 Re D_H`
- Re=200でdelta=10 mm、実際には20 mm entranceを付加
- steady / incompressible / laminar / single-phase
- no-slip
- velocity inlet
- pressure outlet
- constant properties
- thermal case: wall 20°C, inlet 80°C
- FLUENT 3D, tetrahedral mesh

重要指標:

- Di_P
- helix mass-flow fraction
- Nu / Di_T（熱再開後）

Porwalではforwardとreverseでhelix流量割合のRe依存が異なる。これを「正常なGMFの流量分配」の参考にする。

**注意:** 修論中の一部friction-factor相関式には符号/指数の内部不整合が疑われるため、相関式を無批判にgold standardにしない。

---

## 8. Han et al. 2024 benchmark仕様

論文: *A comparative study of enhanced thermal performance in Tesla-type microchannels*, Applied Thermal Engineering 239, 122157.

### Geometry (Table 1)

- heatsink length L = 30 mm
- heatsink width W = 10 mm
- channel width W_ch = 0.45 mm
- channel height H_ch = 0.50 mm
- fin thickness W_t = 0.22 mm
- unit bend length L_t = 1.40 mm
- bend radius R_d = 0.15 mm
- narrow section width W_d = 0.22 mm

### Numerical model

- 3D FEM / COMSOL
- copper + DI water
- laminar incompressible Newtonian
- single-phase
- fluid properties temperature-dependent
- inlet temperature = 303.15 K
- velocity inlet based on Re
- outlet pressure = 0
- bottom heat flux = 16.5 W/cm²
- other surfaces adiabatic
- conjugate solid-fluid interface

### Definitions

`Re = rho u_m D_h / mu`

`D_h = 2 H_ch W_ch / (H_ch + W_ch)`

`Di_p = ΔP_B / ΔP_F`

`Di_t = Nu_B / Nu_F`

`PEC = (Nu/Nu0) / (f/f0)^(1/3)`

### 注意

- AbstractはRe=100–450、numerical BC sectionはRe=50–400。図ごとに対象Reを確認する。
- Abstractの`Di_p ≈ 1.9`は高Re側の値であり、Re=100 targetではない。
- 現在の2D flow-only solverでHanのR_thやPECをいきなり一致させない。まず流動のみ比較。

---

## 9. 既知のバグ・技術リスク

### 修正済み

`np.roll` streaming + geometryが配列端に接触していたため、上下壁が周期境界化していた。

対策: solid paddingを必ず設ける。

### 未解決

1. mass imbalance ≈5.1%
2. 12 cells/mmのgeometry resolution
3. Zou-He / outlet BCによる反射・密度偏りの可能性
4. rasterized斜め壁のstaircase error
5. 2Dではout-of-plane vorticityを表現できない

### 判断ルール

- straightが通らない → solver/BCを直す
- straightが通るがGamboaが通らない → geometry/Re定義/BC/grid/raster/2Dを順に疑う
- Gamboaが通るが7月形状がDi≈1 → geometry起因が有力

---

## 10. 次のVS Codeセッションで最初にすること

1. リポジトリ構造を確認
2. 現在使っているLBM/FVコードを特定
3. `np.roll` wall-padding修正がコードに残っていることを確認
4. straight caseのmass imbalanceとfDReを計算
5. 12→24 cells/mmのgrid check
6. Gamboa Fig.1を見ながらgeometry generatorを実装
7. Re=100 forward/reverseを実行
8. summary JSONを保存

**熱解析のコードは消さなくてよいが、呼ばない。**
