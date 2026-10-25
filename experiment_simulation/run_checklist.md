# 1 構成を最初から最後まで回すチェックリスト

2026-10-24 に symtvm_fwd を実際に通した手順。つまずいた箇所は対処を書き込んである。
`<構成>` は symtvm_fwd / symtvm_rev / tvm_fwd / tvm_rev のどれか。

格子は原典の 1/10 が目標（SYMTVM 約 28 万、TVM 約 18.5 万セル）。

---

## A. Fluent Launcher

| 項目 | 設定 |
|---|---|
| モード | Meshing |
| Double Precision | ON |
| Meshing / Solver Processes | 4 |
| Working Directory | `C:\syuron\experiment_simulation\<構成>`（日本語を含む元のパスは使えない） |

## B. Meshing（Watertight Geometry）

| # | タスク | 設定 |
|---|---|---|
| B1 | Import Geometry | File Format **Mesh**、File Units **mm**、`<構成>\zones` の 10 ファイルを Ctrl+A。警告「Free faces」が出ないこと |
| B2 | Add Local Sizing | no |
| B3 | Generate the Surface Mesh | Min **0.03**、Max **0.25**（SYMTVM）/ **0.3**（TVM、まず試す）、Growth 1.2、Size Functions Curvature & Proximity、Cells Per Gap **2** |
| B4 | Describe Geometry | **both fluid and solid regions and/or voids**、cap openings **No**、Share Topology **Yes**、他は No |
| B5 | Apply Share Topology | 既定のまま。「2 fluid/solid regions and 0 voids」と出ること |
| B6 | Update Boundaries | **sym_left / sym_right を symmetry** に。inlet / outlet はそのまま |
| B7 | Create Regions | Fluid Regions 1 |
| B8 | Update Regions | **bottom_heat = solid、fluid_top = fluid**（既定では両方 fluid になる） |
| B9 | Add Boundary Layers | 5 層、smooth-transition、Add in fluid-regions、Grow on only-walls |
| B10 | Generate the Volume Mesh | poly-hexcore。セル数と最小 Orthogonal Quality（0.1 以上）を記録 |
| B11 | 保存 | File → Write → Mesh、`<構成>\<構成>_1of10.msh.h5` |

セル数が目標から外れたら B3 の Max Size を変えて B3 以降をやり直す（セル数はおよそ Max Size の −3 乗に比例）。

## C. Solution の設定

| # | 場所 | 設定 |
|---|---|---|
| C1 | General | Check で総体積 3.54e-7 m³、最小体積が正。Pressure-Based、Steady、Gravity なし |
| C2 | Models → Energy | ON |
| C3 | Models → Viscous | k-epsilon、**Realizable**、**Enhanced Wall Treatment** |
| C4 | Named Expressions | Tk、rho_w、mu_w、k_w（式は `ansys_setup.md` §4）。**1 ダイアログに 1 式** |
| C5 | Materials | Fluent Database から copper（solid）と water-liquid（fluid）を Copy |
| C6 | water-liquid | Density = expression `rho_w`、Viscosity = `mu_w`、Thermal Conductivity = `k_w`、**Cp は polynomial**（4 係数: 8958.9 / −40.535 / 0.11243 / −0.0001014）。Change/Create |
| C7 | **Cell Zone Conditions** | 流体の領域を開いて **Material Name = water-liquid**、固体は **copper**。**名前ではなく Material Name を見る**（1 回目は air のまま計算した） |
| C8 | inlet | 右クリック → Type → mass-flow-inlet。0.0010918 kg/s、Normal to Boundary、Intensity 5 % + Hydraulic Diameter（SYMTVM 0.000877、TVM 0.000514 m）、Thermal 300 K |
| C9 | outlet | Gauge 0 Pa、Backflow Total Temperature 300 K |
| C10 | bottom_heat:1 | Thermal → Heat Flux **1000000**。（相方の shadow が無い bottom_heat の wall がこれ） |
| C11 | 他の wall | 既定のまま（断熱。界面の対は Coupled） |
| C12 | Methods | SIMPLE、Pressure Standard、Momentum / k / ε / Energy すべて Second Order Upwind |
| C13 | Residuals | continuity 1e-5、速度 1e-6、energy 1e-9、k / ε 1e-5 |
| C14 | Report Definitions | `p_in`（inlet の Area-Weighted Average、Static Pressure）、`tmax_bottom`（bottom_heat:1 の Facet Maximum、Static Temperature）。File / Plot / Console にチェック |
| C15 | Initialization | Hybrid → Initialize |
| C16 | 保存 | File → Write → Case & Data、`<構成>_1of10_1p09gs` |
| C17 | Run Calculation | 3000 回。p_in と tmax_bottom が数百反復動かなくなれば止めてよい |

## D. 計算後の確認（必ずやる）

Results → Reports で

| 確認 | 正しい値 |
|---|---|
| Fluxes → Mass Flow Rate（inlet / outlet） | ±0.0010918 kg/s |
| Fluxes → Total Heat Transfer Rate（bottom_heat:1） | 177 W |
| Surface Integrals、inlet の Density | 約 996（1.225 なら空気のまま） |
| Surface Integrals、inlet の Velocity Magnitude | SYMTVM 約 1.2、TVM 約 2.4 m/s |

最終の p_in、tmax_bottom、セル数を `results_log.md` に書く。図は `ansys_setup.md` §8。

---

## E. 2 構成目以降（symtvm_fwd のケースをひな形にする）

symtvm_fwd の `symtvm_fwd_1of10_1p09gs.cas.h5` には、材料・式・境界条件・解き方・レポート・
断面 `y1mm`・コンター図（fig11_temp / fig12_temp / fig14_vel）・視点 `paper_view` がすべて入っている。
**流体のセル領域名が全構成で `fluid_top` にそろったので、Replace Mesh で材料も引き継がれる。**

| # | 作業 | 内容 |
|---|---|---|
| E1 | メッシュ | Launcher（Meshing、Working Directory = `<構成>`）で B1〜B11。`<構成>\zones` を読む。保存名 `<構成>_1of10.msh.h5` |
| E2 | ひな形を開く | Launcher（**Solution**、Working Directory = `<構成>`）→ File → Read → Case & Data で `..\symtvm_fwd\symtvm_fwd_1of10_1p09gs.cas.h5` |
| E3 | 差し替え | Domain → Replace Mesh → `<構成>_1of10.msh.h5`。Console に「No zone with name ... skipped」が出たら、その領域の設定を確認 |
| E4 | 確認 | Cell Zone: fluid_top = water-liquid、bottom_heat = copper。bottom_heat:1 の Heat Flux 1e6。inlet 0.0010918 kg/s。General → Check で総体積 3.54e-7 |
| E5 | TVM だけ | inlet の Hydraulic Diameter を **0.000514** に。Fig. 14 の色範囲を **0.2〜4.2** に |
| E6 | 計算 | Hybrid → Initialize → 保存（名前を手入力 `<構成>_1of10_1p09gs`）→ 3000 回 Calculate。値が平らになったら Stop して上書き保存 |
| E7 | 確認 | D の 4 項目（流量、177 W、密度 996、流速） |
| E8 | 図 | Results → Graphics → Contours の一覧から fig11_temp / fig12_temp / fig14_vel をそれぞれ Display → View → Views で `paper_view` を Apply → File → Save Picture（`<構成>_fig11` など） |
| E9 | 記録 | p_in、tmax_bottom、セル数を `results_log.md` へ |

**TVM の向きの確認**: TVM は上下非対称。Fig. 12 を描いたら原典 Fig. 12(a)(b) と最初の島の位置を見比べ、
上下逆なら Console で `/display/views/camera/up-vector 0 0 1` → `/display/views/auto-scale`。
