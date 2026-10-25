# ANSYS Fluent 設定手順 — Xia et al. (2025) と同じ条件で 4 構成を回す

出典: Y. Xia et al., *Numerical study on heat transfer and flow characteristics of
symmetric Tesla-type microchannel heat sinks*, Appl. Therm. Eng. 258 (2025) 124611。
以下はすべて原典 2.2〜2.4 節・Table 2〜6 の **[論文記載]**。原典に無い項目は **[推測]** と明記。

---

## 0. 使うファイル

各構成フォルダの `fluid.stl`（水）と `solid_copper.stl`（銅）の 2 つ。蓋は使わない（上面は断熱）。

| 構成 | フォルダ | 入口 | 備考 |
|---|---|---|---|
| SYMTVM-Forward | `symtvm_fwd/` | x = 0 の面 | 対称型・順方向 |
| SYMTVM-Reverse | `symtvm_rev/` | x = 0 の面 | 形状を 180° 回転済み。入口は同じ側でよい |
| TVM-Forward | `tvm_fwd/` | x = 0 の面 | 従来型・順方向 |
| TVM-Reverse | `tvm_rev/` | x = 0 の面 | 同上 |

座標は Y が鉛直。流れ +X（0〜50 mm）、幅 Z（±1.77 mm）、厚み Y（固体 0〜2.0、流体 0.5〜2.0）。

---

## 1. ジオメトリ（SpaceClaim / DesignModeler）

1. `solid_copper.stl` と `fluid.stl` を**同じ座標系で**読み込む（単位 mm、移動・回転しない）。
2. ファセット → ソリッドへ変換（SpaceClaim: Facets → Convert to Solid、Merge faces を有効）。
   島が 30 個あるので面が細かく割れる。変換が荒れる場合は
   `handover/cad/xia2025/*_outline.dxf` から Fusion でソリッドを作り STEP で持ち込む。
3. 2 ボディを 1 コンポーネントに入れ **Share Topology** を有効にする（界面でメッシュを一致させるため）。
4. 境界に名前付き選択（Named Selection）を付ける。

| 名前 | 面 | 領域 |
|---|---|---|
| `inlet` | 流体の x = 0 端面 | fluid |
| `outlet` | 流体の x = 50 端面 | fluid |
| `fluid_top` | 流体の y = 2.0 面（蓋に当たる面） | fluid |
| `bottom_heat` | 固体の y = 0 面 | solid |
| `sym_left` / `sym_right` | 固体の z = ±1.77 面 | solid |
| `solid_top` | 固体の y = 2.0 面 | solid |
| `solid_ends` | 固体の x = 0 / x = 50 端面 | solid |
| 界面 | 流体と固体の接する面全部 | Share Topology なら自動で interior / coupled wall になる |

---

## 2. メッシュ

**[論文記載]** 格子独立性は 1.0918 g/s で確認し、採用格子は

| 構成 | 要素数 | 根拠 |
|---|---:|---|
| TVM | **1,846,186** | 3,261,266 まで増やしても ΔP・Tmax の変化 < 0.75 %（Table 5） |
| SYMTVM | **2,816,560** | 4,828,512 まで増やしても変化 < 0.97 %（Table 6） |

要素の種類・境界層の層数・y+ は原典に書かれていない。**[推測]** Fluent Meshing の
poly-hexcore で壁面に inflation 5 層程度、流路幅 0.31 mm に対して 10 セル以上を目安にし、
上の要素数になるようサイズを調整する。Enhanced Wall Treatment を使うので y+ は 1 付近が望ましい。

---

## 3. ソルバ設定（Fluent 2022R1）

| 項目 | 設定 |
|---|---|
| ソルバ | Pressure-based、**Steady**、重力なし |
| エネルギー | ON |
| 乱流モデル | **Realizable k-ε**、Near-wall: **Enhanced Wall Treatment**（EWF/EWT） |
| 圧力-速度連成 | **SIMPLE** |
| 圧力の離散化 | Standard |
| 運動量・エネルギー | **Second Order Upwind** |
| k, ε | 原典に記載なし。**[推測]** Second Order Upwind |
| 輻射 | なし |

モデル定数は Fluent 既定（C2 = 1.9, σk = 1.0, σε = 1.2）。原典の「σε = 12」は 1.2 の誤植と見る。

### 収束判定

| 方程式 | 残差 |
|---|---|
| 連続 | < 1e-5 |
| 運動量（x, y, z） | < 1e-6 |
| エネルギー | < 1e-9 |

k, ε は記載なし。**[推測]** 1e-5。併せて ΔP と底面最高温度のモニタが頭打ちになることを確認する。

---

## 4. 材料

### 固体: 銅（Table 3）

| 密度 | 比熱 | 熱伝導率 |
|---|---|---|
| 8978 kg/m³ | 381 J/(kg·K) | 387.6 W/(m·K) |

### 流体: 脱イオン水、温度依存（Table 4）

Fluent の各物性を **piecewise-polynomial** または **expression** で入れる。
**T の単位に注意**: 密度だけ °C、残り 3 つは K。

| 物性 | 式 |
|---|---|
| 密度 ρ [kg/m³]（T in °C） | (998.84 + 18.225 T − 7.92e-3 T² − 5.545e-5 T³ + 1.498e-7 T⁴ − 3.933e-10 T⁵) / (1 + 1.816e-2 T) |
| 粘度 μ [Pa·s]（T in K） | 2.414e-5 × 10^(247.8 / (T − 140)) |
| 比熱 cp [J/(kg·K)]（T in K） | 8958.9 − 40.535 T + 0.11243 T² − 1.014e-4 T³ |
| 熱伝導率 k [W/(m·K)]（T in K） | −0.58166 + 6.3556e-3 T − 7.964e-6 T² |

Fluent の expression では温度変数が K なので、密度は `(T - 273.15)` を °C として代入する。

**Named Expressions にそのまま貼れる形**（300 K で ρ 995.9、μ 8.54e-4、cp 4179、k 0.608 を検算済み）

```
Tk    = StaticTemperature/1[K]
rho_w = (998.84 + 18.225*(Tk-273.15) - 7.92e-3*(Tk-273.15)^2 - 5.545e-5*(Tk-273.15)^3 + 1.498e-7*(Tk-273.15)^4 - 3.933e-10*(Tk-273.15)^5) / (1 + 1.816e-2*(Tk-273.15)) * 1[kg/m^3]
mu_w  = 2.414e-5 * 10^(247.8/(Tk-140)) * 1[Pa*s]
cp_w  = (8958.9 - 40.535*Tk + 0.11243*Tk^2 - 1.014e-4*Tk^3) * 1[J/kg/K]
k_w   = (-0.58166 + 6.3556e-3*Tk - 7.964e-6*Tk^2) * 1[W/m/K]
```

water-liquid をデータベースから Copy し、Density / Cp / Thermal Conductivity / Viscosity を expression にして上の名前を入れる。
銅は Fluent Database の copper が Table 3 と同値（8978 / 381 / 387.6）なのでそのまま Copy でよい。

---

## 5. 境界条件

| 境界 | タイプ | 値 |
|---|---|---|
| `inlet` | **Mass-flow-inlet** | 0.5459 〜 1.6377 g/s（下表）、**T = 300 K** |
| `outlet` | Pressure-outlet | **0 Pa（ゲージ）**、逆流温度 300 K |
| `bottom_heat` | Wall、Heat flux | **1.0e6 W/m²** |
| `fluid_top` | Wall、no-slip、断熱（heat flux 0） | 蓋の代わり |
| `sym_left` / `sym_right` | **Symmetry** | 固体側面 |
| `solid_top`, `solid_ends` | Wall、断熱 | |
| 流体–固体界面 | Coupled wall（no-slip、温度連続） | Share Topology なら自動 |

乱流の入口条件は原典に記載なし。**[推測]** Intensity 5 %、Hydraulic diameter = D_h
（TVM 0.514 mm、SYMTVM 0.877 mm）。

### 入口質量流量（Fig. 10 の 11 点、刻み 0.10918 g/s）

計算単位 1 本あたりの値。Fluent には kg/s で入れる。

| # | g/s | kg/s | # | g/s | kg/s |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.5459 | 5.459e-4 | 7 | 1.2010 | 1.2010e-3 |
| 2 | 0.6551 | 6.551e-4 | 8 | 1.3102 | 1.3102e-3 |
| 3 | 0.7643 | 7.643e-4 | 9 | 1.4193 | 1.4193e-3 |
| 4 | 0.8734 | 8.734e-4 | 10 | 1.5285 | 1.5285e-3 |
| 5 | 0.9826 | 9.826e-4 | 11 | 1.6377 | 1.6377e-3 |
| 6 | 1.0918 | 1.0918e-3 | | | |

まず **1.0918 g/s**（原典の格子独立性・等高線図の条件）で 4 構成を回し、
図 11〜17 と見比べてから残りの流量へ広げる。

Re の目安（原典）: 0.5459 g/s で TVM 585、SYMTVM 513。
D_h = 2 Wch Hch / (Wch + Hch)、Wch は TVM 0.31 / SYMTVM 0.62、Hch 1.5 mm。

---

## 6. 取り出す量（2.4 節）

| 量 | 定義 | Fluent での取り方 |
|---|---|---|
| ΔP | 入口と出口の面積平均圧力の差 | Surface Integrals → Area-Weighted Average, Pressure |
| Tmax, Tmin, δT = Tmax − Tmin | `bottom_heat` 面の温度の最大・最小 | Surface Integrals → Facet Max / Min |
| Tw,ave | `bottom_heat` の面積平均温度 | Area-Weighted Average |
| Tin, Tout | 入口・出口の平均温度 | Mass-Weighted Average |
| Rt | (Tmax − Tmin) / q_w | 後処理 |
| h_av | q Ly W / [A_sf (Tw,ave − 0.5 (Tin − Tout))]、A_sf は界面面積 | 界面面積は Report → Surface Integrals → Area |
| Nu | h_av D_h / k_f,ave | 後処理 |
| f | 2 ΔP D_h / (ρ u_in² L) | 後処理 |
| Di | ΔP_reverse / ΔP_forward（同じ流量） | rev と fwd の ΔP の比 |
| α | (Tmax − Tmin) / T_ave | 後処理 |
| PEC | (Nu/Nu0) / (ΔP/ΔP0)^(1/3)、基準 0 は矩形流路 RM | RM は本リポジトリに無い |

### 原典との照合に使える数値（1.0918 g/s）

**訂正（2026-10-24）: Table 5・6（格子独立性）の値は Reverse の値である。** 以前は Forward と推測していたが、
Fig. 13 の棒グラフで Reverse の棒（SYMTVM-Reverse ≈ 14,500、TVM-Reverse ≈ 56,500 Pa）と一致し、
Fig. 10a の Tmax も Reverse の曲線（≈ 355〜356 K）と一致する。

| 構成 | ΔP [Pa] | 底面 Tmax [K] | 出典 |
|---|---:|---:|---|
| SYMTVM-Forward | ≈ 10,500 | ≈ 363 | Fig. 13 / Fig. 10a の目視（±500 Pa / ±2 K） |
| SYMTVM-Reverse | **14,458** | **355.27** | Table 6 |
| TVM-Forward | ≈ 24,500 | ≈ 366 | Fig. 13 / Fig. 10a の目視 |
| TVM-Reverse | **56,352** | **356.55** | Table 5 |

比で見る量（格子の粗さの影響を受けにくい）

| 量 | 原典 | 出典 |
|---|---|---|
| Di = ΔP_rev / ΔP_fwd（SYMTVM） | ≈ 1.38 | Fig. 13 から計算 |
| Di（TVM） | ≈ 2.3（1.6377 g/s で 2.48） | Fig. 13 と本文 |
| Tmax の順位 | SYMTVM-R < TVM-R < SYMTVM-F < TVM-F | Fig. 10a |
| ΔP の順位 | TVM-R > TVM-F > SYMTVM-R > SYMTVM-F | Fig. 13 |
| 熱伝達率 h_rev / h_fwd | SYMTVM 1.106〜1.117、TVM 1.267〜1.285 | 本文 Fig. 18a の説明 |

---|---:|---:|
| TVM-Forward | 56,352 | 356.55 |
| SYMTVM-Forward | 14,458 | 355.27 |

（Table 5・6 の第 4 格子。格子独立性の表は Forward で取られていると読む **[推測]**。）
この 2 値が数 % 以内で合えば設定は原典と同じと判断してよい。

---

## 7. 原典と本 CAD の既知の差

- 本 CAD は Fig.1(b)(e) の **計算単位 1 本**（幅 3.54 mm）。原典も同じ単位で計算している。
- 原典の全長は 50.2 mm 相当になる寸法記載があり、本 CAD は 50.0 に切ってある（`handover/docs/xia2025_geometry.md`）。
- 島の縁の寸法 2.00（原典ラベル）は本 CAD で 2.165。影響は小さいが ΔP に数 % 出る可能性がある。

---

## 8. 原典 Fig. 11・12・14 と同じ図を作る（Fluent の Results）

原典の座標は z が厚み方向（底 z = 0、流路は z = 0.5〜2.0）。**本 CAD は Y が厚み**なので、原典の「z = 1 mm」は本モデルの **y = 1 mm** に当たる。
図の向きは原典が流れ方向を横、幅方向を縦に描いているので、**上から（−Y 方向に）見下ろす**。
4 構成で色の範囲を必ずそろえる（Auto Range を外して下表の Min / Max を入れる）。

| 原典 | 中身 | 断面・面 | 量 | 色の範囲（原典） |
|---|---|---|---|---|
| Fig. 11 | 熱源面の温度分布 | `bottom_heat:1`（y = 0） | Static Temperature | 324〜364 K |
| Fig. 12 | z = 1 mm の温度分布（流体＋銅） | 平面 y = 0.001 m | Static Temperature | 302〜362 K |
| Fig. 14 | 中央面の速度分布（流体のみ） | 平面 y = 0.001 m（流体側だけ表示） | Velocity Magnitude | TVM 0.2〜4.2、SYMTVM 0〜2 m/s |
| Fig. 15 | z = 1 mm の流線 | 同上 | Pathlines（入口から） | — |
| Fig. 17 | 乱流運動エネルギー | 同上 | Turbulent Kinetic Energy | 原典を見て合わせる |

### 手順

1. **断面を作る**: Results タブ → Surface → Create → **Plane**。Method を **XZ Plane**（または Point and Normal で
   点 (0, 0.001, 0)、法線 (0, 1, 0)）、Y = **0.001**。名前 `y1mm`。
2. **Fig. 11**: Results → Graphics → Contours → New。Contours of = Temperature / Static Temperature、
   Surfaces = `bottom_heat:1`、Options の Auto Range を外し Min 324 / Max 364、Filled にチェック。Save/Display。
   底面なので**下から（+Y 方向に）見上げる**。
3. **Fig. 12**: 同じく Contours、Surfaces = `y1mm`、Min 302 / Max 362。上から見下ろす。
4. **Fig. 14**: Contours of = Velocity / Velocity Magnitude、Surfaces = `y1mm`。
   平面は銅も切るが銅の速度は 0 なので、流体だけ見たいときは Contours の「Zones」で流体セル領域だけにする
   （または平面を作るときに Zones で流体だけ選ぶ）。
5. **Fig. 15**: Results → Graphics → Pathlines。Release from Surfaces = `inlet`、Style を line、Step Size を小さめに。
   原典は渦の場所を拡大しているので、2〜3 段目の弁を拡大して保存。
6. **画像の保存**: File → Save Picture。Format PNG、Resolution を 2000 px 以上、Coloring は Color。
   ファイル名は `<構成>_fig11.png` のようにそろえる。

### 定量の図にするとき（任意）

画像だけでなく数値で比べたいときは、Results → Plots → XY Plot で
中心線（z = 0、y = 1 mm）や熱源面の流れ方向の温度分布を出し、File → Write でテキストに書き出す。
原典 Fig. 11 は熱源面の流れ方向の温度上昇（入口付近 ≈ 324〜330 K → 出口付近 ≈ 356〜364 K）が読めるので、
x 方向の温度分布を 4 構成で重ねれば定量比較になる。
