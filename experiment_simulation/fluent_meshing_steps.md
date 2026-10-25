# Fluent Meshing やり直し手順（ゾーン分割済み STL を使う）

`split_zones.py` で `fluid.stl` と `solid_copper.stl` を境界ごとのファイルに分けてある
（`<構成>/zones/`、10 ファイル）。Fluent は STL をファイル名のゾーンとして読むので、
**読み込んだ時点で inlet / outlet / bottom_heat などの名前が付いている。** Separate は要らない。

| ファイル | 面 | 計算での役割 |
|---|---|---|
| `inlet.stl` | 流体の x = 0 端面 | 質量流量入口 |
| `outlet.stl` | 流体の x = 50 端面 | 圧力出口 |
| `fluid_top.stl` | 流体の y = 2.0 面 | 断熱・no-slip 壁（蓋の代わり） |
| `fluid_bottom.stl` | 流体の y = 0.5 面（溝の底） | 固体との界面（**固体側にはこの面を入れない。二重になる**） |
| `fluid_walls.stl` | 島と外壁 | 固体との界面（同上） |
| `bottom_heat.stl` | 固体の y = 0 面 | 熱流束 1e6 W/m² |
| `solid_top.stl` | 固体の y = 2.0 面 | 断熱 |
| `sym_left.stl` / `sym_right.stl` | 固体の z = ∓1.77 面 | 対称 |
| `solid_ends.stl` | 固体の x = 0 / 50 端面 | 断熱 |

---

## 1. Fluent を起動し直す

今の Fluent は File → Exit で閉じる（保存しなくてよい）。Fluent Launcher で

| 項目 | 設定 |
|---|---|
| モード | **Meshing** |
| Dimension | 3D |
| Double Precision | ON |
| Meshing / Solver Processes | 4 |
| Working Directory | `C:\syuron\experiment_simulation\symtvm_fwd` |

Start。

## 2. 読み込み

1. Workflow → **Watertight Geometry**。
2. Import Geometry: File Format **Mesh**、File Units **mm**。
3. File Names の「...」で `C:\syuron\experiment_simulation\symtvm_fwd\zones` を開き、
   **Ctrl+A で 10 ファイルを全部選ぶ**。開く。
4. Import Geometry を押す。

Outline View → Geometry Objects（または Boundary → Manage の一覧）に 10 個のゾーンが
ファイル名で並んでいれば成功（fluid_walls と solid_top は島ごとに自動で細分されるが問題ない）。Bounds で X 0〜50、Y 0〜2、Z ±1.77 を確認。

## 3. （通常は不要）1 つのオブジェクトにまとめ、重複節点を結合する

**2026-10-24 の実測では、10 ファイルは 1 オブジェクトの 10 パーツとして読まれ、節点も自動で結合された**
（ログに「consisting of 10 parts/objects」「0 boundary faces removed」、free faces の警告なし）。
その場合は本節を飛ばして §4 へ。Create Regions で領域が正しく出ないときだけ以下を試す。

1. Outline View → Geometry Objects で全オブジェクトを選び、右クリック → **Merge**。名前は `heatsink`。
2. Console で重複節点（ファイルの継ぎ目で重なっている点）を結合する:

```
/boundary/merge-nodes
```

   聞かれたら対象は全ゾーン（`*` または Enter で既定）、許容値は既定のまま。

3. 自由節点（つながっていない点）が無いことを確認:

```
/boundary/count-free-nodes
```

   0 なら閉じている。0 でなければそのときの表示を見せてほしい。

（コマンド名が通らなければ、Boundary タブ → Manage / Modify のあたりに
Merge Nodes の GUI があるのでそれでもよい。）

## 4. Workflow を上から順に

| タスク | 設定 |
|---|---|
| Add Local Sizing | `no` のまま次へ |
| Generate the Surface Mesh | Min Size **0.03**、Max Size **0.3**（mm）、Growth Rate 1.2。Generate |
| Describe Geometry | **The geometry consists of both fluid and solid regions**。Change all fluid-fluid boundary types は No |
| Enclose Fluid Regions (Capping) | 何もせず次へ（端面はすでにある） |
| Create Regions | Estimated Number of Fluid Regions = 1。Create |
| Update Regions | 一覧で `fluid` 側が fluid、`solid` 側が solid になっているか確認。違っていたら Region Type を直す |
| Add Boundary Layers | Number of Layers **5**、Offset Method smooth-transition、Grow on は fluid 側の壁（fluid_walls, fluid_bottom, fluid_top） |
| Generate the Volume Mesh | Fill With **poly-hexcore**。Generate |

名前に inlet / outlet が入っているので、Update Regions でそれぞれ
velocity-inlet / pressure-outlet に自動で割り当たる（入口の種類は Solution で mass-flow-inlet に変える）。

目標の要素数は論文の SYMTVM 約 280 万、TVM 約 185 万。
Max Size 0.3 で足りなければ 0.2 に下げる。

## 5. Solution へ

Workflow の最後で **Switch to Solution**。以降は `ansys_setup.md` §3〜§5 のとおり。
最初に File → Write → Case で `symtvm_fwd.cas.h5` を保存しておく。
