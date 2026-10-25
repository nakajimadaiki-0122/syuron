# experiment_simulation — ANSYS Fluent 用（Xia 2025 の共役伝熱、4 構成）

`handover/cad/xia2025/` から、共役伝熱に要る 2 ファイルだけを構成ごとに抜き出したもの。
元ファイルと中身は同一（コピー）。寸法の根拠は `handover/docs/xia2025_geometry.md`。

| フォルダ | 構成 | 島の数 |
|---|---|---:|
| `symtvm_fwd/` | 対称型 Tesla（SYMTVM）・順方向 | 30 |
| `symtvm_rev/` | 対称型 Tesla（SYMTVM）・逆方向 | 30 |
| `tvm_fwd/`    | 従来型 Tesla（TVM）・順方向 | 10 |
| `tvm_rev/`    | 従来型 Tesla（TVM）・逆方向 | 10 |

各フォルダの中身

| ファイル | 領域 | 材料 | 範囲 Y [mm] |
|---|---|---|---|
| `fluid.stl`        | 流体領域（流路の空洞） | 水 | 0.5 〜 2.0 |
| `solid_copper.stl` | 固体領域（流路を彫った銅ブロック、蓋なし） | 銅 | 0 〜 2.0 |

2 つを重ねると 50 × 3.54 × 2.0 mm の直方体になる。流体と固体の接する面が共役伝熱の interface。

## 座標（Y が鉛直）

| 方向 | 軸 | 範囲 [mm] |
|---|---|---|
| 流れ方向（入口 x = 0、出口 x = 50） | +X | 0 〜 50.0 |
| 流路幅 | Z | −1.77 〜 +1.77 |
| 厚み | +Y | 0 〜 2.0 |

Reverse は Forward を 180 度回転した形状なので、入口は常に x = 0 側でよい。

## 手順書

- `fluent_meshing_steps.md` — Fluent Meshing でメッシュを作るまで（`zones/` の分割 STL を使う）
- `ansys_setup.md` — ソルバ設定・材料・境界条件・取り出す量（論文の条件）
- `split_zones.py` — `fluid.stl` / `solid_copper.stl` を境界ごとの STL（`zones/`、10 ファイル）に分ける。
  平らな面は三角形を作り直している（元の分割は細長い三角形が多く Fluent が読み込み時に削って穴が開く）。
  固体側の界面は書き出さない（流体側の壁面が共有面を兼ねる。二重にすると free faces になる）。
  `handover/.venv` の Python で実行（shapely, scipy, mapbox_earcut が要る）

## 境界条件の要約（詳細と手順は `ansys_setup.md`）

| 面 | 条件 |
|---|---|
| 底面 y = 0 | 一定熱流束 |
| 上面 y = 2.0 | 断熱（蓋は不要） |
| 側面 z = ±1.77 | 対称 |
| 入口 x = 0 | 質量流量 0.5459 〜 1.6377 g/s（1 本あたり）、脱イオン水 |
| 出口 x = 50 | 圧力出口 |

乱流モデルは Realizable k-ε + Enhanced Wall Treatment、SIMPLE、2 次風上。
