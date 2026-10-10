# cad/xia2025 — Xia et al. (2025) Fig.3 の 4 構成

出典: Y. Xia et al., *Numerical study on heat transfer and flow characteristics of
symmetric Tesla-type microchannel heat sinks*, Appl. Therm. Eng. 258 (2025) 124611。
寸法の来歴は `../../docs/xia2025_geometry.md`、生成は `../../scripts/31_make_xia2025.py`。

## 座標の約束（**Y 軸が鉛直**）

板は **XZ 平面（地面）に寝ている**。2026-10-10 に Z 厚みから変更した。

| 方向 | 軸 | 範囲 [mm] |
|---|---|---|
| 流路の進行方向 | +X | 0 〜 50.0 |
| 流路幅の方向 | Z | −1.77 〜 +1.77 |
| 板の厚み（鉛直） | +Y | 0 〜 2.0（流路の底 0.5、上面 2.0） |

DXF も同じ XZ 平面に載せてある（押し出し方向 (0, −1, 0) の OCS）。
Z 厚みの旧版が要るときは `python scripts/31_make_xia2025.py --up z`。

## ファイル

| 名前 | 構成 | 島 |
|---|---|---:|
| `xia2025_symtvm_fwd_*` | SYMTVM-Forward（対称型・順方向） | 30 |
| `xia2025_symtvm_rev_*` | SYMTVM-Reverse | 30 |
| `xia2025_tvm_fwd_*` | TVM-Forward（従来型・順方向） | 10 |
| `xia2025_tvm_rev_*` | TVM-Reverse | 10 |

| 役割 | 中身 |
|---|---|
| `_fluid.stl` | 流体体積（SYMTVM 78.20 mm³ / TVM 44.18 mm³） |
| `_block.stl` | 流路を彫った銅ブロック 50 × 3.54 × 2.0 mm。入口・出口は端面に開口 |
| `_outline.dxf` | 流路の 2D 輪郭（閉じた LWPOLYLINE。押し出し 1.5 mm 用） |

Reverse は Forward の 180 度回転（鏡像ではない。Fig.3 の島の向きに合わせた）。

## 検査

全モデルで非多様体辺 0、体積は期待値と一致、STL の法線は全面が頂点順と整合。
数値は `../../results/cad/cad_xia2025.json`、図は
`../../figures/geometry/xia2025_four.png`・`xia2025_detail.png`・`xia2025_3d.png`。

**CFD はまだ回していない。** STEP は出せない（cadquery が Python 3.14 に未対応）。
