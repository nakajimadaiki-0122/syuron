# cad/ の一覧

`scripts/27_cad_index.py` が生成。図は `figures/cad_index.png`。

## 命名規則

    {形状}_{こぶの大きさ}_{役割}.{拡張子}

| 形状 | 意味 |
|---|---|
| `valve1` | Gamboa のバルブ 1 個（プレナム込み。CFD の計算領域そのもの） |
| `valve1bare` | 同、プレナムを外したバルブ本体だけ |
| `stair{n}` | 段々（n 段。1 段おきに反転して平均方向を水平にした構成） |
| `side{n}fwd` | 直線流路の側面にこぶ n 個（こぶが順方向を向く） |
| `side{n}rev` | 同、こぶが逆方向を向く（x 反転） |

| 大きさ | ループ外半径 R | こぶ頂部 |
|---|---:|---:|
| `small` | 1.9 w_v | 4.3 w_v |
| `std` | 2.35 w_v（Gamboa Table 2） | 5.2 w_v |
| `large` | 3.0 w_v | 6.5 w_v |

| 役割 | 中身 |
|---|---|
| `_outline.dxf` | 2D 輪郭（外形 + 島）。CAD で押し出す用 |
| `_fluid.stl` | 流体体積。固体から引く用 |
| `_plate.stl` | 溝を彫った板（裏面にポート穴）。実物モデル |
| `_lid.stl` | 蓋（平板） |

## 形状の一覧

| ファイル | 形状 | 大きさ | R [w_v] | 外形 [mm] | 流体体積 [mm³] |
|---|---|---|---:|---|---:|
| `side4fwd_large_fluid.stl` | 直線 + 側面こぶ 4 個（順方向） | large | 3.0 | 44.5 × 8.1 × 1.0 | 122.78 |
| `side4fwd_small_fluid.stl` | 直線 + 側面こぶ 4 個（順方向） | small | 1.9 | 30.2 × 5.9 × 1.0 | 81.13 |
| `side4fwd_std_fluid.stl` | 直線 + 側面こぶ 4 個（順方向） | std | 2.35 | 36.1 × 6.8 × 1.0 | 98.21 |
| `side4rev_large_fluid.stl` | 直線 + 側面こぶ 4 個（逆方向） | large | 3.0 | 44.5 × 8.1 × 1.0 | 122.78 |
| `side4rev_small_fluid.stl` | 直線 + 側面こぶ 4 個（逆方向） | small | 1.9 | 30.2 × 5.9 × 1.0 | 81.13 |
| `side4rev_std_fluid.stl` | 直線 + 側面こぶ 4 個（逆方向） | std | 2.35 | 36.1 × 6.8 × 1.0 | 98.21 |
| `stair2_std_fluid.stl` | 段々 2 段 | std | 2.35 | 21.8 × 10.2 × 1.0 | 58.95 |
| `stair4_large_fluid.stl` | 段々 4 段 | large | 3.0 | 49.0 × 12.1 × 1.0 | 132.74 |
| `stair4_small_fluid.stl` | 段々 4 段 | small | 1.9 | 38.0 × 8.9 × 1.0 | 94.05 |
| `stair4_std_fluid.stl` | 段々 4 段 | std | 2.35 | 42.5 × 10.2 × 1.0 | 109.90 |
| `stair6_std_fluid.stl` | 段々 6 段 | std | 2.35 | 61.3 × 10.2 × 1.0 | 156.03 |
| `stair8_std_fluid.stl` | 段々 8 段 | std | 2.35 | 81.0 × 10.2 × 1.0 | 204.56 |
| `valve1_std_fluid.stl` | バルブ 1 個（プレナム込み、CFD 領域） | std | 2.35 | 21.7 × 14.7 × 1.0 | 102.80 |
| `valve1bare_std_fluid.stl` | バルブ 1 個（プレナム無し） | std | 2.35 | 10.1 × 8.0 × 1.0 | 24.27 |

各形状には同じ名前で `_outline.dxf` と `_plate.stl` がある
（`valve1_std` には `_lid.stl` も）。

寸法は流路幅 w_v = 1.0 mm、深さ 1.0 mm。変えるときは生成し直す。

```
scripts/15_gamboa_full_geom.py --wv-mm 1.0      # valve1 の DXF
scripts/20_make_3d_model.py    --wv-mm 1.0      # valve1 の STL
scripts/22_make_chain.py -n 8 --size std        # 段々（任意の段数）
scripts/26_make_variants.py                     # 変種一式（大きさ 3 通り）
scripts/27_cad_index.py                         # この一覧を作り直す
```
