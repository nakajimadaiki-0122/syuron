# CAD 出力

Gamboa 2005 optimized valve の形状。単位は mm、既定は流路幅 w_v = 1.0 mm、
深さ 1.0 mm（正方形断面）。寸法を変えるときは生成し直す。

```
../.venv/Scripts/python.exe scripts/15_gamboa_full_geom.py --wv-mm 1.0   # DXF
../.venv/Scripts/python.exe scripts/20_make_3d_model.py --wv-mm 1.0 --depth-mm 1.0
```

| ファイル | 内容 |
|---|---|
| `gamboa_fluid_2d.dxf` | 流体の 2D 輪郭（プレナム込み）。外形 + 島の 2 本の閉じたポリライン |
| `gamboa_valve_only_fluid_2d.dxf` | 同上、プレナムを外した版 |
| `gamboa_fluid_3d.stl` | 流体体積（プレナム込み）。CFD 領域そのもの。102.80 mm³ |
| `gamboa_valve_only_fluid_3d.stl` | 流体体積（バルブ本体のみ）。24.27 mm³ |
| `gamboa_plate_3d.stl` | 溝を彫った板。25.7 × 18.7 × 2.0 mm、裏面に φ2 mm のポート穴 2 つ |
| `gamboa_lid_3d.stl` | 蓋（平板 25.7 × 18.7 × 1.0 mm）。板と貼り合わせると流路が閉じる |

## 注意

- **STL は流体（または板）の実体**である。流体 STL は固体ブロックから引く用途を想定
- 島（涙滴形）は 2D CFD では自由な島だが、板モデルでは**溝底から立つ柱**になる
- すべて閉じたメッシュ（各辺がちょうど 2 枚の三角形に共有される）であることを
  生成のたびに確認している。`results/cad_models.json` に三角形数・体積・判定が残る
- 三角形分割は `src/mesh.py` の耳刈り。`shapely` の Delaunay は境界辺を保たず
  非多様体になるので使わない
