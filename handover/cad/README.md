# CAD 出力

**どのファイルがどの形かは [INDEX.md](INDEX.md) と `figures/cad_index.png` を見てください。**
一覧は `scripts/27_cad_index.py` が STL から自動生成します。

## 命名規則

    {形状}_{こぶの大きさ}_{役割}.{拡張子}

| 形状 | 意味 |
|---|---|
| `valve1` | Gamboa のバルブ 1 個（プレナム込み。CFD の計算領域そのもの） |
| `valve1bare` | 同、プレナムを外したバルブ本体だけ |
| `stair{n}` | 段々（n 段。1 段おきに反転して平均方向を水平にした構成） |
| `side{n}fwd` | 直線流路の側面にこぶ n 個（こぶが順方向を向く） |
| `side{n}rev` | 同、こぶが逆方向を向く（x 反転） |

| 大きさ | ループ外半径 R | こぶ頂部 | 備考 |
|---|---:|---:|---|
| `small` | 1.9 w_v | 4.3 w_v | **R < 1.85 では島が自己交差して形状が壊れる**（実測） |
| `std` | 2.35 w_v | 5.2 w_v | Gamboa Table 2 の最適値 |
| `large` | 3.0 w_v | 6.5 w_v | |

| 役割 | 中身 |
|---|---|
| `_outline.dxf` | 2D 輪郭（外形 + 島）。CAD で押し出す用。単位 mm |
| `_fluid.stl` | 流体体積。固体から引く用 |
| `_plate.stl` | 溝を彫った板（裏面にポート穴）。実物モデル |
| `_lid.stl` | 蓋（平板）。板と貼り合わせると流路が閉じる |

## 作り直し方

既定は流路幅 w_v = 1.0 mm、深さ 1.0 mm（正方形断面）、溝の下の肉厚 1.0 mm、
板の縁 2.0 mm、ポート径 2.0 mm。

```
../.venv/Scripts/python.exe scripts/15_gamboa_full_geom.py --wv-mm 1.0   # valve1 の DXF
../.venv/Scripts/python.exe scripts/20_make_3d_model.py    --wv-mm 1.0   # valve1 の STL
../.venv/Scripts/python.exe scripts/22_make_chain.py -n 8 --size std     # 段々（任意の段数）
../.venv/Scripts/python.exe scripts/26_make_variants.py                 # 変種一式（大きさ 3 通り）
../.venv/Scripts/python.exe scripts/27_cad_index.py                     # 一覧を作り直す
```

## 注意

- **STL は流体（または板）の実体**である。流体 STL は固体ブロックから引く用途を想定
- 島（涙滴形）は 2D CFD では自由な島だが、板モデルでは**溝底から立つ柱**になる
- 流路幅 1 mm に φ2 mm ポートは入らないので、**両端に半径 1.6 mm の丸い溜まり**を
  付けてその中心をポートにしている
- すべて閉じたメッシュ（各辺がちょうど 2 枚の三角形に共有される）であることを
  生成のたびに確認している。結果は `results/cad_*.json`
- 三角形分割は `src/mesh.py` の耳刈り。`shapely` の Delaunay は境界辺を保たず
  非多様体になるので使わない
- **`side{n}fwd` / `side{n}rev` は出口区間を閉じ壁に置き換えた周期構成**なので、
  Gamboa の単発バルブとは接合部の位相が違う（逆流がループへ入るのに 132 度曲がる）。
  Gamboa の Di を再現する形状ではなく、「直線流路 + 側面のこぶ」という構成そのものを
  試すための形状である
- 段々（`stair{n}`）は入口・出口面が水平から 24.05 度傾いている。
  CFD にかけるならプレナムを付けるか助走の直管で向きを揃えること
