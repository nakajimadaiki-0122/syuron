# patterned_loop — SYMTVM のループ内壁にパターンを入れた独自形状

元の形状は Xia et al. (2025) の SYMTVM（流路幅 0.30 mm、深さ 1.5 mm、計算単位 50 × 3.54 × 2.0 mm）。
幹の流路、菱形、入口、出口は変えず、ループ部分の外壁と島の壁だけにパターンを入れる。

| フォルダ | パターン | 寸法 |
|---|---|---|
| `pattern1_wavy/` | 1: 壁を正弦波にする（深さ方向は一様） | 片振幅 0.05 mm、波長 約 0.40 mm |
| `pattern2_hemisphere/` | 2: 半球が流路の内側へ出っぱる | 半径 0.1 mm、中心間隔 0.4 mm、深さ方向 3 段 |
| `pattern3_dimple/` | 3: 壁を半球状にへこませる | 半径 0.1 mm、中心間隔 0.4 mm、深さ方向 3 段 |

半球は底面・上面に置かない。深さ方向の 3 段は底（y = 0.5）から 0.35、0.75、1.15 mm。

## 各フォルダの中身

| ファイル | 内容 |
|---|---|
| `symtvm_fwd/`, `symtvm_rev/` | Forward と Reverse。Reverse は Forward の 180 度回転 |
| `<構成>/zones/*.stl` | Fluent Meshing が読むファイル。パターン 1 は 10 個、パターン 2、3 は `hemispheres.stl` を加えた 11 個 |
| `<構成>/fluid.stl`, `solid_copper.stl` | 水と銅の全体形状 |
| `pattern1_wavy/<構成>/outline.dxf` | パターン 1 の平面図（押し出して CAD にできる） |
| `preview_*.png`, `overview_cad.png` | 形状の図。`overview_cad.png` が全体像 |
| `params.json` | 寸法と検査結果 |

パターン 2、3 は深さ方向に形が変わるので、平面図の DXF は作っていない。

## 作り直し方

`handover` の venv の Python で実行する。Fluent は使わない。

```
python make_patterned.py --out pattern1_wavy
python make_hemisphere.py --mode bump
python make_hemisphere.py --mode dimple
python make_overview.py
```

## 検査結果

Python での検査はすべて通っている。

**Fluent の確認結果（2026-11-01、パターン 2 の Forward のみ）: 読み込みから体積メッシュまで通った。**

| 工程 | 結果 |
|---|---|
| 読み込み | 消された面 0、Free faces の警告なし |
| Share Topology | 成功。水と銅の 2 領域。折れた面 4 枚は Fluent が自動で修正し、最大 skewness 0.77 |
| 体積メッシュ | 3,726,718 セル、最小 Orthogonal Quality 0.10（Min 0.03、Max 0.25、Cells Per Gap 2、境界層 5 層） |

### 保存したメッシュ（パターン 2、3）

学生版の上限は 1,048,576 セル（ライセンスのメッセージで確認）。超えると保存時に強制終了する。

共通の設定: Min 0.06、Max 1.0、増加率 1.4、Cells Per Gap 1、曲率の刻み 30°、境界層 3 層、poly-hexcore。

| ファイル | セル数 | 最小 Orthogonal Quality |
|---|---:|---:|
| `pattern2_hemisphere/symtvm_fwd/symtvm_fwd_coarse.msh.h5` | 974,377 | 0.05 |
| `pattern2_hemisphere/symtvm_rev/symtvm_rev_coarse.cas.h5` | 974,377 | 0.05 |
| `pattern3_dimple/symtvm_fwd/symtvm_fwd_coarse.msh.h5` | 1,025,083 | 0.05 |
| `pattern3_dimple/symtvm_rev/symtvm_rev_coarse.cas.h5` | 1,025,083 | 0.05 |

4 つとも Solution モードで読み込み、Mesh Check を通した（総体積 3.54e-7 m³、最小体積は正）。**計算はまだ流していない。**

Reverse は回転した STL から直接メッシュを作ると Share Topology に失敗した（折れた面 13〜19 枚）。
そこで Forward のメッシュを点 (25, 0, 0) mm を通る y 軸まわりに 180 度回し、inlet と outlet の名前を入れ替えた。
入口は x = 0、出口は x = 50 mm にあることを確認済み。順流と逆流でメッシュが同じになる。

- Reverse は `.cas.h5`。メッシュのほかは Fluent の既定の設定しか入っていないので、材料や境界条件を設定する。
- `Replace Mesh` で `.cas.h5` を読めるかは未確認。
- `*.h5` は Git に入れない。4 つで約 590 MB あり、OneDrive に同期される。
- Min 0.06 mm は半球の半径の約 6 割で、半球 1 周に約 10 セル。
- パターンなしと比べるときは、元の形も同じ設定でメッシュを作り直す。
- パターン 1 は Fluent で未確認。

| 項目 | パターン 1 | パターン 2 | パターン 3 |
|---|---|---|---|
| 閉じていない辺 | 0 | 0 | 0 |
| 水 + 銅の体積 | 354.00 mm³ | 354.00 mm³ | 354.00 mm³ |
| 水の体積（元 78.20 mm³） | 78.20 mm³ | 76.19 mm³ | 80.26 mm³ |
| 壁の三角形の最大内角 | 90° | 163.8° | 163.8° |

半球の体積は多面体で近似しているので、理論値より約 5 % 小さい。

最大内角 171.4° の三角形が 78 枚あった版は Fluent に 156 枚の面を消され、Free faces が出た。163.8° まで下げると消されなくなった。

## Fluent で使うときの注意

- 手順は `experiment_simulation/run_checklist.md` と同じで、`zones` の読み込み元だけ変える。
- 読み込み後に Console の `boundary faces removed` がすべて 0 で、`Free faces` の警告が出ないことを確認する。
- 半球を細かく解像するとセル数が上限を超える。今の設定と制約は上の「保存したメッシュ」を参照。
- 流体のセル領域名は `fluid_top` のままなので、ひな形ケースの Replace Mesh が使える。

## 試作の制約

- パターン 1 は上から切削できる形。波の曲率半径は約 0.08 mm で、φ0.1 mm のエンドミルが入る計算。
- パターン 2、3 は垂直な壁の横向きの半球なので、上から入れるエンドミルでは作れない。
- パターン 3 のへこみの近くに残る銅の最小厚さは約 0.09 mm。
