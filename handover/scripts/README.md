# scripts/

番号は歴史的な実行順で、依存関係の順ではありません。**用途で引いてください。**
すべて `handover/` を作業ディレクトリにして実行します。

```
cd handover
./.venv/Scripts/python.exe scripts/01_validate_solver.py     # Windows
./.venv/bin/python         scripts/01_validate_solver.py     # macOS / Linux
```

## 1. まず通すもの（検証）

| # | 内容 | 時間 |
|---|---|---|
| `01_validate_solver.py` | **2 次元ソルバの解析解検証**（f·Re = 96、Nu）。ここが崩れていたら他は全部無意味 | 11 秒 |
| `16_validate_io_bc.py` | 入口・出口 BC の検証（入口一様流速。Gamboa の条件） | 3 分 |
| `18_validate_io_bc_hiRe.py` | 同、放物線入口で高 Re まで。**安定限界の確認** | 5 分 |
| `24_validate_3d.py` | **3 次元ソルバの解析解検証**（矩形ダクト f·Re = 56.91 ほか） | 2 分 |

## 2. 幾何をつくる・確かめる

| # | 内容 |
|---|---|
| `13_fig2_topology.py` | **Fig.2 から接合部の位相を実測**。周期版が再現できなかった原因の根拠 |
| `14_fig2_plenum.py` | Fig.2 からプレナムの円を実測（スケールの独立検証にもなる） |
| `23_digitize_fig7.py` | **Fig.7 の opt CFD を画素実測**。比較の目標値はここから |
| `15_gamboa_full_geom.py` | 単発形状（プレナム込み）の生成・Fig.2 重ね・DXF 出力 |
| `08_validate_gamboa.py` | 生成した周期版形状の自己検証（流路幅・θ） |
| `05`, `06` | 7 月形状の DXF から接合角を実測し、Gamboa の β と比較 |

## 3. 流れを解く

| # | 内容 | 使う場面 |
|---|---|---|
| `17_run_gamboa_full.py` | **案 1。単発モデル（プレナム込み）の順流・逆流** | Gamboa の再現はこれ |
| `25_run_gamboa_3d.py` | 同じ形状を 3 次元で | 実物・Porwal・Han へ進むとき |
| `04_run_unitcell.py` | 周期単位セル（体積力駆動） | 7 月形状の Re スイープ |
| `09_gamboa_case2.py` | 案 2（周期 + 閉じ壁） | **Gamboa 再現には使えないと判明** |
| `02_run_flow.py` | 旧・全長モデル | **Re ≥ 200 で発散する。使わない** |

## 4. 結果を見る

| # | 内容 |
|---|---|
| `19_analyze_full_fields.py` | 単発モデルの場から φ_loop と流れ場の図 |
| `12_gamboa_diagnose.py` | 周期版の順逆の流れ場診断 |
| `07_uncertainty.py` | Di の格子・収束依存性（計算済み JSON を読むだけ） |
| `10_make_figures.py` | 発表用の図（Di 対 Re、接合角、数表） |
| `11_field_plot.py` | 速度場・流線 |
| `03_compare.py` | 順逆比較と作図 |
| `21_make_report.py` | **閲覧用 HTML**（3D 表示つき、1 ファイル完結） |

## 5. CAD を出す

| # | 内容 |
|---|---|
| `20_make_3d_model.py` | 単発バルブの STL（流体・溝板・蓋） |
| `22_make_chain.py` | **段々流路**（任意の段数、こぶの大きさ） |
| `26_make_variants.py` | **試作用の変種一式**（段々 / 側面こぶ × 大きさ 3 通り） |
| `27_cad_index.py` | `cad/INDEX.md` と `figures/geometry/cad_index.png` を作り直す |

出力の命名規則は `cad/README.md` を参照。
