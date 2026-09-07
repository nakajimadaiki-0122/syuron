# results/

計算結果と実測値。**用途ごとにサブディレクトリへ分けています。**
`.npz`（場のデータ）と `ck_*.npy`（チェックポイント）は git 管理外です。

| ディレクトリ | 中身 |
|---|---|
| `validation/` | 解析解による検証。`io_bc_*`（2 次元の入口出口 BC）、`validate_3d`（3 次元ダクト） |
| `geometry/` | 幾何の実測・生成。`fig2_*`（Fig.2 の位相・プレナム）、`fig7_digitized`（Fig.7 の目標値）、`gamboa_full_geom`（生成形状と Fig.2 の重ね）、`loop_junction_angles`（7 月形状の接合角） |
| `single/` | **案 1（プレナム込み単発モデル）の計算結果。**`gamboa_full_Re{Re}_cpm{n}_{fwd,rev}.json` と `_fields.npz`、`gamboa_full_phi_*`（φ_loop） |
| `periodic/` | **案 2（周期セル）の結果。**`cell_*`（7 月形状の Re スイープ）、`gamboa_case2_*`、`di_uncertainty`、`gamboa_diag_*`、`_baseline_2026-08-09/` |
| `cad/` | CAD 出力のメタ情報（三角形数・体積・非多様体辺の検査結果） |

## 案 1 と案 2 を混同しないこと

- `single/` は**プレナム込みの単発バルブ**。Gamboa Fig.7 と比べるのはこちら
- `periodic/` は**周期配列**。7 月形状の性能評価には使えるが、
  Gamboa の Di は再現できない（接合部の位相が違う。`progress/2026-09-06.md` §3）
