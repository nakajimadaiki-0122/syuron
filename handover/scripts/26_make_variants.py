#!/usr/bin/env python3
"""
26_make_variants.py -- 試作する流路の変種をまとめて作る
=======================================================

2 種類 × こぶの大きさ 3 通りを出す。

**流路 1（段々）** `figures/stair8_std.png` の形。
Gamboa の 1 段は入口から alpha だけ折れて出るので、1 段おきに上下反転して
折れを相殺し、入口の向きを +24.05 度から始めて平均方向を水平にしてある。
`src/gamboa_full.py: chain()`

**流路 2（直線 + 側面のこぶ）** `figures/fwd_vs_rev.png` の形。
主流路はまっすぐで、同じループが側面に並ぶ。**こぶの向きが順方向のものと
逆方向のもの**（x 反転）の 2 通りを出す。`src/gamboa.py: straight_with_loops()`

こぶの大きさはループ外半径 R で決まる（こぶの頂部は y = w/2 + 2R）。

| 呼び名 | R [w_v] | こぶ頂部 [w_v] | 備考 |
|---|---:|---:|---|
| small | 1.9 | 4.3 | **R < 1.85 では島が自己交差して形状が壊れる**（実測） |
| std | 2.35 | 5.2 | Gamboa Table 2 の最適値 |
| large | 3.0 | 6.5 | |

出力（`cad/`、単位 mm）は各変種について DXF・流体 STL・溝板 STL。
名前は `{形状}_{大きさ}_{役割}`（例: `side4fwd_large_fluid.stl`）。
全モデルで「各辺がちょうど 2 枚の三角形に共有される」ことを確認している。

注意
----
流路 2 は出口区間を閉じ壁に置き換えた周期構成なので、**Gamboa の単発バルブとは
接合部の位相が違う**（逆流がループへ入るのに 132 度曲がる必要がある。
2026-09-06 の知見）。Gamboa の Di を再現する形状ではなく、
「直線流路 + 側面のこぶ」という構成そのものを試すための形状である。
"""
import argparse, json, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import gamboa as G
import gamboa_full as GF
import mesh as M
import solid3d as S3

CAD = os.path.join(HERE, "..", "cad")
SIZES = {"small": 1.9, "std": 2.35, "large": 3.0}


def emit(tag, poly, depth, wall_t, margin, port_d, pad_r, note):
    """DXF・流体 STL・溝板 STL を書き、検査結果を返す。"""
    # 端に丸い溜まりを付けてポートを収める（流路幅 1 mm に φ2 は入らない）
    b = poly.bounds
    p_in = (b[0] + 0.5 * pad_r, 0.0)
    p_out = (b[2] - 0.5 * pad_r, 0.0)
    ys = [c[1] for c in np.asarray(poly.exterior.coords)]
    poly = S3.port_pads(poly, [p_in, p_out], pad_r)

    rings = M.polygon_rings(M.clean_polygon(poly))
    f_dxf = os.path.join(CAD, f"{tag}_outline.dxf")
    npts = M.write_dxf(f_dxf, rings, note=f"{tag} ({note})")

    ff = S3.fluid_solid(poly, depth)
    bad1, ne1 = M.check_manifold(ff)
    n1 = M.write_stl(os.path.join(CAD, f"{tag}_fluid.stl"), ff)
    vol = M.mesh_volume(ff)

    pf, ext = S3.grooved_plate(poly, depth, wall_t, margin, [p_in, p_out],
                               port_d)
    bad2, ne2 = M.check_manifold(pf)
    n2 = M.write_stl(os.path.join(CAD, f"{tag}_plate.stl"), pf)

    ok = (bad1 == 0 and bad2 == 0)
    print(f"  {tag:26s} 外形 {b[2]-b[0]:6.2f} x {b[3]-b[1]:5.2f} mm  "
          f"流体 {n1:6d} 三角形 {vol:8.2f} mm^3  板 {n2:6d} 三角形  "
          f"{'OK' if ok else f'**非多様体 {bad1}/{bad2}**'}")
    return dict(tag=tag, note=note, dxf_points=npts,
                bounds=[float(v) for v in b],
                length=float(b[2] - b[0]), height=float(b[3] - b[1]),
                fluid_triangles=n1, fluid_volume_mm3=float(vol),
                plate_triangles=n2, plate_extent=[float(v) for v in ext],
                watertight=bool(ok), poly=poly)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", type=int, default=4, help="流路 1 の段数")
    ap.add_argument("--loops", type=int, default=4, help="流路 2 のこぶの数")
    ap.add_argument("--gap", type=float, default=1.0, help="こぶの間隔 [w_v]")
    ap.add_argument("--wv-mm", type=float, default=1.0)
    ap.add_argument("--depth-mm", type=float, default=None)
    ap.add_argument("--wall-mm", type=float, default=1.0)
    ap.add_argument("--margin-mm", type=float, default=2.0)
    ap.add_argument("--port-mm", type=float, default=2.0)
    ap.add_argument("--arc-pts", type=int, default=180)
    a = ap.parse_args()
    wv = a.wv_mm
    depth = a.depth_mm if a.depth_mm else wv
    pad_r = max(0.8 * a.port_mm, 0.75 * wv)
    os.makedirs(CAD, exist_ok=True)

    print("=" * 96)
    print(f" 流路の変種   w_v = {wv} mm, 深さ = {depth} mm, "
          f"ポート径 {a.port_mm} mm（端に半径 {pad_r:.2f} mm の溜まり）")
    print("=" * 96)

    recs, shapes = [], []
    print(f"\n流路 1: 段々（{a.stages} 段、1 段おきに反転、平均方向は水平）")
    for name, R in SIZES.items():
        poly, info = GF.chain(a.stages, w=wv, gap=0.0, arc_pts=a.arc_pts, R=R)
        if info["overlaps"]:
            print(f"  ** {name}: 段どうしが {len(info['overlaps'])} 箇所重なる。"
                  f"流路が短絡するので使わないこと **")
        r = emit(f"stair{a.stages}_{name}", poly, depth, a.wall_mm,
                 a.margin_mm, a.port_mm, pad_r,
                 f"staircase {a.stages} stages, R={R} w_v")
        r.update(kind="staircase", size=name, R=R, stages=a.stages,
                 overlaps=len(info["overlaps"]))
        shapes.append((f"stair {name} (R={R})", r.pop("poly")))
        recs.append(r)

    print(f"\n流路 2: 直線流路の側面にこぶ（{a.loops} 個、gap = {a.gap} w_v）")
    for orient, mirror in (("fwd", False), ("rev", True)):
        for name, R in SIZES.items():
            poly, info = G.straight_with_loops(a.loops, R=R, gap=a.gap, w=wv,
                                               mirror=mirror,
                                               arc_pts=a.arc_pts)
            r = emit(f"side{a.loops}{orient}_{name}", poly, depth, a.wall_mm,
                     a.margin_mm, a.port_mm, pad_r,
                     f"straight+{a.loops} loops, {orient}, R={R} w_v")
            r.update(kind="side_loops", size=name, R=R, orientation=orient,
                     loops=a.loops, period=info["period"],
                     loop_top=info["loop_top"])
            shapes.append((f"side {orient} {name} (R={R})", r.pop("poly")))
            recs.append(r)

    ok = all(r["watertight"] for r in recs)
    print("\n" + ("すべて閉じたメッシュ。CAD 取り込み可。" if ok
                  else "**非多様体のモデルがある**"))
    j = os.path.join(HERE, "..", "results", "cad", "cad_variants.json")
    json.dump(dict(wv_mm=wv, depth_mm=depth, wall_mm=a.wall_mm,
                   port_mm=a.port_mm, pad_r_mm=pad_r, sizes=SIZES,
                   stages=a.stages, loops=a.loops, gap=a.gap,
                   models=recs, all_watertight=bool(ok)),
              open(j, "w"), indent=1)
    print(f"saved {os.path.normpath(j)}")

    # --- 一覧図 ---
    n = len(shapes)
    fig, axes = plt.subplots(n, 1, figsize=(13, 1.9 * n))
    for ax, (lab, poly) in zip(axes, shapes):
        P = np.asarray(poly.exterior.coords)
        ax.fill(P[:, 0], P[:, 1], color="#cfe3e4")
        ax.plot(P[:, 0], P[:, 1], "-", color="#0e7c7b", lw=0.9)
        for r in poly.interiors:
            Q = np.asarray(r.coords)
            ax.fill(Q[:, 0], Q[:, 1], color="white")
            ax.plot(Q[:, 0], Q[:, 1], "-", color="#0e7c7b", lw=0.9)
        ax.set_aspect("equal")
        ax.set_xlim(-3, 95)
        ax.set_ylabel(lab, fontsize=7, rotation=0, ha="right", va="center")
        ax.set_yticks([])
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("x [mm]")
    fig.suptitle("channel variants (top 3: staircase, bottom 6: straight + side bumps)",
                 fontsize=10)
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "geometry", "variants.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
