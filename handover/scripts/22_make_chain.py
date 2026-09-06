#!/usr/bin/env python3
"""
22_make_chain.py -- バルブを直列に繋いだ流路（多段）を作る
==========================================================

Gamboa の 1 段は入口から alpha だけ折れて出るので、そのまま繋ぐと段ごとに
曲がって螺旋になる。**1 段おきに上下反転**して折れを相殺し、さらに入口の向きを
+alpha_axis/2（= 24.05 度）から始めると、向きが ±24.05 度を往復して
**全体の平均方向が水平**になる。Tesla の原特許や 7 月形状と同じジグザグ中心線。

出力（`cad/`、単位 mm）

| ファイル | 内容 |
|---|---|
| `gamboa_chain{n}_2d.dxf` | 多段流路の 2D 輪郭（外形 + 島 n 個） |
| `gamboa_chain{n}_fluid_3d.stl` | 押し出した流体体積 |
| `gamboa_chain{n}_plate_3d.stl` | 溝を彫った板（裏面に入口・出口ポート） |
| `figures/gamboa_chain{n}.png` | 形状と接続の確認図 |

注意
----
- **入口面と出口面は水平から 24.05 度傾いている**（平均方向は水平）。
  CFD にかけるならプレナムを付けるか、助走の直管を足して向きを揃えること
- 段どうしが重なっていないか（流路が短絡しないか）を毎回検査している
- gap を入れると接続部 1 本ぶんだけ出口が横にずれる（gap=1 で 0.41 w_v）。
  gap=0 なら出口の y は厳密に 0
"""
import argparse, json, os, sys

import numpy as np
from shapely.geometry import box, Point
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import gamboa_full as GF
import mesh as M

CAD = os.path.join(HERE, "..", "cad")


def fluid_solid(poly, depth):
    poly = M.clean_polygon(poly)
    faces = _cap(poly, 0.0, False) + _cap(poly, depth, True)
    faces += _wall(_closed(poly.exterior, True), 0.0, depth, True)
    for r in poly.interiors:
        faces += _wall(_closed(r, False), 0.0, depth, True)
    return faces


def _cap(poly, z, up):
    if poly.geom_type == "MultiPolygon":
        out = []
        for g in poly.geoms:
            out += _cap(g, z, up)
        return out
    out = []
    for t in M.triangulate_polygon(poly):
        v = [(t[i][0], t[i][1], z) for i in range(3)]
        a, b, c = t
        cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        out.append(v[::-1] if (cr > 0) != up else v)
    return out


def _wall(P, z0, z1, outward):
    out = []
    for a, b in zip(P[:-1], P[1:]):
        q = [(a[0], a[1], z0), (b[0], b[1], z0),
             (b[0], b[1], z1), (a[0], a[1], z1)]
        t = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
        out += t if outward else [u[::-1] for u in t]
    return out


def _closed(ring, ccw):
    P = np.asarray(ring.coords, float)
    if not np.allclose(P[0], P[-1]):
        P = np.vstack([P, P[:1]])
    a = 0.5 * np.sum(P[:-1, 0] * P[1:, 1] - P[1:, 0] * P[:-1, 1])
    return P[::-1] if (a > 0) != ccw else P


def grooved_plate(poly, depth, wall_t, margin, ports, port_d):
    """流路を溝として彫った板。入口・出口は裏面のポート穴。

    面はすべて `src/mesh.py` の耳刈りで三角形化する（共線点を落としてから
    分割するので、多段の入り組んだ形でも刈り切れる）。短冊に切る必要はない。
    切ると切り口の点が外周に増えて、側壁と辺が食い違う。
    """
    poly = M.clean_polygon(poly)
    x0, y0, x1, y1 = poly.bounds
    rect = box(x0 - margin, y0 - margin, x1 + margin, y1 + margin)
    holes = unary_union([Point(*c).buffer(0.5 * port_d, resolution=32)
                         for c in ports])
    faces = _cap(rect.difference(holes), -wall_t, False)
    faces += _cap(rect.difference(poly), depth, True)
    faces += _cap(poly.difference(holes), 0.0, True)
    faces += _wall(_closed(poly.exterior, True), 0.0, depth, False)
    for r in poly.interiors:
        faces += _wall(_closed(r, False), 0.0, depth, False)
    for h in (holes.geoms if holes.geom_type == "MultiPolygon" else [holes]):
        faces += _wall(_closed(h.exterior, True), -wall_t, 0.0, False)
    faces += _wall(_closed(rect.exterior, True), -wall_t, depth, True)
    return faces, rect.bounds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--stages", type=int, default=4)
    ap.add_argument("--gap", type=float, default=0.0,
                    help="段間の直線区間 [w_v]。Porwal は 1 D_H を使う")
    ap.add_argument("--wv-mm", type=float, default=1.0)
    ap.add_argument("--depth-mm", type=float, default=None)
    ap.add_argument("--wall-mm", type=float, default=1.0)
    ap.add_argument("--margin-mm", type=float, default=2.0)
    ap.add_argument("--port-mm", type=float, default=2.0)
    ap.add_argument("--arc-pts", type=int, default=180)
    ap.add_argument("--no-level", action="store_true",
                    help="平均方向を水平に揃えない（鎖が 24 度下がる）")
    a = ap.parse_args()
    wv = a.wv_mm
    depth = a.depth_mm if a.depth_mm else wv
    os.makedirs(CAD, exist_ok=True)

    poly, info = GF.chain(a.stages, w=wv, gap=a.gap, arc_pts=a.arc_pts,
                          level=not a.no_level)
    b = info["bounds"]
    print("=" * 84)
    print(f" Gamboa optimized を {a.stages} 段連結   w_v = {wv} mm, "
          f"gap = {a.gap} w_v, 深さ = {depth} mm")
    print("=" * 84)
    print(f"  外形 {b[2]-b[0]:.2f} x {b[3]-b[1]:.2f} mm   "
          f"流体断面積 {info['total_area']:.2f} mm^2   島 {info['n_islands']} 個")
    print(f"  入口 ({-0.5*wv:.2f}, 0.00) 向き "
          f"{np.degrees(np.arctan2(*np.array(info['stages'][0]['inlet_dir'])[::-1])):+.2f} 度")
    print(f"  出口 ({info['outlet'][0]:.2f}, {info['outlet'][1]:.2f}) 向き "
          f"{np.degrees(np.arctan2(info['outlet_dir'][1], info['outlet_dir'][0])):+.2f} 度")
    print(f"  段どうしの重なり: {len(info['overlaps'])} 件 "
          f"{'（流路が短絡していない）' if not info['overlaps'] else '**要確認**'}")

    # 流路幅 1 w_v にポート（既定 φ2 mm）は入らないので、両端に丸い溜まりを付け、
    # その中心をポートにする。溜まりが無いと穴が流路からはみ出して
    # 側壁と穴の壁が食い違い、非多様体になる。
    pad_r = max(0.8 * a.port_mm, 0.75 * wv)
    d0 = np.array(info["stages"][0]["inlet_dir"])
    d1 = np.array(info["outlet_dir"])
    p_in = np.array([-0.5 * wv, 0.0]) + 0.5 * pad_r * d0
    p_out = np.array(info["outlet"]) - 0.5 * pad_r * d1
    poly = unary_union([poly,
                        Point(*p_in).buffer(pad_r, resolution=48),
                        Point(*p_out).buffer(pad_r, resolution=48)])
    print(f"  入口・出口に半径 {pad_r:.2f} mm の溜まりを付けた"
          f"（ポート径 {a.port_mm} mm を収めるため）")
    tag = f"chain{a.stages}"
    rings = M.polygon_rings(poly)
    f_dxf = os.path.join(CAD, f"gamboa_{tag}_2d.dxf")
    npts = M.write_dxf(f_dxf, rings, note=f"gamboa {a.stages}-stage, w_v={wv}mm")
    print(f"\n  {os.path.basename(f_dxf):34s} {npts:6d} 点（外形 + 島 "
          f"{len(rings)-1} 個）")

    faces = fluid_solid(poly, depth)
    bad, ne = M.check_manifold(faces)
    n = M.write_stl(os.path.join(CAD, f"gamboa_{tag}_fluid_3d.stl"), faces)
    vol = M.mesh_volume(faces)
    print(f"  {'gamboa_'+tag+'_fluid_3d.stl':34s} {n:6d} 三角形  "
          f"体積 {vol:9.3f} mm^3  非多様体辺 {bad}/{ne}")

    pf, ext = grooved_plate(poly, depth, a.wall_mm, a.margin_mm,
                            [tuple(p_in), tuple(p_out)], a.port_mm)
    bad2, ne2 = M.check_manifold(pf)
    n2 = M.write_stl(os.path.join(CAD, f"gamboa_{tag}_plate_3d.stl"), pf)
    print(f"  {'gamboa_'+tag+'_plate_3d.stl':34s} {n2:6d} 三角形  "
          f"体積 {M.mesh_volume(pf):9.3f} mm^3  非多様体辺 {bad2}/{ne2}")
    print(f"     板 {ext[2]-ext[0]:.1f} x {ext[3]-ext[1]:.1f} x "
          f"{a.wall_mm+depth:.1f} mm、ポート径 {a.port_mm} mm")

    ok = (bad == 0 and bad2 == 0 and not info["overlaps"])
    print("\n  " + ("すべて閉じたメッシュ。CAD 取り込み可。" if ok
                    else "**要確認（非多様体辺または段の重なりがある）**"))

    j = os.path.join(HERE, "..", "results", f"cad_{tag}.json")
    json.dump(dict(stages=a.stages, gap=a.gap, wv_mm=wv, depth_mm=depth,
                   bounds=b, area_mm2=info["total_area"],
                   n_islands=info["n_islands"], outlet=info["outlet"],
                   outlet_dir=info["outlet_dir"], overlaps=info["overlaps"],
                   fluid_triangles=n, fluid_volume_mm3=vol,
                   plate_triangles=n2, watertight=bool(ok)),
              open(j, "w"), indent=1)
    print(f"  saved {os.path.normpath(j)}")

    # --- 図 ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(13, 4.2))
    P = np.asarray(poly.exterior.coords)
    ax.fill(P[:, 0], P[:, 1], color="#cfe3e4", zorder=1)
    ax.plot(P[:, 0], P[:, 1], "-", color="#0e7c7b", lw=1.0, zorder=3)
    for r in poly.interiors:
        Q = np.asarray(r.coords)
        ax.fill(Q[:, 0], Q[:, 1], color="white", zorder=2)
        ax.plot(Q[:, 0], Q[:, 1], "-", color="#0e7c7b", lw=1.0, zorder=3)
    for s in info["stages"]:
        ax.plot(*s["inlet"], "o", color="#c0521c", ms=4, zorder=4)
    ax.plot(*info["outlet"], "s", color="#c0521c", ms=5, zorder=4)
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_xlabel("x [mm]")
    ax.set_title(f"Gamboa optimized x {a.stages} stages (gap = {a.gap} w_v)  "
                 f"— dots: stage inlets, square: outlet")
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", f"gamboa_{tag}.png")
    fig.savefig(png, dpi=130)
    print(f"  saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
