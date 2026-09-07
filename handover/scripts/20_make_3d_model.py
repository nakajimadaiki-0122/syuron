#!/usr/bin/env python3
"""
20_make_3d_model.py -- Gamboa バルブの 3D モデル（CAD / 造形用）を作る
=====================================================================

出力（`cad/`、単位 mm、STL は binary）

| ファイル | 内容 |
|---|---|
| `valve1_std_fluid.stl` | **流体体積**。CFD の領域そのもの。固体から引く用 |
| `valve1_std_plate.stl` | **溝を彫った板**。流路は上面に開き、入口・出口は**板の裏からのポート穴** |
| `valve1_std_lid.stl` | 板と同じ外形の蓋（平板）。貼り合わせると流路が閉じる |
| `valve1bare_std_fluid.stl` | プレナムを外した流体体積 |

寸法の既定値
  w_v = 1.0 mm（流路幅）、深さ = 1.0 mm（正方形断面）、
  溝の下の肉厚 = 1.0 mm、y 方向の縁 = 2.0 mm、蓋の厚み = 1.0 mm

注意
----
- 島（涙滴形）は 2D CFD では自由な島だが、**実物では溝の底から立ち上がる柱**になる。
  板モデルではそうなっている
- 入口・出口はプレナム中心の**裏面ポート**（既定 φ2 mm）。継手を下から差す想定
- 三角形分割は `src/mesh.py` の耳刈り（境界の辺を保つ）。
  **各辺がちょうど 2 枚に共有されること**を毎回確認している
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


def cap(poly, z, up):
    """z 一定の面。up=True なら法線 +z。MultiPolygon も受ける。"""
    if poly.geom_type == "MultiPolygon":
        out = []
        for g in poly.geoms:
            out += cap(g, z, up)
        return out
    out = []
    for t in M.triangulate_polygon(poly):
        v = [(t[i][0], t[i][1], z) for i in range(3)]
        a, b, c = t
        cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if (cr > 0) != up:
            v = v[::-1]
        out.append(v)
    return out


def wall(P, z0, z1, outward, skip=None):
    """閉じた折れ線 P（末尾に始点を含む）を z0..z1 へ押し出す。"""
    out = []
    for a, b in zip(P[:-1], P[1:]):
        if skip is not None and skip(a, b):
            continue
        q = [(a[0], a[1], z0), (b[0], b[1], z0),
             (b[0], b[1], z1), (a[0], a[1], z1)]
        t = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
        if not outward:
            t = [u[::-1] for u in t]
        out += t
    return out


def closed(ring, ccw):
    """リング座標を、指定の向き（ccw=True で反時計回り）の閉じた列にする。"""
    P = np.asarray(ring.coords, float)
    if not np.allclose(P[0], P[-1]):
        P = np.vstack([P, P[:1]])
    a = 0.5 * np.sum(P[:-1, 0] * P[1:, 1] - P[1:, 0] * P[:-1, 1])
    if (a > 0) != ccw:
        P = P[::-1]
    return P


def fluid_solid(poly, depth):
    """流体多角形を深さ方向に押し出した閉じた立体。"""
    poly = M.clean_polygon(poly)
    faces = cap(poly, 0.0, up=False) + cap(poly, depth, up=True)
    faces += wall(closed(poly.exterior, True), 0.0, depth, outward=True)
    for r in poly.interiors:
        faces += wall(closed(r, False), 0.0, depth, outward=True)
    return faces


def grooved_plate(poly, depth, wall_t, margin_xy, ports, port_d):
    """流路を溝として彫った板。

    入口・出口は**板の裏から垂直に開けたポート穴**にする。
    板の端で流路を開口させる案はやめた。プレナムの直線壁は片方が傾いており
    （Fig.2 どおり）、板の外形と一致しないため、端面が T 字接合になって
    非多様体メッシュになる。裏からのポートなら接合が素直で、実物としても
    継手を下から差せるので扱いやすい。

    ports  : [(x, y), ...] ポート中心
    port_d : ポート直径
    """
    poly = M.clean_polygon(poly)
    x0, y0, x1, y1 = poly.bounds
    X0, X1 = x0 - margin_xy, x1 + margin_xy
    Y0, Y1 = y0 - margin_xy, y1 + margin_xy
    rect = box(X0, Y0, X1, Y1)
    holes = [Point(*c).buffer(0.5 * port_d, resolution=32) for c in ports]
    hole_union = unary_union(holes) if holes else None

    faces = []
    # 板の底（法線 -z、ポート穴を抜く）
    faces += cap(rect.difference(hole_union) if hole_union else rect,
                 -wall_t, up=False)
    # 上面（法線 +z、流路を抜く）
    faces += cap(rect.difference(poly), depth, up=True)
    # 溝の底（法線 +z、ポート穴を抜く）
    floor = poly.difference(hole_union) if hole_union else poly
    faces += cap(floor, 0.0, up=True)
    # 溝の側壁（法線は溝の内側 = 流体側）
    faces += wall(closed(poly.exterior, True), 0.0, depth, outward=False)
    for r in poly.interiors:
        faces += wall(closed(r, False), 0.0, depth, outward=False)
    # ポートの内壁（法線は穴の内側）
    for h in holes:
        faces += wall(closed(h.exterior, True), -wall_t, 0.0, outward=False)
    # 板の外側面
    faces += wall(closed(rect.exterior, True), -wall_t, depth, outward=True)
    return faces, (X0, X1, Y0, Y1)


def plain_plate(x0, x1, Y0, Y1, thick):
    r = box(x0, Y0, x1, Y1)
    faces = cap(r, 0.0, up=False) + cap(r, thick, up=True)
    faces += wall(closed(r.exterior, True), 0.0, thick, outward=True)
    return faces


def emit(name, faces, note=""):
    bad, ne = M.check_manifold(faces)
    n = M.write_stl(os.path.join(CAD, name), faces)
    v = M.mesh_volume(faces)
    print(f"  {name:32s} {n:6d} 三角形  体積 {v:9.3f} mm^3  "
          f"非多様体辺 {bad}/{ne}  {note}")
    return dict(file=name, triangles=n, volume_mm3=v, nonmanifold=bad,
                edges=ne)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="optimized")
    ap.add_argument("--wv-mm", type=float, default=1.0)
    ap.add_argument("--depth-mm", type=float, default=None,
                    help="流路の深さ [mm]。既定は w_v（正方形断面）")
    ap.add_argument("--wall-mm", type=float, default=1.0, help="溝の下の肉厚")
    ap.add_argument("--margin-mm", type=float, default=2.0, help="板の縁")
    ap.add_argument("--port-mm", type=float, default=2.0,
                    help="入口・出口ポートの直径（板の裏から開ける）")
    ap.add_argument("--lid-mm", type=float, default=1.0, help="蓋の厚み")
    ap.add_argument("--arc-pts", type=int, default=180)
    a = ap.parse_args()
    wv = a.wv_mm
    depth = a.depth_mm if a.depth_mm else wv
    os.makedirs(CAD, exist_ok=True)

    poly, info = GF.build(a.which, w=wv, plenum_out="faithful",
                          arc_pts=a.arc_pts)
    poly_v, _ = GF.build(a.which, w=wv, plenum_out="faithful",
                         arc_pts=a.arc_pts, include_plenums=False)

    print("=" * 84)
    print(f" Gamboa {a.which} の 3D モデル   w_v = {wv} mm, 深さ = {depth} mm, "
          f"肉厚 = {a.wall_mm} mm, 縁 = {a.margin_mm} mm")
    print("=" * 84)
    print(f"  流体の断面積 {poly.area:.3f} mm^2  "
          f"外形 {poly.bounds[2]-poly.bounds[0]:.2f} x "
          f"{poly.bounds[3]-poly.bounds[1]:.2f} mm")

    rec = []
    rec.append(emit("valve1_std_fluid.stl", fluid_solid(poly, depth),
                    "CFD 領域そのもの"))
    rec.append(emit("valve1bare_std_fluid.stl",
                    fluid_solid(poly_v, depth), "プレナム無し"))
    # ポートはプレナムの**内側**へ寄せる。中心は直線壁の上にあるので、
    # そこへ穴を開けると半分が板の外に出て閉じたメッシュにならない。
    Rp = info["Rp"]
    al = np.radians(info["alpha_deg"])
    e = np.array([np.sin(al), -np.cos(al)])          # 出口区間の向き
    Lc = np.array(info["plenum_left_centre"])
    Rc = np.array(info["plenum_right_centre"])
    ports = [tuple(Lc + np.array([0.5 * Rp, 0.0])), tuple(Rc - 0.5 * Rp * e)]
    fp, ext = grooved_plate(poly, depth, a.wall_mm, a.margin_mm, ports,
                            a.port_mm)
    rec.append(emit("valve1_std_plate.stl", fp,
                    f"板 {ext[1]-ext[0]:.1f} x {ext[3]-ext[2]:.1f} x "
                    f"{a.wall_mm+depth:.1f} mm, ポート径 {a.port_mm} mm"))
    rec.append(emit("valve1_std_lid.stl", plain_plate(*ext, a.lid_mm), "蓋"))

    ok = all(r["nonmanifold"] == 0 for r in rec)
    print("\n  " + ("すべて閉じたメッシュ（非多様体辺なし）。"
                    "3D プリント・CAD 取り込み可。" if ok
                    else "**非多様体辺がある。修正が要る。**"))
    j = os.path.join(HERE, "..", "results", "cad", "cad_models.json")
    json.dump(dict(wv_mm=wv, depth_mm=depth, wall_mm=a.wall_mm,
                   margin_mm=a.margin_mm, port_mm=a.port_mm,
                   ports=[list(p) for p in ports], plate_extent=list(ext),
                   models=rec, all_watertight=bool(ok)),
              open(j, "w"), indent=1)
    print(f"  saved {os.path.normpath(j)}")


if __name__ == "__main__":
    main()
