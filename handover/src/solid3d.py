"""
solid3d.py -- 2 次元の流体多角形から 3 次元の立体（STL）を組む
==============================================================

`20_make_3d_model.py` と `22_make_chain.py` で同じ処理を書いていたのを集約した。

  fluid_solid  : 流体多角形を深さ方向へ押し出した立体（CFD 領域そのもの）
  grooved_plate: 流路を溝として彫った板。入口・出口は**裏面のポート穴**
  plain_plate  : 蓋（平板）

**押し出す前に必ず `mesh.clean_polygon` を通すこと。** 三角形分割は共線点を
落としてから行うので、元の多角形のまま側壁を貼ると上下面と側壁で辺が食い違い、
非多様体メッシュになる（2026-09-07 に実測）。ここでは各関数の冒頭で通している。
"""
from __future__ import annotations
import numpy as np
from shapely.geometry import box, Point
from shapely.ops import unary_union

import mesh as M


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
        out.append(v[::-1] if (cr > 0) != up else v)
    return out


def wall(P, z0, z1, outward):
    """閉じた折れ線 P（末尾に始点を含む）を z0..z1 へ押し出す。"""
    out = []
    for a, b in zip(P[:-1], P[1:]):
        q = [(a[0], a[1], z0), (b[0], b[1], z0),
             (b[0], b[1], z1), (a[0], a[1], z1)]
        t = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
        out += t if outward else [u[::-1] for u in t]
    return out


def closed(ring, ccw):
    """リング座標を、指定の向きの閉じた列にする。"""
    P = np.asarray(ring.coords, float)
    if not np.allclose(P[0], P[-1]):
        P = np.vstack([P, P[:1]])
    a = 0.5 * np.sum(P[:-1, 0] * P[1:, 1] - P[1:, 0] * P[:-1, 1])
    return P[::-1] if (a > 0) != ccw else P


def fluid_solid(poly, depth):
    """流体多角形を深さ方向へ押し出した閉じた立体。"""
    poly = M.clean_polygon(poly)
    faces = cap(poly, 0.0, False) + cap(poly, depth, True)
    faces += wall(closed(poly.exterior, True), 0.0, depth, True)
    for r in poly.interiors:
        faces += wall(closed(r, False), 0.0, depth, True)
    return faces


def grooved_plate(poly, depth, wall_t, margin, ports, port_d):
    """流路を溝として彫った板。入口・出口は裏面のポート穴。

    ports  : [(x, y), ...] ポート中心。**流路の内側**に置くこと
             （流路幅より穴が大きいとはみ出して非多様体になる）
    """
    poly = M.clean_polygon(poly)
    x0, y0, x1, y1 = poly.bounds
    rect = box(x0 - margin, y0 - margin, x1 + margin, y1 + margin)
    holes = unary_union([Point(*c).buffer(0.5 * port_d, resolution=32)
                         for c in ports])
    faces = cap(rect.difference(holes), -wall_t, False)
    faces += cap(rect.difference(poly), depth, True)
    faces += cap(poly.difference(holes), 0.0, True)
    faces += wall(closed(poly.exterior, True), 0.0, depth, False)
    for r in poly.interiors:
        faces += wall(closed(r, False), 0.0, depth, False)
    for h in (holes.geoms if holes.geom_type == "MultiPolygon" else [holes]):
        faces += wall(closed(h.exterior, True), -wall_t, 0.0, False)
    faces += wall(closed(rect.exterior, True), -wall_t, depth, True)
    return faces, rect.bounds


def plain_plate(x0, x1, y0, y1, thick):
    """蓋（平板）。"""
    r = box(x0, y0, x1, y1)
    faces = cap(r, 0.0, False) + cap(r, thick, True)
    faces += wall(closed(r.exterior, True), 0.0, thick, True)
    return faces


def port_pads(poly, pts, pad_r, resolution=48):
    """流路の端に丸い溜まりを付ける（細い流路にポートを収めるため）。"""
    return unary_union([poly] + [Point(*p).buffer(pad_r, resolution=resolution)
                                 for p in pts])
