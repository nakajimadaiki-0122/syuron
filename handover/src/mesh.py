"""
mesh.py -- 穴つき多角形の厳密な三角形分割と STL 出力
====================================================

なぜ自前で書くか
----------------
`shapely.ops.triangulate` は**非拘束** Delaunay なので、境界の辺が三角形の辺として
現れる保証がない。押し出して側壁を貼ると、上下面と側壁で辺が食い違い、
**非多様体（辺を共有する三角形が 2 枚でない）**メッシュになる。
2026-09-06 に実測: 流体形状で 42 / 11228 辺、板モデルで 191 / 10874 辺が非多様体だった。

ここでは

  1. 穴（島）を橋渡しして単一の単純多角形にまとめる
  2. 耳刈り（ear clipping）で三角形化する

の順で処理する。耳刈りは境界の辺を必ず保つので、押し出しても閉じたメッシュになる。

凸性の判定は符号付き面積、内包判定は重心座標で行う（numpy でまとめて評価）。
"""
from __future__ import annotations
import struct
from collections import Counter

import numpy as np
from shapely.geometry import LineString, Polygon


def _ring(coords):
    """閉じた座標列から末尾の重複点を落とした (n,2) を返す。"""
    P = np.asarray(coords, float)
    if len(P) > 1 and np.allclose(P[0], P[-1]):
        P = P[:-1]
    return P


def _signed_area(P):
    return 0.5 * float(np.sum(P[:, 0] * np.roll(P[:, 1], -1)
                              - np.roll(P[:, 0], -1) * P[:, 1]))


def bridge_holes(poly: Polygon) -> np.ndarray:
    """穴つき多角形を、橋渡しで単一の単純多角形（CCW）にする。

    穴の最右点 M から、外側リングの頂点のうち **M から見通せる**最も近いものへ
    橋を架ける（往復するので面積は変わらない）。

    見通し判定は「元の多角形の内部に収まる」だけでは足りない。
    **すでに架けた橋と交差しないこと**も要る（交差すると耳刈りが破綻して
    面積が合わなくなる。2026-09-06 に 4 段連結で +29.6 % ずれた）。
    穴は右のものから順に処理する。
    """
    outer = _ring(poly.exterior.coords)
    if _signed_area(outer) < 0:
        outer = outer[::-1]
    holes = sorted(poly.interiors, key=lambda r: -r.bounds[2])
    for ring in holes:
        hole = _ring(ring.coords)
        if _signed_area(hole) > 0:            # 穴は外側と逆回りにする
            hole = hole[::-1]
        m = int(np.argmax(hole[:, 0]))
        M = hole[m]
        cur = LineString(np.vstack([outer, outer[:1]]))
        d = np.hypot(*(outer - M).T)
        for v in np.argsort(d):
            seg = LineString([tuple(M), tuple(outer[v])])
            if poly.covers(seg) and not cur.crosses(seg):
                break
        else:
            raise ValueError("穴への橋を架けられない")
        hole_cycle = np.vstack([hole[m:], hole[:m], hole[m][None]])
        outer = np.vstack([outer[:v + 1], hole_cycle, outer[v:]])
    return outer


def drop_collinear(P, tol=1e-9):
    """ほぼ一直線に並ぶ中間点を落とす。

    円弧の細分や shapely の交差で共線点が大量に入る。共線点が残ると
    「耳」の判定で辺の上に載る点として扱われ、刈れる耳が無くなって
    耳刈りが止まる（2026-09-06 に 800 点の単純多角形で破綻した）。
    """
    n = len(P)
    keep = np.ones(n, bool)
    for i in range(n):
        a, b, c = P[(i - 1) % n], P[i], P[(i + 1) % n]
        cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        L = max(np.hypot(*(b - a)), np.hypot(*(c - b)), 1e-30)
        if abs(cr) <= tol * L or (np.hypot(*(b - a)) < tol):
            keep[i] = False
    return P[keep] if keep.sum() >= 3 else P


def clean_polygon(poly):
    """全リングから共線点・重複点を落とした多角形を返す。

    **押し出す前に必ず通すこと。** 三角形分割は内部で共線点を落とすので、
    元の多角形のまま側壁を貼ると、上下面と側壁で辺が食い違って非多様体になる。
    """
    def fix(P):
        # 点を落とすと隣が新たに共線になることがあるので、変化しなくなるまで回す
        for _ in range(8):
            Q = drop_collinear(P)
            if len(Q) == len(P):
                return Q
            P = Q
        return P

    ext = fix(_ring(poly.exterior.coords))
    ints = [fix(_ring(r.coords)) for r in poly.interiors]
    return Polygon(ext, [r for r in ints if len(r) >= 3])


def earclip(P: np.ndarray) -> list:
    """単純多角形（CCW、末尾の重複なし）を耳刈りで三角形化する。

    橋渡しで作った多角形は橋の 2 辺が重なるので、内包判定は**厳密な内部**
    （重心座標が正）で行う。境界上の点（橋の重複頂点など）で耳を潰さないため。
    戻り値は頂点添字の 3 つ組のリスト。
    """
    n = len(P)
    idx = list(range(n))
    tris = []
    while len(idx) > 3:
        m = len(idx)
        # 凹頂点（reflex）だけを内包判定の対象にする
        A = P[[idx[i - 1] for i in range(m)]]
        B = P[idx]
        C = P[[idx[(i + 1) % m] for i in range(m)]]
        cr = ((B[:, 0] - A[:, 0]) * (C[:, 1] - A[:, 1])
              - (B[:, 1] - A[:, 1]) * (C[:, 0] - A[:, 0]))
        reflex = [idx[i] for i in range(m) if cr[i] <= 0]
        Q = P[reflex] if reflex else np.zeros((0, 2))
        found = -1
        for i in range(m):
            if cr[i] <= 1e-14:
                continue
            ia, ib, ic = idx[i - 1], idx[i], idx[(i + 1) % m]
            a, b, c = P[ia], P[ib], P[ic]
            if len(Q):
                d = cr[i]
                # 橋渡しで重複した頂点（耳の 3 点と同じ座標）は判定から外す。
                # 外さないと、橋のスリットをまたぐ耳が永久に刈れなくなる。
                dup = ((np.hypot(Q[:, 0] - a[0], Q[:, 1] - a[1]) < 1e-12)
                       | (np.hypot(Q[:, 0] - b[0], Q[:, 1] - b[1]) < 1e-12)
                       | (np.hypot(Q[:, 0] - c[0], Q[:, 1] - c[1]) < 1e-12))
                w0 = ((b[0] - Q[:, 0]) * (c[1] - Q[:, 1])
                      - (b[1] - Q[:, 1]) * (c[0] - Q[:, 0])) / d
                w1 = ((c[0] - Q[:, 0]) * (a[1] - Q[:, 1])
                      - (c[1] - Q[:, 1]) * (a[0] - Q[:, 0])) / d
                w2 = 1.0 - w0 - w1
                inside = (w0 > 1e-12) & (w1 > 1e-12) & (w2 > 1e-12) & ~dup
                if np.any(inside):
                    continue
            found = i
            break
        if found < 0:                     # 耳が見つからない: 凸頂点を強制的に刈る
            cand = np.argmax(cr)
            if cr[cand] <= 0:
                break
            found = int(cand)
        i = found
        tris.append((idx[i - 1], idx[i], idx[(i + 1) % len(idx)]))
        del idx[i]
    if len(idx) == 3:
        tris.append(tuple(idx))
    return tris


def _flatten(g):
    if g.is_empty:
        return []
    if g.geom_type == "Polygon":
        return [g]
    return [p for p in getattr(g, "geoms", []) if p.geom_type == "Polygon"]


def split_multi_hole(poly, depth=0):
    """穴が 2 つ以上ある多角形を、**穴のあいだで縦に切って** 1 穴以下にする。

    橋渡しは穴 1 つなら安定するが、複数の橋が交差すると耳刈りが破綻する
    （2026-09-06 に 4 段連結で面積が +29.6 % ずれた）。
    切り口は両側の断片で同じ座標になるので、三角形分割しても辺は食い違わない。
    """
    from shapely.geometry import box as _box
    if len(poly.interiors) <= 1 or depth > 8:
        return [poly]
    hs = sorted(poly.interiors, key=lambda r: r.bounds[0])
    a, b = hs[0], hs[1]
    x0, y0, x1, y1 = poly.bounds
    if a.bounds[2] < b.bounds[0]:
        cut = 0.5 * (a.bounds[2] + b.bounds[0])
        left = _box(x0 - 1, y0 - 1, cut, y1 + 1)
        right = _box(cut, y0 - 1, x1 + 1, y1 + 1)
    else:                                   # x で分けられないときは y で分ける
        hs = sorted(poly.interiors, key=lambda r: r.bounds[1])
        a, b = hs[0], hs[1]
        cut = 0.5 * (a.bounds[3] + b.bounds[1])
        left = _box(x0 - 1, y0 - 1, x1 + 1, cut)
        right = _box(x0 - 1, cut, x1 + 1, y1 + 1)
    out = []
    for half in (left, right):
        for p in _flatten(poly.intersection(half)):
            out += split_multi_hole(p, depth + 1)
    return out


def triangulate_polygon(poly: Polygon):
    """穴つき多角形を三角形に分ける。戻り値は (3,2) 配列のリスト。

    境界の辺を必ず保つので、押し出しても閉じたメッシュになる。
    まず橋渡し + 耳刈りを試し、面積が合わなければ穴のあいだで切って再試行する
    （切ると外周に新しい頂点が入るので、押し出しの側壁と辺が食い違う。
      **切るのは最後の手段**）。
    """
    def _tri(p):
        # 共線点は**橋渡しの前**に落とす。橋渡し後に落とすと、呼び出し側が
        # 側壁に使うリングと頂点が食い違って非多様体になる。
        P = bridge_holes(clean_polygon(p))
        return [np.array([P[a], P[b], P[c]]) for a, b, c in earclip(P)]

    def _area(ts):
        a = 0.0
        for t in ts:
            u, v = t[1] - t[0], t[2] - t[0]
            a += 0.5 * abs(u[0] * v[1] - u[1] * v[0])
        return a

    tris = _tri(poly)
    if abs(_area(tris) / poly.area - 1) <= 1e-6:
        return tris
    tris = []
    for piece in split_multi_hole(poly):
        tris += _tri(piece)
    err = _area(tris) / poly.area - 1
    if abs(err) > 1e-6:
        raise ValueError(f"三角形分割の面積が {100*err:+.4f} % ずれる")
    return tris


def check_manifold(faces, tol=1e-7):
    """各辺がちょうど 2 枚の三角形に共有されるか。戻り: (非多様体辺の数, 総辺数)"""
    def key(p):
        return (round(p[0] / tol), round(p[1] / tol), round(p[2] / tol))
    c = Counter()
    for f in faces:
        k = [key(p) for p in f]
        for i in range(3):
            a, b = k[i], k[(i + 1) % 3]
            c[(a, b) if a < b else (b, a)] += 1
    return sum(1 for v in c.values() if v != 2), len(c)


def mesh_volume(faces):
    """符号付き体積（発散定理）。閉じていれば正になる。"""
    v = np.asarray(faces, float)
    return float(np.einsum("ij,ij->i", v[:, 0],
                           np.cross(v[:, 1], v[:, 2])).sum() / 6.0)


def write_stl(path, faces):
    with open(path, "wb") as fh:
        fh.write(b"\0" * 80)
        fh.write(struct.pack("<I", len(faces)))
        for tri in faces:
            v = np.asarray(tri, float)
            n = np.cross(v[1] - v[0], v[2] - v[0])
            ln = float(np.linalg.norm(n))
            n = n / ln if ln > 0 else np.zeros(3)
            fh.write(struct.pack("<3f", *n))
            for p in v:
                fh.write(struct.pack("<3f", *p))
            fh.write(struct.pack("<H", 0))
    return len(faces)


def densify_ring(P, step):
    out = []
    for a, b in zip(P[:-1], P[1:]):
        L = float(np.hypot(*(b - a)))
        k = max(int(np.ceil(L / step)), 1)
        out.append(a + (np.arange(k)[:, None] / k) * (b - a))
    out.append(P[-1:])
    return np.vstack(out)


def densify_polygon(poly, step):
    if not step:
        return poly
    rings = [densify_ring(np.asarray(poly.exterior.coords), step)]
    rings += [densify_ring(np.asarray(r.coords), step) for r in poly.interiors]
    return Polygon(rings[0], rings[1:])


def write_dxf(path, rings, note=""):
    """閉じた LWPOLYLINE として書き出す（単位 mm）。依存を増やさない最小実装。"""
    nl = chr(10)
    head = ("999" + nl + note + nl
            + "0" + nl + "SECTION" + nl + "2" + nl + "HEADER" + nl
            + "9" + nl + "$INSUNITS" + nl + "70" + nl + "4" + nl   # 4 = mm
            + "0" + nl + "ENDSEC" + nl
            + "0" + nl + "SECTION" + nl + "2" + nl + "ENTITIES" + nl)
    body = [head]
    for ring in rings:
        P = np.asarray(ring)
        body.append("0" + nl + "LWPOLYLINE" + nl + "8" + nl + "0" + nl
                    + "100" + nl + "AcDbEntity" + nl
                    + "100" + nl + "AcDbPolyline" + nl)
        body.append("90" + nl + str(len(P)) + nl + "70" + nl + "1" + nl)
        for x, y in P:
            body.append("10" + nl + f"{x:.6f}" + nl
                        + "20" + nl + f"{y:.6f}" + nl)
    body.append("0" + nl + "ENDSEC" + nl + "0" + nl + "EOF" + nl)
    open(path, "w").write("".join(body))
    return sum(len(np.asarray(r)) for r in rings)


def polygon_rings(poly):
    """外形と穴の頂点列を [(n,2), ...] で返す。"""
    return ([np.asarray(poly.exterior.coords)]
            + [np.asarray(r.coords) for r in poly.interiors])
