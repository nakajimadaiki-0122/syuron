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
    """
    outer = _ring(poly.exterior.coords)
    if _signed_area(outer) < 0:
        outer = outer[::-1]
    for ring in poly.interiors:
        hole = _ring(ring.coords)
        if _signed_area(hole) > 0:            # 穴は外側と逆回りにする
            hole = hole[::-1]
        m = int(np.argmax(hole[:, 0]))
        M = hole[m]
        # 見通せる外側頂点を近い順に探す
        d = np.hypot(*(outer - M).T)
        for v in np.argsort(d):
            seg = LineString([tuple(M), tuple(outer[v])])
            if poly.covers(seg):
                break
        else:
            raise ValueError("穴への橋を架けられない")
        hole_cycle = np.vstack([hole[m:], hole[:m], hole[m][None]])
        outer = np.vstack([outer[:v + 1], hole_cycle, outer[v:]])
    return outer


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
                w0 = ((b[0] - Q[:, 0]) * (c[1] - Q[:, 1])
                      - (b[1] - Q[:, 1]) * (c[0] - Q[:, 0])) / d
                w1 = ((c[0] - Q[:, 0]) * (a[1] - Q[:, 1])
                      - (c[1] - Q[:, 1]) * (a[0] - Q[:, 0])) / d
                w2 = 1.0 - w0 - w1
                inside = (w0 > 1e-12) & (w1 > 1e-12) & (w2 > 1e-12)
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


def triangulate_polygon(poly: Polygon):
    """穴つき多角形 -> (頂点 (n,2), 三角形 (m,3))。境界の辺を必ず保つ。"""
    P = bridge_holes(poly)
    tris = earclip(P)
    T = np.asarray(tris, int)
    area = 0.0
    for t in T:
        a, b, c = P[t[0]], P[t[1]], P[t[2]]
        area += 0.5 * abs((b[0] - a[0]) * (c[1] - a[1])
                          - (b[1] - a[1]) * (c[0] - a[0]))
    if abs(area / poly.area - 1) > 1e-6:
        raise ValueError(f"三角形分割の面積が {100*(area/poly.area-1):+.4f} % ずれる")
    return P, T


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
