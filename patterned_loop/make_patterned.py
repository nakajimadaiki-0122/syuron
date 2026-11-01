"""make_patterned.py -- SYMTVM のループ内壁に波状パターンを入れた形状を作る

元の形状は handover/src/xia2025.py の symtvm()（Xia et al. 2025）。
ループ = 涙滴形の島の丸い頭を回り込み、島の上縁に沿って戻る幅 0.30 mm の流路。
その両側の壁（外壁と島の壁）を、壁の法線方向へ正弦波状にずらす。

    offset(s) = A * sin(2 pi s / p) * taper(s)

s は壁に沿った距離。波の中心は元の壁なので、流路の平均幅はほぼ変わらない。
両端は 1 波長かけて振幅を 0 にし、元の壁へなめらかにつなぐ。
Reverse は Forward の 180 度回転で、同じ物理形状を逆向きに流す。

使い方: handover の venv の Python で実行
    python make_patterned.py
    python make_patterned.py --amp 0.05 --pitch 0.40 --walls both --out 別フォルダ名
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import shapely
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, cKDTree
from shapely.geometry import Point, Polygon

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "handover" / "src"))
sys.path.insert(0, str(ROOT / "experiment_simulation"))
import xia2025 as X
import mesh as M
import solid3d as S3
import plotstyle
import split_zones as Z

_spec = importlib.util.spec_from_file_location(
    "make_xia", ROOT / "handover" / "scripts" / "31_make_xia2025.py")
MX = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(MX)

POINTS_PER_WAVE = 16
STRAIGHT_STEP = 0.025  # 直線の壁に足す点の間隔 [mm]。波の点間隔とそろえると鋭い角でも分割が壊れない


def _triangulate_cdt(poly):
    """頂点を増やさない三角形分割。波で頂点が増えると元の耳刈りは後半の区間で失敗し、
    穴の間で切って頂点を足すため、側壁と辺が食い違って非多様体になる。"""
    g = M.clean_polygon(poly)
    tri = shapely.constrained_delaunay_triangles(g)
    return [np.asarray(t.exterior.coords)[:3] for t in tri.geoms]


def _drop_duplicates(P, tol=1e-9):
    """元の drop_collinear の代わり。足した点を消さないよう、重複点だけを落とす。"""
    d = np.hypot(*(P - np.roll(P, 1, axis=0)).T)
    return P[d > tol]


M.triangulate_polygon = _triangulate_cdt
M.drop_collinear = _drop_duplicates

# ---- 蓋の三角形分割 ----
# 波の頂点が 3 つ並ぶ所は、輪郭の点だけで分割すると極端に潰れた三角形（最大内角 173 度、
# 面積 1.6e-5 mm2）になる。Fluent は読み込み時にこれを消し、穴が開いて Free faces が出る
# （2026-11-01 に実機で確認）。壁から OFFSET だけ内側に点を足して潰れた三角形をなくす。
OFFSET = 0.03          # 壁から内側へ足す点の距離 [mm]
_triangulate_grid = Z.triangulate_polygon
_retriangulate = Z.retriangulate_flat


def _conforming(poly, spacing):
    rings = [np.asarray(poly.exterior.coords)[:-1]] + [np.asarray(r.coords)[:-1] for r in poly.interiors]
    edges = [(R, np.roll(R, -1, axis=0)) for R in rings]
    shortest = min(np.hypot(*(a_ - b_).T).min() for a_, b_ in edges)
    if shortest >= OFFSET:
        return _triangulate_grid(poly, spacing)
    bpts = np.vstack(rings)
    cands = []
    for dist, step in ((OFFSET, OFFSET), (2.2 * OFFSET, 1.5 * OFFSET)):
        inner = poly.buffer(-dist)
        for g in (inner.geoms if inner.geom_type == "MultiPolygon" else [inner]):
            if g.is_empty:
                continue
            for r in [g.exterior] + list(g.interiors):
                line = shapely.segmentize(r, step)
                cands.append(np.asarray(line.coords)[:-1])
    x0, y0, x1, y1 = poly.bounds
    grid = np.array([(x, y) for x in np.arange(x0, x1, 0.1) for y in np.arange(y0, y1, 0.1)])
    grid = grid[shapely.contains_xy(poly, grid[:, 0], grid[:, 1])]
    grid = grid[~shapely.dwithin(poly.boundary, shapely.points(grid), 3.3 * OFFSET)]
    cands.append(grid)
    C = np.vstack([c for c in cands if len(c)])
    C = C[shapely.contains_xy(poly, C[:, 0], C[:, 1])]
    C = C[~shapely.dwithin(poly.boundary, shapely.points(C), 0.8 * OFFSET)]
    # 点どうしが近すぎるものを間引く
    tree = cKDTree(C)
    keep = np.ones(len(C), bool)
    for i, j in sorted(tree.query_pairs(0.6 * OFFSET)):
        if keep[i] and keep[j]:
            keep[j] = False
    C = C[keep]
    # 輪郭の辺の直径円に入る点は辺を壊すので捨てる
    tree = cKDTree(C)
    drop = set()
    for a_, b_ in edges:
        mids = 0.5 * (a_ + b_)
        rad = 0.5 * np.hypot(*(a_ - b_).T)
        for hits in tree.query_ball_point(mids, rad * 1.001):
            drop.update(hits)
    C = np.delete(C, sorted(drop), axis=0)
    pts = np.vstack([bpts, C])
    T = pts[Delaunay(pts).simplices]
    c = T.mean(axis=1)
    T = T[shapely.contains_xy(poly, c[:, 0], c[:, 1])]
    d1, d2 = T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]
    area = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]).sum()
    got = set()
    for tri in T:
        k = [Z.key(v) for v in tri]
        for e in ((k[0], k[1]), (k[1], k[2]), (k[2], k[0])):
            got.add(tuple(sorted(e)))
    need = set()
    for a_, b_ in edges:
        for u, v in zip(a_, b_):
            need.add(tuple(sorted((Z.key(u), Z.key(v)))))
    if not need <= got or abs(area - poly.area) > 1e-6 * max(poly.area, 1):
        raise ValueError(f"蓋の分割が輪郭と合わない: 欠けた辺 {len(need - got)}、面積 {area:.5f} / {poly.area:.5f}")
    return T, "offset+delaunay"


def _retriangulate_caps(tris, axis=1, step=None, fixed=frozenset(), spacing=None):
    """流路の蓋には内側の点が増える。split_zones の「頂点が増えたら止める」検査は
    step を与えると外れるので、辺を分割しない大きな step を渡す（面積の検査は残る）。"""
    return _retriangulate(tris, axis=axis, step=step or 1e9, fixed=fixed, spacing=spacing)


Z.triangulate_polygon = _conforming
Z.retriangulate_flat = _retriangulate_caps


def densify_straight(poly, step=STRAIGHT_STEP):
    """長い直線の壁に点を足す。入口・出口の端面（x = 0, 50）の辺は端面の切り欠きと
    そろえるため触らない。"""
    def fix(coords):
        P = np.asarray(coords)[:-1]
        out = []
        for a, b in zip(P, np.roll(P, -1, axis=0)):
            out.append(a)
            on_end = (abs(a[0]) < 1e-9 and abs(b[0]) < 1e-9) or                      (abs(a[0] - X.LY) < 1e-9 and abs(b[0] - X.LY) < 1e-9)
            L = float(np.hypot(*(b - a)))
            k = int(np.ceil(L / step))
            if not on_end and k >= 2:
                out.extend(a + (b - a) * (j / k) for j in range(1, k))
        return np.asarray(out)
    return Polygon(fix(poly.exterior.coords), [fix(r.coords) for r in poly.interiors])


def wall_path(C, R, a0_deg, straight_len, amp, pitch):
    """角度 a0 から +75 度まで時計回りの円弧と、続く -15 度の直線に波を載せる。"""
    a1 = 90.0 - 0.5 * X.TIP_DEG
    arc_len = R * np.radians(a0_deg - a1)
    total = arc_len + straight_len
    n_wave = max(1, int(round(total / pitch)))
    p = total / n_wave
    n = int(round(n_wave * POINTS_PER_WAVE))
    s = np.linspace(0.0, total, n + 1)
    taper = np.minimum(1.0, np.minimum(s, total - s) / p)
    taper = taper * taper * (3.0 - 2.0 * taper)
    off = amp * np.sin(2.0 * np.pi * s / p) * taper
    pts = np.zeros((len(s), 2))
    for i, (si, oi) in enumerate(zip(s, off)):
        if si <= arc_len:
            th = np.radians(a0_deg) - si / R
            pts[i] = C + (R + oi) * np.array([np.cos(th), np.sin(th)])
        else:
            T = C + R * X.N_DN
            pts[i] = T + (si - arc_len) * X.D_DN + oi * X.N_DN
    return pts, p


def _ring_open(coords):
    P = [tuple(q) for q in coords]
    return P[:-1] if P[0] == P[-1] else P


def _insert_on_segment(ring, q, tol=1e-6):
    """q を含む辺に頂点として q を挿入する。すでに頂点ならそのまま。"""
    d = np.hypot(*(np.asarray(ring) - q).T)
    if d.min() < tol:
        return ring
    n = len(ring)
    for i in range(n):
        a, b = np.array(ring[i]), np.array(ring[(i + 1) % n])
        if abs(np.linalg.norm(q - a) + np.linalg.norm(q - b) - np.linalg.norm(b - a)) < tol:
            return ring[:i + 1] + [tuple(q)] + ring[i + 1:]
    raise ValueError(f"辺の上に点がない: {q}")


def splice(ring, start, end, via, new_pts):
    """start から via を通って end へ至る区間を new_pts に置き換える。"""
    A = np.asarray(ring)
    i0 = int(np.argmin(np.hypot(*(A - start).T)))
    if np.hypot(*(A[i0] - start)) > 1e-3:
        raise ValueError(f"始点が輪郭にない: {start}")
    ring = ring[i0:] + ring[:i0]
    A = np.asarray(ring)
    i1 = int(np.argmin(np.hypot(*(A - end).T)))
    if np.hypot(*(A[i1] - end)) > 1e-3:
        raise ValueError(f"終点が輪郭にない: {end}")
    forward = A[:i1 + 1]
    if np.hypot(*(forward - via).T).min() < 1e-3:
        return [tuple(q) for q in new_pts] + ring[i1 + 1:]
    backward = np.vstack([A[i1:], A[:1]])
    if np.hypot(*(backward - via).T).min() >= 1e-3:
        raise ValueError(f"経由点が輪郭にない: {via}")
    return [tuple(q) for q in new_pts[::-1]] + ring[1:i1]


def build(amp, pitch, walls):
    poly, info = X.symtvm()
    outer = _ring_open(poly.exterior.coords)
    holes = [_ring_open(r.coords) for r in poly.interiors]
    Cu = [np.array(c) for c in info["centres_upper"]]
    half = np.radians(0.5 * X.TIP_DEG)
    used_pitch = []

    for k, C in enumerate(Cu):
        T_out = C + X.R_OUT * X.N_DN
        T_in = C + X.R_IN * X.N_DN
        tip = C + np.array([X.R_IN / np.sin(half), 0.0])
        straight = float(np.dot(tip - T_in, X.D_DN))
        foot = T_out + straight * X.D_DN
        entry = np.array([info["x_in"], info["inlet_y"][1]]) if k == 0 else X._entry_pt(Cu, k)
        a_entry = np.degrees(np.arctan2(*(entry - C)[::-1])) % 360.0
        if walls in ("both", "outer"):
            pts, p = wall_path(C, X.R_OUT, a_entry, straight, amp, pitch)
            used_pitch.append(p)
            pts[0], pts[-1] = entry, foot
            for sgn in (1.0, -1.0):
                m = np.array([1.0, sgn])
                outer = _insert_on_segment(outer, foot * m)
                outer = splice(outer, entry * m, foot * m, T_out * m, pts * m)
        if walls in ("both", "island"):
            # 島の頭の円弧が弦で切られる点（中心から見て 206.5 度付近）
            a_chord = 75.0 + np.degrees(np.arccos(info["chord_off"] / X.R_IN))
            guess = C + X.R_IN * np.array([np.cos(np.radians(a_chord)), np.sin(np.radians(a_chord))])
            for sgn in (1.0, -1.0):
                m = np.array([1.0, sgn])
                inside = C * m
                j = next(i for i, h in enumerate(holes) if Polygon(h).contains(Point(inside)))
                H = np.asarray(holes[j])
                start = H[int(np.argmin(np.hypot(*(H - guess * m).T)))]
                upper_start = start * m
                a0 = np.degrees(np.arctan2(*(upper_start - C)[::-1])) % 360.0
                pts, p = wall_path(C, X.R_IN, a0, straight, amp, pitch)
                if sgn > 0:
                    used_pitch.append(p)
                pts[0], pts[-1] = upper_start, tip
                holes[j] = splice(holes[j], start, tip * m, T_in * m, pts * m)

    out = Polygon(outer, holes)
    if not out.is_valid:
        raise ValueError("波を入れた多角形が不正（壁が交差している）")
    return poly, out, info, used_pitch


def curvature_radius(ring):
    """輪郭の各頂点の外接円半径。"""
    P = np.asarray(ring)
    a, b, c = np.roll(P, 1, 0), P, np.roll(P, -1, 0)
    ab, bc, ca = (np.linalg.norm(a - b, axis=1), np.linalg.norm(b - c, axis=1),
                  np.linalg.norm(c - a, axis=1))
    cr = np.abs((b - a)[:, 0] * (c - a)[:, 1] - (b - a)[:, 1] * (c - a)[:, 0])
    with np.errstate(divide="ignore", invalid="ignore"):
        return ab * bc * ca / (2.0 * cr)


def edge_hist(tris):
    v = np.round(tris.reshape(-1, 3).astype(np.float32), 6)
    _, inv = np.unique(v, axis=0, return_inverse=True)
    idx = inv.reshape(-1, 3)
    e = np.sort(np.vstack([idx[:, [0, 1]], idx[:, [1, 2]], idx[:, [2, 0]]]), axis=1)
    _, c = np.unique(e, axis=0, return_counts=True)
    vals, cnt = np.unique(c, return_counts=True)
    return {int(a): int(b) for a, b in zip(vals, cnt)}


def verify_zones(d):
    names = ["inlet", "outlet", "fluid_top", "fluid_bottom", "fluid_walls",
             "bottom_heat", "solid_top", "sym_left", "sym_right", "solid_ends"]
    parts = {n: Z.read_stl(d / "zones" / f"{n}.stl") for n in names}
    tris = np.vstack(list(parts.values()))
    hmin = {}
    for n, t in parts.items():
        ab, bc, ca = t[:, 1] - t[:, 0], t[:, 2] - t[:, 1], t[:, 0] - t[:, 2]
        area = 0.5 * np.linalg.norm(np.cross(ab, -ca), axis=1)
        emax = np.max(np.stack([np.linalg.norm(q, axis=1) for q in (ab, bc, ca)]), axis=0)
        hmin[n] = float((2 * area / emax).min())
    return dict(triangles=int(len(tris)), edge_use=edge_hist(tris), min_triangle_height_mm=hmin)


def draw_previews(base, wavy, amp, pitch, walls, out):
    plotstyle.use_jp()
    rev = X.reverse(wavy)
    fig, axes = plt.subplots(2, 1, figsize=(16, 4.6))
    for ax, poly, title in ((axes[0], wavy, "SYMTVM-Forward（波状パターン入り）"),
                            (axes[1], rev, "SYMTVM-Reverse（同じ形を 180 度回転）")):
        MX.draw(ax, poly, title)
    axes[-1].set_xlabel("x [mm]")
    fig.tight_layout()
    fig.savefig(out / "preview_overview.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6.2))
    P = np.asarray(wavy.exterior.coords)
    ax.fill(P[:, 0], P[:, 1], color="#9fc9cb", label="流体（パターンあり）")
    for r in wavy.interiors:
        Q = np.asarray(r.coords)
        ax.fill(Q[:, 0], Q[:, 1], color="white")
        ax.plot(Q[:, 0], Q[:, 1], color="#0e7c7b", lw=0.9)
    ax.plot(P[:, 0], P[:, 1], color="#0e7c7b", lw=0.9)
    for i, r in enumerate([base.exterior] + list(base.interiors)):
        Q = np.asarray(r.coords)
        ax.plot(Q[:, 0], Q[:, 1], color="#c2410c", lw=0.8, ls="--",
                label="元の壁（Xia 2025）" if i == 0 else None)
    x0 = X.X_C1 + X.PITCH
    ax.set_xlim(x0 - 1.3, x0 + 4.2)
    ax.set_ylim(-0.1, 1.72)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]（CAD では z）")
    ax.set_title(f"ループ 1 個の拡大  片振幅 {amp} mm、波長 約 {pitch} mm、対象の壁: {walls}", loc="left")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "preview_loop_zoom.png", dpi=170)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amp", type=float, default=0.05, help="片振幅 [mm]")
    ap.add_argument("--pitch", type=float, default=0.40, help="目標の波長 [mm]")
    ap.add_argument("--walls", choices=["both", "outer", "island"], default="both")
    ap.add_argument("--out", default=".", help="patterned_loop からの出力先")
    a = ap.parse_args()
    out = (HERE / a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    base, wavy, info, used = build(a.amp, a.pitch, a.walls)
    z_floor, z_top = X.H_SOLID - X.H_CH, X.H_SOLID
    h_half = 0.5 * X.W_UNIT
    cuts = [X.X_C1 + X.PITCH * k + 3.0 for k in range(X.N_STAGE - 1)]
    dense = densify_straight(wavy)
    ff = S3.fluid_solid_cut(dense, z_floor, z_top, cuts)
    bf = MX.block_solid(dense, z_floor, z_top, h_half, X.LY, cuts)
    rot = lambda F: [[(X.LY - q[0], -q[1], q[2]) for q in t] for t in F]

    rs = np.concatenate([curvature_radius(r.coords[:-1])
                         for r in [wavy.exterior] + list(wavy.interiors)])
    report = dict(
        source="Xia et al. 2025 SYMTVM に本研究でループ内壁の波状パターンを追加",
        amplitude_mm=a.amp, target_pitch_mm=a.pitch, walls=a.walls,
        pitch_used_mm=[float(min(used)), float(max(used))],
        plan_area_base_mm2=float(base.area), plan_area_pattern_mm2=float(wavy.area),
        wall_perimeter_base_mm=float(base.length), wall_perimeter_pattern_mm=float(wavy.length),
        min_gap_outer_to_island_mm=float(min(wavy.exterior.distance(r) for r in wavy.interiors)),
        min_gap_outer_to_island_base_mm=float(min(base.exterior.distance(r) for r in base.interiors)),
        max_abs_y_mm=float(max(abs(wavy.bounds[1]), abs(wavy.bounds[3]))),
        curvature_radius_percentile_1_mm=float(np.nanpercentile(rs[np.isfinite(rs)], 1)),
        models={},
    )

    configs = {"symtvm_fwd": (wavy, ff, bf),
               "symtvm_rev": (X.reverse(wavy), rot(ff), rot(bf))}
    for cfg, (poly, f, b) in configs.items():
        d = out / cfg
        d.mkdir(exist_ok=True)
        bad_f, _ = M.check_manifold(f)
        bad_b, _ = M.check_manifold(b)
        M.write_stl(d / "fluid.stl", MX.to_yup(f))
        M.write_stl(d / "solid_copper.stl", MX.to_yup(b))
        rings = M.polygon_rings(poly)
        MX.write_dxf_ezdxf(str(d / "outline.dxf"), rings, f"{cfg} wavy loop", "y",
                           X.LY, X.W_UNIT, X.H_SOLID, X.H_CH, 1.0)
        report["models"][cfg] = dict(
            nonmanifold_fluid=int(bad_f), nonmanifold_solid=int(bad_b),
            fluid_volume_mm3=float(M.mesh_volume(f)), solid_volume_mm3=float(M.mesh_volume(b)),
            total_volume_mm3=float(M.mesh_volume(f) + M.mesh_volume(b)))

    Z.HERE = out
    Z.main(list(configs))
    print("分割手法の内訳:", Z.METHODS)
    for cfg in configs:
        report["models"][cfg].update(verify_zones(out / cfg))

    draw_previews(base, wavy, a.amp, a.pitch, a.walls, out)
    (out / "params.json").write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
