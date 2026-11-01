"""make_hemisphere.py -- SYMTVM のループ側壁に半球を等間隔に置いた CAD を作る

パターン 2: --mode bump    半球が流路の内側へ出っぱる（pattern2_hemisphere/）
パターン 3: --mode dimple  半球状に壁をへこませる（pattern3_dimple/）

半径 0.1 mm、中心間隔 0.4 mm（壁に沿う方向も深さ方向も同じ）、深さ方向 3 段。
置くのはループ部分の外壁と島の壁のみ。底面・上面と、幹の流路には置かない。
流路の平面形は元のまま。

半球は壁を展開した (s, z) 座標で作り、壁の法線方向へ高さ sqrt(r^2 - rho^2) だけずらす。
直線部では厳密に半球で、円弧部（R 0.58 / 0.88）では壁の曲率の分だけわずかに変形する。

使い方: handover の venv の Python で実行。Fluent は使わない。
"""
import argparse
import json

import numpy as np
import shapely
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, cKDTree
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union

import make_patterned as P

X, M, S3, Z, MX, HERE = P.X, P.M, P.S3, P.Z, P.MX, P.HERE
# 平面の壁は元の円弧なので、波のような潰れた三角形は出ない。蓋は元の方法で分割する
Z.triangulate_polygon = P._triangulate_grid

R_HEMI = 0.10
SPACING = 0.40
RIM_POINTS = 24
Z_FLOOR, Z_TOP = 0.5, 2.0
ROWS = Z_FLOOR + 0.75 + SPACING * np.array([-1.0, 0.0, 1.0])
MAX_WALL_ANGLE = 150.0   # 壁の三角形の最大内角の目標 [deg]
WALL_STEP_S, WALL_STEP_Z = 0.05, 0.10   # 壁の内側に置く点の間隔 [mm]
SAMPLE_PITCH = 0.48       # wall_path の 1 波分。POINTS_PER_WAVE = 16 で点間隔 約 0.03 mm
MARGIN_START, MARGIN_END_OUTER, MARGIN_END_ISLAND = 0.20, 0.20, 0.25


def _tri_angles(T):
    """三角形 (N, 3, 2 or 3) の 3 つの内角 [deg]。"""
    out = []
    for i in range(3):
        u = T[:, (i + 1) % 3] - T[:, i]
        v = T[:, (i + 2) % 3] - T[:, i]
        c = (u * v).sum(axis=1) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1))
        out.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return np.stack(out, axis=1)


class Panel:
    """ループの壁 1 枚。s は壁に沿う距離、inward_sign は流体側の法線の符号。"""

    def __init__(self, kind, C, R, a0, straight, pts, fluid_sign, margin_end):
        self.kind, self.C, self.R, self.a0 = kind, C, R, a0
        self.arc = R * np.radians(a0 - 75.0)
        self.L = self.arc + straight
        self.pts = pts
        self.s = np.linspace(0.0, self.L, len(pts))
        self.fluid_sign = fluid_sign
        span = self.L - MARGIN_START - margin_end
        n = int(np.floor(span / SPACING)) + 1
        s0 = MARGIN_START + 0.5 * (span - (n - 1) * SPACING)
        self.centres = [(s0 + i * SPACING, z) for i in range(n) for z in ROWS]

    def frame(self, s):
        if s <= self.arc:
            th = np.radians(self.a0) - s / self.R
            n = np.array([np.cos(th), np.sin(th)])
            return self.C + self.R * n, n
        return self.C + self.R * X.N_DN + (s - self.arc) * X.D_DN, X.N_DN.copy()


def build_plan():
    """元の平面形。ループの壁だけ等間隔の点に置き換え、壁の情報を返す。"""
    poly, info = X.symtvm()
    outer = P._ring_open(poly.exterior.coords)
    holes = [P._ring_open(r.coords) for r in poly.interiors]
    Cu = [np.array(c) for c in info["centres_upper"]]
    half = np.radians(0.5 * X.TIP_DEG)
    panels = []
    for k, C in enumerate(Cu):
        T_out, T_in = C + X.R_OUT * X.N_DN, C + X.R_IN * X.N_DN
        tip = C + np.array([X.R_IN / np.sin(half), 0.0])
        straight = float(np.dot(tip - T_in, X.D_DN))
        foot = T_out + straight * X.D_DN
        entry = np.array([info["x_in"], info["inlet_y"][1]]) if k == 0 else X._entry_pt(Cu, k)
        a_out = np.degrees(np.arctan2(*(entry - C)[::-1])) % 360.0
        pts, _ = P.wall_path(C, X.R_OUT, a_out, straight, 0.0, SAMPLE_PITCH)
        pts[0], pts[-1] = entry, foot
        panels.append(Panel("outer", C, X.R_OUT, a_out, straight, pts, -1.0, MARGIN_END_OUTER))
        for sgn in (1.0, -1.0):
            m = np.array([1.0, sgn])
            outer = P._insert_on_segment(outer, foot * m)
            outer = P.splice(outer, entry * m, foot * m, T_out * m, pts * m)
        a_chord = 75.0 + np.degrees(np.arccos(info["chord_off"] / X.R_IN))
        guess = C + X.R_IN * np.array([np.cos(np.radians(a_chord)), np.sin(np.radians(a_chord))])
        j = next(i for i, h in enumerate(holes) if Polygon(h).contains(Point(C)))
        H = np.asarray(holes[j])
        start = H[int(np.argmin(np.hypot(*(H - guess).T)))]
        a_in = np.degrees(np.arctan2(*(start - C)[::-1])) % 360.0
        pts, _ = P.wall_path(C, X.R_IN, a_in, straight, 0.0, SAMPLE_PITCH)
        pts[0], pts[-1] = start, tip
        panels.append(Panel("island", C, X.R_IN, a_in, straight, pts, 1.0, MARGIN_END_ISLAND))
        for sgn in (1.0, -1.0):
            m = np.array([1.0, sgn])
            jj = next(i for i, h in enumerate(holes) if Polygon(h).contains(Point(C * m)))
            holes[jj] = P.splice(holes[jj], start * m, tip * m, T_in * m, pts * m)
    return poly, Polygon(outer, holes), panels, info


def dome_param():
    """半球面を平面へ等距離に展開した点と三角形。縁は RIM_POINTS 点で共通。"""
    nv = 4
    pts = [(0.0, 0.0)]              # (極角 v, 方位 phi)
    for j in range(1, nv):
        v = (np.pi / 2) * j / nv
        n = max(6, int(round(RIM_POINTS * np.sin(v))))
        pts += [(v, 2 * np.pi * i / n) for i in range(n)]
    pts += [(np.pi / 2, 2 * np.pi * i / RIM_POINTS) for i in range(RIM_POINTS)]
    A = np.array(pts)
    flat = np.c_[A[:, 0] * np.cos(A[:, 1]), A[:, 0] * np.sin(A[:, 1])]
    return A, Delaunay(flat).simplices


DOME_A, DOME_T = dome_param()


def panel_faces(panel, height_sign):
    """壁 1 枚の三角形（流体から外向きの法線）。height_sign: +1 で流体側へ出っぱる。"""
    n = len(panel.pts)
    shell = [(panel.s[i], Z_FLOOR) for i in range(n)] + [(panel.s[i], Z_TOP) for i in range(n - 1, -1, -1)]
    rims = []
    for sc, zc in panel.centres:
        phis = 2 * np.pi * np.arange(RIM_POINTS) / RIM_POINTS
        rims.append(np.c_[sc + R_HEMI * np.cos(phis), zc + R_HEMI * np.sin(phis)])
    flat = Polygon(shell, rims)
    if not flat.is_valid:
        raise ValueError("半球の縁が壁からはみ出している")
    # 円弧の壁を平面でつなぐと角張るので、壁の内側にも点を置く。
    # 長方形の外周の辺は凸包の辺なので Delaunay でも必ず残る
    inner = []
    arc_end = panel.arc + 0.1
    s_grid = np.r_[np.arange(0.10, min(arc_end, panel.L - 0.10), WALL_STEP_S),
                   np.arange(arc_end, panel.L - 0.10, 4 * WALL_STEP_S)]
    z_grid = np.arange(Z_FLOOR + 0.06, Z_TOP - 0.05, WALL_STEP_Z)
    # 列ごとに半間隔ずらす（点が一直線に並ぶと潰れた三角形ができる）
    grid = np.array([(a_, b_ + (0.5 * WALL_STEP_Z if i_ % 2 else 0.0))
                     for i_, a_ in enumerate(s_grid) for b_ in z_grid])
    grid = grid[grid[:, 1] < Z_TOP - 0.05]
    holes_c = np.array(panel.centres)
    dist = np.min(np.hypot(grid[:, None, 0] - holes_c[None, :, 0], grid[:, None, 1] - holes_c[None, :, 1]), axis=1)
    inner.append(grid[dist > R_HEMI + 0.06])
    for sc, zc in panel.centres:                    # 縁のすぐ外の輪
        phis = 2 * np.pi * (np.arange(RIM_POINTS) + 0.5) / RIM_POINTS
        inner.append(np.c_[sc + (R_HEMI + 0.025) * np.cos(phis), zc + (R_HEMI + 0.025) * np.sin(phis)])
    inner = np.vstack(inner)
    inner = inner[(inner[:, 0] > 0.09) & (inner[:, 0] < panel.L - 0.09)]
    # 潰れた三角形は Fluent に消されて穴になるので、最大内角が大きい三角形の鈍角側の
    # 内部点を取り除いて分割し直す。外周と半球の縁の点は壁・半球と共有なので残す
    fixed_pts = np.vstack([np.asarray(shell), np.vstack(rims)])
    for _ in range(40):
        all_pts = np.vstack([fixed_pts, inner])
        tri = Delaunay(all_pts)
        cand = all_pts[tri.simplices]
        kp = shapely.contains_xy(flat, cand[:, :, 0].mean(axis=1), cand[:, :, 1].mean(axis=1))
        simp = tri.simplices[kp]
        angles = _tri_angles(all_pts[simp])
        bad = angles.max(axis=1) > MAX_WALL_ANGLE
        if not bad.any():
            break
        drop = set()
        for row, ang in zip(simp[bad], angles[bad]):
            interior = [v for v in row if v >= len(fixed_pts)]
            if not interior:
                continue
            obtuse = row[int(np.argmax(ang))]
            drop.add(int(obtuse if obtuse >= len(fixed_pts) else interior[0]) - len(fixed_pts))
        if not drop:
            break
        inner = np.delete(inner, sorted(drop), axis=0)
    wall2d = all_pts[simp]
    keep = np.ones(len(wall2d), bool)
    wall2d = wall2d[keep]
    area2d = 0.5 * np.abs((wall2d[:, 1, 0] - wall2d[:, 0, 0]) * (wall2d[:, 2, 1] - wall2d[:, 0, 1])
                          - (wall2d[:, 1, 1] - wall2d[:, 0, 1]) * (wall2d[:, 2, 0] - wall2d[:, 0, 0])).sum()
    if abs(area2d - flat.area) > 1e-6 * flat.area:
        raise ValueError(f"壁の分割の面積が合わない: {area2d} / {flat.area}")

    def to3d(s, z, h=0.0):
        i = int(np.searchsorted(panel.s, s - 1e-9))
        if i < n and abs(panel.s[i] - s) < 1e-9 and abs(h) == 0.0 and (z == Z_FLOOR or z == Z_TOP):
            xy = panel.pts[i]                       # 平面形の頂点と厳密に一致させる
            return np.array([xy[0], xy[1], z])
        p, nrm = panel.frame(s)
        q = p + height_sign * panel.fluid_sign * h * nrm
        return np.array([q[0], q[1], z])

    tris = []
    for q in wall2d:
        tris.append((q, [to3d(a, b) for a, b in q]))
    for sc, zc in panel.centres:
        v, phi = DOME_A[:, 0], DOME_A[:, 1]
        rho = R_HEMI * np.sin(v)
        s_ = sc + rho * np.cos(phi)
        z_ = zc + rho * np.sin(phi)
        h_ = R_HEMI * np.cos(v)
        h_[np.isclose(v, np.pi / 2)] = 0.0
        V3 = [to3d(a, b, c) for a, b, c in zip(s_, z_, h_)]
        flat2 = np.c_[v * np.cos(phi), v * np.sin(phi)]
        for t in DOME_T:
            tris.append((flat2[t], [V3[i] for i in t]))
    # 展開面の向きをそろえ、平坦な三角形で流体から外向きになるよう全体の向きを決める
    faces = []
    for q2, q3 in tris:
        area = (q2[1][0] - q2[0][0]) * (q2[2][1] - q2[0][1]) - (q2[1][1] - q2[0][1]) * (q2[2][0] - q2[0][0])
        faces.append(q3 if area > 0 else q3[::-1])
    faces = np.array(faces)
    a, b, c = faces[0]
    nrm = np.cross(b - a, c - a)
    sm = faces[0].mean(axis=0)
    idx = int(np.argmin(np.abs(panel.s - float(np.clip(tris[0][0].mean(axis=0)[0], 0, panel.L)))))
    _, n2 = panel.frame(panel.s[idx])
    outward = -panel.fluid_sign * np.array([n2[0], n2[1], 0.0])
    if np.dot(nrm, outward) < 0:
        faces = faces[:, ::-1]
    return faces, len(wall2d)


def segment_key(a, b):
    return tuple(sorted((tuple(np.round(a, 6)), tuple(np.round(b, 6)))))


def extruded_walls(region, cuts, skip, outward):
    out = []
    for g in S3._pieces(region, cuts):
        for P_ in [S3.closed(g.exterior, True)] + [S3.closed(r, False) for r in g.interiors]:
            for a, b in zip(P_[:-1], P_[1:]):
                if S3._on_cut(a, b, cuts) or segment_key(a, b) in skip:
                    continue
                if outward is False and S3._on_cut(a, b, [0.0, X.LY]):
                    continue
                q = [(a[0], a[1], Z_FLOOR), (b[0], b[1], Z_FLOOR), (b[0], b[1], Z_TOP), (a[0], a[1], Z_TOP)]
                t = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
                out += t if outward else [u[::-1] for u in t]
    return out


def mirror(faces):
    F = np.array(faces, float)
    F[:, :, 1] *= -1.0
    return F[:, ::-1]


def solid_block(plan, cuts, skip, panel_tris):
    h = 0.5 * X.W_UNIT
    rect = box(0.0, -h, X.LY, h)
    faces = S3.caps_cut(rect, 0.0, False, cuts)
    faces += S3.caps_cut(rect.difference(plan), Z_TOP, True, cuts)
    faces += S3.caps_cut(plan, Z_FLOOR, True, cuts)
    faces += extruded_walls(plan, cuts, skip, False)
    xs = [0.0] + sorted(cuts) + [X.LY]
    for xa, xb in zip(xs[:-1], xs[1:]):
        for y in (-h, h):
            q = [(xa, y, 0.0), (xb, y, 0.0), (xb, y, Z_TOP), (xa, y, Z_TOP)]
            t = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
            faces += [u[::-1] for u in t] if y > 0 else t
    xs_, ys_ = np.asarray(plan.exterior.coords).T
    m0, m1 = np.abs(xs_) < 1e-9, np.abs(xs_ - X.LY) < 1e-9
    faces += MX._end_face(0.0, ys_[m0].min(), ys_[m0].max(), Z_FLOOR, Z_TOP, h, -1)
    faces += MX._end_face(X.LY, ys_[m1].min(), ys_[m1].max(), Z_FLOOR, Z_TOP, h, +1)
    return [list(map(tuple, t)) for t in faces] + [list(map(tuple, t)) for t in panel_tris[:, ::-1]]


def thin_walls(panels, base, mode):
    """へこみの裏に残る銅の最小厚さ [mm]。"""
    if mode != "dimple":
        return None
    best = np.inf
    rect = box(0.0, -0.5 * X.W_UNIT, X.LY, 0.5 * X.W_UNIT)
    for pn in panels:
        for sc, zc in pn.centres:
            p, _ = pn.frame(sc)
            disc = Point(*p).buffer(R_HEMI, 32)
            if pn.kind == "island":
                ring = next(r for r in base.interiors if Polygon(r).contains(Point(pn.C)))
                # 同じ島の反対側（幹の流路に面する壁）までの距離
                others = [LineString(np.asarray(ring.coords)[[i, i + 1]])
                          for i in range(len(ring.coords) - 1)
                          if not any(np.hypot(*(np.asarray(ring.coords)[i] - q)) < 0.25 for q in [p])]
                d = min(disc.difference(Point(*p).buffer(0)).distance(o) for o in others
                        if np.dot(np.asarray(o.coords)[0] - p, pn.frame(sc)[1]) < -0.02) if others else np.inf
            else:
                d = rect.exterior.distance(disc)
            best = min(best, d)
    return float(best)


def split_hemispheres(zdir, domes_ztop):
    """半球の面を fluid_walls.stl から hemispheres.stl へ分ける。
    分かれ目が境界線になるので、Fluent が半球の縁の鋭い線を保てる。"""
    d = np.stack([domes_ztop[:, :, 0], domes_ztop[:, :, 2], -domes_ztop[:, :, 1]], axis=2)   # Y-up
    walls = Z.read_stl(zdir / "fluid_walls.stl")
    # STL は単精度なので頂点の丸めでは一致しない。三角形の重心の近さで照合する
    dist, _ = cKDTree(d.mean(axis=1)).query(walls.mean(axis=1))
    mask = dist < 2e-5
    if mask.sum() != len(d):
        raise ValueError(f"半球の三角形が一致しない: {mask.sum()} / {len(d)}")
    Z.write_stl(zdir / "hemispheres.stl", walls[mask])
    Z.write_stl(zdir / "fluid_walls.stl", walls[~mask])
    return int(mask.sum())


def verify(d):
    names = ["inlet", "outlet", "fluid_top", "fluid_bottom", "fluid_walls", "hemispheres",
             "bottom_heat", "solid_top", "sym_left", "sym_right", "solid_ends"]
    parts = [Z.read_stl(d / "zones" / f"{n}.stl") for n in names]
    tris = np.vstack(parts)
    return dict(zone_files=len(names), triangles=int(len(tris)), edge_use=P.edge_hist(tris))


def render_loop(panels, plan, base, mode, out):
    """実際に作った面を描く。ループ 1 個分。"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("preview", HERE / "pattern2_hemisphere" / "make_preview.py")
    PV = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(PV)
    outer_p, isl_p = panels[2], panels[3]
    C = outer_p.C
    h_sign = 1.0 if mode == "bump" else -1.0
    x0, x1, y0, y1 = C[0] - 1.25, C[0] + 2.65, C[1] - 0.65, C[1] + 1.15
    floor = PV.plate(x0, x1, y0, y1, Z_FLOOR)
    chord = X.symtvm()[1]["chord_off"]
    island = X.teardrop(C, True, chord)
    tops = [np.c_[np.asarray(g.exterior.coords)[:3], np.full(3, Z_TOP)]
            for g in shapely.constrained_delaunay_triangles(island).geoms]
    band = LineString(outer_p.pts).buffer(0.2, single_sided=True)
    if band.distance(Point(C)) < 0.5:
        band = LineString(outer_p.pts).buffer(-0.2, single_sided=True)
    tops2 = [np.c_[np.asarray(g.exterior.coords)[:3], np.full(3, Z_TOP)]
             for g in shapely.constrained_delaunay_triangles(band).geoms]
    wall_i, wall_o = panel_faces(isl_p, h_sign)[0], panel_faces(outer_p, h_sign)[0]

    def split(F, col_wall):
        c = F.reshape(-1, 3)
        return [(F, col_wall)]

    a = PV.render([(floor, PV.FLOOR), (np.array(tops), PV.COPPER), (wall_i, PV.COPPER)], (.62, -.62, -.62))
    b = PV.render([(floor, PV.FLOOR), (np.array(tops2), PV.COPPER), (wall_o, PV.COPPER)], (-.45, .75, -.6))

    rows = []
    for pn in (outer_p, isl_p):
        for sc, zc in pn.centres:
            if np.isclose(zc, ROWS[1]):
                rows.append(Point(*pn.frame(sc)[0]).buffer(R_HEMI, 64))
    discs = unary_union(rows)
    fluid = base.difference(discs) if mode == "bump" else base.union(discs)
    P.plotstyle.use_jp()
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1])
    ax = fig.add_subplot(gs[0, 0]); ax.imshow(a); ax.axis("off")
    ax.set_title("(a) 島の壁（外壁を外して見た図）", loc="left", fontsize=12)
    ax = fig.add_subplot(gs[0, 1]); ax.imshow(b); ax.axis("off")
    ax.set_title("(b) 外壁（島を外して内側から見た図）", loc="left", fontsize=12)
    ax = fig.add_subplot(gs[1, :])
    view = box(x0, y0, x1, y1)
    shown = fluid.intersection(view)
    for g in (shown.geoms if shown.geom_type == "MultiPolygon" else [shown]):
        Q = np.asarray(g.exterior.coords)
        ax.fill(Q[:, 0], Q[:, 1], color="#bfe0e2", lw=0)
        for r in g.interiors:
            Q = np.asarray(r.coords)
            ax.fill(Q[:, 0], Q[:, 1], color="white", lw=0)
    for r in [shown.exterior] + list(shown.interiors) if shown.geom_type == "Polygon" else []:
        Q = np.asarray(r.coords)
        ax.plot(Q[:, 0], Q[:, 1], color="#7a4a2b", lw=1.0)
    for r in [base.exterior] + list(base.interiors):
        Q = np.asarray(r.coords)
        ax.plot(Q[:, 0], Q[:, 1], color="#c2410c", lw=0.7, ls="--")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
    ax.set_xlabel("流れ方向 x [mm]"); ax.set_ylabel("幅方向 [mm]")
    ax.set_title("(c) 流路の深さの中央の断面  水色 = 水、破線 = 元の壁", loc="left", fontsize=12)
    name = "パターン 2: 半球の出っぱり" if mode == "bump" else "パターン 3: 半球のへこみ"
    fig.suptitle(f"{name}  半径 {R_HEMI} mm、中心間隔 {SPACING} mm、深さ方向 {len(ROWS)} 段（生成した CAD の面をそのまま描画）", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(out / "preview_cad.png", dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["bump", "dimple"], required=True)
    a = ap.parse_args()
    out = HERE / ("pattern2_hemisphere" if a.mode == "bump" else "pattern3_dimple")
    out.mkdir(exist_ok=True)
    base, plan, panels, info = build_plan()
    plan_dense = P.densify_straight(plan, 0.05)   # ループの壁の点間隔（約 0.03）は分割しない
    h_sign = 1.0 if a.mode == "bump" else -1.0

    skip = set()
    upper, upper_domes = [], []
    for pn in panels:
        for sgn in (1.0, -1.0):
            m = np.array([1.0, sgn])
            for p0, p1 in zip(pn.pts[:-1], pn.pts[1:]):
                skip.add(segment_key(p0 * m, p1 * m))
        F, n_wall = panel_faces(pn, h_sign)
        upper.append(F)
        upper_domes.append(F[n_wall:])
    upper = np.vstack(upper)
    domes = np.vstack(upper_domes)
    domes = np.vstack([domes, mirror(domes)])
    panel_tris = np.vstack([upper, mirror(upper)])

    cuts = [X.X_C1 + X.PITCH * k + 3.0 for k in range(X.N_STAGE - 1)]
    ff = S3.caps_cut(plan_dense, Z_FLOOR, False, cuts) + S3.caps_cut(plan_dense, Z_TOP, True, cuts)
    ff += extruded_walls(plan_dense, cuts, skip, True)
    ff = [list(map(tuple, t)) for t in ff] + [list(map(tuple, t)) for t in panel_tris]
    bf = solid_block(plan_dense, cuts, skip, panel_tris)
    rot = lambda F: [[(X.LY - q[0], -q[1], q[2]) for q in t] for t in F]

    n_hemi = 2 * sum(len(pn.centres) for pn in panels)
    v_hemi = n_hemi * (2.0 / 3.0) * np.pi * R_HEMI ** 3
    base_v = base.area * X.H_CH
    report = dict(mode=a.mode, radius_mm=R_HEMI, spacing_mm=SPACING, rows_z_mm=[float(v) for v in ROWS],
                  hemispheres_total=n_hemi, fluid_volume_base_mm3=float(base_v),
                  fluid_volume_expected_mm3=float(base_v - h_sign * v_hemi),
                  min_copper_thickness_behind_dimple_mm=thin_walls(panels, base, a.mode), models={})
    configs = {"symtvm_fwd": (plan_dense, ff, bf), "symtvm_rev": (X.reverse(plan_dense), rot(ff), rot(bf))}
    for cfg, (poly, f, b) in configs.items():
        d = out / cfg
        d.mkdir(exist_ok=True)
        bad_f, _ = M.check_manifold(f)
        bad_b, _ = M.check_manifold(b)
        M.write_stl(d / "fluid.stl", MX.to_yup(f))
        M.write_stl(d / "solid_copper.stl", MX.to_yup(b))
        report["models"][cfg] = dict(nonmanifold_fluid=int(bad_f), nonmanifold_solid=int(bad_b),
                                     fluid_volume_mm3=float(M.mesh_volume(f)),
                                     solid_volume_mm3=float(M.mesh_volume(b)),
                                     total_volume_mm3=float(M.mesh_volume(f) + M.mesh_volume(b)))
    Z.HERE = out
    Z.main(list(configs))
    print("分割手法の内訳:", Z.METHODS)
    rot_np = lambda F: np.stack([X.LY - F[:, :, 0], -F[:, :, 1], F[:, :, 2]], axis=2)
    for cfg in configs:
        d3 = domes if cfg == "symtvm_fwd" else rot_np(domes)
        n_moved = split_hemispheres(out / cfg / "zones", d3)
        report["models"][cfg]["hemisphere_triangles"] = n_moved
        report["models"][cfg].update(verify(out / cfg))
    render_loop(panels, plan, base, a.mode, out)
    (out / "params.json").write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
