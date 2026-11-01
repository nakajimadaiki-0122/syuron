"""パターン 2 のイメージ図: ループの側壁に半球（r = 0.1 mm）を等間隔に並べる。

既存の流路の形は変えない。半球は流路の内側へ突き出す。底面・上面には置かない。
対象はループ部分の外壁と島の壁のみ。Z バッファで描くので隠れ面が正しく消える。
"""
import json
import sys
from pathlib import Path

import numpy as np
import shapely
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "handover" / "src"))
import xia2025 as X
import plotstyle

R_HEMI = 0.10          # 半球の半径 [mm]
SPACING = 0.40         # 中心間隔 [mm]。壁に沿う方向も深さ方向も同じ
Z_FLOOR, Z_TOP = 0.5, 2.0
ROWS = Z_FLOOR + 0.75 + SPACING * np.array([-1, 0, 1])   # 深さ方向 3 段（中央が流路の中央）
STAGE = 1              # 図にする段（0 始まり）
COPPER, FLOOR, HEMI = np.array([.80, .50, .32]), np.array([.90, .78, .68]), np.array([.05, .62, .62])
LIGHT = np.array([-.35, -.45, .82]) / np.linalg.norm([-.35, -.45, .82])


def wall_point(C, R, a0, s):
    """角度 a0 から +75 度まで時計回りの円弧と続く -15 度の直線。位置と法線を返す。"""
    arc = R * np.radians(a0 - 75.0)
    if s <= arc:
        th = np.radians(a0) - s / R
        n = np.array([np.cos(th), np.sin(th)])
        return C + R * n, n
    return C + R * X.N_DN + (s - arc) * X.D_DN, X.N_DN.copy()


def walls(k):
    """k 段目の上側ループの壁: (中心, 半径, 開始角, 全長, 始端余白, 終端余白)。"""
    _, info = X.symtvm()
    Cu = [np.array(c) for c in info["centres_upper"]]
    C = Cu[k]
    tip = C + np.array([X.R_IN / np.sin(np.radians(X.TIP_DEG / 2)), 0.0])
    straight = float(np.dot(tip - (C + X.R_IN * X.N_DN), X.D_DN))
    entry = np.array([info["x_in"], info["inlet_y"][1]]) if k == 0 else X._entry_pt(Cu, k)
    a_out = np.degrees(np.arctan2(*(entry - C)[::-1])) % 360.0
    a_in = 75.0 + np.degrees(np.arccos(info["chord_off"] / X.R_IN))
    out = (C, X.R_OUT, a_out, X.R_OUT * np.radians(a_out - 75) + straight, 0.15, 0.15)
    isl = (C, X.R_IN, a_in, X.R_IN * np.radians(a_in - 75) + straight, 0.15, 0.25)
    return C, out, isl, info


def centres(wall, inward):
    """壁に沿って SPACING ごとの半球中心。余白の中に左右均等に置く。"""
    C, R, a0, L, m0, m1 = wall
    n = int(np.floor((L - m0 - m1) / SPACING)) + 1
    s0 = m0 + 0.5 * (L - m0 - m1 - (n - 1) * SPACING)
    out = []
    for i in range(n):
        p, nrm = wall_point(C, R, a0, s0 + i * SPACING)
        nrm = -nrm if inward else nrm
        for z in ROWS:
            out.append((p[0], p[1], z, nrm[0], nrm[1]))
    return np.array(out)


def hemi_tris(c, r=R_HEMI, nu=30, nv=14):
    """中心 (x, y, z)、向き (nx, ny) の半球面。"""
    x, y, z, nx, ny = c
    n = np.array([nx, ny, 0.0])
    t = np.array([-ny, nx, 0.0])
    k = np.array([0.0, 0.0, 1.0])
    us, vs = np.linspace(0, 2 * np.pi, nu + 1), np.linspace(0, np.pi / 2, nv + 1)
    P = np.zeros((nv + 1, nu + 1, 3))
    for j, v in enumerate(vs):
        for i, u in enumerate(us):
            P[j, i] = np.array([x, y, z]) + r * (np.sin(v) * (np.cos(u) * t + np.sin(u) * k) + np.cos(v) * n)
    tris = []
    for j in range(nv):
        for i in range(nu):
            tris += [[P[j, i], P[j + 1, i], P[j + 1, i + 1]], [P[j, i], P[j + 1, i + 1], P[j, i + 1]]]
    return np.array(tris)


def prism(poly, z0, z1):
    """多角形を押し出した立体の側面と上面。"""
    tris = []
    for ring in [poly.exterior] + list(poly.interiors):
        Q = np.asarray(shapely.segmentize(ring, 0.03).coords)
        for a, b in zip(Q[:-1], Q[1:]):
            A0, B0, A1, B1 = [*a, z0], [*b, z0], [*a, z1], [*b, z1]
            tris += [[A0, B0, B1], [A0, B1, A1]]
    for g in shapely.constrained_delaunay_triangles(poly).geoms:
        q = np.asarray(g.exterior.coords)[:3]
        tris.append(np.c_[q, np.full(3, z1)])
    return np.array(tris, float)


def plate(x0, x1, y0, y1, z):
    return np.array([[[x0, y0, z], [x1, y0, z], [x1, y1, z]], [[x0, y0, z], [x1, y1, z], [x0, y1, z]]], float)


def render(groups, view, w=1500, h=1050, pad=0.05):
    """groups = [(三角形 N*3*3, 色)]。正射影の Z バッファ描画。"""
    d = np.asarray(view, float) / np.linalg.norm(view)
    right = np.cross(d, [0, 0, 1.0]); right /= np.linalg.norm(right)
    up = np.cross(right, d)
    allv = np.vstack([g[0].reshape(-1, 3) for g in groups])
    u, v = allv @ right, allv @ up
    span = max((u.max() - u.min()) / w, (v.max() - v.min()) / h) * (1 + 2 * pad)
    u0, v0 = 0.5 * (u.max() + u.min()), 0.5 * (v.max() + v.min())
    img = np.ones((h, w, 3))
    zb = np.full((h, w), np.inf)
    for tris, col in groups:
        px = (tris @ right - u0) / span + w / 2
        py = h / 2 - (tris @ up - v0) / span
        dp = tris @ d
        nrm = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
        nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
        nrm *= np.where((nrm @ d) > 0, -1.0, 1.0)[:, None]
        light = -d + 0.55 * up + 0.35 * right
        light = light / np.linalg.norm(light)
        shade = 0.40 + 0.70 * np.clip(nrm @ light, 0, 1)
        for i in range(len(tris)):
            x, y, z = px[i], py[i], dp[i]
            xa, xb = max(int(np.floor(x.min())), 0), min(int(np.ceil(x.max())) + 1, w)
            ya, yb = max(int(np.floor(y.min())), 0), min(int(np.ceil(y.max())) + 1, h)
            if xa >= xb or ya >= yb:
                continue
            X_, Y_ = np.meshgrid(np.arange(xa, xb) + .5, np.arange(ya, yb) + .5)
            den = (y[1] - y[2]) * (x[0] - x[2]) + (x[2] - x[1]) * (y[0] - y[2])
            if abs(den) < 1e-12:
                continue
            l0 = ((y[1] - y[2]) * (X_ - x[2]) + (x[2] - x[1]) * (Y_ - y[2])) / den
            l1 = ((y[2] - y[0]) * (X_ - x[2]) + (x[0] - x[2]) * (Y_ - y[2])) / den
            l2 = 1 - l0 - l1
            m = (l0 >= -1e-3) & (l1 >= -1e-3) & (l2 >= -1e-3)
            if not m.any():
                continue
            zz = l0 * z[0] + l1 * z[1] + l2 * z[2]
            sub = zb[ya:yb, xa:xb]
            m &= zz < sub
            sub[m] = zz[m]
            img[ya:yb, xa:xb][m] = np.clip(col * shade[i], 0, 1)
    return img


def main():
    C, out_w, isl_w, info = walls(STAGE)
    c_out, c_isl = centres(out_w, True), centres(isl_w, False)
    chord = info["chord_off"]
    island = X.teardrop(C, True, chord)
    x0, x1, y0, y1 = C[0] - 1.25, C[0] + 2.65, C[1] - 0.65, C[1] + 1.15

    # 外壁 = 外側の輪郭の外に幅 0.2 の帯
    a_out = out_w[2]
    arc = C + X.R_OUT * np.c_[np.cos(np.radians(np.linspace(a_out, 75, 80))), np.sin(np.radians(np.linspace(a_out, 75, 80)))]
    foot = C + X.R_OUT * X.N_DN + (out_w[3] - X.R_OUT * np.radians(a_out - 75)) * X.D_DN
    from shapely.geometry import LineString
    centre_line = LineString(np.vstack([arc, foot]))
    outer_band = centre_line.buffer(0.2, single_sided=True)
    if outer_band.distance(Point(C)) < 0.5:
        outer_band = centre_line.buffer(-0.2, single_sided=True)

    hemi_i = np.vstack([hemi_tris(c) for c in c_isl])
    hemi_o = np.vstack([hemi_tris(c) for c in c_out])
    floor = plate(x0, x1, y0, y1, Z_FLOOR)
    a = render([(floor, FLOOR), (prism(island, Z_FLOOR, Z_TOP), COPPER), (hemi_i, HEMI)], (.62, -.62, -.62))
    b = render([(floor, FLOOR), (prism(outer_band, Z_FLOOR, Z_TOP), COPPER), (hemi_o, HEMI)], (-.45, .75, -.6))

    # 平面図（流路の深さの中央 z = 1.25 の断面）
    poly, _ = X.symtvm()
    row = np.isclose(c_isl[:, 2], ROWS[1])
    discs = [Point(*q[:2]).buffer(R_HEMI, 64) for q in np.vstack([c_isl[row], c_out[np.isclose(c_out[:, 2], ROWS[1])]])]
    plotstyle.use_jp()
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1])
    ax = fig.add_subplot(gs[0, 0]); ax.imshow(a); ax.axis("off")
    ax.set_title("(a) 島の壁の半球（外壁を外して見た図）", loc="left", fontsize=12)
    ax = fig.add_subplot(gs[0, 1]); ax.imshow(b); ax.axis("off")
    ax.set_title("(b) 外壁の半球（島を外して内側から見た図）", loc="left", fontsize=12)
    ax = fig.add_subplot(gs[1, :])
    view = shapely.box(x0, y0, x1, y1)
    fluid = poly.intersection(view).difference(unary_union(discs))
    for g in (fluid.geoms if fluid.geom_type == "MultiPolygon" else [fluid]):
        P = np.asarray(g.exterior.coords)
        ax.fill(P[:, 0], P[:, 1], color="#bfe0e2", lw=0)
        for r in g.interiors:
            Q = np.asarray(r.coords)
            ax.fill(Q[:, 0], Q[:, 1], color="white", lw=0)
    for d_ in discs:
        P = np.asarray(d_.intersection(poly).exterior.coords)
        ax.fill(P[:, 0], P[:, 1], color=tuple(HEMI), lw=0)
    for r in [poly.exterior] + list(poly.interiors):
        Q = np.asarray(r.coords)
        ax.plot(Q[:, 0], Q[:, 1], color="#7a4a2b", lw=1.0)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
    ax.set_xlabel("流れ方向 x [mm]"); ax.set_ylabel("幅方向 [mm]")
    ax.set_title(f"(c) 平面図（流路の深さの中央の断面）  水色 = 水、濃い青緑 = 半球  半径 {R_HEMI} mm、中心間隔 {SPACING} mm、深さ方向 {len(ROWS)} 段",
                 loc="left", fontsize=12)
    fig.suptitle("パターン 2 のイメージ: ループの側壁に半球を等間隔に配置（ループ 1 個分、形状 CAD はまだ作っていない）", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(HERE / "preview_hemisphere.png", dpi=130)

    # 寸法の確認（全段）
    total = 0
    gaps = []
    for k in range(X.N_STAGE):
        _, o, i_, _ = walls(k)
        co, ci = centres(o, True), centres(i_, False)
        total += 2 * (len(co) + len(ci))
        rowo, rowi = co[np.isclose(co[:, 2], ROWS[1])], ci[np.isclose(ci[:, 2], ROWS[1])]
        dist = np.hypot(rowo[:, None, 0] - rowi[None, :, 0], rowo[:, None, 1] - rowi[None, :, 1])
        gaps.append(float(dist.min() - 2 * R_HEMI))
    report = dict(radius_mm=R_HEMI, spacing_mm=SPACING, rows_z_mm=[float(v) for v in ROWS],
                  per_stage_outer=len(c_out), per_stage_island=len(c_isl),
                  total_hemispheres_all_20_loops=total,
                  min_gap_between_facing_hemispheres_mm=min(gaps),
                  hemisphere_volume_share_of_fluid=float(total * (2 / 3) * np.pi * R_HEMI ** 3 / (poly.area * X.H_CH)))
    (HERE / "preview_numbers.json").write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
