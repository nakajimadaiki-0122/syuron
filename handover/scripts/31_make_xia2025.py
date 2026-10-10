#!/usr/bin/env python3
"""
31_make_xia2025.py -- Xia et al. (2025) Fig. 3 の 4 構成を CAD に出す
=====================================================================

形状は `src/xia2025.py`（寸法の来歴は `docs/xia2025_geometry.md`）。

4 構成: SYMTVM-Forward / SYMTVM-Reverse / TVM-Forward / TVM-Reverse。
Reverse は Forward の 180° 回転（同じ流路を逆向きに流し、入口を左に戻したもの）。

出力（`cad/xia2025/`、単位 mm。Fusion 360 で開ける形式）

| ファイル | 中身 |
|---|---|
| `{名前}_fluid.stl` | 流体体積（流路を深さ 1.5 で押し出し。z = 0.5〜2.0） |
| `{名前}_block.stl` | 計算単位のブロック 50 × 3.54 × 2.0 に流路を深さ 1.5 で彫ったもの。入口・出口は端面に開口 |
| `{名前}_outline.dxf` | 流路の 2D 輪郭（Fusion で「挿入 > DXF を挿入」→ 押し出し 1.5 mm） |

座標の約束（**2026-10-10 変更**）
-----------------------------------

既定は **Y 軸が鉛直（Y-up）**。板は XZ 平面に寝て、厚み方向が +Y になる。
`--up z` を付けると従来どおり Z 厚み（板が XY 平面）で出る。

    Y-up への変換は X 軸まわり −90 度の回転   (x, y, z) -> (x, z, -y)

回転なので三角形の向き（法線）は保たれる。流路の進行方向は +X、
流路幅の方向は Z、彫り込みの深さは −Y（板の上面が y = 2.0、流路の底が y = 0.5）。

DXF も同じ平面に載せる。LWPOLYLINE を押し出し方向 (0, -1, 0) の OCS で書くので、
CAD 上では XZ 平面（地面）に乗る。**押し出し方向を読まない簡易ビューアでは
y が反転した 2 次元図形に見える**（形は正しい）。

全モデルで「各辺がちょうど 2 枚の三角形に共有される」ことを確認する。
STEP は出せない（cadquery が Python 3.14 に未対応）。Fusion で面付きソリッドが
要るときは DXF を押し出す方が速い。
"""
import argparse, json, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, box

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import xia2025 as X
import mesh as M
import solid3d as S3
import plotstyle

def to_yup(faces):
    """Z 厚みの三角形列を Y-up へ回す。X 軸まわり −90 度: (x, y, z) -> (x, z, −y)。

    回転（行列式 +1）なので頂点順＝法線の向きはそのまま使える。
    """
    return [[(p[0], p[2], -p[1]) for p in t] for t in faces]


OUT = os.path.join(HERE, "..", "cad", "xia2025")
UP = "y"                               # main() が --up で上書きする（render3d 用）
FIG = os.path.join(HERE, "..", "figures", "geometry")
TAGS = {"SYMTVM-Forward": "xia2025_symtvm_fwd", "SYMTVM-Reverse": "xia2025_symtvm_rev",
        "TVM-Forward": "xia2025_tvm_fwd", "TVM-Reverse": "xia2025_tvm_rev"}


def write_dxf_ezdxf(path, rings, name, up="y"):
    """閉じた LWPOLYLINE の DXF（R2010、単位 mm）。Fusion 360 の「DXF を挿入」用。

    `mesh.write_dxf` の最小実装は閉じフラグが読まれない環境があるので、ここでは ezdxf で書く。
    外形はレイヤ OUTER、島は ISLAND。

    up="y" のときは押し出し方向を (0, −1, 0) にして **XZ 平面（地面）** に載せる。
    DXF の任意軸アルゴリズムにより OCS の (u, v) が世界座標 (u, 0, v) に対応するので、
    v = −y と書けば STL（`to_yup`）と同じ位置・同じ向きになる。
    """
    import ezdxf
    doc = ezdxf.new("R2010")
    doc.header["$INSUNITS"] = 4
    doc.layers.add("OUTER", color=5)
    doc.layers.add("ISLAND", color=1)
    msp = doc.modelspace()
    sy = -1.0 if up == "y" else 1.0           # OCS 上の v = −y（up="y" のとき）
    ext = {"extrusion": (0.0, -1.0, 0.0)} if up == "y" else {}
    n = 0
    for i, ring in enumerate(rings):
        P = np.asarray(ring)
        if np.allclose(P[0], P[-1]):
            P = P[:-1]
        msp.add_lwpolyline([(float(x), sy * float(y)) for x, y in P], close=True,
                           dxfattribs={"layer": "OUTER" if i == 0 else "ISLAND", **ext})
        n += len(P)
    msp.add_text(f"{name}  Xia et al. 2025 ATE 258 124611 Fig.3/4  unit mm"
                 + ("  [Y-up: XZ plane]" if up == "y" else ""),
                 dxfattribs={"layer": "OUTER", "height": 0.5, **ext}
                 ).set_placement((0.0, sy * 2.2))
    doc.saveas(path)
    return n


def _end_face(x0, y_lo, y_hi, z_floor, z_top, h_half, outward_x):
    """端面（x = x0）。全高の長方形から流路の開口を抜いた面を三角形にする。"""
    prof = Polygon([(-h_half, 0.0), (h_half, 0.0), (h_half, z_top), (y_hi, z_top),
                    (y_hi, z_floor), (y_lo, z_floor), (y_lo, z_top), (-h_half, z_top)])
    out = []
    for t in M.triangulate_polygon(prof):
        v = [(x0, float(p[0]), float(p[1])) for p in t]
        a, b, c = (np.array(q) for q in v)
        n = np.cross(b - a, c - a)
        if (n[0] > 0) != (outward_x > 0):
            v = v[::-1]
        out.append(v)
    return out


def _end_face_hole(x0, y_lo, y_hi, z_floor, z_top, h_half, z_lid, outward_x):
    """端面（x = x0）。蓋を載せた版。流路の開口は面の**内側の穴**になる。

    蓋が無い版（`_end_face`）では開口が上端まで抜けるので外形の切り欠きになるが、
    蓋を載せると周囲が固体で囲まれるため穴として扱う。
    """
    shell = [(-h_half, 0.0), (h_half, 0.0), (h_half, z_lid), (-h_half, z_lid)]
    hole = [(y_lo, z_floor), (y_hi, z_floor), (y_hi, z_top), (y_lo, z_top)]
    prof = Polygon(shell, [hole])
    out = []
    for t in M.triangulate_polygon(prof):
        v = [(x0, float(q[0]), float(q[1])) for q in t]
        a, b, c = (np.array(u) for u in v)
        n = np.cross(b - a, c - a)
        if (n[0] > 0) != (outward_x > 0):
            v = v[::-1]
        out.append(v)
    return out


def lid_solid(length, h_half, z0, t):
    """蓋（平板）。z0..z0+t の直方体。流路を塞ぐカバー。"""
    return _box(0.0, length, -h_half, h_half, z0, z0 + t)


def _box(x0, x1, y0, y1, z0, z1):
    """軸に沿った直方体（法線は外向き）。"""
    P = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    out = []
    for a, b, c, d in quads:
        out += [[P[a], P[b], P[c]], [P[a], P[c], P[d]]]
    return out


def capped_solid(poly, z_floor, z_top, h_half, length, cuts, t_lid):
    """蓋を一体にした固体。流路は内部の空洞になり、端面にだけ開口する。

    `block_solid` との違いは 3 点。
      - 上面が z_top + t_lid の**全面**（流路の外だけ、ではない）
      - 流路の天井（z_top の上に蓋が載る面、法線は下向き）を足す
      - 端面の開口が外形の切り欠きではなく**穴**になる
    """
    poly = M.clean_polygon(poly)
    z_lid = z_top + t_lid
    rect = box(0.0, -h_half, length, h_half)
    faces = S3.caps_cut(rect, 0.0, False, cuts)        # 底面
    faces += S3.caps_cut(rect, z_lid, True, cuts)      # 蓋の上面
    faces += S3.caps_cut(poly, z_floor, True, cuts)    # 流路の底
    faces += S3.caps_cut(poly, z_top, False, cuts)     # 流路の天井（蓋の裏）
    for g in S3._pieces(poly, cuts):                   # 溝の壁
        rings = [S3.closed(g.exterior, True)] + [S3.closed(r, False) for r in g.interiors]
        for P in rings:
            for a, b in zip(P[:-1], P[1:]):
                if S3._on_cut(a, b, cuts) or S3._on_cut(a, b, [0.0, length]):
                    continue
                q = [(a[0], a[1], z_floor), (b[0], b[1], z_floor),
                     (b[0], b[1], z_top), (a[0], a[1], z_top)]
                faces += [[q[0], q[2], q[1]], [q[0], q[3], q[2]]]
    xs = [0.0] + sorted(cuts) + [length]               # 外周 y = ±h
    for xa, xb in zip(xs[:-1], xs[1:]):
        for y in (-h_half, h_half):
            q = [(xa, y, 0.0), (xb, y, 0.0), (xb, y, z_lid), (xa, y, z_lid)]
            tq = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
            faces += [u[::-1] for u in tq] if y > 0 else tq
    xs_, ys_ = np.asarray(poly.exterior.coords).T      # 端面（開口は穴）
    m0 = np.abs(xs_) < 1e-9
    m1 = np.abs(xs_ - length) < 1e-9
    faces += _end_face_hole(0.0, ys_[m0].min(), ys_[m0].max(), z_floor, z_top,
                            h_half, z_lid, -1)
    faces += _end_face_hole(length, ys_[m1].min(), ys_[m1].max(), z_floor, z_top,
                            h_half, z_lid, +1)
    return faces


def block_solid(poly, z_floor, z_top, h_half, length, cuts):
    """流路を彫ったブロック。入口・出口は端面に開口する。

    蓋・溝の壁はすべて同じ x = cuts で切った小片から作るので頂点が一致する。
    """
    poly = M.clean_polygon(poly)
    rect = box(0.0, -h_half, length, h_half)
    faces = S3.caps_cut(rect, 0.0, False, cuts)                     # 底面
    faces += S3.caps_cut(rect.difference(poly), z_top, True, cuts)  # 上面（流路の外）
    faces += S3.caps_cut(poly, z_floor, True, cuts)                 # 流路の底
    # 溝の壁（端面 x=0, x=length 上の辺は除く。法線は溝の内側＝固体の外向き）
    for g in S3._pieces(poly, cuts):
        rings = [S3.closed(g.exterior, True)] + [S3.closed(r, False) for r in g.interiors]
        for P in rings:
            for a, b in zip(P[:-1], P[1:]):
                if S3._on_cut(a, b, cuts) or S3._on_cut(a, b, [0.0, length]):
                    continue
                q = [(a[0], a[1], z_floor), (b[0], b[1], z_floor), (b[0], b[1], z_top), (a[0], a[1], z_top)]
                faces += [[q[0], q[2], q[1]], [q[0], q[3], q[2]]]
    # 外周 y = ±h の面（切断区間ごと）
    xs = [0.0] + sorted(cuts) + [length]
    for xa, xb in zip(xs[:-1], xs[1:]):
        for y in (-h_half, h_half):
            q = [(xa, y, 0.0), (xb, y, 0.0), (xb, y, z_top), (xa, y, z_top)]
            t = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]        # 法線 −y
            faces += [u[::-1] for u in t] if y > 0 else t
    # 端面（開口つき）
    xs_, ys_ = np.asarray(poly.exterior.coords).T
    m0 = np.abs(xs_) < 1e-9
    m1 = np.abs(xs_ - length) < 1e-9
    faces += _end_face(0.0, ys_[m0].min(), ys_[m0].max(), z_floor, z_top, h_half, -1)
    faces += _end_face(length, ys_[m1].min(), ys_[m1].max(), z_floor, z_top, h_half, +1)
    return faces


def _orient_check(faces):
    """全三角形の符号付き体積が正になる向きに揃っているか（閉じていれば体積 > 0）。"""
    return M.mesh_volume(faces)


def draw(ax, poly, title):
    P = np.asarray(poly.exterior.coords)
    ax.fill(P[:, 0], P[:, 1], color="#9fc9cb")
    for r in poly.interiors:
        Q = np.asarray(r.coords)
        ax.fill(Q[:, 0], Q[:, 1], color="white")
        ax.plot(Q[:, 0], Q[:, 1], color="#0e7c7b", lw=0.5)
    ax.plot(P[:, 0], P[:, 1], color="#0e7c7b", lw=0.6)
    ax.add_patch(plt.Rectangle((0, -0.5 * X.W_UNIT), X.LY, X.W_UNIT, fill=False, lw=0.6, ls=":"))
    ax.set_aspect("equal")
    ax.set_xlim(-0.5, X.LY + 0.5)
    ax.set_ylim(-0.5 * X.W_UNIT - 0.1, 0.5 * X.W_UNIT + 0.1)
    ax.set_title(title, fontsize=9, loc="left")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-stl", action="store_true")
    ap.add_argument("--lid-mm", type=float, default=1.0,
                    help="蓋（カバー板）の厚み [mm]。既定 1.0")
    ap.add_argument("--up", choices=["y", "z"], default="y",
                    help="鉛直にする軸。既定 y（板は XZ 平面に寝る）")
    a = ap.parse_args()
    global UP
    UP = a.up
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    plotstyle.use_jp()
    z_floor = X.H_SOLID - X.H_CH          # 0.5
    z_top = X.H_SOLID                     # 2.0
    h_half = 0.5 * X.W_UNIT

    four = X.all_four()
    recs = []
    meshes = {}
    print("=" * 96)
    print(f" Xia et al. (2025) Fig.3 の 4 構成   流路幅 {X.W_CH}, R_in {X.R_IN}, R_out {X.R_OUT}, "
          f"ピッチ {X.PITCH}, cy = {X.CY:.3f}, Ly {X.LY}, W {X.W_UNIT}, Hch {X.H_CH}, H {X.H_SOLID}")
    print("=" * 96)
    for name, (poly, info) in four.items():
        tag = TAGS[name]
        rings = M.polygon_rings(M.clean_polygon(poly))
        n_dxf = write_dxf_ezdxf(os.path.join(OUT, f"{tag}_outline.dxf"), rings, name, a.up)
        rec = dict(name=name, tag=tag, islands=len(poly.interiors), area_mm2=float(poly.area),
                   bounds=[float(v) for v in poly.bounds], dxf_points=n_dxf,
                   **{k: v for k, v in info.items() if k != "name"})
        if not a.no_stl:
            # Reverse は Forward のメッシュを z 軸まわりに 180° 回転する
            # （(x, y) → (Ly − x, −y)）。回転は面の向きを保つので頂点順はそのまま。
            if name.endswith("Reverse"):
                fwd = meshes[name.replace("Reverse", "Forward")]
                rot = lambda F: [[(X.LY - q[0], -q[1], q[2]) for q in t] for t in F]
                ff, bf, cf = (rot(fwd[0]), rot(fwd[1]), rot(fwd[2]))
            else:
                # 段ごとに x = 一定で切る（島の無い接続部）。穴の少ない小片にして
                # 三角形分割を安定させ、蓋と側壁の頂点を一致させる
                cuts = [X.X_C1 + X.PITCH * k + 3.0 for k in range(X.N_STAGE - 1)]
                ff = S3.fluid_solid_cut(poly, z_floor, z_top, cuts)
                bf = block_solid(poly, z_floor, z_top, h_half, X.LY, cuts)
                cf = capped_solid(poly, z_floor, z_top, h_half, X.LY, cuts, a.lid_mm)
            meshes[name] = (ff, bf, cf)      # 以降の Reverse 生成は Z 厚みのまま使う
            bad_f, _ = M.check_manifold(ff)
            bad_b, _ = M.check_manifold(bf)
            bad_c, _ = M.check_manifold(cf)
            vol_f = M.mesh_volume(ff)
            vol_b = M.mesh_volume(bf)
            vol_c = M.mesh_volume(cf)
            # 書き出す直前に Y-up へ回す（検査は回転不変なので順序はどちらでもよい）
            tr = to_yup if a.up == "y" else (lambda F: F)
            ff_w, bf_w, cf_w = tr(ff), tr(bf), tr(cf)
            nf = M.write_stl(os.path.join(OUT, f"{tag}_fluid.stl"), ff_w)
            nb = M.write_stl(os.path.join(OUT, f"{tag}_block.stl"), bf_w)
            nc = M.write_stl(os.path.join(OUT, f"{tag}_capped.stl"), cf_w)
            vol_expect = X.LY * X.W_UNIT * X.H_SOLID - vol_f
            rec.update(fluid_triangles=nf, fluid_volume_mm3=float(vol_f),
                       block_triangles=nb, block_volume_mm3=float(vol_b),
                       block_volume_expected=float(vol_expect),
                       capped_triangles=nc, capped_volume_mm3=float(vol_c),
                       capped_volume_expected=float(vol_expect + X.LY * X.W_UNIT * a.lid_mm),
                       lid_mm=a.lid_mm, up_axis=a.up,
                       stl_bounds=[float(v) for v in
                                   np.asarray(ff_w).reshape(-1, 3).min(0)]
                                  + [float(v) for v in
                                     np.asarray(ff_w).reshape(-1, 3).max(0)],
                       watertight=bool(bad_f == 0 and bad_b == 0 and bad_c == 0))
            vc_exp = vol_expect + X.LY * X.W_UNIT * a.lid_mm
            print(f"  {name:16s} 島 {len(poly.interiors):2d}  "
                  f"流体 {nf:6d}面 {vol_f:7.2f}  "
                  f"ブロック {nb:6d}面 {vol_b:7.2f}（期待 {vol_expect:.2f}）  "
                  f"蓋付き {nc:6d}面 {vol_c:7.2f}（期待 {vc_exp:.2f}）  "
                  f"{'OK' if rec['watertight'] else f'**非多様体 {bad_f}/{bad_b}/{bad_c}**'}")
        recs.append(rec)

    if not a.no_stl:
        lf = lid_solid(X.LY, h_half, z_top, a.lid_mm)
        bad_l, _ = M.check_manifold(lf)
        lf_w = to_yup(lf) if a.up == "y" else lf
        nl = M.write_stl(os.path.join(OUT, "xia2025_lid.stl"), lf_w)
        print(f"  {'蓋（4 構成で共通）':16s} {nl:6d}面 "
              f"{M.mesh_volume(lf):7.2f} mm^3  "
              f"{X.LY} x {X.W_UNIT} x {a.lid_mm} mm  "
              f"{'OK' if bad_l == 0 else f'**非多様体 {bad_l}**'}")

    j = os.path.join(HERE, "..", "results", "cad", "cad_xia2025.json")
    os.makedirs(os.path.dirname(j), exist_ok=True)
    json.dump(dict(source="Xia et al., Appl. Therm. Eng. 258 (2025) 124611, Fig.3/Fig.4/Table 2",
                   params=dict(w_ch=X.W_CH, R_in=X.R_IN, R_out=X.R_OUT, tip_deg=X.TIP_DEG,
                               pitch=X.PITCH, x_c1=X.X_C1, cy=float(X.CY), band_deg=X.ANG,
                               Ly=X.LY, W=X.W_UNIT, H_ch=X.H_CH, H=X.H_SOLID, n_stage=X.N_STAGE,
                               lid_mm=a.lid_mm, up_axis=a.up),
                   models=recs), open(j, "w"), indent=1, ensure_ascii=False)
    print(f"saved {os.path.normpath(j)}")

    fig, axes = plt.subplots(4, 1, figsize=(16, 9))
    for ax, (name, (poly, _)) in zip(axes, four.items()):
        draw(ax, poly, f"{name}   ({TAGS[name]})")
    axes[-1].set_xlabel("x [mm]")
    fig.suptitle("Xia et al. (2025) Fig.3 の 4 構成（点線は計算単位 50 × 3.54 mm）", fontsize=11)
    fig.tight_layout()
    png = os.path.join(FIG, "xia2025_four.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")

    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    for ax, name in zip(axes, ("SYMTVM-Forward", "TVM-Forward")):
        poly, _ = four[name]
        P = np.asarray(poly.exterior.coords)
        ax.fill(P[:, 0], P[:, 1], color="#9fc9cb")
        for r in poly.interiors:
            Q = np.asarray(r.coords)
            ax.fill(Q[:, 0], Q[:, 1], color="#8a8a8a")
            ax.plot(Q[:, 0], Q[:, 1], color="k", lw=0.6)
        ax.plot(P[:, 0], P[:, 1], color="k", lw=0.6)
        ax.axhline(0, color="r", lw=0.5, ls="--")
        ax.set_aspect("equal")
        ax.set_xlim(-0.2, 12.5)
        ax.set_ylim(-1.77, 1.77)
        ax.set_xticks(np.arange(0, 13))
        ax.grid(alpha=0.3)
        ax.set_title(f"{name}  先頭 2〜3 段（Fig.4 と比較する用）", fontsize=10)
    fig.tight_layout()
    png = os.path.join(FIG, "xia2025_detail.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")


def _read_faces(path):
    import struct
    b = open(path, "rb").read()
    n = struct.unpack("<I", b[80:84])[0]
    a = np.frombuffer(b[84:], dtype=np.uint8).reshape(n, 50)
    return a[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype(float)


def render3d():
    """fluid / block の等角投影（深さ方向 3 倍）。figures/geometry/xia2025_3d.png"""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.colors as mc
    fig = plt.figure(figsize=(14, 3.1 * len(TAGS)))
    L = np.array([0.3, -0.5, 0.8]); L /= np.linalg.norm(L)
    for i, (name, tag) in enumerate(TAGS.items()):
        for j, (role, col, zs) in enumerate((("block", "#c9c2b4", 3.0), ("fluid", "#5fa8a6", 3.0))):
            tri = _read_faces(os.path.join(OUT, f"{tag}_{role}.stl"))
            if UP == "y":           # 描画用に Z-up へ戻す（matplotlib の 3D は Z 上）
                tri = np.stack([tri[:, :, 0], -tri[:, :, 2], tri[:, :, 1]], axis=-1)
            ax = fig.add_subplot(len(TAGS), 2, 2 * i + j + 1, projection="3d")
            n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-15
            k = 0.45 + 0.55 * np.clip(n @ L, 0, 1)
            rgb = np.clip(np.asarray(mc.to_rgb(col))[None, :] * k[:, None], 0, 1)
            t = tri.copy(); t[:, :, 2] *= zs
            ax.add_collection3d(Poly3DCollection(t, facecolors=rgb, edgecolors="none"))
            lo, hi = t.reshape(-1, 3).min(0), t.reshape(-1, 3).max(0); ext = hi - lo
            ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
            ax.set_box_aspect(tuple(ext / ext.max()), zoom=1.3)
            ax.view_init(elev=45, azim=-62); ax.set_axis_off()
            ax.set_title(f"{name}  {role}   {ext[0]:.1f} × {ext[1]:.2f} × {ext[2]/zs:.1f} mm",
                         fontsize=9, loc="left", pad=0)
    fig.suptitle("Xia et al. (2025) Fig.3 の 4 構成 — 左: 流路を彫ったブロック、右: 流体体積（深さ方向 3 倍）", fontsize=11)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01, hspace=0.05, wspace=0.02)
    png = os.path.join(FIG, "xia2025_3d.png")
    fig.savefig(png, dpi=120)
    print(f"saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
    if "--no-stl" not in sys.argv:
        render3d()
