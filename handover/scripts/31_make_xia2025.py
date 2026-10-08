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

OUT = os.path.join(HERE, "..", "cad", "xia2025")
FIG = os.path.join(HERE, "..", "figures", "geometry")
TAGS = {"SYMTVM-Forward": "xia2025_symtvm_fwd", "SYMTVM-Reverse": "xia2025_symtvm_rev",
        "TVM-Forward": "xia2025_tvm_fwd", "TVM-Reverse": "xia2025_tvm_rev"}


def write_dxf_ezdxf(path, rings, name):
    """閉じた LWPOLYLINE の DXF（R2010、単位 mm）。Fusion 360 の「DXF を挿入」用。

    `mesh.write_dxf` の最小実装は閉じフラグが読まれない環境があるので、ここでは ezdxf で書く。
    外形はレイヤ OUTER、島は ISLAND。
    """
    import ezdxf
    doc = ezdxf.new("R2010")
    doc.header["$INSUNITS"] = 4
    doc.layers.add("OUTER", color=5)
    doc.layers.add("ISLAND", color=1)
    msp = doc.modelspace()
    n = 0
    for i, ring in enumerate(rings):
        P = np.asarray(ring)
        if np.allclose(P[0], P[-1]):
            P = P[:-1]
        msp.add_lwpolyline([(float(x), float(y)) for x, y in P], close=True,
                           dxfattribs={"layer": "OUTER" if i == 0 else "ISLAND"})
        n += len(P)
    msp.add_text(f"{name}  Xia et al. 2025 ATE 258 124611 Fig.3/4  unit mm",
                 dxfattribs={"layer": "OUTER", "height": 0.5}).set_placement((0, 2.2))
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
    a = ap.parse_args()
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
        n_dxf = write_dxf_ezdxf(os.path.join(OUT, f"{tag}_outline.dxf"), rings, name)
        rec = dict(name=name, tag=tag, islands=len(poly.interiors), area_mm2=float(poly.area),
                   bounds=[float(v) for v in poly.bounds], dxf_points=n_dxf,
                   **{k: v for k, v in info.items() if k != "name"})
        if not a.no_stl:
            # Reverse は Forward のメッシュを z 軸まわりに 180° 回転する
            # （(x, y) → (Ly − x, −y)）。回転は面の向きを保つので頂点順はそのまま。
            if name.endswith("Reverse"):
                fwd = meshes[name.replace("Reverse", "Forward")]
                ff = [[(X.LY - p[0], -p[1], p[2]) for p in t] for t in fwd[0]]
                bf = [[(X.LY - p[0], -p[1], p[2]) for p in t] for t in fwd[1]]
            else:
                # 段ごとに x = 一定で切る（島の無い接続部）。穴の少ない小片にして
                # 三角形分割を安定させ、蓋と側壁の頂点を一致させる
                cuts = [X.X_C1 + X.PITCH * k + 3.0 for k in range(X.N_STAGE - 1)]
                ff = S3.fluid_solid_cut(poly, z_floor, z_top, cuts)
                bf = block_solid(poly, z_floor, z_top, h_half, X.LY, cuts)
            meshes[name] = (ff, bf)
            bad_f, _ = M.check_manifold(ff)
            nf = M.write_stl(os.path.join(OUT, f"{tag}_fluid.stl"), ff)
            vol_f = M.mesh_volume(ff)
            bad_b, _ = M.check_manifold(bf)
            nb = M.write_stl(os.path.join(OUT, f"{tag}_block.stl"), bf)
            vol_b = M.mesh_volume(bf)
            vol_expect = X.LY * X.W_UNIT * X.H_SOLID - vol_f
            rec.update(fluid_triangles=nf, fluid_volume_mm3=float(vol_f),
                       block_triangles=nb, block_volume_mm3=float(vol_b),
                       block_volume_expected=float(vol_expect),
                       watertight=bool(bad_f == 0 and bad_b == 0))
            print(f"  {name:16s} 島 {len(poly.interiors):2d}  流体 {nf:6d} 三角形 {vol_f:7.2f} mm^3  "
                  f"ブロック {nb:6d} 三角形 {vol_b:7.2f} mm^3（期待 {vol_expect:.2f}）  "
                  f"{'OK' if rec['watertight'] else f'**非多様体 {bad_f}/{bad_b}**'}")
        recs.append(rec)

    j = os.path.join(HERE, "..", "results", "cad", "cad_xia2025.json")
    os.makedirs(os.path.dirname(j), exist_ok=True)
    json.dump(dict(source="Xia et al., Appl. Therm. Eng. 258 (2025) 124611, Fig.3/Fig.4/Table 2",
                   params=dict(w_ch=X.W_CH, R_in=X.R_IN, R_out=X.R_OUT, tip_deg=X.TIP_DEG,
                               pitch=X.PITCH, x_c1=X.X_C1, cy=float(X.CY), band_deg=X.ANG,
                               Ly=X.LY, W=X.W_UNIT, H_ch=X.H_CH, H=X.H_SOLID, n_stage=X.N_STAGE),
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
