#!/usr/bin/env python3
"""
15_gamboa_full_geom.py -- 単発バルブ全体（プレナム込み）の生成・検証・CAD 出力
==============================================================================

CLAUDE.md の方針「幾何を近似で作らない」に従い、CFD の前に必ず通す。

検証
  1. 流路幅が全域で w か（骨格点から直接確認できる量を表示）
  2. **Fig.2 との重ね合わせ**（IoU、はみ出しの内訳）。plenum_out="faithful" で比較
  3. 接合部が同一直線か（生成形状から測り直す。Fig.2 実測と突き合わせる）

出力
  figures/gamboa_full.png        生成形状と Fig.2 の重ね
  results/gamboa_full_geom.json  骨格点・検証値
  cad/valve1_std_outline.dxf      2D 輪郭（プレナム込み、単位 mm）
  cad/valve1bare_std_outline.dxf  プレナムを外した版

命名規則は `{形状}_{こぶの大きさ}_{役割}`。cad/README.md の一覧表を参照。

（--which reference で作ると同じ名前で上書きされる。使い分けるなら退避すること）

STL（3D）は `20_make_3d_model.py` が作る。三角形分割は耳刈りで境界辺を保つので、
閉じたメッシュになる（ここでの Delaunay 版は非多様体辺が出ていた）。
"""
import argparse, json, os, struct, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from shapely.geometry import Polygon
from shapely.ops import triangulate
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import gamboa_full as GF

FIG2 = os.path.join(HERE, "..", "docs", "gamboa2005_fig2.png")


# ----------------------------------------------------------------- Fig.2 比較
def fig2_mask():
    a = mpimg.imread(FIG2)
    a = a if a.ndim == 2 else a[..., 0]
    lab, n = ndi.label(a < 0.5, structure=np.ones((3, 3)))
    sz = ndi.sum(a < 0.5, lab, range(1, n + 1))
    oid = int(np.argmax(sz)) + 1
    outer = ndi.binary_fill_holes(lab == oid)
    for s, k in sorted(((sz[k - 1], k) for k in range(1, n + 1) if k != oid),
                       reverse=True):
        ys, xs = np.where(lab == k)
        if outer[ys, xs].all():
            island = ndi.binary_fill_holes(lab == k)
            break
    return outer & ~island


def overlay(poly, scale):
    """生成形状を Fig.2 の画素格子へ載せて IoU を出す。"""
    from matplotlib.path import Path
    ref = fig2_mask()
    H, W = ref.shape
    w_px, x0, yc = scale
    X, Y = np.meshgrid(np.arange(W) + 0.5, np.arange(H) + 0.5)
    # 画素座標 -> w_v 座標
    U = (X - x0) / w_px
    V = -(Y - yc) / w_px
    pts = np.column_stack([U.ravel(), V.ravel()])
    gen = Path(np.asarray(poly.exterior.coords)).contains_points(pts)
    for ring in poly.interiors:
        gen &= ~Path(np.asarray(ring.coords)).contains_points(pts)
    gen = gen.reshape(H, W)
    inter = (gen & ref).sum()
    union = (gen | ref).sum()
    return dict(iou=float(inter / union), only_fig2=int((ref & ~gen).sum()),
                only_gen=int((gen & ~ref).sum()), inter=int(inter)), gen, ref


# --------------------------------------------------------------------- CAD 出力
def write_dxf(path, rings, unit_note=""):
    """閉じた LWPOLYLINE として書き出す（依存を増やさない最小実装）。"""
    def hdr():
        return ("999\n" + f"gamboa valve {unit_note}\n"
                "0\nSECTION\n2\nHEADER\n"
                "9\n$INSUNITS\n70\n4\n"          # 4 = millimeters
                "0\nENDSEC\n"
                "0\nSECTION\n2\nENTITIES\n")
    body = [hdr()]
    for ring in rings:
        P = np.asarray(ring)
        body.append("0\nLWPOLYLINE\n8\n0\n100\nAcDbEntity\n100\nAcDbPolyline\n")
        body.append(f"90\n{len(P)}\n70\n1\n")     # 70=1: closed
        for x, y in P:
            body.append(f"10\n{x:.6f}\n20\n{y:.6f}\n")
    body.append("0\nENDSEC\n0\nEOF\n")
    open(path, "w").write("".join(body))


def densify_ring(P, step):
    out = []
    for a, b in zip(P[:-1], P[1:]):
        L = float(np.hypot(*(b - a)))
        k = max(int(np.ceil(L / step)), 1)
        out.append(a + (np.arange(k)[:, None] / k) * (b - a))
    out.append(P[-1:])
    return np.vstack(out)


def triangulate_polygon(poly, step=None):
    """穴つき多角形を三角形分割する（境界を細分した Delaunay + 内部判定）。

    非拘束 Delaunay なので、境界が粗いと凹部で三角形が落ちて面積が欠ける。
    境界を step 間隔に細分すると拘束 Delaunay に十分近づく。
    """
    if step:
        rings = [densify_ring(np.asarray(poly.exterior.coords), step)]
        rings += [densify_ring(np.asarray(r.coords), step) for r in poly.interiors]
        poly = Polygon(rings[0], rings[1:])
    tris = [t for t in triangulate(poly) if poly.contains(t.representative_point())]
    area = sum(t.area for t in tris)
    return tris, area


def write_stl(path, poly, depth, name="gamboa", step=None):
    """2D 断面を深さ方向へ押し出して binary STL を書く。"""
    tris, area = triangulate_polygon(poly, step=step)
    if abs(area / poly.area - 1) > 1e-3:
        print(f"   警告: 三角形分割の面積が {100*(area/poly.area-1):+.3f} % ずれる")
    faces = []
    for t in tris:
        c = np.asarray(t.exterior.coords)[:3]
        faces.append([(x, y, 0.0) for x, y in c][::-1])          # 下面（法線 -z）
        faces.append([(x, y, depth) for x, y in c])              # 上面（法線 +z）
    for ring in [poly.exterior] + list(poly.interiors):
        P = np.asarray(ring.coords)
        for a, b in zip(P[:-1], P[1:]):
            faces.append([(a[0], a[1], 0.0), (b[0], b[1], 0.0),
                          (b[0], b[1], depth)])
            faces.append([(a[0], a[1], 0.0), (b[0], b[1], depth),
                          (a[0], a[1], depth)])
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(struct.pack("<I", len(faces)))
        for tri in faces:
            v = np.array(tri, float)
            nvec = np.cross(v[1] - v[0], v[2] - v[0])
            ln = np.linalg.norm(nvec)
            nvec = nvec / ln if ln > 0 else np.zeros(3)
            f.write(struct.pack("<3f", *nvec))
            for p in v:
                f.write(struct.pack("<3f", *p))
            f.write(struct.pack("<H", 0))
    return len(faces), area


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="optimized")
    ap.add_argument("--wv-mm", type=float, default=1.0,
                    help="CAD 出力での流路幅 w_v [mm]")
    ap.add_argument("--depth-mm", type=float, default=None,
                    help="押し出し深さ [mm]。既定は w_v（正方形断面）")
    ap.add_argument("--arc-pts", type=int, default=240, help="CAD 用の円弧分割数")
    a = ap.parse_args()

    # --- 無次元（w_v = 1）で組んで検証 ---
    poly_f, info_f = GF.build(a.which, w=1.0, plenum_out="faithful")
    poly_v, info_v = GF.build(a.which, w=1.0, plenum_out="vertical")
    print("=" * 78)
    print(f" Gamboa {a.which} 単発バルブ（プレナム込み）")
    print("=" * 78)
    for k in ("theta_deg", "beta_deg", "alpha_deg", "splitter_tip", "centre",
              "outlet_tangent", "E", "Q", "island_tip", "J", "outlet_entry",
              "plenum_left_centre", "plenum_right_centre", "bounds"):
        v = info_f[k]
        print(f"  {k:22s} " + (f"{v:.4f}" if isinstance(v, float)
                               else str([round(x, 4) for x in v]) if isinstance(v, list)
                               else str(v)))
    print(f"  流体面積 (faithful) = {poly_f.area:.4f} w_v^2, "
          f"(vertical) = {poly_v.area:.4f} w_v^2")

    # --- 1. 接合部の同一直線性（生成形状から）---
    e = np.array([np.sin(np.radians(info_f['alpha_deg'])),
                  -np.cos(np.radians(info_f['alpha_deg']))])
    ang = np.degrees(np.arctan2(e[1], e[0]))
    print(f"\n1. 出口区間の向き = 主流路軸から {abs(ang):.3f} 度 "
          f"（Fig.2 実測 48.03 度、Table 2 の alpha = {info_f['alpha_deg']} は法線から）")
    print(f"   ループ下流枝と出口区間は**構成上同一直線**（同じ直線 L_out 上に "
          f"T2, E, Q, 円弧接点がある）")
    print(f"   主流路下壁の終端 Q = x {info_f['Q'][0]:.3f} w_v "
          f"（Fig.2 実測: 下壁の水平は 6.94 まで、帯は 7.90 で消える）")

    # --- 2. Fig.2 との重ね合わせ ---
    scale = (44.76, 246.38, 271.38)      # 14_fig2_plenum.py の実測
    ov, gen, ref = overlay(poly_f, scale)
    print(f"\n2. Fig.2 との重ね合わせ（scale: w_v = {scale[0]} px）")
    print(f"   IoU = {ov['iou']:.4f}   共通 {ov['inter']} px   "
          f"Fig.2 のみ {ov['only_fig2']} px   生成のみ {ov['only_gen']} px")
    print(f"   （周期版 gamboa.py の IoU は 0.878〜0.898。"
          f"プレナムを含むぶん面積が大きいので直接比較はできない）")

    # --- 3. CAD 出力 ---
    wv = a.wv_mm
    depth = a.depth_mm if a.depth_mm else wv
    cad = os.path.join(HERE, "..", "cad")
    os.makedirs(cad, exist_ok=True)
    poly_mm, info_mm = GF.build(a.which, w=wv, plenum_out="faithful",
                                arc_pts=a.arc_pts)
    rings = [np.asarray(poly_mm.exterior.coords)] + \
            [np.asarray(r.coords) for r in poly_mm.interiors]
    f_dxf = os.path.join(cad, "valve1_std_outline.dxf")
    write_dxf(f_dxf, rings, unit_note=f"w_v = {wv} mm")
    print(f"\n3. CAD 出力（w_v = {wv} mm, 押し出し {depth} mm）")
    print(f"   {os.path.normpath(f_dxf)}  "
          f"（外形 {len(rings[0])} 点 + 島 {len(rings)-1} 個）")

    # プレナムを外した版（バルブ本体だけ）
    poly_valve, _ = GF.build(a.which, w=wv, plenum_out="faithful",
                             arc_pts=a.arc_pts, include_plenums=False)
    rings_v = [np.asarray(poly_valve.exterior.coords)] + \
              [np.asarray(r.coords) for r in poly_valve.interiors]
    f_dxf2 = os.path.join(cad, "valve1bare_std_outline.dxf")
    write_dxf(f_dxf2, rings_v, unit_note=f"valve only, w_v = {wv} mm")
    print(f"   {os.path.normpath(f_dxf2)}")
    print("   STL（3D モデル）は scripts/20_make_3d_model.py で作る")

    # --- 図 ---
    fig, axes = plt.subplots(2, 1, figsize=(11, 12))
    ax = axes[0]
    for pol, col, lab in ((poly_f, "C0", "faithful (Fig.2)"),
                          (poly_v, "C1", "vertical outlet plane (CFD)")):
        P = np.asarray(pol.exterior.coords)
        ax.plot(P[:, 0], P[:, 1], col + "-", lw=1.2, label=lab)
        for r in pol.interiors:
            Q = np.asarray(r.coords)
            ax.plot(Q[:, 0], Q[:, 1], col + "-", lw=1.2)
    ax.plot(*np.array([info_f["splitter_tip"], info_f["island_tip"],
                       info_f["Q"], info_f["J"], info_f["outlet_entry"]]).T,
            "k.", ms=6)
    ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title("generated geometry [w_v]")
    ax = axes[1]
    rgb = np.zeros(ref.shape + (3,))
    rgb[..., 0] = ref            # 赤 = Fig.2
    rgb[..., 1] = gen            # 緑 = 生成
    ax.imshow(1 - rgb)
    ax.set_title(f"overlay on Fig.2   IoU = {ov['iou']:.4f}  "
                 "(red only = Fig.2, green only = generated)")
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "geometry", "gamboa_full.png")
    fig.savefig(png, dpi=130)
    print(f"\nsaved {os.path.normpath(png)}")

    out = dict(info=info_f, overlay=ov, area_faithful=poly_f.area,
               area_vertical=poly_v.area, cad=dict(wv_mm=wv, depth_mm=depth))
    j = os.path.join(HERE, "..", "results", "geometry", "gamboa_full_geom.json")
    json.dump(out, open(j, "w"), indent=1, default=str)
    print(f"saved {os.path.normpath(j)}")


if __name__ == "__main__":
    main()
