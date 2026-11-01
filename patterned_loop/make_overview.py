"""パターン 2、3 の全体像を描く。生成した fluid.stl をそのまま描画する（水の形）。

図の水の形では、出っぱり（パターン 2）は水のへこみ、へこみ（パターン 3）は水のふくらみに見える。
半球の面だけ色を変える。Fluent は使わない。
"""
import importlib.util

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.ops import unary_union

import make_hemisphere as H

Z, X, P, HERE = H.Z, H.X, H.P, H.HERE
spec = importlib.util.spec_from_file_location("preview", HERE / "pattern2_hemisphere" / "make_preview.py")
PV = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PV)

WATER, DOME = np.array([.45, .75, .80]), np.array([.90, .42, .15])


def load(path):
    """Y-up の STL を描画用の座標 (x, 平面の y, 高さ) へ戻す。"""
    t = Z.read_stl(path)
    return np.stack([t[:, :, 0], -t[:, :, 2], t[:, :, 1]], axis=2)


def groups(t, xlim=None):
    if xlim:
        c = t[:, :, 0].mean(axis=1)
        t = t[(c > xlim[0]) & (c < xlim[1])]
    n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)
    dome = (np.abs(n[:, 2]) > 0.12) & (np.abs(n[:, 2]) < 0.999)
    return [(t[~dome], WATER), (t[dome], DOME)]


def main():
    base, plan, panels, _ = H.build_plan()
    for folder, mode, name in (("pattern2_hemisphere", "bump", "パターン 2: 半球の出っぱり"),
                               ("pattern3_dimple", "dimple", "パターン 3: 半球のへこみ")):
        out = HERE / folder
        fwd = load(out / "symtvm_fwd" / "fluid.stl")
        rev = load(out / "symtvm_rev" / "fluid.stl")
        view = (0.18, 0.72, -0.67)
        a = PV.render(groups(fwd), view, w=2600, h=470, pad=0.01)
        b = PV.render(groups(rev), view, w=2600, h=470, pad=0.01)
        c = PV.render(groups(fwd, (5.0, 17.5)), view, w=2600, h=760, pad=0.02)

        discs = []
        for pn in panels:
            for sc, zc in pn.centres:
                if np.isclose(zc, H.ROWS[1]):
                    p = pn.frame(sc)[0]
                    discs += [Point(p[0], p[1]).buffer(H.R_HEMI, 32), Point(p[0], -p[1]).buffer(H.R_HEMI, 32)]
        discs = unary_union(discs)
        section = base.difference(discs) if mode == "bump" else base.union(discs)

        P.plotstyle.use_jp()
        fig = plt.figure(figsize=(16, 12.5))
        gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1.55, 0.72])
        for i, (img, title) in enumerate(((a, "(a) Forward の全体（長さ 50 mm）。入口は左"),
                                          (b, "(b) Reverse の全体（Forward を 180 度回転）。入口は左"),
                                          (c, "(c) Forward の拡大（x = 5〜17.5 mm）。オレンジが半球の面"))):
            ax = fig.add_subplot(gs[i])
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(title, loc="left", fontsize=12)
        ax = fig.add_subplot(gs[3])
        geoms = section.geoms if section.geom_type == "MultiPolygon" else [section]
        for g in geoms:
            Q = np.asarray(g.exterior.coords)
            ax.fill(Q[:, 0], Q[:, 1], color="#73bfcc", lw=0)
            for r in g.interiors:
                Q = np.asarray(r.coords)
                ax.fill(Q[:, 0], Q[:, 1], color="white", lw=0)
        ax.add_patch(plt.Rectangle((0, -0.5 * X.W_UNIT), X.LY, X.W_UNIT, fill=False, lw=0.8, ls=":"))
        ax.set_xlim(-0.3, X.LY + 0.3)
        ax.set_ylim(-0.5 * X.W_UNIT - 0.1, 0.5 * X.W_UNIT + 0.1)
        ax.set_aspect("equal")
        ax.set_xlabel("流れ方向 x [mm]")
        ax.set_title("(d) Forward の平面図（流路の深さの中央の断面）。点線は計算単位 50 × 3.54 mm", loc="left", fontsize=12)
        fig.suptitle(f"{name} の全体像（水の形）  半径 {H.R_HEMI} mm、中心間隔 {H.SPACING} mm、深さ方向 {len(H.ROWS)} 段", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(out / "overview_cad.png", dpi=130)
        plt.close(fig)
        print("saved", out / "overview_cad.png")


if __name__ == "__main__":
    main()
