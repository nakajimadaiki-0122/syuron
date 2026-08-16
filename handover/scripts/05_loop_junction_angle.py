#!/usr/bin/env python3
"""
05_loop_junction_angle.py -- ループ流路中心線と主流路の接合角を DXF から直接測る
================================================================================

方法
----
1. 外形輪郭のうち、主流路上壁 y = +0.5 から離れて再び戻るまでの区間を
   「ループ外側輪郭」とする
2. 島（内側輪郭）から、y = +0.5 上に乗る底辺を除いた部分を「ループ内側輪郭」とする
3. 両輪郭を細かく標本化し、内側の各点に対する外側の最近点との中点を取る
   → これがループ流路の中心線
4. 中心線の両端（= ループが主流路へ開く側）で接線方向を求める
5. 主流路の +x 方向を 0 度として角度を報告する

DXF の輪郭は bulge を持たない厳密な区分線形なので、線分上の標本化は
補間ではなく厳密である。中心線が y = +0.5 に到達しない場合のみ、
終端の直線部を延長する（延長の可否は直線性の残差で判定し、報告する）。
"""
import json, os, sys
import numpy as np
import ezdxf
from shapely.geometry import Polygon

STEP = 0.001          # 輪郭標本化の刻み [mm]
MAIN_HALF = 0.5       # 主流路の半幅 [mm]


def load_rings(path):
    doc = ezdxf.readfile(path)
    rings = []
    for e in doc.modelspace():
        P = np.array([(q[0], q[1]) for q in e.get_points()])
        if np.allclose(P[0], P[-1]):
            P = P[:-1]
        bul = max(abs(q[4]) for q in e.get_points())
        rings.append(dict(pts=P, area=Polygon(P).area, max_bulge=bul))
    rings.sort(key=lambda r: -r["area"])
    return rings[0], rings[1:]


def densify(P, closed=False, step=STEP):
    """折れ線を等間隔で標本化する。区分線形なので厳密。"""
    Q = np.vstack([P, P[:1]]) if closed else P
    out = []
    for a, b in zip(Q[:-1], Q[1:]):
        L = np.hypot(*(b - a))
        n = max(int(np.ceil(L / step)), 1)
        t = np.arange(n)[:, None] / n
        out.append(a + t * (b - a))
    out.append(Q[-1:])
    return np.vstack(out)


def outer_loop_arc(outer_pts, x_lo, x_hi):
    """外形輪郭のうち、主流路上壁を離れてから戻るまでの区間を抜き出す。

    y > 0.5 の連続区間を列挙し、島の x 範囲を含むものを選ぶ。
    """
    above = outer_pts[:, 1] > MAIN_HALF + 1e-9
    runs, i = [], 0
    n = len(outer_pts)
    while i < n:
        if above[i]:
            j = i
            while j + 1 < n and above[j + 1]:
                j += 1
            runs.append((i, j))
            i = j + 1
        else:
            i += 1
    hit = [(i0, i1) for (i0, i1) in runs
           if outer_pts[i0:i1 + 1, 0].min() < x_lo
           and outer_pts[i0:i1 + 1, 0].max() > x_hi]
    if len(hit) != 1:
        raise ValueError(f"島 x=({x_lo:.3f},{x_hi:.3f}) を含む外側区間が "
                         f"{len(hit)} 個。区間: "
                         f"{[(round(outer_pts[a,0],2), round(outer_pts[b,0],2)) for a,b in runs]}")
    i0, i1 = hit[0][0] - 1, hit[0][1] + 1      # 両端の y=0.5 上の点を含める
    seg = outer_pts[i0:i1 + 1]
    if abs(seg[0, 1] - MAIN_HALF) > 1e-9 or abs(seg[-1, 1] - MAIN_HALF) > 1e-9:
        raise ValueError(f"外側輪郭の端点が y=0.5 に乗っていない: "
                         f"{seg[0]} .. {seg[-1]}")
    return seg


def island_flank(isl_pts):
    """島の輪郭から y = +0.5 上の底辺を除き、下流端 → 上流端の順に並べ替える。"""
    # 4 番目のループは作図誤差で下流端の角が y=0.5 から 4.5 um 浮いている。
    # ラスタ化（12 cells/mm = 83 um）では見えない量なので許容する。
    on = np.isclose(isl_pts[:, 1], MAIN_HALF, atol=0.01)
    if on.sum() != 2:
        raise ValueError(f"島が y=0.5 に接する点が {on.sum()} 個。2 個を想定")
    k = np.flatnonzero(on)
    P = np.vstack([isl_pts, isl_pts])           # 巡回させて底辺以外を取り出す
    a, b = k[0], k[1]
    seg1 = P[a:b + 1]                            # a -> b
    seg2 = P[b:a + len(isl_pts) + 1]             # b -> a(次周)
    # 底辺は 2 点だけの直線。長い方が側面
    flank = seg1 if len(seg1) > len(seg2) else seg2
    # 下流端（x 大）が先頭になるよう向きを揃える
    if flank[0, 0] < flank[-1, 0]:
        flank = flank[::-1]
    return flank


def centerline(inner, outer):
    """内側輪郭の各標本点に対し、外側輪郭上の最近点との中点を取る。"""
    I = densify(inner)
    O = densify(outer)
    d2 = ((I[:, None, 0] - O[None, :, 0]) ** 2
          + (I[:, None, 1] - O[None, :, 1]) ** 2)
    j = d2.argmin(1)
    Q = O[j]
    width = np.hypot(*(I - Q).T)
    return 0.5 * (I + Q), width, I, Q


def terminal_tangent(C, n_fit, label):
    """中心線の端点側 n_fit 点に直線を当て、接線方向と直線性の残差を返す。

    C[0] が端点。方向は「端点から内側へ向かう向き」の逆、
    すなわち **ループから主流路へ出ていく向き** に揃える。
    """
    S = C[:n_fit]
    c = S.mean(0)
    u, s, vt = np.linalg.svd(S - c)
    d = vt[0]
    resid = float(s[1] / np.sqrt(len(S)))        # 直線からの RMS ずれ [mm]
    # 端点 C[0] から外向き（= C[n_fit-1] から C[0] へ向かう向き）
    if np.dot(d, C[0] - S[-1]) < 0:
        d = -d
    ang = float(np.degrees(np.arctan2(d[1], d[0])))
    return dict(label=label, dir=d.tolist(), angle_from_plusx_deg=ang,
                straightness_rms_mm=resid, fit_points=int(n_fit),
                fit_length_mm=float(np.hypot(*(S[-1] - S[0]))),
                endpoint=C[0].tolist())


def extend_to_wall(endpoint, direction):
    """端点から接線方向へ延長して y = +0.5 と交わる点を返す。"""
    if abs(direction[1]) < 1e-12:
        return None
    t = (MAIN_HALF - endpoint[1]) / direction[1]
    return (np.array(endpoint) + t * np.array(direction)).tolist(), float(t)


def main():
    dxf = sys.argv[1] if len(sys.argv) > 1 else "data/tesla_channel.dxf"
    outer, islands = load_rings(dxf)
    print(f"DXF: {dxf}")
    print(f"  外形: {len(outer['pts'])} 点, max|bulge| = {outer['max_bulge']:.3e}")
    for i, r in enumerate(islands):
        b = Polygon(r["pts"]).bounds
        print(f"  島 {i}: {len(r['pts'])} 点, area = {r['area']:.4f}, "
              f"x = {b[0]:.3f}..{b[2]:.3f}, max|bulge| = {r['max_bulge']:.3e}")
    if max(r["max_bulge"] for r in [outer] + islands) > 1e-12:
        print("  警告: bulge が非零。円弧を折れ線として扱っている")

    islands = sorted(islands, key=lambda r: Polygon(r["pts"]).bounds[0])
    results = []
    for li, isl in enumerate(islands):
        b = Polygon(isl["pts"]).bounds
        arc = outer_loop_arc(outer["pts"], b[0], b[2])
        flank = island_flank(isl["pts"])
        C, width, I, Q = centerline(flank, arc)

        print(f"\n===== ループ {li + 1}  (島 x = {b[0]:.3f} .. {b[2]:.3f}) =====")
        print(f"  外側輪郭: x = {arc[0,0]:.5f} .. {arc[-1,0]:.5f} で y=0.5 を離脱/復帰")
        print(f"  島の接地点: x = {flank[0,0]:.5f}（下流端）, "
              f"{flank[-1,0]:.5f}（上流端）")
        core = width[(C[:, 1] > 1.2)]
        print(f"  流路幅（中心線上、y>1.2 の区間）: "
              f"{core.min():.4f} .. {core.max():.4f} mm  平均 {core.mean():.4f}")
        print(f"  中心線 端点: 下流側 ({C[0,0]:.4f}, {C[0,1]:.4f}), "
              f"上流側 ({C[-1,0]:.4f}, {C[-1,1]:.4f})")

        # 中心線が有効なのは、内外輪郭が実際に向かい合っている区間だけ。
        # 島の底の角では最近点が主流路上壁に飛ぶため幅が跳ね、中心線が退化する。
        w0 = float(np.median(width))
        valid = np.abs(width - w0) < 0.05 * w0
        k = np.flatnonzero(valid)
        if k.size == 0:
            raise ValueError("有効な中心線区間がない")
        k0, k1 = k.min(), k.max()
        print(f"  有効区間: 標本 {k0}..{k1} / {len(C)}  "
              f"（幅 {w0:.4f} mm の ±5 % 以内）")
        print(f"  退化して切り捨てた長さ: 下流側 "
              f"{np.hypot(*(C[k0]-C[0])):.4f} mm, 上流側 "
              f"{np.hypot(*(C[k1]-C[-1])):.4f} mm")
        Cv = C[k0:k1 + 1]
        print(f"  有効中心線 端点: 下流側 ({Cv[0,0]:.4f}, {Cv[0,1]:.4f}), "
              f"上流側 ({Cv[-1,0]:.4f}, {Cv[-1,1]:.4f})")

        r = dict(loop=li + 1, island_x=[b[0], b[2]],
                 outer_leave_x=float(arc[0, 0]), outer_return_x=float(arc[-1, 0]),
                 width_mean_mm=float(core.mean()),
                 width_min_mm=float(core.min()), width_max_mm=float(core.max()),
                 valid_span=[int(k0), int(k1)], n_samples=int(len(C)))
        for side, Cs in (("downstream", Cv), ("upstream", Cv[::-1])):
            print(f"\n  --- {side} 側接合 ---")
            rows = []
            for n_fit in [100, 200, 400, 800]:
                if n_fit > len(Cs):
                    continue
                t = terminal_tangent(Cs, n_fit, side)
                ext = extend_to_wall(t["endpoint"], t["dir"])
                t["cross_y05"] = None if ext is None else ext[0]
                t["extend_len_mm"] = None if ext is None else abs(ext[1])
                rows.append(t)
                cx = "--" if ext is None else f"{ext[0][0]:.4f}"
                el = "--" if ext is None else f"{abs(ext[1]):.4f}"
                print(f"    fit {n_fit:4d} 点 ({t['fit_length_mm']:.3f} mm): "
                      f"角度 {t['angle_from_plusx_deg']:+8.3f} deg, "
                      f"直線残差 {t['straightness_rms_mm']:.2e} mm, "
                      f"y=0.5 交点 x = {cx}, 延長 {el} mm")
            r[side] = rows
        results.append(r)
        if li == 0:
            np.savez('results/loop1_centerline.npz', C=C, width=width,
                     I=I, Q=Q, arc=arc, flank=flank, valid=valid)

    os.makedirs("results", exist_ok=True)
    json.dump(results, open("results/loop_junction_angles.json", "w"),
              indent=1, default=str)
    print("\nsaved results/loop_junction_angles.json")


if __name__ == "__main__":
    main()
