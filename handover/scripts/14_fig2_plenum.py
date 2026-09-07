#!/usr/bin/env python3
"""
14_fig2_plenum.py -- Fig.2 からプレナムの円と接続位置を実測する
================================================================

案 1（プレナム込みの単発モデル）を組むために必要な、まだ未確定だった量を測る。

  1. 左右プレナムの円の中心と半径（論文記載は R_p = 5 w_v）
     -> 半径が 5 w_v になるかは **スケール（w_v [px]）の独立検証**になる
  2. 出口区間が右プレナムの円弧に入る点（= LENOUT の終点）
  3. LENOUT の起点をどこに取れば Table 2 の 2.94 w_v になるか

方法
----
流体マスクの行走査。左プレナムは各行の最左区間の右端が円弧、
右プレナムは各行の最右区間の左端が円弧に載る。直線壁の点は入らない。
円は最小二乗（Kasa 法）で当て、残差 3 px 超を 1 回除いて再当てはめする。

w_v の基準は **線の中心**で測る（塗り潰し領域の端で測ると線幅ぶん太くなり、
2026-08-09 の実測 44.50 px に対して 8 % ずれる）。
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "..", "docs", "gamboa2005_fig2.png")
OUT = os.path.join(HERE, "..", "results", "geometry", "fig2_plenum.json")


def masks():
    a = mpimg.imread(IMG)
    a = a if a.ndim == 2 else a[..., 0]
    dark = a < 0.5
    lab, n = ndi.label(dark, structure=np.ones((3, 3)))
    sz = ndi.sum(dark, lab, range(1, n + 1))
    oid = int(np.argmax(sz)) + 1
    outer = ndi.binary_fill_holes(lab == oid)
    iid = None
    for s, k in sorted(((sz[k - 1], k) for k in range(1, n + 1) if k != oid),
                       reverse=True):
        ys, xs = np.where(lab == k)
        if outer[ys, xs].all():
            iid = k
            break
    island = ndi.binary_fill_holes(lab == iid)
    return outer & ~island, (lab == oid), island


def line_centre_scale(line, ylo=200, yhi=340, xlo=250, xhi=305):
    """分流板の下面線と主流路下壁の**線の中心**間隔で w_v を測る。"""
    ws, ycs = [], []
    for x in range(xlo, xhi):
        ys = np.where(line[:, x])[0]
        ys = ys[(ys > ylo) & (ys < yhi)]
        if len(ys) < 2:
            continue
        runs = np.split(ys, np.where(np.diff(ys) > 1)[0] + 1)
        if len(runs) < 2:
            continue
        # 分流板は厚みがあり上下 2 本の線になるので、下から 2 本
        # （分流板の下面線と主流路下壁）を使う
        c0 = runs[-2].mean(); c1 = runs[-1].mean()
        ws.append(c1 - c0); ycs.append(0.5 * (c0 + c1))
    ws = np.array(ws); ycs = np.array(ycs)
    keep = np.abs(ws - np.median(ws)) <= 1.0
    return float(ws[keep].mean()), float(ycs[keep].mean()), int(keep.sum())


def fit_circle(P):
    """Kasa 法。戻り: (cx, cy), r, 残差 RMS。"""
    x, y = P[:, 0], P[:, 1]
    A = np.c_[2 * x, 2 * y, np.ones(len(x))]
    b = x ** 2 + y ** 2
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy = sol[0], sol[1]
    r = np.sqrt(sol[2] + cx ** 2 + cy ** 2)
    res = np.hypot(x - cx, y - cy) - r
    return (float(cx), float(cy)), float(r), float(np.sqrt((res ** 2).mean())), res


def fit_circle_robust(P, tol=3.0):
    c, r, rms, res = fit_circle(P)
    keep = np.abs(res) <= max(tol, 2 * rms)
    c, r, rms, _ = fit_circle(P[keep])
    return c, r, rms, int(keep.sum()), int(len(P))


def main():
    fluid, outer_line, island = masks()
    H, W = fluid.shape
    w_px, yc, nfit = line_centre_scale(outer_line)
    x0 = 224 + 0.5 * w_px          # 入口流路の左端 x=224 px から 0.5 w_v
    print(f"w_v = {w_px:.2f} px（線の中心間、n={nfit}）   "
          f"主流路中心線 y = {yc:.2f} px   原点 x = {x0:.2f} px")
    print("（2026-08-09 の実測: w_v = 44.50, yc = 271.25, x0 = 246.25）")

    def wv(x, y):
        return ((x - x0) / w_px, -(y - yc) / w_px)

    # --- 左プレナム: 各行の最左区間の右端が円弧 ---
    PL = []
    for y in range(5, H - 5):
        xs = np.where(fluid[y])[0]
        if len(xs) == 0:
            continue
        segs = np.split(xs, np.where(np.diff(xs) > 1)[0] + 1)
        s = segs[0]
        if s[0] > 40:                    # 左端に接していない行は左プレナムでない
            continue
        if s[-1] > 300:                  # 主流路と繋がる行は除く
            continue
        PL.append((s[-1], y))
    PL = np.array(PL, float)
    cL, rL, rmsL, kL, nL = fit_circle_robust(PL)
    print(f"\n左プレナム: 中心 ({cL[0]:.1f}, {cL[1]:.1f}) px, 半径 {rL:.1f} px "
          f"= {rL/w_px:.3f} w_v   残差 {rmsL:.2f} px  ({kL}/{nL} 点)")
    print(f"   中心 = ({wv(*cL)[0]:.3f}, {wv(*cL)[1]:.3f}) w_v   "
          f"（論文 R_p = 5、幾何から予想される中心 x = -5.475, y = 0）")

    # --- 右プレナム: 下半分の各行の最右区間の左端が円弧 ---
    PR = []
    for y in range(430, H - 5):
        xs = np.where(fluid[y])[0]
        if len(xs) == 0:
            continue
        segs = np.split(xs, np.where(np.diff(xs) > 1)[0] + 1)
        s = segs[-1]
        if s[-1] - s[0] < 30:
            continue
        PR.append((s[0], y))
    PR = np.array(PR, float)
    cR, rR, rmsR, kR, nR = fit_circle_robust(PR)
    print(f"\n右プレナム: 中心 ({cR[0]:.1f}, {cR[1]:.1f}) px, 半径 {rR:.1f} px "
          f"= {rR/w_px:.3f} w_v   残差 {rmsR:.2f} px  ({kR}/{nR} 点)")
    print(f"   中心 = ({wv(*cR)[0]:.3f}, {wv(*cR)[1]:.3f}) w_v")

    # --- 出口区間の中心線と右プレナム円の交点 = LENOUT の終点 ---
    # 出口区間の傾き（13_fig2_topology.py の実測）
    alpha_axis = np.radians(48.03)
    e = np.array([np.cos(alpha_axis), np.sin(alpha_axis)])      # 画像座標で下向き
    # 中心線は右プレナムの中心を通ると仮定できるか確認する
    cRv = np.array(wv(*cR))
    print(f"\n出口区間の軸（w_v 系、下向き -48.03 度）と右プレナム中心の関係")
    # 出口中心線: 点 J0 を通り方向 (cos(-48.03), sin(-48.03))
    ev = np.array([np.cos(-alpha_axis), np.sin(-alpha_axis)])
    # 13 の実測から出口区間の左右壁の中点直線（w_v 系）を作る
    top = json.load(open(os.path.join(HERE, "..", "results", "geometry", "fig2_topology.json")))
    cmid = np.array(top["bands"]["outlet_segment"]["mid"]["c"])   # px
    cmidv = np.array(wv(*cmid))
    t = np.dot(cRv - cmidv, ev)
    foot = cmidv + t * ev
    print(f"   出口中心線から右プレナム中心までの法線距離 = "
          f"{np.linalg.norm(cRv - foot):.3f} w_v "
          f"-> {'中心線上にある' if np.linalg.norm(cRv-foot) < 0.3 else '中心線から外れる'}")
    entry = cRv - rR / w_px * ev          # 円弧との交点（上流側）
    print(f"   円弧との交点（= LENOUT の終点） = "
          f"({entry[0]:.3f}, {entry[1]:.3f}) w_v")

    res = dict(scale=dict(w_px=w_px, yc=yc, x0=x0),
               left_plenum=dict(centre_px=cL, r_px=rL, r_wv=rL / w_px,
                                centre_wv=list(wv(*cL)), rms_px=rmsL,
                                n_used=kL, n=nL),
               right_plenum=dict(centre_px=cR, r_px=rR, r_wv=rR / w_px,
                                 centre_wv=list(wv(*cR)), rms_px=rmsR,
                                 n_used=kR, n=nR),
               outlet_entry_wv=list(map(float, entry)))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"\nsaved {os.path.normpath(OUT)}")

    fig, ax = plt.subplots(figsize=(11, 8))
    ax.imshow(fluid, cmap="gray_r")
    th = np.linspace(0, 2 * np.pi, 400)
    for c, r, col in ((cL, rL, "C3"), (cR, rR, "C0")):
        ax.plot(c[0] + r * np.cos(th), c[1] + r * np.sin(th), col + "-", lw=1.0)
        ax.plot([c[0]], [c[1]], col + "o", ms=4)
    ax.plot(PL[:, 0], PL[:, 1], "C3.", ms=1)
    ax.plot(PR[:, 0], PR[:, 1], "C0.", ms=1)
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_title("Fig.2: plenum circle fits")
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "geometry", "fig2_plenum.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
