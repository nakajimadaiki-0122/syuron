#!/usr/bin/env python3
"""
13_fig2_topology.py -- Fig.2（Gamboa optimized の実形状）から接合部の位相を実測する
====================================================================================

なぜ必要か
----------
`src/gamboa.py` は周期化のために Gamboa の出口区間を「閉じ壁」に置き換えている。
その結果、生成形状では

    主流路がループ下流開口の先へまっすぐ続き、
    ループ下流枝は主流路へ 48 度で降りてくる側枝になる

一方 Fig.1 / Fig.2 では、実形状は

    主流路は接合部で終わり、流れは alpha だけ折れて出口区間へ入る。
    ループ下流枝は出口区間と同一直線（折れ角ゼロ）である

これが正しければ、逆流時に「出口区間からまっすぐループへ入る」という
Tesla バルブの中核機構が周期化で失われている。目視ではなく実測して確定する。

方法
----
docs/gamboa2005_fig2.png（PDF 埋め込み画像 997x714、再標本化なし）を連結成分に分け、
optimum の外形（実線なので単一成分）を塗って島を抜き、流体領域を得る。
reference は破線なので多数の小片に分かれ、自動的に除かれる。

斜め帯（ループ下流枝 + 出口区間）は水平線で切ると 1 本の区間になるので、
**行ごとの区間の左端・右端・中点**を集めて直線を当てる。接合部の上下で別々に当て、
角度差と法線オフセットを比較する。

出力
  1. 主流路の下壁が水平を保つ x 範囲（= 主流路の終端）
  2. 斜め帯の壁の直線当てはめ（接合部の上 = ループ下流枝、下 = 出口区間）
  3. 同一直線性の判定
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
OUT = os.path.join(HERE, "..", "results", "geometry", "fig2_topology.json")

# 2026-08-09 の実測値（docs/gamboa2005_geometry.md 5.2）
W_PX_REF = 44.50
YC_REF = 271.25
X0_REF = 246.25

# 帯を切り出す行の範囲（画像座標、y は下向き）
ROWS_LOOP = (95, 235)      # 接合部より上: ループ下流枝
ROWS_OUT = (300, 370)      # 接合部より下: 出口区間（右プレナムに入る手前まで）


def fluid_mask():
    a = mpimg.imread(IMG)
    a = a if a.ndim == 2 else a[..., 0]
    lab, n = ndi.label(a < 0.5, structure=np.ones((3, 3)))
    sz = ndi.sum(a < 0.5, lab, range(1, n + 1))
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
    return outer & ~island, island, oid, iid


def row_segments(fluid, y):
    xs = np.where(fluid[y])[0]
    if len(xs) == 0:
        return []
    return np.split(xs, np.where(np.diff(xs) > 1)[0] + 1)


def band_points(fluid, rows, x_ref):
    """各行で x_ref に最も近い区間を斜め帯とみなし、左端・右端・中点を返す。"""
    L, R, M = [], [], []
    for y in range(*rows):
        segs = row_segments(fluid, y)
        if not segs:
            continue
        seg = min(segs, key=lambda s: abs(0.5 * (s[0] + s[-1]) - x_ref(y)))
        if not (20 <= len(seg) <= 90):   # 破線の切れ端・プレナムとの合体を除く
            continue
        L.append((seg[0], y)); R.append((seg[-1], y))
        M.append((0.5 * (seg[0] + seg[-1]), y))
    return (np.array(L, float), np.array(R, float), np.array(M, float))


def fit_line(pts):
    """点群に直線を当てる。戻り: 方向 d（画像座標）、通過点 c、残差 RMS。"""
    c = pts.mean(0)
    _, s, vt = np.linalg.svd(pts - c)
    return vt[0], c, float(s[1] / np.sqrt(len(pts)))


def ang_from_x(d):
    """画像座標（y 下向き）の方向 -> 主流路軸からの角度 [deg]、-90..90。"""
    a = np.degrees(np.arctan2(-d[1], d[0]))
    return (a + 90) % 180 - 90


def normal_offset(a, b):
    """直線 a（角度・通過点）に対する直線 b の通過点の法線距離 [px]。"""
    th = np.radians(a["angle_deg"])
    d = np.array([np.cos(th), -np.sin(th)])
    nrm = np.array([-d[1], d[0]])
    return float(np.dot(np.array(b["c"]) - np.array(a["c"]), nrm))


def main():
    fluid, island, oid, iid = fluid_mask()
    H, W = fluid.shape
    res = {"image": dict(w=W, h=H, outer_component=oid, island_component=iid,
                         fluid_px=int(fluid.sum()))}
    print(f"画像 {W} x {H} px   外形成分 #{oid}   島成分 #{iid}   "
          f"流体 {int(fluid.sum())} px")

    # --- 0. 基準量（主流路の幅と中心線）: 戻り流路の開口より上流で測る ---
    tops, bots = [], []
    for x in range(250, 305):        # 分流板（x=224..321）の下、先端の手前
        segs = row_segments(fluid.T, x)      # 転置して列走査
        for s in segs:
            if s[0] <= YC_REF <= s[-1]:
                tops.append(s[0]); bots.append(s[-1])
    tops = np.array(tops); bots = np.array(bots)
    # 分流板が途切れる列（ループと繋がって帯が伸びる）を中央値で除く
    wid = bots - tops
    keep = np.abs(wid - np.median(wid)) <= 2
    w_px = float(np.mean(wid[keep]) + 1)
    yc = float(np.mean(bots[keep] + tops[keep]) / 2)
    print(f"主流路: 幅 {w_px:.2f} px（08-09 実測 {W_PX_REF}）, "
          f"中心線 y = {yc:.2f} px（同 {YC_REF}）")
    res["scale"] = dict(w_px=w_px, yc=yc, x0=X0_REF)

    def to_wv(x, y=None):
        return ((x - X0_REF) / w_px if y is None
                else ((x - X0_REF) / w_px, -(y - yc) / w_px))

    # --- 1. 主流路はどこで終わるか ---
    xs_scan, bot_y = [], []
    for x in range(340, W):
        segs = row_segments(fluid.T, x)
        seg = next((s for s in segs if s[0] - 1 <= yc <= s[-1] + 1), None)
        if seg is None:
            break
        xs_scan.append(x); bot_y.append(seg[-1])
    xs_scan = np.array(xs_scan); bot_y = np.array(bot_y)
    flat = np.abs(bot_y - bot_y[0]) <= 1
    x_flat_end = int(xs_scan[flat][-1])
    print(f"\n1. 主流路の下壁 y = {bot_y[0]} px が水平な範囲: "
          f"x = {xs_scan[0]}..{x_flat_end} px "
          f"= {to_wv(xs_scan[0]):.2f}..{to_wv(x_flat_end):.2f} w_v")
    print(f"   中心線を含む帯が消える x = {xs_scan[-1]} px "
          f"= {to_wv(xs_scan[-1]):.2f} w_v"
          "   -> 主流路はここで終わる（まっすぐ先へ続かない）")
    res["main_channel"] = dict(bottom_wall_y=int(bot_y[0]),
                               flat_x=[int(xs_scan[0]), x_flat_end],
                               flat_x_wv=[to_wv(xs_scan[0]), to_wv(x_flat_end)],
                               last_x=int(xs_scan[-1]),
                               last_x_wv=to_wv(xs_scan[-1]))

    # --- 2. 斜め帯（ループ下流枝 / 出口区間）の直線当てはめ ---
    ys_i, xs_i = np.where(island)
    tip_x = int(xs_i.max()); tip_y = int(ys_i[xs_i == tip_x].mean())
    print(f"\n2. 島の下流端（尖点） = ({tip_x}, {tip_y}) px "
          f"= ({to_wv(tip_x, tip_y)[0]:.3f}, {to_wv(tip_x, tip_y)[1]:.3f}) w_v")

    # 帯の概略位置（行 y での中心 x）。行 100 と 340 の実測から線形に見積もる。
    def x_ref(y):
        return 413.0 + (626.5 - 413.0) * (y - 100) / (340 - 100)

    parts = {}
    for name, rows in (("loop_branch", ROWS_LOOP), ("outlet_segment", ROWS_OUT)):
        L, R, M = band_points(fluid, rows, x_ref)
        entry = {}
        for side, P in (("left", L), ("right", R), ("mid", M)):
            d, c, rms = fit_line(P)
            entry[side] = dict(n=int(len(P)), angle_deg=float(ang_from_x(d)),
                               rms_px=rms, c=[float(c[0]), float(c[1])])
        # 帯の幅（法線方向）: 中心線の角度で左右端の距離を射影
        th = np.radians(entry["mid"]["angle_deg"])
        width = float(np.mean(R[:, 0] - L[:, 0]) * abs(np.sin(th)))
        entry["width_px"] = width
        entry["width_wv"] = width / w_px
        parts[name] = entry
        print(f"   {name:15s} rows {rows[0]}..{rows[1]}  "
              f"n={entry['mid']['n']:3d}  "
              f"角度 = {entry['mid']['angle_deg']:+7.3f} deg  "
              f"残差 {entry['mid']['rms_px']:.2f} px  "
              f"幅 = {width:.1f} px = {width/w_px:.3f} w_v")
    res["bands"] = parts

    # --- 3. 同一直線性 ---
    print("\n3. 接合部の上下で同一直線か（ループ下流枝 対 出口区間）")
    coll = {}
    for side in ("left", "right", "mid"):
        a = parts["loop_branch"][side]
        b = parts["outlet_segment"][side]
        dang = b["angle_deg"] - a["angle_deg"]
        off = normal_offset(a, b)
        coll[side] = dict(dangle_deg=float(dang), offset_px=float(off),
                          offset_wv=float(off / w_px))
        print(f"   {side:5s}壁: 角度差 = {dang:+6.3f} deg,  "
              f"法線オフセット = {off:+6.2f} px = {off/w_px:+.3f} w_v")
    ok = (abs(coll["mid"]["dangle_deg"]) < 2.0
          and abs(coll["mid"]["offset_wv"]) < 0.15)
    coll["collinear"] = bool(ok)
    res["collinearity"] = coll
    print("   -> 出口区間とループ下流枝は同一直線（接合部に折れ角がない）" if ok
          else "   -> 同一直線ではない")

    # alpha の基準（主流路軸から測るか法線から測るか）
    a_axis = abs(parts["outlet_segment"]["mid"]["angle_deg"])
    print(f"\n4. 出口区間の傾き: 主流路軸から {a_axis:.2f} deg, "
          f"法線から {90 - a_axis:.2f} deg   (Table 2 の alpha = 41.9)")
    res["alpha_check"] = dict(from_axis_deg=a_axis, from_normal_deg=90 - a_axis,
                              table2_alpha=41.9)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"\nsaved {os.path.normpath(OUT)}")

    # --- 図 ---
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.imshow(fluid, cmap="gray_r")
    for name, col in (("loop_branch", "C3"), ("outlet_segment", "C0")):
        for side in ("left", "right"):
            a = parts[name][side]
            th = np.radians(a["angle_deg"])
            d = np.array([np.cos(th), -np.sin(th)])
            t = np.linspace(-330, 330, 2)
            L = np.array(a["c"]) + t[:, None] * d
            ax.plot(L[:, 0], L[:, 1], col + "-", lw=1.0)
    ax.plot([tip_x], [tip_y], "C2o", ms=5)
    ax.axvline(x_flat_end, color="C2", lw=0.8, ls="--")
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_title("Fig.2 optimized: red = loop branch fit, blue = outlet fit")
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "geometry", "fig2_topology.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
