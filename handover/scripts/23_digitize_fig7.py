#!/usr/bin/env python3
"""
23_digitize_fig7.py -- Gamboa Fig.7 の「opt CFD」曲線を数値化する
=================================================================

なぜやり直すか
--------------
Re = 100 / 300 / 500 の目標値（1.02 / 1.12 / 1.37）は 2026-08-09 に目視で読んだもので、
スクリプトが残っていない。**再現の合否がこの読み取りに直接かかっている**ので、
手順を残せる形で測り直し、不確かさも出す。

方法
----
PDF 埋め込み画像（997 x 713、再標本化なし）を使う。Fig.7 には 4 系列がある。

  opt CFD  破線     <- これを取る
  opt EXP  黒丸 + 点線のフィット
  ref CFD  実線（1 つの巨大な連結成分）
  ref EXP  x 印 + 点線のフィット

2 段構えで取る。

  1. **破線の一片を連結成分として拾う**（面積 30〜120 px）。実線は 1 個の巨大成分、
     丸印は約 500 px、点線は 20 px 未満なので分離できる。
     曲線が混み合わない Re >= 350 ではこれできれいに取れる
  2. Re < 350 は他の曲線と融合するので、**実線（1 個の巨大成分）を差し引いてから**
     列ごとに左へ追跡する。傾きは直近 60 px の点列から見積もる

**Re < 約 290 では破線と実線が重なって分離できない。**（実測: Re=300 では
破線 1.161 と実線 1.124 が別々の区間として見えるが、Re=250 以下では 1 本になる）
その領域の値は「実線の位置 ± 線幅」としか言えないので None を返す。
2026-08-09 の目視値（Re=100 で 1.02、200 で 1.07）は実線の位置とほぼ同じで、
**低 Re では両曲線が区別できていなかった**ことを意味する。

軸の較正は枠線から（左下が Re=0, Di=1、右上が Re=2000, Di=2）。
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import ndimage as ndi
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "..", "docs", "gamboa2005_fig7.png")
OUT = os.path.join(HERE, "..", "results", "fig7_digitized.json")

RE_RANGE = (0.0, 2000.0)
DI_RANGE = (1.0, 2.0)
DASH_AREA = (30, 120)          # 破線の一片とみなす面積 [px]


def main():
    a = mpimg.imread(IMG)
    a = a if a.ndim == 2 else a[..., 0]
    dark = a < 0.5
    H, W = dark.shape
    vc = np.flatnonzero(dark.sum(0) > 0.8 * H)
    hr = np.flatnonzero(dark.sum(1) > 0.8 * W)
    x0, x1, y0, y1 = int(vc.min()), int(vc.max()), int(hr.min()), int(hr.max())
    print(f"画像 {W} x {H} px   枠 x {x0}..{x1}, y {y0}..{y1}")

    inner = np.zeros_like(dark)
    inner[y0 + 2:y1 - 1, x0 + 2:x1 - 1] = dark[y0 + 2:y1 - 1, x0 + 2:x1 - 1]
    # 凡例（枠の左上）を外す。破線がここへ入るのは Di > 1.7 かつ Re < 940 だが、
    # 実際に破線が 1.7 に届くのは Re = 950 付近なので曲線は落ちない。
    inner[:int(y0 + 0.30 * (y1 - y0)), :int(x0 + 0.47 * (x1 - x0))] = False

    def to_re(px):
        return (px - x0) / (x1 - x0) * (RE_RANGE[1] - RE_RANGE[0]) + RE_RANGE[0]

    def to_di(py):
        return DI_RANGE[1] - (py - y0) / (y1 - y0) * (DI_RANGE[1] - DI_RANGE[0])

    # --- 1. 破線の一片を連結成分として拾う ---
    lab, n = ndi.label(inner, structure=np.ones((3, 3)))
    objs = ndi.find_objects(lab)
    sz = ndi.sum(inner, lab, range(1, n + 1))
    seg = []
    for k, sl in enumerate(objs):
        if not (DASH_AREA[0] <= sz[k] <= DASH_AREA[1]):
            continue
        ys, xs = np.where(lab[sl] == k + 1)
        seg.append((float(sl[1].start + xs.mean()),
                    float(sl[0].start + ys.mean())))
    seg.sort()
    # 右端に点線の切れ端が混ざる。画像では上へ単調に進むはずなので外れ値を落とす
    clean = []
    for cx, cy in seg:
        if clean and cy > clean[-1][1] + 3:
            continue
        clean.append((cx, cy))
    seg = clean
    print(f"破線の一片 {len(seg)} 個   "
          f"Re {to_re(seg[0][0]):.0f}..{to_re(seg[-1][0]):.0f}")

    # --- 2. 実線（ref CFD）を分離してから低 Re 側を追跡する ---
    big = int(np.argmax(sz)) + 1
    solid = (lab == big)
    rest = inner & ~solid
    print(f"実線（ref CFD）= 成分 #{big}、面積 {sz[big-1]:.0f} px")

    def runs(col, arr):
        ys = np.flatnonzero(arr[:, col])
        if len(ys) == 0:
            return []
        segs = np.split(ys, np.flatnonzero(np.diff(ys) > 1) + 1)
        return [float(np.mean(s)) for s in segs if len(s) <= 14]

    pts = list(seg)
    miss = 0
    for c in range(int(seg[0][0]) - 1, x0 + 2, -1):
        tail = [q for q in pts if 0 <= q[0] - c <= 60][:40]
        if len(tail) >= 5 and (tail[-1][0] - tail[0][0]) >= 20:
            A = np.polyfit([q[0] for q in tail], [q[1] for q in tail], 1)
            slope = float(np.clip(A[0], -3.0, 3.0))
        else:
            slope = 0.0
        last = min(pts, key=lambda q: q[0])
        pred = last[1] + slope * (c - last[0])
        gate = 5.0 + 1.5 * abs(slope) * max(1, last[0] - c)
        cand = [(abs(y - pred), y) for y in runs(c, rest)
                if abs(y - pred) <= gate]
        if cand:
            pts.append((float(c), min(cand)[1]))
            miss = 0
        else:
            miss += 1
            if miss > 25:
                break

    # 実線の位置も記録する（低 Re でどれだけ近いかを示すため）
    solid_pts = []
    for c in range(x0 + 3, x1 - 3):
        r = runs(c, solid)
        if r:
            solid_pts.append((float(c), float(np.mean(r))))

    pts.sort()
    cx = np.array([q[0] for q in pts])
    cy = np.array([q[1] for q in pts])
    re_pts, di_pts = to_re(cx), to_di(cy)
    print(f"追跡後の点 {len(pts)} 個   Re {re_pts.min():.0f}..{re_pts.max():.0f}")

    # --- 3. 低 Re 側の補足測定 ---
    # 追跡は破線の切れ目と交差で不安定なので、目標 Re の列を直接見て、
    # **実線を除いた区間**のうち外挿値にいちばん近いものを取る。
    # 別々の区間として見えている間だけ採用する（重なったら None）。
    # 外挿は左端 2 点の傾きによる直線で行う（PCHIP の外挿は暴れる）
    _p = sorted(pts)[:6]
    _r = to_re(np.array([q[0] for q in _p]))
    _d = to_di(np.array([q[1] for q in _p]))
    _a = np.polyfit(_r, _d, 1)

    def f_ext(t):
        return _a[0] * t + _a[1]
    extra = {}
    for t in (250, 275, 300, 320):
        c = int(round(x0 + t / RE_RANGE[1] * (x1 - x0)))
        pred = float(f_ext(t))
        cands = [(abs(to_di(y) - pred), to_di(y)) for y in runs(c, rest)]
        cands = [c_ for c_ in cands if c_[0] < 0.03]
        if cands:
            v = min(cands)[1]
            extra[t] = float(v)
            pts.append((float(c), float(y0 + (DI_RANGE[1] - v)
                                        * (y1 - y0) / (DI_RANGE[1] - DI_RANGE[0]))))
    print(f"低 Re の補足測定: " + ", ".join(f"Re={k} -> {v:.4f}"
                                            for k, v in sorted(extra.items()))
          + ("（それ以外は実線と重なって分離できない）" if extra else "なし"))
    pts.sort()
    cx = np.array([q[0] for q in pts])
    cy = np.array([q[1] for q in pts])
    re_pts, di_pts = to_re(cx), to_di(cy)

    keep = np.concatenate([[True], np.diff(re_pts) > 1e-9])
    f = PchipInterpolator(re_pts[keep], di_pts[keep], extrapolate=False)

    px_err = (DI_RANGE[1] - DI_RANGE[0]) / (y1 - y0)
    prev = {100: 1.02, 200: 1.07, 300: 1.12, 500: 1.37, 600: 1.48,
            1000: 1.73, 2000: 1.92}
    print(f"\n1 px = {px_err:.4f} Di。破線の太さは 2〜3 px なので "
          f"読み取り誤差は ±{1.5*px_err:.3f} 程度\n")
    print(f"{'Re':>6}{'Di（実測）':>12}{'08-09 の目視':>14}{'差':>9}")
    out = {}
    for t in (100, 200, 300, 400, 500, 600, 750, 1000, 1250, 1500, 2000):
        inside = re_pts.min() <= t <= re_pts.max()
        v = float(f(t)) if inside else float("nan")
        out[str(t)] = None if not inside else round(v, 4)
        p = prev.get(t)
        print(f"{t:>6}{(f'{v:.4f}' if inside else '範囲外'):>12}"
              f"{(f'{p:.2f}' if p else '—'):>14}"
              f"{(f'{v-p:+.3f}' if (p and inside) else '—'):>9}")

    # 実線との間隔（低 Re で分離できているかの目安）
    sol = {int(round(to_re(q[0]))): to_di(q[1]) for q in solid_pts}
    print()
    print("実線（ref CFD）との差 = opt CFD − ref CFD")
    gaps = {}
    for t in (300, 350, 400, 500, 1000):
        if out.get(str(t)) and t in sol:
            gaps[t] = round(out[str(t)] - sol[t], 4)
            print(f"  Re={t:4d}: {out[str(t)]:.4f} − {sol[t]:.4f} = {gaps[t]:+.4f}")
    print("  Re <= 250 では 2 本が 1 つの区間に融合し、分離できない")

    res = dict(source="Gamboa 2005 Fig.7, opt CFD（破線）",
               ref_cfd_solid={str(k): round(v, 4) for k, v in sol.items()
                              if k in (100, 150, 200, 250, 300, 350, 400, 500)},
               gap_to_ref=gaps,
               note="Re < 約 290 では破線と実線が重なり分離できない",
               image=os.path.basename(IMG), frame=[x0, x1, y0, y1],
               calibration=dict(Re=list(RE_RANGE), Di=list(DI_RANGE)),
               read_error_Di=float(1.5 * px_err), n_points=len(pts),
               Re_covered=[float(re_pts.min()), float(re_pts.max())],
               values=out, low_re_extra=extra,
               points=[[round(float(r), 2), round(float(d), 4)]
                       for r, d in zip(re_pts, di_pts)])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"\nsaved {os.path.normpath(OUT)}")

    fig, ax = plt.subplots(figsize=(9, 6.4))
    ax.imshow(dark, cmap="gray_r")
    ax.plot(cx, cy, "-", color="#0e7c7b", lw=1.6, label="tracked opt CFD")
    ax.plot([q[0] for q in seg], [q[1] for q in seg], ".", color="#c0521c",
            ms=4, label="dash centroids")
    ax.set_title("Fig.7 digitization: opt CFD (dashed)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "fig7_digitized.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
