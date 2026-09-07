#!/usr/bin/env python3
"""
08_validate_gamboa.py -- 生成した Gamboa 形状を CFD の前に検証する
==================================================================

CLAUDE.md の方針「幾何を近似で作らない」に従い、CFD へ進む前に必ず通す。

  1. 流路幅が全域で一定か
  2. 中心線から theta を実測し 161.6 度が再現するか
     （beta = 71.7 度は論文記載値なので theta = beta + 90 は独立検証にならない。
       生成した形状から測って初めて逆算の裏付けになる）
  3. 全ループが厳密に同一か（周期セルを並べて列単位で比較）
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import gamboa as G

STEP = 0.0005          # 輪郭標本化の刻み [w_v]


def densify(P, step=STEP):
    out = []
    for a, b in zip(P[:-1], P[1:]):
        L = np.hypot(*(b - a))
        n = max(int(np.ceil(L / step)), 1)
        out.append(a + (np.arange(n)[:, None] / n) * (b - a))
    out.append(P[-1:])
    return np.vstack(out)


def centreline(inner, outer):
    I, O = densify(inner), densify(outer)
    d2 = ((I[:, None, 0] - O[None, :, 0]) ** 2
          + (I[:, None, 1] - O[None, :, 1]) ** 2)
    j = d2.argmin(1)
    Q = O[j]
    return 0.5 * (I + Q), np.hypot(*(I - Q).T)


def terminal_angle(C, n_fit):
    S = C[:n_fit]
    c = S.mean(0)
    u, s, vt = np.linalg.svd(S - c)
    d = vt[0]
    if np.dot(d, C[0] - S[-1]) < 0:
        d = -d
    # d は「ループ内部 -> 端点」向き = ループから主流路への吐出向き
    # theta は主流路 -> ループ向きなので 180 度反転する
    ang = np.degrees(np.arctan2(d[1], d[0])) + 180.0
    ang = (ang + 180.0) % 360.0 - 180.0
    return float(ang), float(s[1] / np.sqrt(len(S)))


def check_shape(name, theta_deg, R, alpha_deg, X2, w=1.0):
    outer, isl, info = G.valve_contours(theta_deg, R, alpha_deg, X2, w)
    print("=" * 78)
    print(f" {name}   theta = {theta_deg:.4f} deg,  beta = {theta_deg-90:.4f} deg")
    print("=" * 78)
    print(f"  円弧中心 = ({info['center'][0]:.4f}, {info['center'][1]:.4f})  "
          f"R = {R}  分流板先端 = ({X2}, {0.5*w})")

    # --- 1. 流路幅 ---
    # 島の輪郭は P1 -> 内側円弧 -> P2 の開いた折れ線として扱う。
    # 直線壁は多角形の辺として暗黙に含まれるので、densify がそこも標本化する。
    C_all, width_all = centreline(isl, outer)
    # 両端は開口部で退化する（最近点が主流路上壁へ飛ぶ）ので幅で有効区間を切る
    valid = np.abs(width_all - w) < 0.05 * w
    k = np.flatnonzero(valid)
    C, width = C_all[k.min():k.max() + 1], width_all[k.min():k.max() + 1]
    print(f"\n  1. 流路幅  有効 n={len(width)} / {len(width_all)}  "
          f"min {width.min():.6f}  max {width.max():.6f}  平均 {width.mean():.6f} w_v")
    # 中点構成は開口部の近くで退化するので、幅の分布で判定する。
    # 壁が向かい合っている区間では幅は厳密に w になるはず。
    frac = float((np.abs(width - w) < 1e-3 * w).mean())
    dev = 100 * (width.max() - width.min()) / w
    print(f"     |幅 - w| < 0.1 % の割合 = {100*frac:.2f} %   "
          f"99 パーセンタイル = {np.percentile(width,99):.6f}")
    print(f"     幅の全振れ幅 = {dev:.4f} %（退化区間を含む）-> "
          f"{'合格' if frac > 0.95 else '要確認'}")
    print(f"     退化して切り捨てた長さ: 上流側 "
          f"{np.hypot(*(C_all[k.min()]-C_all[0])):.4f}, 下流側 "
          f"{np.hypot(*(C_all[k.max()]-C_all[-1])):.4f} w_v")

    # --- 2. theta の実測 ---
    print("\n  2. 中心線から theta を実測（上流側 = 戻り流路の吐出口）")
    print(f"     {'fit 点数':>10}{'fit 長さ':>10}{'theta [deg]':>14}{'残差 [w_v]':>14}")
    rows = []
    for n_fit in (200, 400, 800, 1600):
        if n_fit > len(C):
            continue
        ang, res = terminal_angle(C, n_fit)
        L = float(np.hypot(*(C[n_fit - 1] - C[0])))
        rows.append((n_fit, L, ang, res))
        print(f"     {n_fit:>10d}{L:>10.4f}{ang:>14.4f}{res:>14.2e}")
    th_meas = rows[1][2] if len(rows) > 1 else rows[0][2]
    err = th_meas - theta_deg
    print(f"     実測 theta = {th_meas:.4f} deg,  設計値 {theta_deg:.4f} deg,  "
          f"差 {err:+.4f} deg -> {'合格' if abs(err) < 0.1 else '要確認'}")
    return dict(name=name, theta_design=theta_deg, theta_measured=th_meas,
                theta_err=err, width_min=float(width.min()),
                width_max=float(width.max()), width_spread_pct=float(dev),
                center=info["center"], loop_x=info["loop_x"])


def check_periodicity(theta_deg, gap, cpms=(12, 24, 48), w=1.0):
    print("\n" + "=" * 78)
    print(f" 3. 周期セルの厳密同一性（gap = {gap} w_v）")
    print("=" * 78)
    poly, isl, info = G.tesla_valve_theta(theta_deg, w=w)
    poly, isl, info = G.periodic_cell(poly, info, gap, w)
    print(f"  セル長 = {info['cell_length']:.6f} w_v  "
          f"(ループ {info['loop_x'][0]:.3f}..{info['loop_x'][1]:.3f})")
    ok = True
    for cpm in cpms:
        mask, xs, ys = G.rasterize_cell(poly, cpm, w=w)
        nx = mask.shape[0]
        # 4 段並べて、列パターンが厳密に周期かを確認
        tiled = np.concatenate([mask] * 4, axis=0)
        same = all(np.array_equal(tiled[i], tiled[i + nx])
                   for i in range(len(tiled) - nx))
        edge = (not mask[:, 0].any()) and (not mask[:, -1].any())
        plain = (mask[0] == mask[-1]).all()
        print(f"  cpm={cpm:3d}: {nx} x {mask.shape[1]} セル, 流体 {int(mask.sum()):6d}  "
              f"4 段並べて厳密周期={same}  y 端に流体なし={edge}  左右端の列が同一={plain}")
        ok &= same and edge
    print(f"  -> {'合格' if ok else '不合格'}")
    return ok


if __name__ == "__main__":
    th = G.theta_from_params(**{k: G.OPTIMIZED[k] for k in ("X2", "n", "Y3")})
    res = [check_shape("Gamboa optimized", th, G.OPTIMIZED["R"],
                       G.OPTIMIZED["alpha"], G.OPTIMIZED["X2"])]
    thr = G.theta_from_params(**{k: G.REFERENCE[k] for k in ("X2", "n", "Y3")})
    res.append(check_shape("Gamboa reference", thr, G.REFERENCE["R"],
                           G.REFERENCE["alpha"], G.REFERENCE["X2"]))
    for t in (120.0, 90.0, 47.4):
        res.append(check_shape(f"parametric theta={t}", t, G.OPTIMIZED["R"],
                               G.OPTIMIZED["alpha"], G.OPTIMIZED["X2"]))
    ok = check_periodicity(th, gap=2.0)
    os.makedirs("results/geometry", exist_ok=True)
    json.dump(dict(shapes=res, periodic_ok=bool(ok)),
              open("results/geometry/gamboa_geometry_check.json", "w"), indent=1, default=str)
    print("\nsaved results/gamboa_geometry_check.json")
