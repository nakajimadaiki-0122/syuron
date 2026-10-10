#!/usr/bin/env python3
"""
32_run_xia2025.py -- Xia 2025 の TVM / SYMTVM を解く
====================================================

計画と制約は `docs/xia2025_sim_plan.md` にある。**先に読むこと。**

要点だけ再掲する。

- 原典の Re は **513〜1755**（k-e 乱流モデル）。**本ソルバでは届かない。**
  安定条件 Re <= 10 * N_Dh（N_Dh = D_h あたりのセル数）と計算時間で律速される
- できるのは**低 Re での流動の比較**。Di の大小関係と立ち上がりを見る
- Re は原典の定義（**D_h = 2 W_ch H_ch / (W_ch + H_ch)**）。
  TVM 0.5138 mm、SYMTVM 0.8774 mm

2 つの解き方
------------

`--mode io`（既定、2 次元向き）
    全長 50 mm をそのまま解く。入口一様流速・出口定圧。
    流路の両端は水平な直線区間なので、プレナムを足さなくても BC 面が立つ。

`--mode periodic`（3 次元向き）
    流路はピッチ 4.72 mm で厳密に周期なので、1 ピッチを x 周期で解く。
    体積力駆動なので入口・出口 BC の不安定が無い。
    **切断面の断面形が前後で一致することを確認してから回す**（README B-12）。

向き
----

Xia の Forward / Reverse は**同じ流路を逆向きに流したもの**（Fig.3 は入口を
左に揃えて描いてあるので 180 度回転に見える）。したがって Forward 形状の
マスクを x 反転すれば Reverse になる。形状ファイルを 2 つ作る必要はない。

    python scripts/32_run_xia2025.py --shape symtvm --dir fwd --dim 2 --Re 100
    python scripts/32_run_xia2025.py --shape tvm --dir rev --dim 3 --mode periodic
"""
import argparse, json, os, sys, time

import numpy as np
from shapely.geometry import box

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import xia2025 as X
import gamboa_full as GF          # rasterize を借りる
import lbm
import lbm3d as L3
import postproc as pp

OUT = os.path.join(HERE, "..", "results", "xia")

# 原典 Table 2 の入口幅（mm）。D_h = 2 Wch Hch / (Wch + Hch)
W_IN = {"tvm": 0.31, "symtvm": 0.62}
SHAPE_NAME = {"tvm": "TVM-Forward", "symtvm": "SYMTVM-Forward"}


def dh_mm(shape):
    w, h = W_IN[shape], X.H_CH
    return 2.0 * w * h / (w + h)


def geometry(shape, mode, pitch_index=0, pitches=1):
    """流体の多角形と説明を返す。

    mode="periodic" のときは島の無い接続部で 1 ピッチ切り出す。
    切断位置は STL の `cuts` と同じ x = X_C1 + PITCH*k + 3.0。
    """
    poly, info = X.all_four()[SHAPE_NAME[shape]]
    note = f"{SHAPE_NAME[shape]} 全長 {X.LY} mm"
    extra = {}
    if mode == "periodic":
        x0 = X.X_C1 + X.PITCH * pitch_index + 3.0
        x1 = x0 + X.PITCH * pitches
        cut = poly.intersection(box(x0, -X.W_UNIT, x1, X.W_UNIT))
        if cut.geom_type != "Polygon":
            raise SystemExit(f"切り出しが {cut.geom_type} になった。位置を変えること")
        note = (f"{SHAPE_NAME[shape]} {pitches} ピッチ "
                f"x = {x0:.2f}〜{x1:.2f} mm")
        # **多角形は切らずに返す。** 周期性の確認で x1 の先の断面を見るため
        extra = dict(x0=x0, x1=x1, pitches=pitches,
                     cut_area_mm2=float(cut.area), cut_islands=len(cut.interiors),
                     cut_bounds=[float(v) for v in cut.bounds])
        return poly, dict(note=note, bounds=extra["cut_bounds"],
                          area_mm2=extra["cut_area_mm2"],
                          islands=extra["cut_islands"], **extra)
    return poly, dict(note=note, bounds=[float(v) for v in poly.bounds],
                      area_mm2=float(poly.area), islands=len(poly.interiors))


def rasterize_exact(poly, x0, x1, cpm, pad=1):
    """x = x0〜x1 をちょうど nx セルで刻む。**列を一切捨てない。**

    `gamboa_full.rasterize` は空の端列を落とすので、周期セルでは
    切断面がずれて位相が壊れる。周期長がセルの整数倍になる cpm を使い、
    こちらを使うこと。
    """
    from matplotlib.path import Path
    span = x1 - x0
    nx = int(round(span * cpm))
    if abs(span * cpm - nx) > 1e-6:
        raise SystemExit(f"周期長 {span} mm x {cpm} cells/mm = {span*cpm} は整数でない。"
                         f"cpm を 25 の倍数にすること（4.72 x 25 = 118）")
    ymin, ymax = poly.bounds[1], poly.bounds[3]
    ny = int(np.ceil((ymax - ymin) * cpm))
    # nx + 1 列を取る。最後の 1 列は周期性の確認にだけ使って捨てる
    # （周期配列では列 0 と列 nx が同じ位置にあたる。列 nx-1 ではない）
    xs = x0 + (np.arange(nx + 1) + 0.5) / cpm
    ys = ymin + (np.arange(ny) + 0.5) / cpm
    Xg, Yg = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([Xg.ravel(), Yg.ravel()])
    inside = Path(np.asarray(poly.exterior.coords)).contains_points(pts)
    for ring in poly.interiors:
        inside &= ~Path(np.asarray(ring.coords)).contains_points(pts)
    core = inside.reshape(nx + 1, ny)
    m = np.zeros((nx + 1, ny + 2 * pad), bool)  # y の端は必ず固体（B-1）
    m[:, pad:pad + ny] = core
    return m[:nx], xs[:nx], ys, m[nx]           # 最後の列は確認用に別途返す


def check_periodic(first, wrap):
    """周期配列として繋がるか。

    列 0 が、周期長だけ進んだ位置の断面（`wrap`）と一致するかを見る。
    `mask[-1]` と比べるのは 1 セルぶんずれていて誤り。
    """
    same = bool(np.array_equal(first, wrap))
    return same, int(first.sum()), int(wrap.sum()), int((first ^ wrap).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", required=True, choices=["tvm", "symtvm"])
    ap.add_argument("--dir", default="both", choices=["both", "fwd", "rev"])
    ap.add_argument("--dim", type=int, default=2, choices=[2, 3])
    ap.add_argument("--mode", default="io", choices=["io", "periodic"])
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--cpm", type=int, default=48, help="cells/mm")
    ap.add_argument("--U", type=float, default=0.05, help="入口平均流速（格子単位）")
    ap.add_argument("--pitch-index", type=int, default=0)
    ap.add_argument("--pitches", type=int, default=1,
                    help="周期セルに入れるピッチ数。TVM は島が上下交互なので 2 が要る")
    ap.add_argument("--iters", type=int, default=600000)
    ap.add_argument("--dp-window", type=int, default=10000)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    Dh = dh_mm(a.shape)
    n_dh = Dh * a.cpm                       # D_h あたりのセル数
    re_max = 10.0 * n_dh                    # 安定限界の見込み（B-9 から）
    poly, info = geometry(a.shape, a.mode, a.pitch_index, a.pitches)
    wrap2 = None
    if a.mode == "periodic":
        m2, xs, ys, wrap2 = rasterize_exact(poly, info["x0"], info["x1"], a.cpm)
    else:
        m2, xs, ys = GF.rasterize(poly, a.cpm, w=1.0)
    depth = max(2, int(round(X.H_CH * a.cpm)))
    mask = m2 if a.dim == 2 else L3.extrude(m2, depth)
    nu = L3.nu_lattice(a.U, n_dh, a.Re) if a.dim == 3 else lbm.nu_lattice(a.U, n_dh, a.Re)

    tag = f"xia_{a.shape}_{a.mode}_{a.dim}d_Re{a.Re:.0f}_cpm{a.cpm}"
    print("=" * 92)
    print(f" {info['note']}   {a.dim} 次元 / {a.mode}   Re = {a.Re}, "
          f"cpm = {a.cpm}" + (f", 深さ {depth} セル" if a.dim == 3 else ""))
    print("=" * 92)
    b = info["bounds"]
    print(f"  外形 {b[2]-b[0]:.2f} x {b[3]-b[1]:.2f} mm、島 {info['islands']} 個")
    print(f"  D_h = {Dh:.4f} mm = {n_dh:.1f} セル、流路幅 "
          f"{W_IN[a.shape]:.2f} mm = {W_IN[a.shape]*a.cpm:.1f} セル")
    print(f"  格子 {' x '.join(str(v) for v in mask.shape)} = "
          f"{np.prod(mask.shape)/1e6:.2f} M、流体 {int(mask.sum())/1e3:.0f} k")
    print(f"  nu = {nu:.5f}, tau_p = {lbm.tau_from_nu(nu):.4f}")
    print(f"  安定限界の見込み Re <= {re_max:.0f}"
          + ("" if a.Re <= re_max else "   ** 超えている。発散する見込み **"))
    if a.Re > re_max:
        print("  条件を出さない方針（docs/xia2025_sim_plan.md §4）。中止する。")
        return

    if a.mode == "periodic":
        wrap = wrap2 if a.dim == 2 else L3.extrude(wrap2[None], depth)[0]
        same, na, nb, nd = check_periodic(mask[0], wrap)
        print(f"  周期性の確認: 両端の断面 {'一致' if same else '**不一致**'} "
              f"（流体セル {na} / {nb}、差 {nd}）", flush=True)
        if not same:
            print("  切断面が一致しない。位相を壊すので中止する（README B-12）。")
            return

    out = dict(tag=tag, shape=a.shape, name=SHAPE_NAME[a.shape], mode=a.mode,
               dim=a.dim, Re=a.Re, cpm=a.cpm, U=a.U, Dh_mm=Dh, n_Dh=n_dh,
               re_max_est=re_max, depth=(depth if a.dim == 3 else None),
               nu=nu, geometry=info, grid=[int(v) for v in mask.shape],
               fluid=int(mask.sum()))

    todo = [("fwd", mask), ("rev", mask[::-1].copy())]
    if a.dir != "both":
        todo = [t for t in todo if t[0] == a.dir]

    # 格子 -> 物理のスケール。長さの基準は D_h（原典の定義）
    scale = pp.lattice_to_physical_scale(a.Re, Dh * 1e-3, a.U)

    for lab, m in todo:
        t0 = time.time()
        if a.mode == "periodic":
            A = float(m.sum() / m.shape[0])          # 断面の平均流体セル数
            mdot = a.U * A
            print(f"  [{lab}] 体積力駆動、断面の流体 {A:.0f} セル", flush=True)

            def repp(it, d, g, err):
                # lbm3d.solve_duct_periodic は (it, 残差, G, 流量誤差) を渡す
                if it % 5000 == 0:
                    print(f"    it={it:7d} res={d:.2e} G={g:.4e} "
                          f"流量誤差={err:+.2e}", flush=True)

            if a.dim == 2:
                raise SystemExit("periodic は 3 次元用（2 次元は io を使う）")
            ux, uy, uz, rho, inf = L3.solve_duct_periodic(
                m, nu, mdot, iters=a.iters, report=repp)
            dp_lat = inf["dp_lattice"]
            dp_std = 0.0
            conv, unsteady = inf["converged"], False
            spread = float("nan")
            resid = float("nan")
            iters = inf["iters"]
            out.setdefault("periodic_info", {})[lab] = dict(
                G=float(inf["G"]), mdot=float(inf["mdot"]),
                mdot_target=float(mdot), tau_p=float(inf["tau_p"]))
        else:
            n_in = int(m[0].sum())
            U_in = a.U
            print(f"  [{lab}] 入口 {n_in} セル, U_in = {U_in:.5f}", flush=True)

            def rep2(it, d, mi, mo, dp=None):
                if it % 20000 == 0:
                    print(f"    it={it:7d} res={d:.2e} dp={dp:.6e}", flush=True)

            def rep3(it, d, dp, sp):
                if it % 5000 == 0:
                    print(f"    it={it:7d} res={d:.2e} dp={dp:.6e} "
                          f"質量ばらつき={sp:.2e}", flush=True)

            if a.dim == 2:
                ux, uy, p, rho, f, inf = lbm.solve_flow_io(
                    m, nu, U_in, iters=a.iters, dp_window=a.dp_window,
                    report=rep2)
            else:
                ux, uy, uz, rho, p, f, inf = L3.solve_flow_io(
                    m, nu, U_in, iters=a.iters, dp_window=a.dp_window,
                    report=rep3)
            dp_lat = inf["dp_lattice_mean"]
            dp_std = inf["dp_lattice_std"]
            conv, unsteady = inf["converged"], inf["unsteady"]
            spread = inf.get("mass_spread_final", float("nan"))
            resid = inf["residual"]
            iters = inf["iters"]

        r = dict(label=lab, dp_lattice=float(dp_lat),
                 dp_Pa=pp.to_pascal(dp_lat, scale),
                 dp_Pa_std=pp.to_pascal(dp_std, scale),
                 converged=bool(conv), unsteady=bool(unsteady),
                 iters=int(iters), residual=float(resid),
                 mass_spread=float(spread),
                 wall_s=round(time.time() - t0, 1))
        out[lab] = r
        print(f"  [{lab}] dp = {r['dp_Pa']:.4f} ± {r['dp_Pa_std']:.4f} Pa, "
              f"収束 {conv} it={iters} {r['wall_s']/60:.0f} 分", flush=True)

    if "fwd" in out and "rev" in out:
        out["Di"] = out["rev"]["dp_Pa"] / out["fwd"]["dp_Pa"]
        print(f"  Di = {out['Di']:.4f}")

    j = os.path.join(OUT, f"{tag}_{a.dir}.json")
    json.dump(out, open(j, "w"), indent=1, default=float)
    print(f"saved {os.path.normpath(j)}")


if __name__ == "__main__":
    main()
