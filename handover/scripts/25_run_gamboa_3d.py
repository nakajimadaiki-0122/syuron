#!/usr/bin/env python3
"""
25_run_gamboa_3d.py -- Gamboa 単発モデルを 3 次元で解く
======================================================

2 次元（`17_run_gamboa_full.py`）と**同じ形状・同じ流速・同じ流体**で、
深さ方向だけを有限にする。違いは上下壁があることだけなので、
差はそのまま「3 次元にした効果」になる。

Re の定義について
-----------------
Gamboa は 2 次元なので D_h = 2 w_v。3 次元の正方形断面では D_h = w_v。
ここでは**2 次元と同じ条件で解く**ため、粘性を 2 次元と同じ式

    nu = U * (2 * cpm) / Re

で決める（= Gamboa の Re 定義）。正方形断面としての Reynolds 数は
その半分になるので、両方を出力する。

比較の要点
----------
平行平板の f*Re = 96 に対し、正方形ダクトは 56.91。
**同じ流量なら 3 次元の方が圧損は大きく**（壁が増える）、
ループへの流量配分も変わる。Di がどう動くかがこの計算の目的。
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import gamboa_full as GF
import lbm3d as L3
import postproc as pp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=100.0,
                    help="Gamboa の定義（D_h = 2 w_v）")
    ap.add_argument("--cpm", type=int, default=8, help="流路幅のセル数")
    ap.add_argument("--depth", type=int, default=None,
                    help="流路の深さ [セル]。既定は cpm（正方形断面）")
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--iters", type=int, default=600000)
    ap.add_argument("--dir", default="both", choices=["both", "fwd", "rev"])
    ap.add_argument("--dp-window", type=int, default=10000)
    ap.add_argument("--out", default=None)
    ap.add_argument("--save-fields", action="store_true")
    a = ap.parse_args()
    depth = a.depth if a.depth else a.cpm
    w_mm = 1.0

    poly, info = GF.build("optimized", w=1.0, plenum_out="vertical")
    m2, xs, ys = GF.rasterize(poly, a.cpm, w=1.0)
    mask3 = L3.extrude(m2, depth)
    nx, ny, nz = mask3.shape
    nu = L3.nu_lattice(a.U, 2.0 * a.cpm, a.Re)
    print("=" * 84)
    print(f" Gamboa optimized 3D   Re = {a.Re}（D_h = 2w 基準）, "
          f"cpm = {a.cpm}, 深さ = {depth} セル")
    print("=" * 84)
    print(f"  領域 {nx} x {ny} x {nz} = {nx*ny*nz/1e6:.2f} M セル、"
          f"流体 {int(mask3.sum())/1e3:.0f} k セル（{100*mask3.mean():.1f} %）")
    print(f"  nu = {nu:.5f}, tau_p = {L3.tau_from_nu(nu):.4f}, "
            f"正方形断面としての Re = {a.Re/2:.0f}")
    print(f"  f は {19*nx*ny*nz*8/1e6:.0f} MB", flush=True)

    out = dict(Re=a.Re, Re_square=a.Re / 2, cpm=a.cpm, depth=depth, U=a.U,
               nu=nu, nx=nx, ny=ny, nz=nz, fluid=int(mask3.sum()),
               dim=3, which="optimized")
    todo = [("fwd", mask3), ("rev", mask3[::-1].copy())]
    if a.dir != "both":
        todo = [t for t in todo if t[0] == a.dir]

    for lab, m in todo:
        n_in = int(m[0].sum())
        Q = a.U * a.cpm * depth               # 目標体積流量（流路断面）
        U_in = Q / n_in
        print(f"  [{lab}] 入口 {n_in} セル, U_in = {U_in:.5f}", flush=True)

        def rep(it, d, dp, spread):
            if it % 5000 == 0:
                print(f"    it={it:7d} res={d:.2e} dp={dp:.6e} "
                      f"質量ばらつき={spread:.2e}", flush=True)

        t0 = time.time()
        ck = os.path.join(HERE, "..", "results", "single", f"ck3d_{lab}_Re{int(a.Re)}_cpm{a.cpm}.npy")
        ux, uy, uz, rho, p, f, inf = L3.solve_flow_io(
            m, nu, U_in, iters=a.iters, dp_window=a.dp_window, report=rep,
            ckpt=ck)
        dp_lat = inf["dp_lattice_mean"]
        scale = pp.lattice_to_physical_scale(a.Re, 2.0 * w_mm * 1e-3, a.U)
        dp_Pa = pp.to_pascal(dp_lat, scale)
        r = dict(label=lab, n_inlet=n_in, U_in=U_in,
                 dp_lattice=dp_lat, dp_Pa=dp_Pa,
                 dp_Pa_std=pp.to_pascal(inf["dp_lattice_std"], scale),
                 converged=inf["converged"], unsteady=inf["unsteady"],
                 iters=inf["iters"], residual=inf["residual"],
                 mass_spread=inf["mass_spread_final"],
                 wall_s=round(time.time() - t0, 1))
        out[lab] = r
        print(f"  [{lab}] dp = {dp_lat:.6e} = {dp_Pa:.4f} Pa, "
              f"収束 {inf['converged']}{'（非定常）' if inf['unsteady'] else ''} "
              f"it={inf['iters']} 質量ばらつき={inf['mass_spread_final']:.1e} "
              f"{r['wall_s']/60:.0f} 分", flush=True)
        if os.path.exists(ck):
            os.remove(ck)
        if a.save_fields:
            npz = os.path.join(HERE, "..", "results", "single", f"gamboa3d_Re{int(a.Re)}_cpm{a.cpm}_{lab}.npz")
            np.savez_compressed(npz, ux=ux.astype(np.float32),
                                uy=uy.astype(np.float32),
                                uz=uz.astype(np.float32),
                                p=p.astype(np.float32), mask=m, xs=xs, ys=ys)
            print(f"  saved {os.path.normpath(npz)}", flush=True)

    o = a.out or os.path.join(HERE, "..", "results", "single", f"gamboa3d_Re{int(a.Re)}_cpm{a.cpm}"
                              + ("" if a.dir == "both" else f"_{a.dir}")
                              + ".json")
    if a.dir != "both":
        other = "rev" if a.dir == "fwd" else "fwd"
        op = o.replace(f"_{a.dir}.json", f"_{other}.json")
        if os.path.exists(op):
            j = json.load(open(op))
            if other in j:
                out[other] = j[other]
    os.makedirs(os.path.dirname(o), exist_ok=True)
    if "fwd" in out and "rev" in out:
        out["Di"] = out["rev"]["dp_Pa"] / out["fwd"]["dp_Pa"]
        print(f"\n  Di(3D) = {out['rev']['dp_Pa']:.4f} / "
              f"{out['fwd']['dp_Pa']:.4f} = {out['Di']:.4f}")
        two = os.path.join(HERE, "..", "results", "single", f"gamboa_full_Re{int(a.Re)}_cpm16_fwd.json")
        if os.path.exists(two):
            j = json.load(open(two))
            if "Di" in j:
                print(f"  2 次元（cpm=16）の Di = {j['Di']:.4f}  "
                      f"-> 差 {out['Di']-j['Di']:+.4f}")
                out["Di_2d"] = j["Di"]
    json.dump(out, open(o, "w"), indent=1, default=float)
    print(f"saved {os.path.normpath(o)}")


if __name__ == "__main__":
    main()
