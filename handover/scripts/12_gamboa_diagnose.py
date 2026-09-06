#!/usr/bin/env python3
"""
12_gamboa_diagnose.py -- Gamboa 形状の順逆の流れ場を直接見る（診断用）

背景
----
案 2（周期セル + 直線補正）で得た Di_valve は Re = 100 で 0.9925、Re = 500 で 0.9678。
Gamboa Fig.7 の 2D CFD（1.02, 1.37）と符号ごと食い違う。
幾何からの期待（逆流時はループ出口が主流に対向 -> 逆流の損失が大きい）とも矛盾する。

そこでまず**流れ場そのもの**を見る。出力は

  - phi_loop（各 x 断面でループ側を通る質量流量の割合、rho*u_x の積分）
  - ループ枝の質量流量（左枝・右枝の断面を横切る流量）
  - 速度場・流線の図

順流・逆流とも**同じ向き（+x）に駆動した格子**で計算するので、
逆流ケースは形状を x 反転したものである（09 と同じ）。
"""
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import gamboa as G
import lbm
import postproc as pp


def run(mask, Re, U, cpm, iters, tol, label):
    nu = lbm.nu_lattice(U, 2.0 * cpm, Re)
    mdot_target = U * cpm
    t0 = time.time()
    ux, uy, p, rho, f, info = lbm.solve_flow_periodic(
        mask, nu, mdot_target, iters=iters, tol=tol)
    print(f"  {label}: G={info['G']:.6e} dp_lat={info['dp_lattice']:.6e} "
          f"conv={info['converged']} it={info['iters']} "
          f"dG={info['dG_final']:.2e} err={info['mdot_err']:.2e} "
          f"{time.time()-t0:.0f}s", flush=True)
    return ux, uy, p, rho, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=500.0)
    ap.add_argument("--cpm", type=int, default=12)
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--gap", type=float, default=1.0)
    ap.add_argument("--iters", type=int, default=3_000_000)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    w = 1.0
    th = G.theta_from_params(G.OPTIMIZED["X2"], G.OPTIMIZED["n"], G.OPTIMIZED["Y3"])
    poly, isl, info = G.gamboa_valve("optimized", w=w)
    poly, isl, info = G.periodic_cell(poly, info, a.gap * w, w=w)
    mask, xs, ys = G.rasterize_cell(poly, a.cpm, w=w)
    print(f"Re={a.Re} cpm={a.cpm} gap={a.gap} mask={mask.shape} "
          f"fluid={mask.sum()} theta={th:.3f}", flush=True)

    res = {}
    fields = {}
    for lab, m, X in (("fwd", mask, xs),
                      ("rev", mask[::-1, :].copy(), -xs[::-1])):
        ux, uy, p, rho, inf = run(m, a.Re, a.U, a.cpm, a.iters, a.tol, lab)
        # phi_loop: 主流路帯 |y| <= 0.5 の外を通る質量流量の割合
        lp = pp.loop_mass_fraction(ux, m, ys, main_half_width_mm=0.5 * w, rho=rho)
        gflux = rho * ux
        mtot = np.array([gflux[i, m[i]].sum() for i in range(m.shape[0])])
        main = np.abs(ys) <= 0.5 * w + 1e-9
        mloop = np.array([gflux[i, m[i] & ~main].sum() for i in range(m.shape[0])])
        res[lab] = dict(G=inf["G"], dp_lattice=inf["dp_lattice"],
                        converged=inf["converged"], iters=inf["iters"],
                        mdot_err=inf["mdot_err"],
                        phi_loop_mean=lp["phi_loop_mean"],
                        phi_loop_max=lp["phi_loop_max"],
                        phi_loop_min=float(np.nanmin(np.where(
                            [(m[i] & ~main).sum() > 0 for i in range(m.shape[0])],
                            1.0 - np.array([gflux[i, m[i] & main].sum() for i in range(m.shape[0])]) / np.maximum(mtot, 1e-30),
                            np.nan))),
                        mdot_total=float(mtot.mean()))
        fields[lab] = dict(ux=ux, uy=uy, p=p, rho=rho, mask=m, x=X,
                           mloop=mloop, mtot=mtot)
        print(f"    phi_loop mean={lp['phi_loop_mean']:+.4f} "
              f"max={lp['phi_loop_max']:+.4f}", flush=True)

    Di = res["rev"]["dp_lattice"] / res["fwd"]["dp_lattice"]
    print(f"\n  Di_cell = {Di:.5f}   (gap={a.gap}, 直線補正なし)")
    res["Di_cell"] = Di
    res["meta"] = dict(Re=a.Re, cpm=a.cpm, gap=a.gap, U=a.U, theta_deg=th,
                       nx=int(mask.shape[0]), ny=int(mask.shape[1]),
                       fluid=int(mask.sum()))

    out = a.out or f"results/gamboa_diag_Re{int(a.Re)}_cpm{a.cpm}.json"
    os.makedirs("results", exist_ok=True)
    json.dump(res, open(out, "w"), indent=1, default=float)
    print(f"saved {out}")

    # ---- 図 ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(11, 11))
    for ax, lab in zip(axes[:2], ("fwd", "rev")):
        d = fields[lab]
        m = d["mask"]
        sp = np.hypot(d["ux"], d["uy"]) / a.U
        sp = np.where(m, sp, np.nan)
        ax.imshow(sp.T, origin="lower", cmap="viridis",
                  extent=[0, m.shape[0], ys[0] * a.cpm, ys[-1] * a.cpm],
                  aspect="equal", vmin=0, vmax=np.nanmax(sp))
        Xg, Yg = np.meshgrid(np.arange(m.shape[0]),
                             np.arange(m.shape[1]), indexing="ij")
        u = np.where(m, d["ux"], np.nan); v = np.where(m, d["uy"], np.nan)
        ax.streamplot(np.arange(m.shape[0]),
                      (np.arange(m.shape[1]) + ys[0] * a.cpm),
                      u.T, v.T, color="w", density=1.6, linewidth=0.5,
                      arrowsize=0.7)
        ax.set_title(f"{lab}  |u|/U   Re={a.Re}  cpm={a.cpm}  "
                     f"phi_loop={res[lab]['phi_loop_mean']:+.3f}")
    ax = axes[2]
    for lab, c in (("fwd", "C0"), ("rev", "C3")):
        d = fields[lab]
        ax.plot(d["mloop"] / np.maximum(d["mtot"], 1e-30), c, label=lab)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("i (格子, 駆動方向 +x)"); ax.set_ylabel("phi_loop")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    png = out.replace("results/", "figures/").replace(".json", ".png")
    os.makedirs("figures", exist_ok=True)
    fig.savefig(png, dpi=130)
    print(f"saved {png}")


if __name__ == "__main__":
    main()
