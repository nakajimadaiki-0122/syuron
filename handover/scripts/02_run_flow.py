#!/usr/bin/env python3
"""
02_run_flow.py -- 任意形状の順流／逆流を計算する
================================================

使い方
------
    python scripts/02_run_flow.py --dxf data/tesla_channel.dxf --case fwd --Re 100
    python scripts/02_run_flow.py --dxf data/tesla_channel.dxf --case rev --Re 100
    python scripts/02_run_flow.py --dxf data/straight_channel.dxf --case fwd --Re 100

`--case rev` は幾何を x 方向に鏡像反転する。入口・出口の数値処理が
順流と完全に同一になるため、Δp の比が純粋に幾何の効果を表す。

計算時間の目安（1 コア, 12 cells/mm, Re=100）
    直線流路  約  50 s（16,000 反復で収束）
    Tesla    約 150 s（16,000-20,000 反復で収束）
収束しない場合は --iters を増やして再実行すれば
results/ck_<case>.npy から再開する。
"""
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import geometry as geo
import lbm
import postproc as pp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dxf", required=True)
    ap.add_argument("--case", default="fwd", choices=["fwd", "rev"])
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--cpm", type=int, default=12, help="cells per mm")
    ap.add_argument("--U", type=float, default=0.05, help="lattice velocity")
    ap.add_argument("--W", type=float, default=1.0, help="main channel width [mm]")
    ap.add_argument("--iters", type=int, default=40000)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()

    tag = a.tag or f"{os.path.splitext(os.path.basename(a.dxf))[0]}_{a.case}_Re{int(a.Re)}"
    os.makedirs(a.outdir, exist_ok=True)

    poly = geo.load_dxf(a.dxf)
    mask, xs, ys = geo.rasterize(poly, a.cpm)
    if a.case == "rev":
        mask = geo.mirror_x(mask)
    chk = geo.check_rasterization(poly, mask, a.cpm)
    print("geometry:", json.dumps(chk), flush=True)
    if abs(chk["error_pct"]) > 0.5:
        print("  WARNING: rasterization area error > 0.5 %. Increase --cpm.", flush=True)

    # Dh = 2W （2 次元平面視 = 平行平板とみなす）
    Dh_cells = 2.0 * a.W * a.cpm
    nu = lbm.nu_lattice(a.U, Dh_cells, a.Re)
    tau = lbm.tau_from_nu(nu)
    Ma = a.U / np.sqrt(lbm.CS2)
    print(f"nu_lat={nu:.5f}  tau_p={tau:.4f}  Ma={Ma:.4f}", flush=True)
    if tau < 0.52:
        print("  WARNING: tau_p < 0.52. Raise --U or --cpm.", flush=True)
    if Ma > 0.1:
        print("  WARNING: Ma > 0.1. Compressibility error will be large.", flush=True)

    ck = os.path.join(a.outdir, f"ck_{tag}.npy")
    f0 = np.load(ck) if os.path.exists(ck) else None

    t0 = time.time()
    ux, uy, p, rho, f = lbm.solve_flow(
        mask, nu, a.U, iters=a.iters, tol=2e-6, f0=f0,
        report=lambda i, d: print(f"  it={i} d={d:.3e}", flush=True))
    wall = time.time() - t0
    np.save(ck, f)

    Dh_m = 2.0 * a.W * 1e-3
    scale = pp.lattice_to_physical_scale(a.Re, Dh_m, a.U)
    dp_lat = pp.pressure_drop(p, ux, mask)
    res = dict(tag=tag, dxf=a.dxf, case=a.case, Re=a.Re, cells_per_mm=a.cpm,
               U_lattice=a.U, nu_lattice=nu, tau_p=tau, Mach=Ma,
               wall_time_s=round(wall, 1),
               dp_lattice=dp_lat, dp_Pa=pp.to_pascal(dp_lat, scale),
               **chk)
    lf = pp.loop_mass_fraction(ux, mask, ys, rho=rho)
    res.update({k: v for k, v in lf.items() if k != "phi"})
    try:
        res["u_loop_over_u_main"] = pp.velocity_ratio(ux, mask, ys)
    except Exception:
        res["u_loop_over_u_main"] = None

    print(json.dumps(res, indent=1, default=str), flush=True)
    json.dump(res, open(os.path.join(a.outdir, f"flow_{tag}.json"), "w"),
              indent=1, default=str)
    np.savez_compressed(os.path.join(a.outdir, f"flow_{tag}.npz"),
                        ux=ux, uy=uy, p=p, rho=rho, mask=mask, xs=xs, ys=ys,
                        phi=lf["phi"])
    print("DONE", tag, flush=True)


if __name__ == "__main__":
    main()
