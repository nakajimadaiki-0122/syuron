#!/usr/bin/env python3
"""
04_run_unitcell.py -- 周期単位セルで順流・逆流を解く
====================================================

`02_run_flow.py` の入口・出口境界条件は Re >= 200 で発散する
（2026-08-09 に直線流路でも再現。README §8 B-6）。本スクリプトは
入口・出口を持たず、流れ方向に周期境界＋体積力駆動で 1 段分だけを解く。
使う数値機構は `01_validate_solver.py` で解析解に対し検証済みのものと同一。

使い方
------
    python scripts/04_run_unitcell.py --dxf data/tesla_channel.dxf --case fwd --Re 100
    python scripts/04_run_unitcell.py --dxf data/tesla_channel.dxf --case rev --Re 100
    # 切り出し位置の感度を見る
    python scripts/04_run_unitcell.py --dxf data/tesla_channel.dxf --case fwd --Re 100 --scan-offsets 5

得られる Di_p は**1 段あたり**の値であり、Porwal・Gamboa の定義と一致する。
装置全体の Di（02_run_flow.py の出力）とは直接比較できない（README §9）。

駆動は流量一定。体積力 G を制御して質量流量を目標値に合わせ、
収束後の G から Δp = G * L_cell を得る。したがって Di_p = G_逆 / G_順。
"""
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import geometry as geo
import lbm
import postproc as pp


def run_one(sub, ys, nu, U, W_cells, iters, report_prefix="",
            tol=1e-9, mdot_tol=1e-6):
    """単位セル 1 本を解く。mdot_target は主流路幅 W と平均流速 U から決める。"""
    mdot_target = U * W_cells          # rho0 = 1, 断面 W_cells の平均流速 U
    t0 = time.time()
    ux, uy, p, rho, f, info = lbm.solve_flow_periodic(
        sub, nu, mdot_target, iters=iters, tol=tol, mdot_tol=mdot_tol,
        report=lambda i, d, G, e: print(
            f"  {report_prefix}it={i} d={d:.2e} G={G:.3e} mdot_err={e:+.2e}",
            flush=True))
    info["wall_time_s"] = round(time.time() - t0, 1)
    info["mdot_target"] = mdot_target
    return ux, uy, p, rho, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dxf", required=True)
    ap.add_argument("--case", default="fwd", choices=["fwd", "rev"])
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--cpm", type=int, default=12, help="cells per mm")
    ap.add_argument("--U", type=float, default=0.05, help="lattice mean velocity")
    ap.add_argument("--W", type=float, default=1.0, help="main channel width [mm]")
    ap.add_argument("--pitch", type=float, default=6.0, help="unit cell length [mm]")
    ap.add_argument("--offset", type=int, default=None,
                    help="切り出し開始インデックス。既定は候補の中央")
    ap.add_argument("--scan-offsets", type=int, default=0,
                    help="切り出し位置を N 通り試して感度を見る")
    ap.add_argument("--iters", type=int, default=400000)
    ap.add_argument("--tol", type=float, default=1e-9,
                    help="収束判定 max|du|。既定 1e-9")
    ap.add_argument("--mdot-tol", type=float, default=1e-6,
                    help="流量制御の許容誤差。既定 1e-6")
    ap.add_argument("--suffix", default="", help="タグ末尾に付ける識別子")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    poly = geo.load_dxf(a.dxf)
    mask, xs, ys = geo.rasterize(poly, a.cpm)
    period_cells = int(round(a.pitch * a.cpm))

    exact = geo.periodic_windows(mask, period_cells)
    seam = geo.seam_windows(mask, period_cells)
    if not seam:
        raise SystemExit(f"ピッチ {a.pitch} mm で継ぎ目が一致する窓がない。"
                         f"--pitch を確認すること。")
    i_def, n_mis, note = geo.choose_unit_cell(mask, ys, period_cells)
    print(f"厳密周期の窓: {len(exact)} 通り / 継ぎ目一致の窓: {len(seam)} 通り", flush=True)
    print(f"既定の切り出し: i0={i_def} (x={xs[i_def]:.3f} mm), {note}, "
          f"窓内部の非周期列 {n_mis}", flush=True)
    if n_mis > 0:
        print(f"  注意: 単位セルは隣の段と厳密に同一ではない（{n_mis} 列が相違）。"
              f"継ぎ目は一致しているので壁は連続する。"
              f"原因は DXF のループごとの折れ線近似の差（README §7.2）。", flush=True)

    seam_i = [i for i, _ in seam]
    if a.scan_offsets > 0:
        pool = [i for i, n in seam if n == 0] or seam_i
        offsets = [int(pool[k]) for k in
                   np.linspace(0, len(pool) - 1, a.scan_offsets).astype(int)]
    else:
        offsets = [a.offset if a.offset is not None else i_def]
    for o in offsets:
        if o not in set(seam_i):
            raise SystemExit(f"offset {o} は継ぎ目が一致しない。"
                             f"候補例: {seam_i[:5]} ...")

    Dh_cells = 2.0 * a.W * a.cpm
    nu = lbm.nu_lattice(a.U, Dh_cells, a.Re)
    W_cells = a.W * a.cpm
    print(f"nu_lat={nu:.5f}  tau_p={lbm.tau_from_nu(nu):.4f}  "
          f"Ma={a.U/np.sqrt(lbm.CS2):.4f}  unit={period_cells} cells", flush=True)

    Dh_m = 2.0 * a.W * 1e-3
    scale = pp.lattice_to_physical_scale(a.Re, Dh_m, a.U)

    out = []
    for o in offsets:
        sub, xsub = geo.cut_unit_cell(mask, xs, o, period_cells)
        if a.case == "rev":
            sub = geo.mirror_x(sub)
        tag = (f"{os.path.splitext(os.path.basename(a.dxf))[0]}_cell"
               f"_{a.case}_Re{int(a.Re)}_cpm{a.cpm}_o{o}{a.suffix}")
        print(f"\n=== {tag}  (x = {xsub[0]:.3f} .. {xsub[-1]:.3f} mm, "
              f"fluid {sub.sum()} cells) ===", flush=True)
        ux, uy, p, rho, info = run_one(sub, ys, nu, a.U, W_cells, a.iters,
                                       tol=a.tol, mdot_tol=a.mdot_tol)

        lf = pp.loop_mass_fraction(ux, sub, ys, rho=rho)
        res = dict(tag=tag, dxf=a.dxf, case=a.case, Re=a.Re, cells_per_mm=a.cpm,
                   U_lattice=a.U, nu_lattice=nu, Mach=float(a.U/np.sqrt(lbm.CS2)),
                   pitch_mm=a.pitch, offset=o,
                   x_start_mm=float(xsub[0]), x_end_mm=float(xsub[-1]),
                   fluid_cells=int(sub.sum()),
                   tol=a.tol, mdot_tol=a.mdot_tol,
                   G=info["G"], converged=info["converged"], iters=info["iters"],
                   dG_final=info["dG_final"], u_residual=info["u_residual"],
                   mdot=info["mdot"], mdot_target=info["mdot_target"],
                   mdot_err=info["mdot_err"], tau_p=info["tau_p"],
                   wall_time_s=info["wall_time_s"],
                   dp_lattice=info["dp_lattice"],
                   dp_Pa=pp.to_pascal(info["dp_lattice"], scale))
        res.update({k: v for k, v in lf.items() if k != "phi"})
        try:
            res["u_loop_over_u_main"] = pp.velocity_ratio(ux, sub, ys)
        except Exception:
            res["u_loop_over_u_main"] = None

        print(json.dumps({k: v for k, v in res.items()}, indent=1, default=str),
              flush=True)
        json.dump(res, open(os.path.join(a.outdir, f"cell_{tag}.json"), "w"),
                  indent=1, default=str)
        np.savez_compressed(os.path.join(a.outdir, f"cell_{tag}.npz"),
                            ux=ux, uy=uy, p=p, rho=rho, mask=sub,
                            xs=xsub, ys=ys, phi=lf["phi"])
        out.append(res)

    if len(out) > 1:
        dp = np.array([r["dp_Pa"] for r in out])
        phi = np.array([r["phi_loop_mean"] for r in out])
        print(f"\n--- 切り出し位置の感度 ({len(out)} 通り) ---")
        print(f"  dp_Pa      : {dp.min():.4f} .. {dp.max():.4f}  "
              f"(spread {100*(dp.max()-dp.min())/dp.mean():.3f} %)")
        print(f"  phi_loop   : {phi.min():.5f} .. {phi.max():.5f}  "
              f"(spread {100*(phi.max()-phi.min())/phi.mean():.3f} %)")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
