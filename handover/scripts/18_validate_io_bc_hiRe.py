#!/usr/bin/env python3
"""
18_validate_io_bc_hiRe.py -- 入口・出口 BC を高 Re で検証する（放物線入口）
==========================================================================

`16_validate_io_bc.py` は入口を**一様流速**にしている（Gamboa の条件）。
これは壁ぎわで速度が不連続なので、高 Re では入口直後の剪断が強く、
Re = 500・L = 120 w_v の直線流路では it = 33000 で発散した（2026-09-06）。

ここでは入口を**放物線分布**にして、流れを最初から発達させた状態で解く。
こうすると助走区間が要らないので短い流路で済み、BC そのものの精度と安定性を
高 Re で分離して確認できる。

判定
  1. u_max/u_avg が離散理論値（NY=16 で約 1.489）に一致するか
  2. f_Darcy * Re = 96（弱圧縮性による加速ぶんを補正した値）
  3. Re = 500, 1000 で発散しないか
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import lbm


def run(Re, cpm, L_wv, U=0.05, iters=400000, tol=1e-8, pad=2, ramp=5000):
    nx = int(round(L_wv * cpm))
    mask = np.zeros((nx, cpm + 2 * pad), bool)
    mask[:, pad:-pad] = True
    Dh = 2.0 * cpm
    nu = lbm.nu_lattice(U, Dh, Re)
    yy = (np.arange(cpm) + 0.5) / cpm
    prof = 6.0 * U * yy * (1.0 - yy)          # 断面平均が U になる放物線

    ux, uy, p, rho, f, info = lbm.solve_flow_io(
        mask, nu, prof, iters=iters, tol=tol, ramp=ramp)

    i0, i1 = int(0.30 * nx), int(0.80 * nx)
    im = (i0 + i1) // 2

    def pm(i):
        w = np.maximum(ux[i, pad:-pad], 1e-12)
        return float(np.average(p[i, pad:-pad], weights=w))

    u_avg = float(ux[im, pad:-pad].mean())
    rho_m = float(rho[im, pad:-pad].mean())
    A = float(cpm)
    mom = (1.2 * (rho * ux ** 2)[i0, pad:-pad].sum() / A
           - 1.2 * (rho * ux ** 2)[i1, pad:-pad].sum() / A)
    dpdx = (pm(i0) - pm(i1)) / (i1 - i0)
    dpdx_f = (pm(i0) - pm(i1) + mom) / (i1 - i0)
    Re_m = u_avg * Dh / nu
    fl = np.array([(rho * ux)[i, pad:-pad].sum() for i in range(1, nx - 1)])
    return dict(Re=Re, cpm=cpm, L_wv=L_wv, nx=nx, nu=nu, tau_p=info["tau_p"],
                converged=info["converged"], iters=info["iters"],
                residual=info["residual"], wall_s=info["wall_s"],
                u_max_over_avg=float(ux[i1, pad:-pad].max()
                                     / ux[i1, pad:-pad].mean()),
                Re_measured=Re_m,
                fRe=float(dpdx * Dh / (0.5 * rho_m * u_avg ** 2) * Re_m),
                fRe_corrected=float(dpdx_f * Dh / (0.5 * rho_m * u_avg ** 2)
                                    * Re_m),
                mdot_interior_spread=float((fl.max() - fl.min()) / fl.mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpm", type=int, default=16)
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--L", type=float, default=20.0)
    ap.add_argument("--Re", default="100,300,500,1000")
    ap.add_argument("--iters", type=int, default=400000)
    ap.add_argument("--out", default="results/io_bc_validation_hiRe.json")
    a = ap.parse_args()

    print("=" * 84)
    print(" 入口放物線 + 出口定圧（Zou-He）の高 Re 検証   平行平板, Dh = 2w, "
          f"L = {a.L} w_v")
    print("=" * 84)
    print(f"{'Re':>6}{'tau_p':>9}{'umax/uavg':>11}{'f*Re':>9}{'補正後':>9}"
          f"{'err%':>8}{'質量ばらつき':>13}{'収束':>7}{'it':>9}{'s':>7}")
    rows = []
    for r_ in a.Re.split(","):
        try:
            r = run(float(r_), a.cpm, a.L, U=a.U, iters=a.iters)
        except FloatingPointError as e:
            print(f"{float(r_):>6.0f}   発散: {e}")
            rows.append(dict(Re=float(r_), diverged=True, message=str(e)))
            continue
        rows.append(r)
        print(f"{r['Re']:>6.0f}{r['tau_p']:>9.4f}{r['u_max_over_avg']:>11.5f}"
              f"{r['fRe']:>9.3f}{r['fRe_corrected']:>9.3f}"
              f"{100*(r['fRe_corrected']/96-1):>8.3f}"
              f"{r['mdot_interior_spread']:>13.2e}{str(r['converged']):>7}"
              f"{r['iters']:>9d}{r['wall_s']:>7.0f}", flush=True)

    ok = all(not r.get("diverged") and abs(r["fRe_corrected"] / 96 - 1) < 0.02
             for r in rows)
    print("\n判定: " + ("合格" if ok else "要確認"))
    print("  f*Re の離散値: NY=16 で約 95.6（NY=12: 95.34, NY=20: 95.76）")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(dict(cpm=a.cpm, U=a.U, L=a.L, rows=rows, passed=bool(ok)),
              open(a.out, "w"), indent=1)
    print(f"saved {a.out}")


if __name__ == "__main__":
    main()
