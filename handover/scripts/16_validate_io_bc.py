#!/usr/bin/env python3
"""
16_validate_io_bc.py -- 入口・出口 BC（`lbm.solve_flow_io`）の解析解検証
=======================================================================

Gamboa の単発モデルは周期化できないので入口・出口 BC が要る。
既存の `solve_flow`（非平衡外挿）は Re >= 200 で直線流路でも発散した（README §8 B-6）。
Zou-He に置き換えた `solve_flow_io` について、**Tesla 形状へ進む前に**
平行平板の解析解で検証する（CLAUDE.md の方針）。

検証項目（発達領域で測る）
  1. u_max / u_avg -> 1.5（離散理論値は NY 依存。NY=16 なら約 1.489）
  2. f_Darcy * Re = 96
  3. 質量保存（**BC 列を除いた内部**の断面流量のばらつき）
  4. **Re = 500 で発散しないこと**

注意（2026-09-06 に確認）
------------------------
- BC 列（i=0, i=nx-1）の断面流量は内部と一致しない。入口列は一様流速を課すので
  局所的に外れる。**質量保存は内部で見ること。** 内部は 1e-5 で一定だった
- LBM は弱圧縮性なので出口で rho を固定すると下流ほど rho が下がり u が増える。
  この加速ぶん（運動量流束の変化）が Δp に乗る。生の f*Re は約 +1 % 高く出るので、
  運動量補正した値も併記する
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import lbm


def run(Re, cpm, L_wv, U=0.05, iters=600000, tol=1e-8, pad=2):
    ny_core = cpm
    nx = int(round(L_wv * cpm))
    mask = np.zeros((nx, ny_core + 2 * pad), bool)
    mask[:, pad:-pad] = True
    Dh = 2.0 * cpm
    nu = lbm.nu_lattice(U, Dh, Re)

    def rep(it, d, mi, mo, dp=None):
        if it % 25000 == 0:
            print(f"    it={it:7d}  res={d:.3e}  mdot in/out = "
                  f"{mi:.6e} / {mo:.6e}  dp={dp:.6e}", flush=True)

    ux, uy, p, rho, f, info = lbm.solve_flow_io(
        mask, nu, U, iters=iters, tol=tol, report=rep)

    # --- 発達領域（下流側 65 %〜92 %）。助走 L_e ~ 0.05 Re Dh に注意 ---
    i0, i1 = int(0.65 * nx), int(0.92 * nx)
    prof = ux[i1, pad:-pad]
    umax_uavg = float(prof.max() / prof.mean())

    def pm(i):
        w = np.maximum(ux[i, pad:-pad], 1e-12)
        return float(np.average(p[i, pad:-pad], weights=w))

    im = (i0 + i1) // 2
    u_avg = float(ux[im, pad:-pad].mean())
    rho_m = float(rho[im, pad:-pad].mean())
    dpdx = (pm(i0) - pm(i1)) / (i1 - i0)
    fD = dpdx * Dh / (0.5 * rho_m * u_avg ** 2)
    Re_meas = u_avg * Dh / nu

    # 加速ぶんを差し引いた摩擦だけの勾配（beta = 1.2: 放物線分布の運動量係数）
    A = float(ny_core)
    mom = (1.2 * (rho * ux ** 2)[i0, pad:-pad].sum() / A
           - 1.2 * (rho * ux ** 2)[i1, pad:-pad].sum() / A)
    dpdx_fric = (pm(i0) - pm(i1) + mom) / (i1 - i0)
    fRe_corr = dpdx_fric * Dh / (0.5 * rho_m * u_avg ** 2) * Re_meas

    # 質量保存: BC 列を除いた内部の断面流量のばらつき
    fl = np.array([(rho * ux)[i, pad:-pad].sum() for i in range(1, nx - 1)])
    spread = float((fl.max() - fl.min()) / fl.mean())

    return dict(Re=Re, cpm=cpm, L_wv=L_wv, nx=nx, nu=nu, tau_p=info["tau_p"],
                converged=info["converged"], iters=info["iters"],
                residual=info["residual"], wall_s=info["wall_s"],
                mdot_in_col=info["mdot_in"], mdot_out_col=info["mdot_out"],
                mdot_interior=float(fl.mean()), mdot_interior_spread=spread,
                u_max_over_avg=umax_uavg, u_avg=u_avg, rho_mid=rho_m,
                Re_measured=Re_meas, fRe=float(fD * Re_meas),
                fRe_corrected=float(fRe_corr), station=[i0, i1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpm", type=int, default=16)
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--cases", default="100:30,500:120",
                    help="Re:長さ[w_v] をカンマ区切り")
    ap.add_argument("--iters", type=int, default=800000)
    ap.add_argument("--out", default="results/io_bc_validation.json")
    a = ap.parse_args()

    print("=" * 86)
    print(" 入口一様流速 + 出口定圧（Zou-He）の解析解検証   平行平板, Dh = 2w")
    print("=" * 86)
    hdr = (f"{'Re':>6}{'L[w_v]':>8}{'nx':>7}{'tau_p':>8}{'umax/uavg':>11}"
           f"{'f*Re':>9}{'補正後':>9}{'内部ばらつき':>13}{'収束':>7}{'s':>7}")
    print(hdr)
    rows = []
    for c in a.cases.split(","):
        Re, L = c.split(":")
        r = run(float(Re), a.cpm, float(L), U=a.U, iters=a.iters)
        rows.append(r)
        print(f"{r['Re']:>6.0f}{r['L_wv']:>8.0f}{r['nx']:>7d}{r['tau_p']:>8.4f}"
              f"{r['u_max_over_avg']:>11.5f}{r['fRe']:>9.3f}"
              f"{r['fRe_corrected']:>9.3f}{r['mdot_interior_spread']:>13.2e}"
              f"{str(r['converged']):>7}{r['wall_s']:>7.0f}", flush=True)

    ok = all(abs(r["fRe_corrected"] / 96 - 1) < 0.02
             and r["mdot_interior_spread"] < 5e-3 and r["converged"]
             for r in rows)
    print("\n判定: " + ("合格（補正後 f*Re が 2 % 以内、内部の質量ばらつき 0.5 % 以内）"
                        if ok else "要確認"))
    print("  離散理論値の目安: u_max/u_avg は NY=12 で 1.48276、NY=20 で 1.49377")
    print("  f*Re の離散値: NY=12 で 95.34、NY=20 で 95.76（NY=16 なら約 95.6）")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(dict(cpm=a.cpm, U=a.U, rows=rows, passed=bool(ok)),
              open(a.out, "w"), indent=1)
    print(f"saved {a.out}")


if __name__ == "__main__":
    main()
