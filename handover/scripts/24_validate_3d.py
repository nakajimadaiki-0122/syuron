#!/usr/bin/env python3
"""
24_validate_3d.py -- 3 次元ソルバ（D3Q19 TRT）の解析解検証
==========================================================

**Tesla 形状を 3 次元で解く前に必ず通す。** 2 次元でやったのと同じ手順。

矩形ダクトの完全発達層流には解析解がある（Shah & London）。

| 縦横比 b/a | f * Re | u_max / u_avg |
|---|---|---|
| 1（正方形） | 56.91 | 2.0962 |
| 2 | 62.19 | 2.0435 |
| 4 | 72.93 | 1.9912 |
| 無限（平行平板） | 96 | 1.5 |

**正方形は平行平板の 96 に対して 56.91。** 2 次元の結果をそのまま
3 次元の実物に当てはめられない理由がこれ。

f * Re は水力直径基準: f = (dp/dx) * D_h / (rho u^2 / 2)、D_h = 4A/P。
周期境界 + 体積力なので助走区間は要らず、Δp = G * nx で厳密に取れる。
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import lbm3d as L3

# Shah & London の値（縦横比 -> f*Re, u_max/u_avg）
EXACT = {1.0: (56.91, 2.0962), 2.0: (62.19, 2.0435), 4.0: (72.93, 1.9912)}


def run(ny, nz, Re, U=0.05, nx=6, iters=400000, tol=1e-10):
    """周期ダクト（断面 ny x nz セル）を体積力で駆動する。"""
    mask = L3.duct_mask(nx, ny, nz)
    A = ny * nz
    P = 2.0 * (ny + nz)
    Dh = 4.0 * A / P                      # セル単位の水力直径
    nu = L3.nu_lattice(U, Dh, Re)
    mdot = U * A
    ux, uy, uz, rho, info = L3.solve_duct_periodic(mask, nu, mdot, iters=iters,
                                                   tol=tol)
    core = mask[nx // 2]
    prof = ux[nx // 2][core]
    u_avg = float(prof.mean())
    u_max = float(prof.max())
    dpdx = info["G"]                      # 体積力 = -dp/dx
    rho_m = float(rho[nx // 2][core].mean())
    f_darcy = dpdx * Dh / (0.5 * rho_m * u_avg ** 2)
    Re_meas = u_avg * Dh / nu
    return dict(ny=ny, nz=nz, aspect=max(ny, nz) / min(ny, nz), Dh=Dh,
                Re=Re, nu=nu, tau_p=info["tau_p"], iters=info["iters"],
                converged=info["converged"], u_avg=u_avg, u_max=u_max,
                u_max_over_avg=u_max / u_avg, fRe=float(f_darcy * Re_meas),
                wall_s=info["wall_s"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--cells", type=int, default=20, help="短辺のセル数")
    ap.add_argument("--aspects", default="1,2,4")
    ap.add_argument("--out", default="results/validate_3d.json")
    a = ap.parse_args()

    print("=" * 78)
    print(f" D3Q19 TRT の解析解検証   矩形ダクト、Re = {a.Re}、"
          f"短辺 {a.cells} セル")
    print("=" * 78)
    print(f"{'縦横比':>7}{'断面':>10}{'tau_p':>8}{'u_max/u_avg':>13}{'解析解':>9}"
          f"{'f*Re':>9}{'解析解':>9}{'誤差%':>8}{'s':>7}")
    rows = []
    for s in a.aspects.split(","):
        asp = float(s)
        ny, nz = a.cells, int(round(a.cells * asp))
        r = run(ny, nz, a.Re)
        ex = EXACT.get(asp)
        r["exact_fRe"], r["exact_umax"] = (ex if ex else (None, None))
        err = 100 * (r["fRe"] / ex[0] - 1) if ex else float("nan")
        r["fRe_err_pct"] = err
        rows.append(r)
        print(f"{asp:>7.0f}{f'{ny}x{nz}':>10}{r['tau_p']:>8.4f}"
              f"{r['u_max_over_avg']:>13.4f}{ex[1] if ex else 0:>9.4f}"
              f"{r['fRe']:>9.3f}{ex[0] if ex else 0:>9.2f}{err:>8.2f}"
              f"{r['wall_s']:>7.0f}", flush=True)

    ok = all(abs(r["fRe_err_pct"]) < 3.0 for r in rows)
    print("\n判定: " + ("合格（f*Re が解析解の 3 % 以内）" if ok else "要確認"))
    print("  離散誤差は 1/N^2 で減る。短辺 20 セルなら数 % 残る。"
          "格子を上げて誤差が下がることを確認すること")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(dict(Re=a.Re, cells=a.cells, rows=rows, passed=bool(ok)),
              open(a.out, "w"), indent=1)
    print(f"saved {a.out}")


if __name__ == "__main__":
    main()
