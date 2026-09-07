"""
lbm3d.py -- 3 次元 D3Q19 / TRT 格子ボルツマン法
===============================================

なぜ 3 次元が要るか
-------------------
2 次元（`lbm.py`）は上下壁の無い「無限に深い流路」を解いている。
実物の正方形断面では

    f * Re = 56.91（正方形ダクト）  対  96（平行平板）

と全く違う。上下壁は損失を増やすだけでなく、ループへの流量配分も変える。
Gamboa の再現自体は 2D CFD が相手なので 2 次元でよいが、
**実物・Porwal（正方形断面）・Han（共役伝熱）へ進むには 3 次元が必要**。

構成は 2 次元版と同じ考え方で揃えてある。

  - TRT 衝突（magic = 1/4）。外力項は**奇数成分を omega_minus で緩和**する
    （2 次元で 33 % の圧損誤差を出したのと同じ罠。README §8 B-2）
  - 壁はハーフウェイ・バウンスバック
  - `solve_duct_periodic`: 流れ方向に周期、体積力駆動。流量一定制御。
    Δp = G * nx。検証用
  - `solve_flow_io`: 入口一様流速・出口定圧。非平衡バウンスバック
    （Zou-He の 3 次元版として広く使われる形）
  - 収束判定は Δp の安定性 + **内部の断面流量のばらつき**

メモリ
------
f は (19, nx, ny, nz) の float64。cpm = 8 の Gamboa 領域（144x133x12）で 35 MB、
cpm = 12（216x200x16）で 105 MB、cpm = 16（288x266x20）で 233 MB。
一時配列を含めて 4〜5 倍を見込むこと。
"""
from __future__ import annotations
import time

import numpy as np

# D3Q19
_C = [(0, 0, 0),
      (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
      (1, 1, 0), (-1, -1, 0), (1, -1, 0), (-1, 1, 0),
      (1, 0, 1), (-1, 0, -1), (1, 0, -1), (-1, 0, 1),
      (0, 1, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1)]
C = np.array(_C, dtype=int)
EX, EY, EZ = C[:, 0], C[:, 1], C[:, 2]
WT = np.array([1 / 3] + [1 / 18] * 6 + [1 / 36] * 12)
OPP = np.array([_C.index((-c[0], -c[1], -c[2])) for c in _C])
CS2 = 1.0 / 3.0
MAGIC = 0.25

Q = len(_C)


def tau_from_nu(nu):
    return nu / CS2 + 0.5


def nu_lattice(U, Dh_cells, Re):
    return U * Dh_cells / Re


def _feq(rho, ux, uy, uz):
    usq = ux ** 2 + uy ** 2 + uz ** 2
    out = np.empty((Q,) + rho.shape)
    for k in range(Q):
        cu = EX[k] * ux + EY[k] * uy + EZ[k] * uz
        out[k] = WT[k] * rho * (1 + cu / CS2 + cu ** 2 / (2 * CS2 ** 2)
                                - usq / (2 * CS2))
    return out


def _moments(f, mask, gx=0.0):
    rho = np.where(mask, np.maximum(f.sum(0), 1e-8), 1.0)
    ux = ((f * EX[:, None, None, None]).sum(0) + 0.5 * gx) / rho * mask
    uy = (f * EY[:, None, None, None]).sum(0) / rho * mask
    uz = (f * EZ[:, None, None, None]).sum(0) / rho * mask
    return rho, ux, uy, uz


def _bb_masks(solid):
    return [np.roll(np.roll(np.roll(solid, EX[k], 0), EY[k], 1), EZ[k], 2)
            for k in range(Q)]


def _stream(fpost, bb):
    fnew = np.empty_like(fpost)
    for k in range(Q):
        fnew[k] = np.roll(np.roll(np.roll(fpost[k], EX[k], 0), EY[k], 1),
                          EZ[k], 2)
        fnew[k][bb[k]] = fpost[OPP[k]][bb[k]]
    return fnew


def _collide(f, feq, om_p, om_m, F=None):
    fp = 0.5 * (f + f[OPP]); fm = 0.5 * (f - f[OPP])
    ep = 0.5 * (feq + feq[OPP]); em = 0.5 * (feq - feq[OPP])
    out = f - om_p * (fp - ep) - om_m * (fm - em)
    if F is not None:
        Fp = 0.5 * (F + F[OPP]); Fm = 0.5 * (F - F[OPP])
        out += (1 - 0.5 * om_p) * Fp + (1 - 0.5 * om_m) * Fm
    return out


def solve_duct_periodic(mask, nu, mdot_target, iters=200000, tol=1e-9,
                        mdot_tol=1e-6, ctrl_every=200, relax=0.5, report=None):
    """流れ方向（x）に周期、体積力駆動。流量一定制御。検証用。

    mask : (nx, ny, nz) bool。**y, z の端は固体**であること（np.roll 対策）
    戻り値: ux, uy, uz, rho, info（info["dp_lattice"] = G * nx）
    """
    nx, ny, nz = mask.shape
    if mask[:, 0].any() or mask[:, -1].any() or mask[:, :, 0].any() \
            or mask[:, :, -1].any():
        raise ValueError("y, z の端に流体セルがある。pad を使うこと")
    solid = ~mask
    tau_p = tau_from_nu(nu)
    om_p, om_m = 1.0 / tau_p, 1.0 / (0.5 + MAGIC / (tau_p - 0.5))
    bb = _bb_masks(solid)

    f = np.zeros((Q, nx, ny, nz))
    for k in range(Q):
        f[k] = WT[k]
    f *= mask[None]
    A = float(mask.sum() / nx)                    # 断面の流体セル数
    G = 12.0 * nu * (mdot_target / max(A, 1.0)) / max(A, 1.0)
    ux = np.zeros_like(mask, float)
    hist, converged = [], False
    t0 = time.time()
    for it in range(iters):
        rho, ux_n, uy, uz = _moments(f, mask, G)
        if it % ctrl_every == 0 and it > 0:
            if not np.isfinite(ux_n).all():
                raise FloatingPointError(f"diverged at it={it}")
            d = float(np.max(np.abs(ux_n - ux)))
            mdot = float((rho * ux_n)[mask].sum() / nx)
            err = mdot / mdot_target - 1.0
            hist.append((it, d, G, mdot, err))
            if report:
                report(it, d, G, err)
            dG = float("inf")
            if len(hist) >= 3:
                gs = [h[2] for h in hist[-3:]]
                dG = max(abs(g / gs[-1] - 1.0) for g in gs[:-1])
            if dG < tol and abs(err) < mdot_tol:
                ux = ux_n
                converged = True
                break
            G *= float(np.clip(1.0 - relax * err, 0.5, 2.0))
        ux = ux_n
        feq = _feq(rho, ux, uy, uz)
        F = np.empty_like(f)
        for k in range(Q):
            cu = EX[k] * ux + EY[k] * uy + EZ[k] * uz
            F[k] = WT[k] * ((EX[k] - ux) / CS2 + cu * EX[k] / CS2 ** 2) * G
        fpost = np.where(mask[None], _collide(f, feq, om_p, om_m, F), f)
        f = _stream(fpost, bb) * mask[None]

    rho, ux, uy, uz = _moments(f, mask, G)
    info = dict(G=G, converged=bool(converged), iters=it + 1,
                dp_lattice=G * nx, tau_p=tau_p, nu=nu,
                mdot=float((rho * ux)[mask].sum() / nx),
                wall_s=round(time.time() - t0, 1), history=hist)
    return ux, uy, uz, rho, info


def solve_flow_io(mask, nu, U_in, iters=400000, tol=1e-7, check_every=200,
                  ramp=5000, rho_out=1.0, dp_tol=1e-5, dp_window=10000,
                  avg_tol=1e-3, mass_tol=2e-3, report=None, ckpt=None,
                  ckpt_every=50000):
    """入口一様流速・出口定圧（x 方向）。

    境界は**非平衡バウンスバック**（Zou-He の 3 次元版として広く使われる形）。

        未知成分 i について  f_i = f_opp(i) + feq_i(rho,u) - feq_opp(i)(rho,u)

    速度 BC では rho を質量収支から、圧力 BC では ux を質量収支から決める。

    収束判定は 2 次元版と揃えてある。Δp の窓平均が動かなくなり、かつ
    **内部の断面流量のばらつきが mass_tol 未満**であること。
    """
    nx, ny, nz = mask.shape
    if mask[:, 0].any() or mask[:, -1].any() or mask[:, :, 0].any() \
            or mask[:, :, -1].any():
        raise ValueError("y, z の端に流体セルがある。pad を使うこと")
    solid = ~mask
    tau_p = tau_from_nu(nu)
    om_p, om_m = 1.0 / tau_p, 1.0 / (0.5 + MAGIC / (tau_p - 0.5))
    bb = _bb_masks(solid)
    inl, out = mask[0], mask[-1]
    if not inl.any() or not out.any():
        raise ValueError("入口列または出口列に流体セルがない")
    kin = [k for k in range(Q) if EX[k] > 0]        # 入口の未知成分
    kout = [k for k in range(Q) if EX[k] < 0]       # 出口の未知成分
    k0 = [k for k in range(Q) if EX[k] == 0]

    f = np.zeros((Q, nx, ny, nz))
    for k in range(Q):
        f[k] = WT[k]
    f *= mask[None]
    ux = np.zeros_like(mask, float)
    hist, converged, unsteady = [], False, False
    spread = float("nan")
    t0 = time.time()
    for it in range(iters):
        rho, ux_n, uy, uz = _moments(f, mask)
        if it % check_every == 0 and it > 0:
            if not np.isfinite(ux_n).all():
                raise FloatingPointError(f"diverged at it={it}")
            d = float(np.max(np.abs(ux_n - ux))) / max(U_in, 1e-12)
            p = rho * CS2
            w1 = np.maximum(ux_n[1][mask[1]], 1e-12)
            w2 = np.maximum(ux_n[-2][mask[-2]], 1e-12)
            dp_now = float(np.average(p[1][mask[1]], weights=w1)
                           - np.average(p[-2][mask[-2]], weights=w2))
            g = rho * ux_n
            fl = np.array([g[i][mask[i]].sum() for i in range(1, nx - 1, 4)])
            spread = float((fl.max() - fl.min()) / max(abs(fl.mean()), 1e-30))
            hist.append((it, d, dp_now, spread))
            if report:
                report(it, d, dp_now, spread)
            nback = max(int(dp_window // check_every), 2)
            ddp = (abs(dp_now / hist[-1 - nback][2] - 1.0)
                   if len(hist) > nback else float("inf"))
            mass_ok = spread < mass_tol
            steady = it > 2 * ramp and mass_ok and (d < tol or ddp < dp_tol)
            stat = False
            if (not steady and mass_ok and it > max(3 * ramp, 3 * dp_window)
                    and len(hist) > 2 * nback):
                m1 = np.mean([h[2] for h in hist[-nback:]])
                m0 = np.mean([h[2] for h in hist[-2 * nback:-nback]])
                stat = abs(m1 / m0 - 1.0) < avg_tol
            if steady or stat:
                ux = ux_n
                converged, unsteady = True, bool(stat and not steady)
                break
        if ckpt and it % ckpt_every == 0 and it > 0:
            np.save(ckpt, f)
        ux = ux_n
        feq = _feq(rho, ux, uy, uz)
        fpost = np.where(mask[None], _collide(f, feq, om_p, om_m), f)
        f = _stream(fpost, bb) * mask[None]

        # ---- 入口: 一様流速、非平衡バウンスバック ----
        u0 = U_in * min(1.0, (it + 1) / max(ramp, 1))
        s0 = sum(f[k, 0] for k in k0)
        sneg = sum(f[k, 0] for k in kout)
        r_in = np.where(inl, (s0 + 2.0 * sneg) / (1.0 - u0), 1.0)
        zero = np.zeros_like(r_in)
        fe = _feq(r_in, np.where(inl, u0, 0.0), zero, zero)
        for k in kin:
            f[k, 0] = np.where(inl, f[OPP[k], 0] + fe[k] - fe[OPP[k]], f[k, 0])

        # ---- 出口: 定圧、非平衡バウンスバック ----
        s0 = sum(f[k, -1] for k in k0)
        spos = sum(f[k, -1] for k in kin)
        ue = np.where(out, -1.0 + (s0 + 2.0 * spos) / rho_out, 0.0)
        fe = _feq(np.full_like(ue, rho_out), ue, np.zeros_like(ue),
                  np.zeros_like(ue))
        for k in kout:
            f[k, -1] = np.where(out, f[OPP[k], -1] + fe[k] - fe[OPP[k]],
                                f[k, -1])

    rho, ux, uy, uz = _moments(f, mask)
    p = rho * CS2
    nback = max(int(dp_window // check_every), 2)
    win = [h[2] for h in hist[-nback:]] if hist else [float("nan")]
    info = dict(converged=bool(converged), unsteady=bool(unsteady),
                iters=it + 1, tau_p=tau_p, nu=nu, U_in=U_in,
                dp_lattice_mean=float(np.mean(win)),
                dp_lattice_std=float(np.std(win)),
                dp_samples=len(win), mass_spread_final=spread,
                residual=float(hist[-1][1]) if hist else float("nan"),
                wall_s=round(time.time() - t0, 1), history=hist)
    return ux, uy, uz, rho, p, f, info


def duct_mask(nx, ny_core, nz_core, pad=1):
    """矩形ダクトのマスク（y, z に固体パディング）。"""
    m = np.zeros((nx, ny_core + 2 * pad, nz_core + 2 * pad), bool)
    m[:, pad:-pad, pad:-pad] = True
    return m


def extrude(mask2d, depth_cells, pad=1):
    """2 次元マスクを z 方向へ押し出して 3 次元にする。

    mask2d : (nx, ny) bool（`gamboa_full.rasterize` の出力）
    depth  : 流路の深さ（セル数）。z の両端に pad 個の固体を付ける
    """
    nx, ny = mask2d.shape
    m = np.zeros((nx, ny, depth_cells + 2 * pad), bool)
    m[:, :, pad:-pad] = mask2d[:, :, None]
    return m
