"""
lbm.py -- 格子ボルツマン法 D2Q9 / TRT 衝突演算子
================================================

検証済み（scripts/01_validate_solver.py, 2026-08-09）
-----------------------------------------------------
平行平板・完全発達層流に対して
    f_Darcy * Re = 96      → NY=40 で 95.940（誤差 -0.062 %）
    u_max / u_avg = 1.5    → NY=40 で 1.49844
いずれも 2 次収束を確認（NY=12→20→30→40 で収束次数 2.00 / 1.99 / 2.02）。

過去に踏んだバグ
----------------
1. TRT の外力項。Guo 外力の **奇数成分は omega_minus** で緩和させる必要がある。
   BGK と同じく omega_plus を一律適用すると f*Re = 128 になり、
   **圧力損失が一律 33 % 過大**になる。速度分布の形（u_max/u_avg = 1.5）は
   正しいままなので、分布図を見ても気づけない。
2. 流体セルが配列端に接すると np.roll で壁が消える（geometry.py 冒頭を参照）。
"""
from __future__ import annotations
import numpy as np

# D2Q9 格子
EX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
EY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
WT = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])
CS2 = 1.0 / 3.0
MAGIC = 0.25          # TRT の magic parameter。1/4 で壁位置が厳密に中間


def nu_lattice(U: float, Dh_cells: float, Re: float) -> float:
    """目標 Re から格子動粘性係数を決める。Dh_cells = 水力直径のセル数。"""
    return U * Dh_cells / Re


def tau_from_nu(nu: float) -> float:
    return nu / CS2 + 0.5


def solve_flow_periodic(mask, nu, mdot_target, iters=200000, tol=1e-9,
                        mdot_tol=1e-6, G0=None, report=None, f0=None,
                        ctrl_every=200, relax=0.5):
    """流れ方向に周期境界、体積力駆動で単位セルを解く。

    `02_run_flow.py` の入口・出口境界条件は Re >= 200 で発散する
    （2026-08-09 に直線流路でも再現。README §8 B-6）。本関数は
    `scripts/01_validate_solver.py` で解析解に対し検証済みの機構
    （周期境界 + Guo 外力の TRT 分割）をそのまま 2 次元マスクへ拡張したもので、
    入口・出口処理を一切持たない。

    駆動は**流量一定**で行う。体積力 G を制御して質量流量を mdot_target に
    合わせ、収束後の G から単位セルあたりの圧力損失を得る。

        dp/dx = -G   →   Δp（単位セル）= G * nx

    順流・逆流を同一流量で比較するので、Di_p = G_逆 / G_順 になる。

    Parameters
    ----------
    mask         : (nx, ny) bool, True = 流体。**y 方向の両端は固体であること**
                   （np.roll による周期化を防ぐ。geometry.rasterize の pad を使う）
                   x 方向は周期境界なので端に流体があってよい。
    nu           : 格子動粘性係数
    mdot_target  : 目標質量流量（格子単位、断面あたり sum(rho*ux)）
    tol          : 収束判定。直近 3 回の検査で体積力 G の相対変化がこれ未満、
                   かつ質量流量誤差が mdot_tol 未満なら収束とみなす。
                   **速度残差ではなく G を見る**のは、Δp = G * nx であり
                   Di_p が依存するのは G だからである。速度残差は微小な
                   音響モードのせいで 1e-8 付近で下げ止まることがある。
    G0           : 体積力の初期値。None なら平行平板の解析解から見積もる
    ctrl_every   : 体積力を更新する間隔
    relax        : 体積力更新の緩和係数

    Returns
    -------
    ux, uy, p, rho, f, info
    """
    nx, ny = mask.shape
    if mask[:, 0].any() or mask[:, -1].any():
        raise ValueError("y 方向の端に流体セルがある。np.roll で上下が周期化する。"
                         "geometry.rasterize の pad を使うこと。")
    solid = ~mask
    tau_p = tau_from_nu(nu)
    om_p = 1.0 / tau_p
    tau_m = 0.5 + MAGIC / (tau_p - 0.5)
    om_m = 1.0 / tau_m

    if G0 is None:
        # 平行平板 u_avg = G H^2 / (12 nu) からの粗い見積もり
        H = float(mask.sum() / nx)
        G0 = 12.0 * nu * (mdot_target / max(H, 1.0)) / max(H, 1.0) ** 2
    G = float(G0)

    if f0 is not None:
        f = f0.copy()
    else:
        f = np.zeros((9, nx, ny))
        for k in range(9):
            f[k] = WT[k]
        f *= mask[None, :, :]

    bb = [np.roll(np.roll(solid, EX[k], 0), EY[k], 1) for k in range(9)]

    ux = np.zeros((nx, ny))
    uy = np.zeros((nx, ny))
    converged = False
    hist = []
    for it in range(iters):
        rho = np.where(mask, np.maximum(f.sum(0), 1e-8), 1.0)
        # Guo: u = (sum f e + F/2) / rho
        ux_n = ((f * EX[:, None, None]).sum(0) + 0.5 * G) / rho * mask
        uy_n = (f * EY[:, None, None]).sum(0) / rho * mask

        if it % ctrl_every == 0 and it > 0:
            if not np.isfinite(ux_n).all():
                raise FloatingPointError(f"diverged at it={it}")
            d = float(np.max(np.abs(ux_n - ux)))
            mdot = float((rho * ux_n)[mask].sum() / nx)
            err = mdot / mdot_target - 1.0
            hist.append((it, d, G, mdot, err))
            if report:
                report(it, d, G, err)
            if len(hist) >= 3:
                _Gs = [h[2] for h in hist[-3:]]
                dG_now = max(abs(g / _Gs[-1] - 1.0) for g in _Gs[:-1])
            else:
                dG_now = float("inf")
            # 収束判定は **G の安定性** で行う。Δp = G * nx なので、Di_p が
            # 依存するのは G であって速度場の残差ではない。
            # 速度残差 d は微小な音響モードのせいで 1e-8 程度で下げ止まり、
            # 絶対値基準では到達できないことがある（2026-08-09、cpm=24 で確認）。
            dG = dG_now
            if dG < tol and abs(err) < mdot_tol:
                ux, uy = ux_n, uy_n
                converged = True
                break
            # 流量一定制御（G と mdot はほぼ線形なので比例更新でよい）
            G *= float(np.clip(1.0 - relax * err, 0.5, 2.0))
        ux, uy = ux_n, uy_n

        usq = ux ** 2 + uy ** 2
        feq = np.empty_like(f)
        Fraw = np.empty_like(f)
        for k in range(9):
            cu = EX[k] * ux + EY[k] * uy
            feq[k] = WT[k] * rho * (1 + cu / CS2 + cu ** 2 / (2 * CS2 ** 2)
                                    - usq / (2 * CS2))
            Fraw[k] = WT[k] * ((EX[k] - ux) / CS2
                               + cu * EX[k] / CS2 ** 2) * G

        fp = 0.5 * (f + f[OPP]); fm = 0.5 * (f - f[OPP])
        ep = 0.5 * (feq + feq[OPP]); em = 0.5 * (feq - feq[OPP])
        # 外力項も TRT 分割する。奇数成分を om_p で緩和すると
        # f*Re が 33 % 過大になる（README §8 B-2）
        Fp = 0.5 * (Fraw + Fraw[OPP]); Fm = 0.5 * (Fraw - Fraw[OPP])
        Fk = (1 - 0.5 * om_p) * Fp + (1 - 0.5 * om_m) * Fm

        fpost = f - om_p * (fp - ep) - om_m * (fm - em) + Fk
        fpost = np.where(mask[None], fpost, f)

        fnew = np.empty_like(fpost)
        for k in range(9):
            fnew[k] = np.roll(np.roll(fpost[k], EX[k], 0), EY[k], 1)
            fnew[k][bb[k]] = fpost[OPP[k]][bb[k]]
        f = fnew * mask[None]

    rho = np.where(mask, np.maximum(f.sum(0), 1e-8), 1.0)
    p = rho * CS2
    mdot = float((rho * ux)[mask].sum() / nx)
    dG_final = float("nan")
    if len(hist) >= 3:
        _Gs = [h[2] for h in hist[-3:]]
        dG_final = max(abs(g / _Gs[-1] - 1.0) for g in _Gs[:-1])
    info = dict(G=G, converged=bool(converged), iters=it + 1,
                mdot=mdot, mdot_err=mdot / mdot_target - 1.0,
                dG_final=dG_final,
                u_residual=float(hist[-1][1]) if hist else float("nan"),
                dp_lattice=G * nx, tau_p=tau_p, history=hist)
    return ux, uy, p, rho, f, info


def solve_flow(mask, nu, U, iters=40000, tol=2e-6, report=None, f0=None):
    """定常流れを解く。

    境界条件
      入口 (i=0)   : 放物線速度分布、非平衡外挿
      出口 (i=nx-1): 密度 1.0 固定、速度は外挿
      壁            : ハーフウェイ・バウンスバック

    Parameters
    ----------
    mask  : (nx, ny) bool, True = 流体（geometry.rasterize の出力）
    nu    : 格子動粘性係数
    U     : 入口の断面平均速度（格子単位）。Ma = U/sqrt(1/3) を 0.1 以下に保つ
    tol   : 収束判定。2000 反復あたりの max|Δux| < tol * U で停止
    f0    : 分布関数の初期値。中断・再開に使う

    Returns
    -------
    ux, uy, p, rho, f
    """
    nx, ny = mask.shape
    solid = ~mask
    tau_p = tau_from_nu(nu)
    om_p = 1.0 / tau_p
    tau_m = 0.5 + MAGIC / (tau_p - 0.5)
    om_m = 1.0 / tau_m

    # 入口の放物線分布（開口部にのみ与える）
    j = np.where(mask[0])[0]
    if len(j) == 0:
        raise ValueError("inlet column contains no fluid cell")
    j0, j1 = j.min(), j.max() + 1
    h = j1 - j0
    yy = (np.arange(j0, j1) - j0 + 0.5) / h
    uin = np.zeros(ny)
    uin[j0:j1] = 6.0 * U * yy * (1.0 - yy)

    if f0 is not None:
        f = f0.copy()
    else:
        f = np.zeros((9, nx, ny))
        for k in range(9):
            f[k] = WT[k]
        f *= mask[None, :, :]

    # バウンスバック元セルを事前計算（毎反復の np.roll を避ける）
    bb = [np.roll(np.roll(solid, EX[k], 0), EY[k], 1) for k in range(9)]

    ux = np.zeros((nx, ny))
    uy = np.zeros((nx, ny))
    converged = False
    for it in range(iters):
        rho = np.where(mask, np.maximum(f.sum(0), 1e-8), 1.0)
        ux_n = (f * EX[:, None, None]).sum(0) / rho * mask
        uy_n = (f * EY[:, None, None]).sum(0) / rho * mask

        if it % 2000 == 0 and it > 0:
            d = float(np.max(np.abs(ux_n - ux)))
            if report:
                report(it, d)
            if d < tol * U:
                ux, uy = ux_n, uy_n
                converged = True
                break
        ux, uy = ux_n, uy_n

        usq = ux ** 2 + uy ** 2
        feq = np.empty_like(f)
        for k in range(9):
            cu = EX[k] * ux + EY[k] * uy
            feq[k] = WT[k] * rho * (1 + cu / CS2 + cu ** 2 / (2 * CS2 ** 2)
                                    - usq / (2 * CS2))

        fp = 0.5 * (f + f[OPP]); fm = 0.5 * (f - f[OPP])
        ep = 0.5 * (feq + feq[OPP]); em = 0.5 * (feq - feq[OPP])
        fpost = f - om_p * (fp - ep) - om_m * (fm - em)
        fpost = np.where(mask[None], fpost, f)

        fnew = np.empty_like(fpost)
        for k in range(9):
            fnew[k] = np.roll(np.roll(fpost[k], EX[k], 0), EY[k], 1)
            fnew[k][bb[k]] = fpost[OPP[k]][bb[k]]

        # 入口
        rin = rho[1]
        for k in range(9):
            cu = EX[k] * uin
            feq_in = WT[k] * rin * (1 + cu / CS2 + cu ** 2 / (2 * CS2 ** 2)
                                    - uin ** 2 / (2 * CS2))
            cu2 = EX[k] * ux[1] + EY[k] * uy[1]
            feq_1 = WT[k] * rho[1] * (1 + cu2 / CS2 + cu2 ** 2 / (2 * CS2 ** 2)
                                      - (ux[1] ** 2 + uy[1] ** 2) / (2 * CS2))
            fnew[k, 0] = np.where(mask[0], feq_in + (f[k, 1] - feq_1), fnew[k, 0])
        # 出口
        uo, vo = ux[-2], uy[-2]
        for k in range(9):
            cu = EX[k] * uo + EY[k] * vo
            feq_o = WT[k] * (1 + cu / CS2 + cu ** 2 / (2 * CS2 ** 2)
                             - (uo ** 2 + vo ** 2) / (2 * CS2))
            cu2 = EX[k] * ux[-2] + EY[k] * uy[-2]
            feq_2 = WT[k] * rho[-2] * (1 + cu2 / CS2 + cu2 ** 2 / (2 * CS2 ** 2)
                                       - (ux[-2] ** 2 + uy[-2] ** 2) / (2 * CS2))
            fnew[k, -1] = np.where(mask[-1], feq_o + (f[k, -2] - feq_2), fnew[k, -1])

        f = fnew * mask[None]

    rho = np.where(mask, np.maximum(f.sum(0), 1e-8), 1.0)
    p = rho * CS2
    if not converged and report:
        report(-1, float("nan"))
    return ux, uy, p, rho, f
