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


def _regularized(fcol, j, rho, ux, uy, unknown):
    """正則化境界（Latt et al. 2008）で境界節点の 9 成分を再構成する。

    Zou-He は tau が 0.5 に近いと不安定になる（2026-09-06 に確認:
    tau_p = 0.5096 の Re = 500 では、放物線入口の完全発達流でも
    it = 11500 で発散した）。正則化 BC は非平衡部を応力テンソルへ射影して
    捨てるので、低 tau でも安定である。

    手順
      1. 既知成分から f_neq = f - f_eq(rho, u) を作る
      2. 未知成分の f_neq は**対向成分から取る**（非平衡部のバウンスバック）
      3. Pi_ab = sum_i c_ia c_ib f_neq_i
      4. f_i = f_eq_i + w_i/(2 cs^4) * (c_ia c_ib - cs^2 delta_ab) Pi_ab

    fcol     : (9, n) 境界列の分布関数（この場で書き換える）
    unknown  : 未知成分の添字。西面なら [1,5,8]、東面なら [3,6,7]
    """
    usq = ux ** 2 + uy ** 2
    feq = np.empty((9, len(j)))
    for k in range(9):
        cu = EX[k] * ux + EY[k] * uy
        feq[k] = WT[k] * rho * (1 + cu / CS2 + cu ** 2 / (2 * CS2 ** 2)
                                - usq / (2 * CS2))
    fneq = fcol[:, j] - feq
    for k in unknown:                      # 非平衡部のバウンスバック
        fneq[k] = fneq[OPP[k]]
    Pxx = (fneq * (EX ** 2)[:, None]).sum(0)
    Pyy = (fneq * (EY ** 2)[:, None]).sum(0)
    Pxy = (fneq * (EX * EY)[:, None]).sum(0)
    for k in range(9):
        Q = ((EX[k] ** 2 - CS2) * Pxx + 2.0 * EX[k] * EY[k] * Pxy
             + (EY[k] ** 2 - CS2) * Pyy)
        fcol[k, j] = feq[k] + WT[k] / (2.0 * CS2 ** 2) * Q


def solve_flow_io(mask, nu, U_in, iters=400000, tol=1e-7, check_every=500,
                  ramp=5000, report=None, f0=None, rho_out=1.0,
                  ckpt=None, ckpt_every=50000, dp_tol=1e-5, bc="zouhe",
                  dp_window=10000, avg_tol=1e-3, avg_burn=0.3,
                  mass_tol=2e-3):
    """入口一様流速・出口定圧で解く（Gamboa の単発モデル用）。

    既存の `solve_flow` との違い
    ---------------------------
    1. 入口分布を**一様流速**にできる（Gamboa 2.2 節の条件）
    2. 入口・出口とも **Zou-He** に置き換えた。`solve_flow` の非平衡外挿は
       Re >= 200 で発散する（README §8 B-6）
    3. 入口速度を `ramp` 反復かけて 0 から立ち上げる（初期の音響衝撃を避ける）
    4. チェックポイントを書ける（長時間計算の中断・再開用）

    前提
    ----
    - i = 0 の列が入口面、i = nx-1 の列が出口面。どちらも**平らな鉛直面**であること
      （プレナムの直線壁）。傾いた面には使えない
    - y 方向の両端は固体（np.roll による周期化を防ぐ。README §8 B-1）

    Parameters
    ----------
    mask   : (nx, ny) bool
    nu     : 格子動粘性係数
    U_in   : 入口面での流速（格子単位）。スカラなら一様、入口流体セル数と同じ長さの
             配列なら任意分布（放物線など）を課せる。**流路の平均流速ではない**。
             プレナム面が流路より広ければ U_in = U_channel * (流路幅 / 入口面の高さ)
    tol    : 収束判定。check_every 反復あたりの max|Δux| / U_in がこれ未満
    dp_tol : Δp の相対変化が **dp_window 反復のあいだ** これ未満なら収束とみなす。
             直近 3 点だけで見ると一時的な平坦部で止まるので窓を取る
    avg_tol: 非定常な流れ用の判定。Δp が振動して定常判定に入らない場合、
             **直近 dp_window の平均**と**その 1 つ前の dp_window の平均**が
             この相対差以内なら「統計的に定常」とみなし、平均値を採用する。
             Re が高いと流れ自体が非定常になる（2026-09-07、Re=500 で確認）
    mass_tol: 内部（BC 列を除く）の断面流量のばらつきがこれ未満でなければ
             収束と判定しない。**Δp が落ち着いても質量が溜まり込んでいることがある**
    avg_burn: 平均判定を始めるまでの助走（iters に対する割合ではなく、
             ramp*3 と dp_window*3 の大きい方を下限とする）
    rho_out: 出口の密度（= 圧力）。1.0 が「ゼロゲージ圧」
    bc     : "zouhe"（既定）| "regularized"。
             **正則化 BC は本ソルバ（TRT, Lambda=1/4）では逆に不安定**
             （2026-09-06 実測: Re=300, cpm=16 で Zou-He は安定、正則化は it=7000 で発散）。
             TRT では非平衡部の奇数成分が omega_minus でゆっくりしか緩和せず、
             正則化が応力テンソルへ射影する際にその成分を捨ててしまうためと考えられる

    Returns
    -------
    ux, uy, p, rho, f, info
    """
    import os, time
    nx, ny = mask.shape
    if mask[:, 0].any() or mask[:, -1].any():
        raise ValueError("y 方向の端に流体セルがある。pad を使うこと。")
    solid = ~mask
    tau_p = tau_from_nu(nu)
    om_p = 1.0 / tau_p
    tau_m = 0.5 + MAGIC / (tau_p - 0.5)
    om_m = 1.0 / tau_m

    U_in_arr = np.atleast_1d(np.asarray(U_in, float))
    U_ref = float(np.mean(U_in_arr))
    jin = np.where(mask[0])[0]
    jout = np.where(mask[-1])[0]
    if len(jin) == 0 or len(jout) == 0:
        raise ValueError("入口列または出口列に流体セルがない")
    # 角のセル（上下どちらかが固体）は Zou-He が成立しないので feq で埋める
    if U_in_arr.size not in (1, len(jin)):
        raise ValueError("U_in は定数か入口流体セル数と同じ長さの配列")
    u_prof = (np.full(len(jin), U_in_arr[0]) if U_in_arr.size == 1
              else U_in_arr.copy())
    in_core = jin[(mask[0][jin - 1]) & (mask[0][np.minimum(jin + 1, ny - 1)])]
    core_sel = np.isin(jin, in_core)
    out_core = jout[(mask[-1][jout - 1]) & (mask[-1][np.minimum(jout + 1, ny - 1)])]

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
    unsteady = False
    hist = []
    t0 = time.time()
    it = 0
    for it in range(iters):
        rho = np.where(mask, np.maximum(f.sum(0), 1e-8), 1.0)
        ux_n = (f * EX[:, None, None]).sum(0) / rho * mask
        uy_n = (f * EY[:, None, None]).sum(0) / rho * mask

        if it % check_every == 0 and it > 0:
            if not np.isfinite(ux_n).all():
                raise FloatingPointError(f"diverged at it={it}")
            d = float(np.max(np.abs(ux_n - ux))) / max(U_ref, 1e-12)
            mi = (rho * ux_n)[0, jin].sum()
            mo = (rho * ux_n)[-1, jout].sum()
            # Δp（入口面と出口面の流量重み平均圧力の差）。収束判定はこれで行う。
            # 周期版が体積力 G の安定性で判定するのと同じ考え方（README §8 B-7）。
            pp_ = rho * CS2
            w1 = np.maximum(ux_n[1, mask[1]], 1e-12)
            w2 = np.maximum(ux_n[-2, mask[-2]], 1e-12)
            dp_now = float(np.average(pp_[1, mask[1]], weights=w1)
                           - np.average(pp_[-2, mask[-2]], weights=w2))
            hist.append((it, d, float(mi), float(mo), dp_now))
            if report:
                report(it, d, float(mi), float(mo), dp_now)
            # **窓付き**で判定する。直近 3 点だけで見ると、Δp が一時的に
            # 平らになったところで止まってしまう（2026-09-06 に実測。
            # Re=100 は it=20000 -> 30000 で Δp がまだ 3.3 % 動いていた）。
            ddp = float("inf")
            nback = max(int(dp_window // check_every), 2)
            if len(hist) > nback:
                ddp = abs(dp_now / hist[-1 - nback][4] - 1.0)
            # 内部の断面流量が揃っていない = まだ質量が溜まり込んでいる。
            # これを満たさないと Δp の平均が動かなくても「収束」とは言えない
            # （2026-09-07、Re=500 で内部ばらつき 20 % のまま統計判定が通った）
            g = rho * ux_n
            fl = np.array([g[i, mask[i]].sum() for i in range(1, nx - 1, 4)])
            spread = float((fl.max() - fl.min()) / max(abs(fl.mean()), 1e-30))
            mass_ok = spread < mass_tol
            steady = it > ramp * 2 and mass_ok and (d < tol or ddp < dp_tol)
            # 非定常でも「時間平均が動かなくなった」ら止める
            stat = False
            burn = max(3 * ramp, 3 * dp_window)
            if not steady and mass_ok and it > burn and len(hist) > 2 * nback:
                m1 = np.mean([h[4] for h in hist[-nback:]])
                m0 = np.mean([h[4] for h in hist[-2 * nback:-nback]])
                stat = abs(m1 / m0 - 1.0) < avg_tol
            if steady or stat:
                ux, uy = ux_n, uy_n
                converged = True
                unsteady = bool(stat and not steady)
                break
        if ckpt and it % ckpt_every == 0 and it > 0:
            np.save(ckpt, f)
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
        f = fnew * mask[None]

        # ---- 入口: 速度 BC（西面、uy = 0）----
        ramp_f = min(1.0, (it + 1) / max(ramp, 1))
        u0 = u_prof * ramp_f
        j = jin
        f0_, f2, f4 = f[0, 0, j], f[2, 0, j], f[4, 0, j]
        f3, f6, f7 = f[3, 0, j], f[6, 0, j], f[7, 0, j]
        r_in = (f0_ + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - u0)
        if bc == "zouhe":
            f[1, 0, j] = f3 + (2.0 / 3.0) * r_in * u0
            f[5, 0, j] = f7 - 0.5 * (f2 - f4) + (1.0 / 6.0) * r_in * u0
            f[8, 0, j] = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * r_in * u0
        else:
            _regularized(f[:, 0, :], j, r_in, u0, np.zeros_like(u0), (1, 5, 8))

        # ---- 出口: 圧力 BC（東面、uy = 0）----
        j = jout
        f0_, f2, f4 = f[0, -1, j], f[2, -1, j], f[4, -1, j]
        f1, f5, f8 = f[1, -1, j], f[5, -1, j], f[8, -1, j]
        ue = -1.0 + (f0_ + f2 + f4 + 2.0 * (f1 + f5 + f8)) / rho_out
        if bc == "zouhe":
            f[3, -1, j] = f1 - (2.0 / 3.0) * rho_out * ue
            f[7, -1, j] = f5 + 0.5 * (f2 - f4) - (1.0 / 6.0) * rho_out * ue
            f[6, -1, j] = f8 - 0.5 * (f2 - f4) - (1.0 / 6.0) * rho_out * ue
        else:
            _regularized(f[:, -1, :], j, np.full(len(j), rho_out), ue,
                         np.zeros_like(ue), (3, 6, 7))

    rho = np.where(mask, np.maximum(f.sum(0), 1e-8), 1.0)
    ux = (f * EX[:, None, None]).sum(0) / rho * mask
    uy = (f * EY[:, None, None]).sum(0) / rho * mask
    p = rho * CS2
    mdot_in = float((rho * ux)[0, jin].sum())
    mdot_out = float((rho * ux)[-1, jout].sum())
    ddp_final = float("nan")
    nback = max(int(dp_window // check_every), 2)
    if len(hist) > nback:
        ddp_final = abs(hist[-1][4] / hist[-1 - nback][4] - 1.0)
    # 直近 dp_window ぶんの Δp の統計（非定常なときはこちらを使う）
    nback = max(int(dp_window // check_every), 2)
    win = [h[4] for h in hist[-nback:]] if hist else [float("nan")]
    info = dict(dp_lattice_mean=float(np.mean(win)),
                mass_spread_final=float(spread) if hist else float("nan"),
                dp_lattice_std=float(np.std(win)),
                dp_lattice_min=float(np.min(win)),
                dp_lattice_max=float(np.max(win)),
                dp_samples=len(win), unsteady=bool(unsteady),
                converged=bool(converged), iters=it + 1, tau_p=tau_p, nu=nu,
                U_in=U_ref, mdot_in=mdot_in, mdot_out=mdot_out,
                mdot_err=(mdot_out - mdot_in) / max(abs(mdot_in), 1e-30),
                residual=float(hist[-1][1]) if hist else float("nan"),
                ddp_final=ddp_final,
                wall_s=round(time.time() - t0, 1), history=hist)
    return ux, uy, p, rho, f, info
