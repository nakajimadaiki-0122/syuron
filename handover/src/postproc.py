"""
postproc.py -- 圧力損失・ダイオード性・ループ流量比の算出
========================================================

定義を一箇所に集約している。他論文と比較するときは
**必ずここの定義が相手と一致しているか確認すること**。
7 月の解析との不一致の一因が定義差である可能性が残っている。
"""
from __future__ import annotations
import numpy as np

# 物性（水, 300 K）
NU_WATER = 8.576e-7      # m^2/s
RHO_WATER = 997.0        # kg/m^3
CP_WATER = 4180.0        # J/(kg K)
PR_WATER = 5.83


def lattice_to_physical_scale(Re: float, Dh_m: float, U_lattice: float,
                              nu_phys: float = NU_WATER) -> float:
    """格子速度 → 物理速度 [m/s] の換算係数。"""
    U_phys = Re * nu_phys / Dh_m
    return U_phys / U_lattice


def pressure_drop(p, ux, mask, i0=2, i1=None):
    """流速で重み付けした断面平均圧力の差（格子単位）。

    単純平均ではなく流量重み平均にしているのは、ループを含む断面で
    停滞領域の圧力が過大に効くのを避けるため。
    """
    nx = mask.shape[0]
    if i1 is None:
        i1 = nx - 3

    def pm(i):
        j = np.where(mask[i])[0]
        w = np.maximum(ux[i, j], 1e-12)
        return float(np.average(p[i, j], weights=w))
    return pm(i0) - pm(i1)


def to_pascal(dp_lattice: float, scale: float, rho: float = RHO_WATER) -> float:
    """格子圧力差 → Pa。p_lat は rho 単位なので rho * scale^2 を掛ける。"""
    return rho * scale ** 2 * dp_lattice


def diodicity(dp_reverse: float, dp_forward: float) -> float:
    """Di_p = Δp_逆流 / Δp_順流。

    注意: この値が「1 段あたり」か「装置全体」かは幾何依存。
    Porwal はバルブ 1 段ごとに算出しており、直線区間を含む全体値と
    直接比較してはいけない。
    """
    return dp_reverse / dp_forward


def loop_mass_fraction(ux, mask, ys, main_half_width_mm=0.5, rho=None):
    """各 x 断面で、主流路帯の外（= ループ）を通る流量の割合 phi_loop。

    main_half_width_mm : 主流路の半幅。本プロジェクトの形状では 0.5 mm
                         （主流路は y = -0.5 ... +0.5 mm）
    rho : 密度場。**必ず渡すこと。** 質量流量は rho*u_n の断面積分で評価する。
          LBM では圧力が密度に比例するため、上流ほど rho が大きい。速度だけを
          積分すると下流の流量が過大に出て、見かけ上 5-6 % の質量不整合が
          現れる（2026-08-09 に判明。実際の保存誤差は 0.02 % 以下）。
          None の場合は rho=1 とみなすが、その値を質量保存の指標に使わないこと。
    """
    nx = mask.shape[0]
    main = np.abs(ys) <= main_half_width_mm + 1e-9
    g = ux if rho is None else rho * ux          # 質量流束 rho*u_n
    mtot = np.array([g[i, mask[i]].sum() for i in range(nx)])
    mmain = np.array([g[i, mask[i] & main].sum() for i in range(nx)])
    phi = 1.0 - mmain / np.maximum(mtot, 1e-30)
    has_loop = np.array([(mask[i] & ~main).sum() > 0 for i in range(nx)])
    phi_at_loop = np.where(has_loop, phi, np.nan)

    # 参考: 密度を無視した体積流束の不整合（過去の報告値との対応確認用）
    vtot = np.array([ux[i, mask[i]].sum() for i in range(nx)])
    return dict(phi=phi,
                phi_loop_mean=float(np.nanmean(phi_at_loop)),
                phi_loop_max=float(np.nanmax(phi_at_loop)),
                mdot_in=float(mtot[2]), mdot_out=float(mtot[-3]),
                mass_imbalance_pct=float(100 * (mtot[-3] - mtot[2]) / mtot[2]),
                volflux_imbalance_pct=float(100 * (vtot[-3] - vtot[2]) / vtot[2]))


def velocity_ratio(ux, mask, ys, main_half_width_mm=0.5):
    """u_loop / u_main。7 月が報告した量と同じ定義。

    注意: これは **流速比** であって流量比ではない。ループ断面が広ければ
    流速比が小さくても流量は多い。7 月に「ループが死水域」と解釈した
    根拠がこの量で、phi_loop を見ると 21 % 流れていた。
    """
    nx = mask.shape[0]
    main = np.abs(ys) <= main_half_width_mm + 1e-9
    um = np.mean([np.abs(ux[i, mask[i] & main]).mean() for i in range(nx)])
    ul = [np.abs(ux[i, mask[i] & ~main]).mean()
          for i in range(nx) if (mask[i] & ~main).sum() > 0]
    return float(np.mean(ul) / um)
