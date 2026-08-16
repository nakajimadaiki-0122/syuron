"""
Phase 0A: Solver validation against analytic solutions (2D parallel plates)
==========================================================================
Targets:
  Hydrodynamic : u_max/u_avg = 1.5 exactly
                 f_Darcy * Re = 96   (Dh = 2H)
  Thermal      : Nu = 7.541  (both walls isothermal, T condition)
                 Nu = 8.235  (both walls uniform heat flux, H1 condition)

Flow  : LBM D2Q9, TRT collision, half-way bounce-back, body-force driven,
        periodic streamwise -> exact fully-developed Poiseuille
Heat  : finite-volume steady advection-diffusion, 2nd-order central diffusion,
        upwind advection, scipy.sparse direct solve
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ----------------------------------------------------------------------
# 1. LBM D2Q9 TRT  -- fully developed channel flow
# ----------------------------------------------------------------------
def lbm_poiseuille(NY, tau_plus=0.8, Lambda=0.25, U_target=0.02,
                   nx=8, iters=60000, tol=1e-12):
    """
    NY  : number of FLUID cells across the channel (wall at NY+0.5 offsets)
    tau_plus : symmetric relaxation time  -> nu = cs^2 (tau_plus - 0.5)
    Lambda   : magic parameter, 1/4 puts the bounce-back wall exactly midway
    """
    cs2 = 1.0 / 3.0
    nu = cs2 * (tau_plus - 0.5)
    tau_minus = 0.5 + Lambda / (tau_plus - 0.5)
    omega_p, omega_m = 1.0 / tau_plus, 1.0 / tau_minus

    ex = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
    ey = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
    w = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
    opp = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

    H = float(NY)                      # channel height in lattice units
    # body force chosen so that analytic u_avg = U_target
    # u_avg = G*H^2/(12 nu)  (G = force per unit volume / rho, rho=1)
    G = 12.0 * nu * U_target / H**2

    f = np.zeros((9, nx, NY))
    for k in range(9):
        f[k] = w[k]

    for it in range(iters):
        rho = f.sum(axis=0)
        ux = (f * ex[:, None, None]).sum(axis=0) / rho + 0.5 * G / rho
        uy = (f * ey[:, None, None]).sum(axis=0) / rho

        usq = ux**2 + uy**2
        feq = np.empty_like(f)
        for k in range(9):
            cu = ex[k] * ux + ey[k] * uy
            feq[k] = w[k] * rho * (1 + cu/cs2 + cu**2/(2*cs2**2) - usq/(2*cs2))

        # TRT split
        fp = 0.5 * (f + f[opp])
        fm = 0.5 * (f - f[opp])
        fep = 0.5 * (feq + feq[opp])
        fem = 0.5 * (feq - feq[opp])

        # Guo forcing term, TRT-consistent split
        # (even part relaxes with omega_plus, odd part with omega_minus)
        Fraw = np.empty_like(f)
        for k in range(9):
            cu = ex[k] * ux + ey[k] * uy
            Fraw[k] = w[k] * ((ex[k] - ux) / cs2 + cu * ex[k] / cs2**2) * G
        Fp = 0.5 * (Fraw + Fraw[opp])
        Fm = 0.5 * (Fraw - Fraw[opp])
        Fk = (1 - 0.5*omega_p) * Fp + (1 - 0.5*omega_m) * Fm

        fpost = f - omega_p * (fp - fep) - omega_m * (fm - fem) + Fk

        # stream (periodic in x, bounce-back at y walls)
        fnew = np.empty_like(fpost)
        for k in range(9):
            fnew[k] = np.roll(np.roll(fpost[k], ex[k], axis=0), ey[k], axis=1)

        # half-way bounce-back on bottom (j=0) and top (j=NY-1)
        for k in range(9):
            if ey[k] == 1:
                fnew[k, :, 0] = fpost[opp[k], :, 0]
            if ey[k] == -1:
                fnew[k, :, -1] = fpost[opp[k], :, -1]

        if it % 500 == 0 and it > 0:
            du = np.max(np.abs(ux - ux_old))
            if du < tol:
                break
        ux_old = ux.copy()
        f = fnew

    rho = f.sum(axis=0)
    ux = (f * ex[:, None, None]).sum(axis=0) / rho + 0.5 * G / rho
    prof = ux.mean(axis=0)             # streamwise-averaged profile
    return prof, nu, G, H


def check_hydro(NY):
    prof, nu, G, H = lbm_poiseuille(NY)
    u_avg = prof.mean()
    u_max = prof.max()
    # analytic parabola sampled at the same cell centres
    y = (np.arange(NY) + 0.5)          # cell centres, wall at 0 and H
    u_an = G / (2*nu) * y * (H - y)
    u_max_an_cont = G * H**2 / (8*nu)  # true continuum peak

    Dh = 2.0 * H
    Re = u_avg * Dh / nu
    # force balance: dp/dx = -G  (rho = 1)
    f_darcy = G * Dh / (0.5 * u_avg**2)
    return dict(NY=NY, u_avg=u_avg, ratio=u_max_an_cont/u_avg if False else
                (u_max/u_avg), Re=Re, fRe=f_darcy*Re,
                L2=np.sqrt(np.mean((prof-u_an)**2))/u_avg)


# ----------------------------------------------------------------------
# 2. FV thermal solver -- developing then fully developed
# ----------------------------------------------------------------------
def fv_thermal(NY, NX, Pr, bc='T', u_avg=1.0, nu=1.0, Re=100.0):
    """
    2D steady advection-diffusion in a parallel-plate channel.
    bc = 'T'  : both walls isothermal  Tw = 0, inlet T = 1
    bc = 'H1' : both walls uniform heat flux q'' (into fluid)
    Returns Nu(x) using Dh = 2H.
    """
    H = float(NY)
    dy = 1.0
    dx = 1.0
    y = (np.arange(NY) + 0.5)
    # exact parabolic profile with prescribed u_avg
    u = 6.0 * u_avg * (y/H) * (1 - y/H)
    alpha = nu / Pr                     # thermal diffusivity

    N = NX * NY
    idx = lambda i, j: i*NY + j
    rows, cols, vals = [], [], []
    b = np.zeros(N)

    qflux = 1.0                          # for H1
    Tin = 1.0                            # for T

    for i in range(NX):
        for j in range(NY):
            p = idx(i, j)
            if i == 0:                   # inlet: Dirichlet
                rows.append(p); cols.append(p); vals.append(1.0)
                b[p] = Tin if bc == 'T' else 0.0
                continue

            ap = 0.0
            # --- streamwise advection: first-order upwind (u>0) ---
            ap += u[j] / dx
            rows.append(p); cols.append(idx(i-1, j)); vals.append(-u[j]/dx)

            # --- streamwise diffusion (outflow: zero-gradient at i=NX-1) ---
            if i < NX-1:
                ap += alpha/dx**2
                rows.append(p); cols.append(idx(i+1, j)); vals.append(-alpha/dx**2)
            ap += alpha/dx**2
            rows.append(p); cols.append(idx(i-1, j)); vals.append(-alpha/dx**2)

            # --- transverse diffusion ---
            # south face
            if j == 0:
                if bc == 'T':
                    ap += 2*alpha/dy**2          # wall at Tw=0, half-cell
                else:
                    b[p] += qflux*alpha/dy       # prescribed flux
            else:
                ap += alpha/dy**2
                rows.append(p); cols.append(idx(i, j-1)); vals.append(-alpha/dy**2)
            # north face
            if j == NY-1:
                if bc == 'T':
                    ap += 2*alpha/dy**2
                else:
                    b[p] += qflux*alpha/dy
            else:
                ap += alpha/dy**2
                rows.append(p); cols.append(idx(i, j+1)); vals.append(-alpha/dy**2)

            rows.append(p); cols.append(p); vals.append(ap)

    A = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    T = spla.spsolve(A, b).reshape(NX, NY)

    # ---- Nusselt number ----
    Dh = 2*H
    Nu = np.zeros(NX)
    for i in range(NX):
        Tb = np.sum(u*T[i]) / np.sum(u)          # bulk (velocity-weighted)
        if bc == 'T':
            Tw = 0.0
            # wall gradient from half-cell
            q_s = alpha*(T[i, 0] - Tw)/(0.5*dy)
            q_n = alpha*(T[i, -1] - Tw)/(0.5*dy)
            q = 0.5*(q_s + q_n)
            Nu[i] = q*Dh/(alpha*(Tb - Tw)) if abs(Tb-Tw) > 1e-14 else np.nan
        else:
            Tw_s = T[i, 0] + qflux*0.5*dy
            Tw_n = T[i, -1] + qflux*0.5*dy
            Tw = 0.5*(Tw_s + Tw_n)
            # qflux is prescribed as dT/dy at the wall, so Nu = (dT/dy)*Dh/(Tw-Tb)
            Nu[i] = qflux*Dh/(Tw - Tb) if abs(Tw-Tb) > 1e-14 else np.nan
    return Nu, T


if __name__ == "__main__":
    print("="*74)
    print(" PHASE 0A  SOLVER VALIDATION  (2D parallel plates, Dh = 2H)")
    print("="*74)

    print("\n--- (1) Hydrodynamics: LBM D2Q9 TRT, fully developed ---")
    print(f"{'NY':>5} {'u_max/u_avg':>12} {'Re':>9} {'f*Re':>10} {'err%':>8} {'L2(prof)':>10}")
    for NY in [12, 20, 30, 40]:
        r = check_hydro(NY)
        err = (r['fRe']-96.0)/96.0*100
        print(f"{r['NY']:>5} {r['ratio']:>12.5f} {r['Re']:>9.3f} "
              f"{r['fRe']:>10.4f} {err:>8.3f} {r['L2']:>10.2e}")
    print("  target: u_max/u_avg = 1.5 (continuum),  f*Re = 96")

    print("\n--- (2) Thermal: FV, isothermal walls (T condition) ---")
    print(f"{'NY':>5} {'Nu_fd':>10} {'target':>9} {'err%':>8}")
    for NY in [12, 20, 30, 40]:
        NX = 60*NY
        Nu, _ = fv_thermal(NY, NX, Pr=5.83, bc='T', u_avg=0.02,
                           nu=0.1, Re=None)
        Nu_fd = Nu[int(0.9*NX)]
        print(f"{NY:>5} {Nu_fd:>10.4f} {7.541:>9.3f} {(Nu_fd-7.541)/7.541*100:>8.3f}")

    print("\n--- (3) Thermal: FV, uniform wall heat flux (H1 condition) ---")
    print(f"{'NY':>5} {'Nu_fd':>10} {'target':>9} {'err%':>8}")
    for NY in [12, 20, 30, 40]:
        NX = 60*NY
        Nu, _ = fv_thermal(NY, NX, Pr=5.83, bc='H1', u_avg=0.02,
                           nu=0.1, Re=None)
        Nu_fd = Nu[int(0.9*NX)]
        print(f"{NY:>5} {Nu_fd:>10.4f} {8.235:>9.3f} {(Nu_fd-8.235)/8.235*100:>8.3f}")
