#!/usr/bin/env python3
"""
03_compare.py -- 順流・逆流を並べて比較し、図を出力する
=======================================================

使い方
------
    python scripts/03_compare.py --fwd results/flow_tesla_channel_fwd_Re100.npz \
                                 --rev results/flow_tesla_channel_rev_Re100.npz \
                                 --out figures/fwd_vs_rev.png

判定基準（2026-08-09 改訂）
--------------------------
**Di_p を単独の閾値で判定してはいけない。** Gamboa 2005 の最適化形状ですら
2D CFD で Re = 100 では Di = 1.02 であり、方向依存性は Re が 400 を超えてから
立ち上がる（docs/gamboa2005_geometry.md §4）。

判定は Di(Re) 曲線の立ち上がりで行う。比較対象（Gamboa optimized, 2D CFD）:

    Re   100   200   300   500   600  1000  2000
    Di  1.02  1.07  1.12  1.37  1.48  1.73  1.92
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import postproc as pp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fwd", required=True)
    ap.add_argument("--rev", required=True)
    ap.add_argument("--out", default="figures/fwd_vs_rev.png")
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--W", type=float, default=1.0)
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--ref-phi", type=float, default=18.0,
                    help="reference loop fraction to draw, e.g. Porwal forward max")
    a = ap.parse_args()

    df, dr = np.load(a.fwd), np.load(a.rev)
    scale = pp.lattice_to_physical_scale(a.Re, 2 * a.W * 1e-3, a.U)

    jf = json.load(open(a.fwd.replace(".npz", ".json")))
    jr = json.load(open(a.rev.replace(".npz", ".json")))
    Di = pp.diodicity(jr["dp_Pa"], jf["dp_Pa"])

    print(f"{'':14s}{'dp [Pa]':>10s}{'phi_mean':>10s}{'phi_max':>10s}{'massimb %':>11s}")
    for nm, j in (("forward", jf), ("reverse", jr)):
        print(f"{nm:14s}{j['dp_Pa']:10.3f}{j['phi_loop_mean']:10.4f}"
              f"{j['phi_loop_max']:10.4f}{j['mass_imbalance_pct']:11.2f}")
    # Gamboa 2005 Fig.7, optimized valve, 2D CFD（画素計測による数値化）
    GAMBOA_CFD = {100: 1.02, 200: 1.07, 300: 1.12, 500: 1.37,
                  600: 1.48, 1000: 1.73, 2000: 1.92}
    print(f"\n  Di_p = dp_rev / dp_fwd = {Di:.4f}")
    ref = GAMBOA_CFD.get(int(a.Re))
    if ref is not None:
        print(f"  Gamboa optimized 2D CFD at Re={int(a.Re)}: Di = {ref:.2f}")
    if a.Re <= 200:
        print("  NOTE: Re <= 200 では Di ~ 1 は正常。この Re で形状の良否は判定できない。")
    elif ref is not None and Di < 1 + 0.5 * (ref - 1):
        print("  --> 方向依存性が Gamboa の半分未満。幾何 / BC / 格子を疑う。")

    fig, ax = plt.subplots(3, 1, figsize=(11, 7.2),
                           gridspec_kw={"height_ratios": [1, 1, 0.85]})
    vmax = 0.0
    for d in (df, dr):
        U = np.sqrt(d["ux"] ** 2 + d["uy"] ** 2) * scale * 1000
        vmax = max(vmax, np.nanpercentile(np.where(d["mask"], U, np.nan), 99.5))

    ext = [df["xs"].min(), df["xs"].max(), df["ys"].min(), df["ys"].max()]
    for k, (d, lab) in enumerate([(df, "Forward"), (dr, "Reverse (mirrored)")]):
        U = np.sqrt(d["ux"] ** 2 + d["uy"] ** 2) * scale * 1000
        im = ax[k].imshow(np.where(d["mask"], U, np.nan).T, origin="lower",
                          aspect="equal", extent=ext, vmin=0, vmax=vmax,
                          cmap="viridis")
        ax[k].set_title(f"{lab}   —   flow →", fontsize=10, loc="left")
        ax[k].set_ylabel("y [mm]")
        plt.colorbar(im, ax=ax[k], label="|u| [mm/s]", fraction=0.025, pad=0.01)
    ax[0].set_xticklabels([])

    for d, lab, c in [(df, "Forward", "C0"), (dr, "Reverse", "C3")]:
        phi = d["phi"]
        x = np.linspace(ext[0], ext[1], len(phi))
        ax[2].plot(x, 100 * np.where(np.isfinite(phi), phi, np.nan), c, lw=1.4, label=lab)
    if a.ref_phi:
        ax[2].axhline(a.ref_phi, color="k", ls="--", lw=1,
                      label=f"reference ({a.ref_phi:.0f} %)")
    ax[2].set_xlim(ext[0], ext[1])
    ax[2].set_xlabel("x [mm]"); ax[2].set_ylabel(r"$\phi_{loop}$ [%]")
    ax[2].legend(fontsize=8, ncol=3); ax[2].grid(alpha=.3)
    ax[2].set_title("Mass-flow fraction through the loops", fontsize=10, loc="left")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    plt.tight_layout(); plt.savefig(a.out, dpi=140)
    print(f"\nsaved {a.out}")


if __name__ == "__main__":
    main()
