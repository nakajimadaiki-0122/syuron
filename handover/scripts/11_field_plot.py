#!/usr/bin/env python3
"""
11_field_plot.py -- 保存済みの場データから速度分布・流線・流量配分を作図する
============================================================================

`04_run_unitcell.py` が保存する npz（ux, uy, p, rho, mask, xs, ys, phi）を読む。

    python scripts/11_field_plot.py \
        --fwd results/cell_tesla_channel_cell_fwd_Re100_cpm24_o251.npz \
        --rev results/cell_tesla_channel_cell_rev_Re100_cpm24_o251.npz \
        --out figures/field_Re100.png --Re 100

**温度分布は出せない。** 熱は意図的に外しており、2 次元形状用の熱ソルバは存在しない
（`01_validate_solver.py` の熱パートは平行平板の解析解照合専用）。README §3 を参照。
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import postproc as pp


def load(path):
    d = np.load(path)
    return {k: d[k] for k in d.files}


def panel(ax, d, scale, vmax, title, stream=True):
    ux, uy, mask = d["ux"], d["uy"], d["mask"]
    xs, ys = d["xs"], d["ys"]
    U = np.sqrt(ux ** 2 + uy ** 2) * scale * 1000.0        # mm/s
    ext = [xs.min(), xs.max(), ys.min(), ys.max()]
    im = ax.imshow(np.where(mask, U, np.nan).T, origin="lower", aspect="equal",
                   extent=ext, vmin=0, vmax=vmax, cmap="viridis")
    if stream:
        # streamplot は等間隔格子を要求するので mask 外を 0 にして描く
        X = np.linspace(ext[0], ext[1], ux.shape[0])
        Y = np.linspace(ext[2], ext[3], ux.shape[1])
        ax.streamplot(X, Y, np.where(mask, ux, 0).T, np.where(mask, uy, 0).T,
                      color="w", linewidth=0.5, density=1.5, arrowsize=0.7)
    ax.contour(np.linspace(ext[0], ext[1], mask.shape[0]),
               np.linspace(ext[2], ext[3], mask.shape[1]),
               mask.T.astype(float), levels=[0.5], colors="k", linewidths=1.0)
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_ylabel("y [mm]")
    return im


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fwd", required=True)
    ap.add_argument("--rev", required=True)
    ap.add_argument("--out", default="figures/field.png")
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--W", type=float, default=1.0)
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--no-stream", action="store_true")
    a = ap.parse_args()

    df, dr = load(a.fwd), load(a.rev)
    scale = pp.lattice_to_physical_scale(a.Re, 2 * a.W * 1e-3, a.U)
    vmax = 0.0
    for d in (df, dr):
        U = np.sqrt(d["ux"] ** 2 + d["uy"] ** 2) * scale * 1000.0
        vmax = max(vmax, float(np.nanpercentile(np.where(d["mask"], U, np.nan), 99.5)))

    fig, ax = plt.subplots(3, 1, figsize=(10, 8.2),
                           gridspec_kw={"height_ratios": [1, 1, 0.7]})
    im = panel(ax[0], df, scale, vmax, f"Forward  (Re = {a.Re:.0f})",
               not a.no_stream)
    panel(ax[1], dr, scale, vmax, f"Reverse (mirrored)  (Re = {a.Re:.0f})",
          not a.no_stream)
    for k in (0, 1):
        plt.colorbar(im, ax=ax[k], label="|u| [mm/s]", fraction=0.02, pad=0.01)
    ax[0].set_xticklabels([])

    for d, lab, c in ((df, "Forward", "C0"), (dr, "Reverse", "C3")):
        phi = d["phi"]
        x = np.linspace(d["xs"].min(), d["xs"].max(), len(phi))
        ax[2].plot(x, 100 * np.where(np.isfinite(phi), phi, np.nan), c, lw=1.5,
                   label=lab)
    ax[2].set_xlim(df["xs"].min(), df["xs"].max())
    ax[2].set_xlabel("x [mm]"); ax[2].set_ylabel(r"$\phi_{loop}$ [%]")
    ax[2].set_title("Mass-flow fraction through the loop", fontsize=10, loc="left")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    plt.tight_layout(); plt.savefig(a.out, dpi=160)
    print(f"saved {a.out}")
    print("注意: 温度分布は出せない（熱ソルバは 2 次元形状に未実装。README §3）")


if __name__ == "__main__":
    main()
