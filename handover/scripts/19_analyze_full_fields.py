#!/usr/bin/env python3
"""
19_analyze_full_fields.py -- 案 1 の場からループ流量比と機構を確認する
=====================================================================

`17_run_gamboa_full.py --save-fields` が保存した npz を読んで、

  1. φ_loop(x) = 主流路帯（|y| <= w/2）の外を通る質量流量の割合
     **主流路が存在する x 区間だけ**で評価する。プレナムまで含めると
     帯の外が領域の大半になり φ_loop = 1 になってしまう（17 の値はこの意味で無効）
  2. ループ枝を横切る質量流量（戻り流路と下流枝でそれぞれ）
  3. 速度場・流線の図

を出す。**逆流でループへ十分な流量が入っているか**が Tesla バルブとして
機能しているかの直接の指標になる（Gamboa 3.1 節）。

順流の npz と逆流の npz は同じ格子（逆流は x 反転した形状）なので、
逆流側は x を反転して物理座標に戻してから比べる。
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import gamboa_full as GF


def load(Re, cpm, which):
    p = os.path.join(HERE, "..", "results",
                     f"gamboa_full_Re{Re}_cpm{cpm}_{which}_fields.npz")
    if not os.path.exists(p):
        return None
    d = np.load(p)
    key = which
    return dict(xs=d["xs"], ys=d["ys"], ux=d[f"{key}_ux"], uy=d[f"{key}_uy"],
                p=d[f"{key}_p"], rho=d[f"{key}_rho"], mask=d[f"{key}_mask"])


def phi_loop(fd, xs, ys, x_lo, x_hi, half=0.5):
    """主流路帯の外を通る質量流量の割合を x 断面ごとに返す。"""
    ux, rho, m = fd["ux"], fd["rho"], fd["mask"]
    g = rho * ux
    main = np.abs(ys) <= half + 1e-9
    out = []
    for i in range(len(xs)):
        if not (x_lo <= xs[i] <= x_hi) or not m[i].any():
            out.append(np.nan)
            continue
        tot = g[i, m[i]].sum()
        mm = g[i, m[i] & main].sum()
        out.append(np.nan if abs(tot) < 1e-30 else 1.0 - mm / tot)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=int, default=100)
    ap.add_argument("--cpm", type=int, default=16)
    a = ap.parse_args()

    poly, info = GF.build("optimized", w=1.0, plenum_out="vertical")
    x_lo = -0.5                      # 入口流路の左端
    # 主流路が幅 w を保つのは外壁が y=+0.5 を横切る点 E まで。
    # その先は接合部の楔なので「主流路帯」の定義が壊れる。
    x_hi = float(info["E"][0])
    print(f"主流路の区間: x = {x_lo:.2f} .. {x_hi:.2f} w_v")

    res = {}
    fig, axes = plt.subplots(3, 1, figsize=(12, 13),
                             height_ratios=[1, 1, 0.6])
    for k, (lab, ttl) in enumerate((("fwd", "forward"), ("rev", "reverse"))):
        fd = load(a.Re, a.cpm, lab)
        if fd is None:
            print(f"  {lab}: npz が無い（--save-fields を付けて計算する）")
            continue
        xs, ys = fd["xs"], fd["ys"]
        # 逆流は形状を x 反転して解いている。配列を戻して ux の符号を反転すれば
        # **元の格子座標 xs のまま**物理場になる（xs を反転してはいけない）。
        if lab == "rev":
            for key in ("ux", "uy", "p", "rho", "mask"):
                fd[key] = fd[key][::-1]
            fd["ux"] = -fd["ux"]
        phi = phi_loop(fd, xs, ys, x_lo, x_hi)
        good = ~np.isnan(phi)
        res[lab] = dict(phi_mean=float(np.nanmean(phi)),
                        phi_max=float(np.nanmax(phi)),
                        phi_at_loop_span=[float(xs[good][0]),
                                          float(xs[good][-1])])
        print(f"  {lab}: φ_loop 平均 {res[lab]['phi_mean']:+.4f}  "
              f"最大 {res[lab]['phi_max']:+.4f}")

        ax = axes[k]
        m = fd["mask"]
        sp = np.where(m, np.hypot(fd["ux"], fd["uy"]), np.nan)
        ext = [xs[0], xs[-1], ys[0], ys[-1]]
        im = ax.imshow(sp.T / 0.05, origin="lower", extent=ext, aspect="equal",
                       cmap="viridis", vmin=0, vmax=np.nanpercentile(sp, 99.8) / 0.05)
        ax.streamplot(xs, ys, np.where(m, fd["ux"], np.nan).T,
                      np.where(m, fd["uy"], np.nan).T, color="w",
                      density=2.2, linewidth=0.4, arrowsize=0.6)
        ax.set_title(f"{ttl}  |u|/U_ch   Re={a.Re}, cpm={a.cpm}   "
                     f"phi_loop mean {res[lab]['phi_mean']:+.3f}")
        ax.set_xlim(-1.5, 10.5); ax.set_ylim(-3.0, 5.5)
        plt.colorbar(im, ax=ax, shrink=0.8)

        axes[2].plot(xs, phi, label=f"{lab} (mean {res[lab]['phi_mean']:+.3f})")

    axes[2].axhline(0, color="k", lw=0.5)
    axes[2].set_xlabel("x [w_v]"); axes[2].set_ylabel("phi_loop")
    axes[2].set_xlim(x_lo, x_hi); axes[2].grid(alpha=0.3); axes[2].legend()
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures",
                       f"gamboa_full_fields_Re{a.Re}_cpm{a.cpm}.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")

    if len(res) == 2:
        print(f"\n  φ_loop(逆) / φ_loop(順) = "
              f"{res['rev']['phi_mean']/res['fwd']['phi_mean']:.3f}")
        print("  Gamboa 3.1 節: 高 Di には**逆流時にループへ十分な流量が入る**ことが要る。"
              "\n  この比が 1 を超えていれば機構が働いている。"
              "案 2（閉じ壁）では 0.72 で、逆になっていた。")
        j = os.path.join(HERE, "..", "results",
                         f"gamboa_full_phi_Re{a.Re}_cpm{a.cpm}.json")
        json.dump(res, open(j, "w"), indent=1)
        print(f"saved {os.path.normpath(j)}")


if __name__ == "__main__":
    main()
