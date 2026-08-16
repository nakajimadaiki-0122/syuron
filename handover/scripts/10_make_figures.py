#!/usr/bin/env python3
"""
10_make_figures.py -- 発表用の図をまとめて作る
==============================================

出力
  figures/di_vs_re.png        Di_p 対 Re（7 月形状 3 格子 + Gamboa Fig.7 + 誤差棒）
  figures/junction_angle.png  接合角の比較（7 月形状 vs Gamboa）
  figures/summary_table.md    数表

Gamboa の値は Fig.7 を画素計測で数値化したもの（docs/gamboa2005_geometry.md §4）。
低 Re 側は最適形状と参照形状の曲線が分離できないため ±0.02 の不確かさがある。
"""
import glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Gamboa 2005 Fig.7 optimized 2D CFD（画素計測）
GAMBOA = {100: 1.02, 200: 1.07, 300: 1.12, 500: 1.37, 600: 1.48,
          1000: 1.73, 2000: 1.92}
GAMBOA_ERR = {100: 0.02, 200: 0.02, 300: 0.02}

CUT = {12: 125, 24: 251, 48: 502}       # 直線区間で切った窓


def july_di():
    """7 月形状の Di_p。cpm ごとに Re をまとめる。"""
    out = {}
    for f in glob.glob("results/cell_tesla_channel_cell_fwd_Re*_cpm*_o*.json"):
        if "_tight" in f:
            continue
        j = json.load(open(f))
        cpm, Re = j["cells_per_mm"], int(j["Re"])
        r = f.replace("_fwd_", "_rev_")
        if not os.path.exists(r):
            continue
        out.setdefault(cpm, {})[Re] = json.load(open(r))["dp_Pa"] / j["dp_Pa"]
    # 旧タグ（cpm 無し）も拾う
    for f in glob.glob("results/cell_tesla_channel_cell_fwd_Re*_o145.json"):
        j = json.load(open(f))
        Re = int(j["Re"])
        r = f.replace("_fwd_", "_rev_")
        if os.path.exists(r) and Re not in out.get(12, {}):
            out.setdefault(12, {})[Re] = json.load(open(r))["dp_Pa"] / j["dp_Pa"]
    return out


def fig_di_vs_re(july, path="figures/di_vs_re.png"):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4),
                           gridspec_kw={"width_ratios": [1, 1]})

    # (a) Gamboa と同じ縦軸スケール
    re = sorted(GAMBOA)
    ax[0].plot(re, [GAMBOA[r] for r in re], "k-o", ms=4, lw=1.6,
               label="Gamboa 2005 optimized (2D CFD)")
    for cpm, mk in ((12, "^"), (24, "s"), (48, "o")):
        if cpm not in july:
            continue
        rr = sorted(july[cpm])
        ax[0].plot(rr, [july[cpm][r] for r in rr], mk + "--", ms=4, lw=1.2,
                   label=f"this work, July shape ({cpm} cells/mm)")
    ax[0].axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax[0].set_xscale("log"); ax[0].set_xlabel("Re"); ax[0].set_ylabel("$Di_p$")
    ax[0].set_title("(a) full range", fontsize=10, loc="left")
    ax[0].legend(fontsize=7.5); ax[0].grid(alpha=.3)

    # (b) 1 近傍の拡大
    for cpm, mk in ((12, "^"), (24, "s"), (48, "o")):
        if cpm not in july:
            continue
        rr = sorted(july[cpm])
        ax[1].plot(rr, [july[cpm][r] for r in rr], mk + "--", ms=5, lw=1.2,
                   label=f"{cpm} cells/mm")
    rr = [r for r in sorted(GAMBOA) if r <= 500]
    ax[1].errorbar(rr, [GAMBOA[r] for r in rr],
                   yerr=[GAMBOA_ERR.get(r, 0.0) for r in rr],
                   fmt="k-o", ms=4, lw=1.6, capsize=3,
                   label="Gamboa optimized")
    ax[1].axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax[1].set_xlabel("Re"); ax[1].set_ylabel("$Di_p$")
    ax[1].set_ylim(0.985, 1.15)
    ax[1].set_title("(b) close-up near $Di_p = 1$", fontsize=10, loc="left")
    ax[1].legend(fontsize=7.5); ax[1].grid(alpha=.3)

    os.makedirs("figures", exist_ok=True)
    plt.tight_layout(); plt.savefig(path, dpi=160)
    print(f"saved {path}")


def fig_angle(path="figures/junction_angle.png"):
    """接合角の比較。上流側接合での「主流路 -> ループ」向きを矢印で示す。"""
    shapes = [("July shape\n(this work, measured)", 47.4, "C3"),
              ("Gamboa reference", 98.5, "C0"),
              ("Gamboa optimized", 161.6, "C2")]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, (name, th, c) in zip(axs, shapes):
        ax.axhspan(-0.5, 0.5, color="0.92")
        ax.axhline(0.5, color="k", lw=1.2); ax.axhline(-0.5, color="k", lw=1.2)
        # 逆流は -x 向き
        ax.annotate("", xy=(-1.6, 0), xytext=(1.6, 0),
                    arrowprops=dict(arrowstyle="-|>", color="0.45", lw=2))
        ax.text(-0.9, 0.13, "reverse main flow  ($-x$)", color="0.45", fontsize=8)
        d = np.array([np.cos(np.radians(th)), np.sin(np.radians(th))])
        ax.annotate("", xy=(1.6 * d[0], 0.5 + 1.6 * d[1]), xytext=(0, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=2.4))
        ax.annotate("", xy=(-1.3 * d[0], 0.5 - 1.3 * d[1]), xytext=(0, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=2.4, ls="--"))
        ax.plot([0, 0], [0.5, 2.4], color="0.7", lw=0.8, ls=":")
        ax.plot([-2.2, 2.2], [0.5, 0.5], color="0.7", lw=0.8, ls=":")
        jet = abs(((th + 180) - 180 + 180) % 360 - 180)
        ax.set_title(f"{name}\n"
                     r"$\theta$ = " + f"{th:.1f}°,  " + r"$\beta$ = "
                     + f"{th-90:+.1f}°", fontsize=9)
        ax.text(0.02, 0.03, f"jet vs reverse flow: {jet:.1f}°\n"
                + ("opposing" if jet > 90 else "NOT opposing"),
                transform=ax.transAxes, fontsize=8.5, color=c,
                bbox=dict(fc="white", ec=c, alpha=.85))
        ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.0, 2.6); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Loop axis direction at the upstream junction "
                 "(solid: into loop, dashed: discharge)", fontsize=10)
    plt.tight_layout(); plt.savefig(path, dpi=160)
    print(f"saved {path}")


def table(july, path="figures/summary_table.md"):
    lines = ["# 発表用 数表", "",
             "## 7 月形状の Di_p（周期単位セル、1 段あたり）", ""]
    cpms = sorted(july)
    res = sorted({r for c in cpms for r in july[c]})
    lines.append("| Re | " + " | ".join(f"{c} cells/mm" for c in cpms)
                 + " | Gamboa optimized |")
    lines.append("|---|" + "---|" * (len(cpms) + 1))
    for r in res:
        row = [f"{july[c].get(r, float('nan')):.5f}"
               if r in july[c] else "—" for c in cpms]
        g = f"{GAMBOA[r]:.2f}" if r in GAMBOA else "—"
        lines.append(f"| {r} | " + " | ".join(row) + f" | {g} |")
    lines += ["", "## 接合角", "",
              "| 形状 | θ [deg] | β [deg] | 逆流時の噴出角 | 判定 |",
              "|---|---|---|---|---|",
              "| 7 月形状（DXF 実測） | 47.4 | −42.6 | 47.4° | 並列バイパス |",
              "| Gamboa reference | 98.5 | +8.5 | 98.5° | 対向 |",
              "| Gamboa optimized | 161.6 | +71.6 | 161.6° | 対向 |", ""]
    open(path, "w").write("\n".join(lines))
    print(f"saved {path}")


if __name__ == "__main__":
    j = july_di()
    print("July shape Di_p:", {c: {r: round(v, 5) for r, v in d.items()}
                               for c, d in j.items()})
    fig_di_vs_re(j)
    fig_angle()
    table(j)
