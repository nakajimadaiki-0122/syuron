#!/usr/bin/env python3
"""
06_compare_junction_gamboa.py -- 接合角を Gamboa の beta と同一基準で比較し作図する
==================================================================================

基準の定義（Gamboa 2005 Fig.1 / Table 2 と整合）
-----------------------------------------------
上流側接合において、主流路からループへ入る向きのループ軸方向を theta とする
（+x を 0 度、反時計回りを正）。Gamboa の beta は主流路の法線から測るので

    beta = theta - 90 [deg]

beta > 0 はループ軸が上流側へ傾くことを意味する。このときループから主流路へ
吐き出される流れは +x 成分を持ち、逆流（主流 -x）に対して逆向きに噴き出す。
beta < 0 なら吐き出しは主流と同方向になり、方向依存性は生じない。

この換算が Gamboa Table 2 の 2 行（最適形状 71.7 deg、参照形状 8.53 deg）を
同時に再現することを本スクリプトで確認する。
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Gamboa 2005 Table 2（論文記載値）
GAMBOA = {"Optimized": dict(X2=1.60, n=0.797, Y3=0.608, beta_table=71.7),
          "Reference": dict(X2=1.50, n=0.990, Y3=0.600, beta_table=8.53)}
WV_HALF = 0.5


def gamboa_theta(X2, n, Y3):
    """A=(X2, +0.5) から B=(n*X2, Y3) へ向かう向き（主流路 -> ループ）。"""
    d = np.array([n * X2 - X2, Y3 - WV_HALF])
    return float(np.degrees(np.arctan2(d[1], d[0])))


def main():
    res = json.load(open("results/geometry/loop_junction_angles.json"))

    print("=" * 78)
    print(" 換算式の検証: beta = theta - 90 が Gamboa Table 2 を再現するか")
    print("=" * 78)
    print(f"{'Valve':>12}{'theta [deg]':>14}{'beta = theta-90':>18}{'Table 2':>12}{'diff':>9}")
    for k, g in GAMBOA.items():
        th = gamboa_theta(g["X2"], g["n"], g["Y3"])
        g["theta"] = th
        g["beta"] = th - 90.0
        print(f"{k:>12}{th:>14.3f}{g['beta']:>18.3f}{g['beta_table']:>12.2f}"
              f"{g['beta']-g['beta_table']:>9.3f}")

    # 7 月形状: 上流側接合の「主流路 -> ループ」向きは、測定値（ループ -> 主流路）の逆
    print("\n" + "=" * 78)
    print(" 7 月形状（DXF から実測）")
    print("=" * 78)
    print(f"{'loop':>5}{'discharge [deg]':>18}{'theta [deg]':>14}"
          f"{'beta [deg]':>13}{'残差 [mm]':>12}")
    thetas = []
    for r in res:
        u = r["upstream"][1]          # fit 200 点（0.2 mm）を代表値にする
        disc = u["angle_from_plusx_deg"]
        th = disc + 180.0
        th = th - 360.0 if th > 180 else th
        thetas.append(th)
        print(f"{r['loop']:>5}{disc:>18.3f}{th:>14.3f}{th-90:>13.3f}"
              f"{u['straightness_rms_mm']:>12.1e}")
    th_j = float(np.mean(thetas))
    beta_j = th_j - 90.0

    print("\n" + "=" * 78)
    print(" 同一基準での比較")
    print("=" * 78)
    print(f"{'形状':>22}{'theta [deg]':>14}{'beta [deg]':>13}"
          f"{'逆流時の噴出角':>18}{'判定':>16}")
    rows = [("7 月形状（実測）", th_j, beta_j),
            ("Gamboa Reference", GAMBOA["Reference"]["theta"],
             GAMBOA["Reference"]["beta"]),
            ("Gamboa Optimized", GAMBOA["Optimized"]["theta"],
             GAMBOA["Optimized"]["beta"])]
    for name, th, be in rows:
        # 逆流（主流 -x = 180 deg）に対する噴出方向（theta+180）の成す角
        disc = th + 180.0
        ang = abs(((disc - 180.0) + 180) % 360 - 180)
        verdict = "対向（Tesla 動作）" if ang > 90 else "順方向（並列バイパス）"
        print(f"{name:>22}{th:>14.2f}{be:>13.2f}{ang:>18.2f}{verdict:>16}")

    # ---- 作図 ----
    d = np.load("results/geometry/loop1_centerline.npz")
    C, arc, flank, valid = d["C"], d["arc"], d["flank"], d["valid"]
    k = np.flatnonzero(valid)
    Cv = C[k.min():k.max() + 1]
    u = res[0]["upstream"][1]
    dn = res[0]["downstream"][0]

    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.plot(arc[:, 0], arc[:, 1], "k-", lw=1.6, label="wall (outer / island)")
    ax.plot(flank[:, 0], flank[:, 1], "k-", lw=1.6)
    ax.plot([4.0, 11.5], [0.5, 0.5], color="0.6", lw=1.0, ls="-")
    ax.plot([4.0, 11.5], [-0.5, -0.5], "k-", lw=1.6)
    ax.plot(C[:, 0], C[:, 1], color="0.75", lw=1.0, ls=":",
            label="centreline (incl. degenerate ends)")
    ax.plot(Cv[:, 0], Cv[:, 1], "C0-", lw=2.0, label="centreline (valid span)")
    ax.annotate("", xy=(5.6, -0.05), xytext=(4.4, -0.05),
                arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.6))
    ax.text(4.4, 0.05, "forward (+x)", color="0.4", fontsize=9)

    def wrap(a):
        return (a + 180) % 360 - 180

    for t, c, lab in ((u, "C3", "upstream junction"),
                      (dn, "C1", "downstream junction")):
        p = np.array(t["endpoint"])
        q = np.array(t["cross_y05"])
        ax.plot([p[0], q[0]], [p[1], q[1]], c, lw=2.2, ls="--")
        ax.plot(*q, "o", color=c, ms=7)
        th = wrap(t["angle_from_plusx_deg"] + 180.0)
        txt = (f"{lab}\ntheta = {th:+.1f} deg\nbeta = {th-90:+.1f} deg"
               if lab.startswith("upstream")
               else f"{lab}\ntheta undetermined\n(polyline too coarse)")
        ax.annotate(txt, xy=q, xytext=(q[0] - 1.5, q[1] - 1.55), color=c,
                    fontsize=9, arrowprops=dict(arrowstyle="->", color=c, lw=1.2))

    ax.set_aspect("equal"); ax.set_xlim(4.0, 11.5); ax.set_ylim(-2.1, 3.0)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title("July geometry, loop 1: loop centreline / main-channel junction "
                 "(measured from DXF)", fontsize=11, loc="left")
    ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=.25)
    os.makedirs("figures/geometry", exist_ok=True)
    plt.tight_layout(); plt.savefig("figures/geometry/loop_junction_angle.png", dpi=150)
    print("\nsaved figures/loop_junction_angle.png")


if __name__ == "__main__":
    main()
