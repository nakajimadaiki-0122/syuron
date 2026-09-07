#!/usr/bin/env python3
"""
17_run_gamboa_full.py -- Gamboa 単発モデル（プレナム込み）の順流・逆流を解く
============================================================================

案 1 の本体。`src/gamboa_full.py` の形状を `lbm.solve_flow_io`（入口一様流速・
出口定圧）で解き、Gamboa Fig.7 の 2D CFD 曲線と比べる。

Gamboa の条件（2.2 節）に合わせている点
  - 2D 定常非圧縮層流
  - 入口 BC = **プレナムの直線壁に一様流速**
  - 出口 BC = **プレナムの直線壁で圧力ゼロ**
  - 特性速度 U = 入口流路の断面平均流速、特性長 = w_v、D_h = 2 w_v
  - Re = rho U D_h / mu
  - 格子は流路横断方向 16 要素（Gamboa と同じ。--cpm で変更可）

逆流は形状を x 反転して同じ +x 方向に解く。左右のプレナムはどちらも直線壁が
鉛直で高さ 2 R_p なので、反転しても BC の形が変わらない（`plenum_out="vertical"`）。

流量は順逆で同じにする。入口面のセル数が順逆でわずかに違うので、
一様流速ではなく**体積流量**を揃える（U_in = Q_target / n_inlet）。

Di = Δp_逆 / Δp_順、Δp は入口面と出口面の流量重み平均圧力の差。
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import gamboa_full as GF
import lbm
import postproc as pp


def solve_one(mask, Re, U_ch, cpm, label, iters, tol, ckpt=None, w_mm=1.0,
              dp_tol=1e-5, dp_window=10000, avg_tol=1e-3):
    nx, ny = mask.shape
    n_in = int(mask[0].sum())
    n_out = int(mask[-1].sum())
    Q = U_ch * cpm                      # 目標体積流量（流路幅 = cpm セル）
    U_in = Q / n_in
    nu = lbm.nu_lattice(U_ch, 2.0 * cpm, Re)
    print(f"  [{label}] {nx} x {ny}, 流体 {int(mask.sum())} セル, "
          f"入口 {n_in} / 出口 {n_out} セル, U_in = {U_in:.5f}, "
          f"nu = {nu:.5f}, tau_p = {lbm.tau_from_nu(nu):.4f}", flush=True)

    state = {}

    def rep(it, d, mi, mo, dp=None):
        if it % 20000 == 0:
            print(f"    it={it:7d} res={d:.2e} mdot in/out={mi:.5e}/{mo:.5e} "
                  f"dp={dp:.6e}", flush=True)

    t0 = time.time()
    ux, uy, p, rho, f, info = lbm.solve_flow_io(
        mask, nu, U_in, iters=iters, tol=tol, report=rep, ckpt=ckpt,
        dp_tol=dp_tol, dp_window=dp_window, avg_tol=avg_tol)

    # 非定常なときは瞬時値ではなく**時間平均**を使う（README §8 B-13）
    dp_inst = pp.pressure_drop(p, ux, mask, i0=1, i1=nx - 2)
    dp_lat = info["dp_lattice_mean"]
    scale = pp.lattice_to_physical_scale(Re, 2.0 * w_mm * 1e-3, U_ch)
    dp_Pa = pp.to_pascal(dp_lat, scale)
    dp_Pa_std = pp.to_pascal(info["dp_lattice_std"], scale)
    # 内部の断面流量のばらつき（BC 列を除く）で質量保存を見る
    g = rho * ux
    fl = np.array([g[i, mask[i]].sum() for i in range(1, nx - 1)])
    res = dict(label=label, nx=nx, ny=ny, fluid=int(mask.sum()),
               n_inlet=n_in, n_outlet=n_out, U_in=U_in, U_ch=U_ch, nu=nu,
               tau_p=info["tau_p"], converged=info["converged"],
               iters=info["iters"], residual=info["residual"],
               dp_lattice=dp_lat, dp_Pa=dp_Pa, dp_Pa_std=dp_Pa_std,
               dp_lattice_inst=dp_inst, dp_lattice_std=info["dp_lattice_std"],
               dp_lattice_min=info["dp_lattice_min"],
               dp_lattice_max=info["dp_lattice_max"],
               unsteady=info["unsteady"], dp_samples=info["dp_samples"],
               mdot_interior_mean=float(fl.mean()),
               mdot_interior_spread=float((fl.max() - fl.min()) / fl.mean()),
               wall_s=round(time.time() - t0, 1))
    print(f"  [{label}] dp = {dp_lat:.6e} +- {info['dp_lattice_std']:.1e} "
          f"(格子, 直近 {info['dp_samples']} 点の平均) = {dp_Pa:.4f} +- "
          f"{dp_Pa_std:.4f} Pa, 収束 {info['converged']}"
          f"{'（非定常・統計的）' if info['unsteady'] else ''} "
          f"it={info['iters']} res={info['residual']:.2e} "
          f"{res['wall_s']:.0f}s", flush=True)
    return res, ux, uy, p, rho


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--cpm", type=int, default=16, help="流路幅のセル数")
    ap.add_argument("--U", type=float, default=0.05, help="流路の平均流速（格子）")
    ap.add_argument("--iters", type=int, default=600000)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--which", default="optimized")
    ap.add_argument("--out", default=None)
    ap.add_argument("--save-fields", action="store_true")
    ap.add_argument("--dir", default="both", choices=["both", "fwd", "rev"],
                    help="向きを分けて別プロセスで走らせると並列化できる")
    ap.add_argument("--dp-window", type=int, default=10000,
                    help="Δp の変化を見る窓［反復］")
    ap.add_argument("--avg-tol", type=float, default=1e-3,
                    help="非定常時: 時間平均が動かなくなったと判定する相対差")
    ap.add_argument("--dp-tol", type=float, default=1e-5,
                    help="Δp の相対変化がこれ未満で収束とみなす")
    a = ap.parse_args()
    w_mm = 1.0

    poly, info = GF.build(a.which, w=1.0, plenum_out="vertical")
    mask, xs, ys = GF.rasterize(poly, a.cpm, w=1.0)
    print("=" * 80)
    print(f" Gamboa {a.which} 単発モデル（プレナム込み）  Re = {a.Re}, "
          f"cpm = {a.cpm}")
    print("=" * 80)
    print(f"  領域 {mask.shape[0]} x {mask.shape[1]} セル "
          f"= {xs[-1]-xs[0]:.2f} x {ys[-1]-ys[0]:.2f} w_v, "
          f"流体率 {100*mask.mean():.1f} %")

    out = dict(Re=a.Re, cpm=a.cpm, U=a.U, which=a.which,
               geometry=dict(theta_deg=info["theta_deg"],
                             beta_deg=info["beta_deg"],
                             alpha_deg=info["alpha_deg"],
                             lenout_axis=info["lenout_axis"],
                             Rp=info["Rp"], bounds=info["bounds"]))
    fields = {}
    todo = [("fwd", mask), ("rev", mask[::-1, :].copy())]
    if a.dir != "both":
        todo = [t for t in todo if t[0] == a.dir]
    for lab, m in todo:
        ck = os.path.join(HERE, "..", "results", "single", f"ck_full_{lab}_Re{int(a.Re)}_cpm{a.cpm}.npy")
        r, ux, uy, p, rho = solve_one(m, a.Re, a.U, a.cpm, lab, a.iters,
                                      a.tol, ckpt=ck, w_mm=w_mm,
                                      dp_tol=a.dp_tol, dp_window=a.dp_window,
                                      avg_tol=a.avg_tol)
        # phi_loop（主流路帯の外を通る質量流量の割合）は主流路が水平な区間のみ意味を持つ
        lp = pp.loop_mass_fraction(ux, m, ys, main_half_width_mm=0.5, rho=rho)
        r["phi_loop_mean"] = lp["phi_loop_mean"]
        r["phi_loop_max"] = lp["phi_loop_max"]
        out[lab] = r
        fields[lab] = (ux, uy, p, rho, m)
        if os.path.exists(ck):
            os.remove(ck)

    o = a.out or os.path.join(HERE, "..", "results", "single", f"gamboa_full_Re{int(a.Re)}_cpm{a.cpm}"
                              + ("" if a.dir == "both" else f"_{a.dir}")
                              + ".json")
    if a.dir != "both":
        # 相手の向きの結果があれば読み込んで Di を出す
        other = "rev" if a.dir == "fwd" else "fwd"
        op = o.replace(f"_{a.dir}.json", f"_{other}.json")
        if os.path.exists(op):
            out[other] = json.load(open(op))[other]
    if a.save_fields:
        npz = o.replace(".json", "_fields.npz")
        np.savez_compressed(npz, xs=xs, ys=ys,
                            **{f"{k}_{n}": v for k, (ux, uy, p, rho, m)
                               in fields.items()
                               for n, v in (("ux", ux), ("uy", uy),
                                            ("p", p), ("rho", rho),
                                            ("mask", m))})
        print(f"saved {os.path.normpath(npz)}")

    if "fwd" not in out or "rev" not in out:
        os.makedirs(os.path.dirname(o), exist_ok=True)
        json.dump(out, open(o, "w"), indent=1, default=float)
        print("saved " + os.path.normpath(o) + "  （片方向のみ。Di は未算出）")
        return
    Di = out["rev"]["dp_Pa"] / out["fwd"]["dp_Pa"]
    out["Di"] = Di
    # Fig.7 の opt CFD（破線）を画素実測した値（scripts/23_digitize_fig7.py）。
    # Re <= 250 は実線（ref CFD）と重なって分離できないので入れない。
    # 2026-08-09 の目視値（Re=300 で 1.12）は**実線を読んでいた**誤りだった。
    ref = {300: 1.1622, 400: 1.2523, 500: 1.3667, 600: 1.4808,
           750: 1.5982, 1000: 1.7332, 1250: 1.8079, 1500: 1.8557}
    tgt = ref.get(int(a.Re))
    print(f"\n  Di = Δp_逆 / Δp_順 = {out['rev']['dp_Pa']:.4f} / "
          f"{out['fwd']['dp_Pa']:.4f} = {Di:.4f}")
    if tgt:
        print(f"  Gamboa Fig.7 の 2D CFD = {tgt:.2f}  -> 差 {Di-tgt:+.3f} "
              f"({100*(Di/tgt-1):+.1f} %)")
        out["gamboa_fig7"] = tgt
    if not tgt:
        print("  Fig.7 の実測値なし（Re <= 250 は実線と重なって分離できない）")
    print(f"  参考: 周期版（案 2、閉じ壁）の Di_valve は "
          f"{'0.9925' if a.Re==100 else '0.9678' if a.Re==500 else 'n/a'}")

    os.makedirs(os.path.dirname(o), exist_ok=True)
    json.dump(out, open(o, "w"), indent=1, default=float)
    print(f"\nsaved {os.path.normpath(o)}")

    # --- 図 ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if len(fields) < 2:
        return
    fig, axes = plt.subplots(2, 1, figsize=(11, 12))
    for ax, lab in zip(axes, ("fwd", "rev")):
        ux, uy, p, rho, m = fields[lab]
        sp = np.where(m, np.hypot(ux, uy) / a.U, np.nan)
        ax.imshow(sp.T, origin="lower", cmap="viridis", aspect="equal",
                  vmin=0, vmax=np.nanpercentile(sp, 99.5))
        ax.streamplot(np.arange(m.shape[0]), np.arange(m.shape[1]),
                      np.where(m, ux, np.nan).T, np.where(m, uy, np.nan).T,
                      color="w", density=2.0, linewidth=0.4, arrowsize=0.6)
        ax.set_title(f"{lab}  |u|/U_ch   Re={a.Re}  cpm={a.cpm}  "
                     f"dp={out[lab]['dp_Pa']:.3f} Pa"
                     f"  phi_loop={out[lab]['phi_loop_mean']:+.3f}")
    fig.suptitle(f"Gamboa full model  Di = {Di:.4f}"
                 + (f"  (Fig.7: {tgt})" if tgt else ""))
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "flow", f"gamboa_full_Re{int(a.Re)}_cpm{a.cpm}.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
