#!/usr/bin/env python3
"""
09_gamboa_case2.py -- 案 2（直線接続区間を挟んだ周期配列）の成立確認
====================================================================

Gamboa のバルブは周期形状ではない（出口区間が主流路軸から 48.1 度傾く）。
`src/gamboa.py` はこれを「閉じ壁」に置き換えて周期化する。
バルブ前後に長さ gap の直線区間を足し、直線区間の寄与を差し引いて

    Di_valve = (dp_R - dp_str) / (dp_F - dp_str)

を得る。dp_str は同じ長さ・同じ格子・同じ流量の直線流路セルの圧力損失。

**接続部の形状について**
本構成の接続部は**直線のみで、折れ角を含まない**。閉じ壁は主流路上壁 y=+w/2 に
戻して終端するので、接続部は幅 w の平行平板そのものである。したがって
順逆対称は幾何から自明であり、検証 4 は**形式的確認**である
（鏡像反転しても同一のマスクになるため、計算は同一になる）。

検証項目
  4. 接続部のみのセルで dp_F と dp_R が一致するか（形式的）
  5. gap を変えて Di_valve が一致するか（本質的）
  6. dp_str / dp_F の比（大きいと差し引きで不確かさが増幅する）
"""
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import gamboa as G
import lbm
import postproc as pp


def solve_cell(mask, Re, U, cpm, w_mm, tol=1e-9, iters=3_000_000, label=""):
    """周期セルを流量一定で解き、Δp [Pa] と体積力を返す。"""
    Dh_cells = 2.0 * w_mm * cpm / w_mm * w_mm      # = 2 * cells_per_w
    Dh_cells = 2.0 * cpm
    nu = lbm.nu_lattice(U, Dh_cells, Re)
    mdot_target = U * cpm                           # 主流路幅 = cpm セル
    t0 = time.time()
    ux, uy, p, rho, f, info = lbm.solve_flow_periodic(
        mask, nu, mdot_target, iters=iters, tol=tol)
    scale = pp.lattice_to_physical_scale(Re, 2.0 * w_mm * 1e-3, U)
    dp_Pa = pp.to_pascal(info["dp_lattice"], scale)
    return dict(label=label, nx=int(mask.shape[0]), ny=int(mask.shape[1]),
                fluid=int(mask.sum()), nu=nu, tau_p=info["tau_p"],
                G=info["G"], dp_lattice=info["dp_lattice"], dp_Pa=dp_Pa,
                converged=info["converged"], iters=info["iters"],
                dG=info["dG_final"], mdot_err=info["mdot_err"],
                wall_s=round(time.time() - t0, 1)), ux, rho


def build(theta, gap, cpm, w_mm=1.0):
    poly, isl, info = G.gamboa_valve("optimized", w=w_mm)
    poly, isl, info = G.periodic_cell(poly, info, gap * w_mm, w=w_mm)
    mask, xs, ys = G.rasterize_cell(poly, cpm, w=w_mm)
    return mask, xs, ys, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--cpm", type=int, default=24, help="cells per w_v")
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--gaps", type=float, nargs="+", default=[1.0, 3.0])
    ap.add_argument("--iters", type=int, default=3_000_000)
    ap.add_argument("--out", default="results/gamboa_case2.json")
    a = ap.parse_args()
    w_mm = 1.0
    th = G.theta_from_params(G.OPTIMIZED["X2"], G.OPTIMIZED["n"], G.OPTIMIZED["Y3"])

    print("=" * 84)
    print(f" 案 2 の成立確認   Re = {a.Re},  cpm = {a.cpm},  theta = {th:.4f} deg")
    print("=" * 84)
    print("接続部は**直線のみ・折れ角なし**。検証 4 は形式的確認。\n")

    out = dict(Re=a.Re, cpm=a.cpm, U=a.U, theta_deg=th, gaps=a.gaps, cases=[])

    # ---- 検証 4: 接続部のみ（直線）の順逆対称 ----
    print("--- 検証 4: 接続部のみのセル（直線流路）---")
    Lstr = 4.0
    sp = G.straight_cell(Lstr * w_mm, w=w_mm)
    m_str, _, _ = G.rasterize_cell(sp, a.cpm, w=w_mm)
    r_f, _, _ = solve_cell(m_str, a.Re, a.U, a.cpm, w_mm, iters=a.iters,
                           label="straight fwd")
    r_r, _, _ = solve_cell(m_str[::-1, :].copy(), a.Re, a.U, a.cpm, w_mm,
                           iters=a.iters, label="straight rev")
    d = abs(r_r["dp_Pa"] / r_f["dp_Pa"] - 1.0)
    print(f"  L = {Lstr} w_v, {m_str.shape[0]} x {m_str.shape[1]} セル")
    print(f"  dp_F = {r_f['dp_Pa']:.9f} Pa,  dp_R = {r_r['dp_Pa']:.9f} Pa,  "
          f"相対差 = {d:.2e}")
    # 解析解との照合（f*Re = 96, Dh = 2w）
    H = a.cpm
    fRe_expected = 96.0
    u_avg = a.U
    Dh = 2.0 * a.cpm
    fRe = (r_f["G"] * Dh / (0.5 * u_avg ** 2)) * (u_avg * Dh / r_f["nu"])
    print(f"  f*Re = {fRe:.4f}（解析解 96、離散理論値は格子依存）  "
          f"-> 順逆対称: {'合格（自明）' if d < 1e-12 else '要確認'}")
    out["check4"] = dict(L=Lstr, dp_F=r_f["dp_Pa"], dp_R=r_r["dp_Pa"],
                         rel_diff=d, fRe=fRe, connector="straight only, no bend",
                         nature="formal (mirror-symmetric by construction)")

    # ---- 検証 5 & 6: gap を変えて Di_valve が一致するか ----
    print("\n--- 検証 5, 6: gap 依存性 ---")
    print(f"{'gap':>6}{'セル長':>9}{'dp_F [Pa]':>13}{'dp_R [Pa]':>13}"
          f"{'dp_str [Pa]':>13}{'dp_str/dp_F':>13}{'Di_cell':>10}{'Di_valve':>10}")
    for gap in a.gaps:
        mask, xs, ys, info = build(th, gap, a.cpm, w_mm)
        rf, _, _ = solve_cell(mask, a.Re, a.U, a.cpm, w_mm, iters=a.iters,
                              label=f"valve fwd gap={gap}")
        rr, _, _ = solve_cell(mask[::-1, :].copy(), a.Re, a.U, a.cpm, w_mm,
                              iters=a.iters, label=f"valve rev gap={gap}")
        s = G.straight_cell(gap * w_mm, w=w_mm)
        ms, _, _ = G.rasterize_cell(s, a.cpm, w=w_mm)
        rs, _, _ = solve_cell(ms, a.Re, a.U, a.cpm, w_mm, iters=a.iters,
                              label=f"straight gap={gap}")
        Di_cell = rr["dp_Pa"] / rf["dp_Pa"]
        Di_valve = (rr["dp_Pa"] - rs["dp_Pa"]) / (rf["dp_Pa"] - rs["dp_Pa"])
        ratio = rs["dp_Pa"] / rf["dp_Pa"]
        print(f"{gap:>6.1f}{info['cell_length']:>9.3f}{rf['dp_Pa']:>13.6f}"
              f"{rr['dp_Pa']:>13.6f}{rs['dp_Pa']:>13.6f}{ratio:>13.4f}"
              f"{Di_cell:>10.5f}{Di_valve:>10.5f}")
        out["cases"].append(dict(gap=gap, cell_length=info["cell_length"],
                                 fwd=rf, rev=rr, straight=rs,
                                 Di_cell=Di_cell, Di_valve=Di_valve,
                                 dp_str_over_dp_F=ratio))
    if len(out["cases"]) > 1:
        dv = [c["Di_valve"] for c in out["cases"]]
        dc = [c["Di_cell"] for c in out["cases"]]
        print(f"\n  Di_valve の gap 依存: {min(dv):.5f} .. {max(dv):.5f}  "
              f"差 {max(dv)-min(dv):+.5f}")
        print(f"  Di_cell  の gap 依存: {min(dc):.5f} .. {max(dc):.5f}  "
              f"差 {max(dc)-min(dc):+.5f}（希釈されるので gap 依存で当然）")
        print(f"  -> 補正の成立: "
              f"{'合格' if max(dv)-min(dv) < 0.005 else '要確認'}")
        out["Di_valve_spread"] = float(max(dv) - min(dv))

    os.makedirs("results", exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1, default=str)
    print(f"\nsaved {a.out}")


if __name__ == "__main__":
    main()
