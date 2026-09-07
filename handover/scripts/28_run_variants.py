#!/usr/bin/env python3
"""
28_run_variants.py -- 試作用の変種を解く（2 次元 / 3 次元）
===========================================================

`26_make_variants.py` で作った形状の性能を計算する。

| 形状 | 意味 |
|---|---|
| `stair4_{small,std,large}` | 段々 4 段 |
| `side4_{small,std,large}` | 直線流路の側面にこぶ 4 個 |

**`side4rev` は `side4fwd` の鏡像**なので、`side4` を順流・逆流の両方で解けば
両方の向きを覆える（`side4fwd` の逆流 = `side4rev` の順流）。

両端に**直線壁が鉛直な半円プレナム**（R_p = 5 w_v）を付けてから解く。
段々流路は端の向きが 24 度傾いているので、これが無いと鉛直な BC 面を作れない。
単発モデル（`17_run_gamboa_full.py`）と同じ構成なので、Δp は直接比べられる。

Re は Gamboa の定義（D_h = 2 w_v）。3 次元は同じ条件で深さを有限にする
（正方形断面としての Re はその半分）。

    python scripts/28_run_variants.py --shape side4_std --dim 2 --dir fwd
"""
import argparse, json, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import gamboa as G
import gamboa_full as GF
import lbm
import lbm3d as L3
import postproc as pp

SIZES = {"small": 1.9, "std": 2.35, "large": 3.0}
OUT = os.path.join(HERE, "..", "results", "single")


def build(shape, w=1.0, gap=1.0, arc_pts=720):
    """変種の流体多角形（プレナム込み）と説明を返す。"""
    kind, size = shape.rsplit("_", 1)
    R = SIZES[size]
    if kind.startswith("stair"):
        n = int(kind[5:])
        poly, info = GF.chain(n, w=w, gap=0.0, arc_pts=arc_pts, R=R)
        d0 = np.array(info["stages"][0]["inlet_dir"])
        d1 = np.array(info["outlet_dir"])
        ends = [((-0.5 * w, 0.0), -d0), (info["outlet"], d1)]
        note = f"段々 {n} 段, R={R} w_v"
    elif kind.startswith("side"):
        n = int(kind[4:])
        poly, info = G.straight_with_loops(n, R=R, gap=gap, w=w,
                                           arc_pts=arc_pts)
        b = poly.bounds
        ends = [((b[0], 0.0), (-1.0, 0.0)), ((b[2], 0.0), (1.0, 0.0))]
        note = f"直線 + 側面こぶ {n} 個, R={R} w_v, gap={gap} w_v"
    else:
        raise ValueError(shape)
    full, planes = GF.attach_plenums(poly, ends, w=w)
    return full, dict(shape=shape, kind=kind, size=size, R=R, note=note,
                      bare_bounds=[float(v) for v in poly.bounds],
                      bounds=[float(v) for v in full.bounds], planes=planes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", required=True,
                    help="stair4_std / side4_large など")
    ap.add_argument("--dim", type=int, default=2, choices=[2, 3])
    ap.add_argument("--Re", type=float, default=100.0)
    ap.add_argument("--cpm", type=int, default=12)
    ap.add_argument("--depth", type=int, default=None,
                    help="3 次元の深さ [セル]。既定は cpm（正方形断面）")
    ap.add_argument("--U", type=float, default=0.05)
    ap.add_argument("--dir", default="both", choices=["both", "fwd", "rev"])
    ap.add_argument("--iters", type=int, default=800000)
    ap.add_argument("--dp-window", type=int, default=10000)
    ap.add_argument("--gap", type=float, default=1.0)
    a = ap.parse_args()
    depth = a.depth if a.depth else a.cpm
    w_mm = 1.0

    poly, info = build(a.shape, gap=a.gap)
    m2, xs, ys = GF.rasterize(poly, a.cpm, w=1.0)
    mask = m2 if a.dim == 2 else L3.extrude(m2, depth)
    nu = lbm.nu_lattice(a.U, 2.0 * a.cpm, a.Re)

    print("=" * 88)
    print(f" {a.shape}（{info['note']}）  {a.dim} 次元  Re = {a.Re}, "
          f"cpm = {a.cpm}" + (f", 深さ {depth} セル" if a.dim == 3 else ""))
    print("=" * 88)
    print(f"  外形（プレナム込み） "
          f"{info['bounds'][2]-info['bounds'][0]:.2f} x "
          f"{info['bounds'][3]-info['bounds'][1]:.2f} mm")
    print(f"  格子 {' x '.join(str(v) for v in mask.shape)} = "
          f"{np.prod(mask.shape)/1e6:.2f} M、流体 {int(mask.sum())/1e3:.0f} k")
    print(f"  nu = {nu:.5f}, tau_p = {lbm.tau_from_nu(nu):.4f}", flush=True)

    out = dict(shape=a.shape, dim=a.dim, Re=a.Re, cpm=a.cpm, U=a.U,
               depth=(depth if a.dim == 3 else None), nu=nu, geometry=info,
               grid=[int(v) for v in mask.shape], fluid=int(mask.sum()))

    todo = [("fwd", mask), ("rev", mask[::-1].copy())]
    if a.dir != "both":
        todo = [t for t in todo if t[0] == a.dir]

    for lab, m in todo:
        n_in = int(m[0].sum())
        Q = a.U * a.cpm * (depth if a.dim == 3 else 1)
        U_in = Q / n_in
        print(f"  [{lab}] 入口 {n_in} セル, U_in = {U_in:.5f}", flush=True)

        def rep2(it, d, mi, mo, dp=None):
            if it % 20000 == 0:
                print(f"    it={it:7d} res={d:.2e} dp={dp:.6e}", flush=True)

        def rep3(it, d, dp, spread):
            if it % 5000 == 0:
                print(f"    it={it:7d} res={d:.2e} dp={dp:.6e} "
                      f"質量ばらつき={spread:.2e}", flush=True)

        t0 = time.time()
        if a.dim == 2:
            ux, uy, p, rho, f, inf = lbm.solve_flow_io(
                m, nu, U_in, iters=a.iters, dp_window=a.dp_window, report=rep2)
        else:
            ux, uy, uz, rho, p, f, inf = L3.solve_flow_io(
                m, nu, U_in, iters=a.iters, dp_window=a.dp_window, report=rep3)
        dp_lat = inf["dp_lattice_mean"]
        scale = pp.lattice_to_physical_scale(a.Re, 2.0 * w_mm * 1e-3, a.U)
        r = dict(label=lab, n_inlet=n_in, U_in=U_in, dp_lattice=dp_lat,
                 dp_Pa=pp.to_pascal(dp_lat, scale),
                 dp_Pa_std=pp.to_pascal(inf["dp_lattice_std"], scale),
                 converged=inf["converged"], unsteady=inf["unsteady"],
                 iters=inf["iters"], residual=inf["residual"],
                 mass_spread=inf.get("mass_spread_final", float("nan")),
                 wall_s=round(time.time() - t0, 1))
        out[lab] = r
        print(f"  [{lab}] dp = {r['dp_Pa']:.4f} ± {r['dp_Pa_std']:.4f} Pa, "
              f"収束 {r['converged']}"
              f"{'（非定常）' if r['unsteady'] else ''} it={r['iters']} "
              f"質量ばらつき={r['mass_spread']:.1e} "
              f"{r['wall_s']/60:.0f} 分", flush=True)

    base = f"variant_{a.shape}_{a.dim}d_Re{int(a.Re)}_cpm{a.cpm}"
    o = os.path.join(OUT, base + ("" if a.dir == "both" else f"_{a.dir}")
                     + ".json")
    if a.dir != "both":
        other = "rev" if a.dir == "fwd" else "fwd"
        op = os.path.join(OUT, base + f"_{other}.json")
        if os.path.exists(op):
            j = json.load(open(op))
            if other in j:
                out[other] = j[other]
    os.makedirs(OUT, exist_ok=True)
    if "fwd" in out and "rev" in out:
        out["Di"] = out["rev"]["dp_Pa"] / out["fwd"]["dp_Pa"]
        print(f"\n  Di = {out['rev']['dp_Pa']:.4f} / "
              f"{out['fwd']['dp_Pa']:.4f} = {out['Di']:.4f}")
    json.dump(out, open(o, "w"), indent=1, default=float)
    print(f"saved {os.path.normpath(o)}")


if __name__ == "__main__":
    main()
