#!/usr/bin/env python3
"""
07_uncertainty.py -- Di_p の不確かさを格子と収束判定から評価する
================================================================

比較する 3 構成
  A: cpm = 12, tol = 1e-9   （基準）
  B: cpm = 24, tol = 1e-9   （格子を 2 倍）
  C: cpm = 12, tol = 1e-10  （収束判定を 1 桁厳しく）

Di_p = Δp_逆 / Δp_順。周期単位セルでは Δp = G * L なので Di_p = G_逆 / G_順。
セル長 L は順逆で同じなので、Di_p は L に依存しない。
"""
import glob, json, os
import numpy as np


def load(Re, cpm, suffix=""):
    out = {}
    for case in ("fwd", "rev"):
        pat = f"results/periodic/cell_tesla_channel_cell_{case}_Re{Re}_cpm{cpm}_o*{suffix}.json"
        hits = [f for f in glob.glob(pat)
                if f.endswith(f"{suffix}.json") and
                ("_tight" in f) == (suffix == "_tight")]
        if not hits:
            return None
        out[case] = json.load(open(sorted(hits)[0]))
    return out


def main():
    CONF = [("A", 12, "", "cpm=12, tol=1e-9  (基準)"),
            ("B", 24, "", "cpm=24, tol=1e-9  (格子 2 倍)"),
            ("C", 12, "_tight", "cpm=12, tol=1e-10 (判定 1 桁厳格)")]
    table = {}
    for Re in (100, 500):
        print("=" * 92)
        print(f" Re = {Re}")
        print("=" * 92)
        print(f"{'':4}{'構成':<32}{'dp_F [Pa]':>12}{'dp_R [Pa]':>12}"
              f"{'Di_p':>10}{'iters_F':>9}{'dG_F':>9}")
        for key, cpm, suf, label in CONF:
            d = load(Re, cpm, suf)
            if d is None:
                print(f"{key:<4}{label:<32}{'（未完了）':>12}")
                continue
            F, R = d["fwd"], d["rev"]
            Di = R["dp_Pa"] / F["dp_Pa"]
            table[(Re, key)] = dict(Di=Di, dpF=F["dp_Pa"], dpR=R["dp_Pa"],
                                    label=label)
            print(f"{key:<4}{label:<32}{F['dp_Pa']:>12.6f}{R['dp_Pa']:>12.6f}"
                  f"{Di:>10.6f}{F['iters']:>9d}{F['dG_final']:>9.1e}")

        have = [k for k in "ABC" if (Re, k) in table]
        if "A" in have:
            print()
            base = table[(Re, "A")]["Di"]
            for k in have[1:]:
                d = table[(Re, k)]["Di"] - base
                print(f"    Di_p({k}) - Di_p(A) = {d:+.6f}"
                      f"   （Δp_F は {100*(table[(Re,k)]['dpF']/table[(Re,'A')]['dpF']-1):+.2f} % 変化）")
        print()

    # ---- 不確かさの合成 ----
    print("=" * 92)
    print(" Di_p の不確かさ")
    print("=" * 92)
    for Re in (100, 500):
        if (Re, "A") not in table:
            continue
        base = table[(Re, "A")]["Di"]
        u_grid = abs(table[(Re, "B")]["Di"] - base) if (Re, "B") in table else np.nan
        u_conv = abs(table[(Re, "C")]["Di"] - base) if (Re, "C") in table else np.nan
        u = float(np.hypot(u_grid, u_conv))
        print(f"  Re = {Re:4d}:  Di_p = {base:.5f}"
              f"   格子由来 {u_grid:.5f}   収束由来 {u_conv:.5f}"
              f"   合成 ±{u:.5f}")
        table[(Re, "u")] = u

    if (100, "u") in table and (500, "u") in table:
        d100, d500 = table[(100, "A")]["Di"], table[(500, "A")]["Di"]
        u100, u500 = table[(100, "u")], table[(500, "u")]
        span = abs(d500 - d100)
        us = float(np.hypot(u100, u500))
        print()
        print(f"  Re=100 から Re=500 への Di_p の変化量: {d500-d100:+.5f}")
        print(f"  その不確かさ（合成）:                  ±{us:.5f}")
        print(f"  信号対不確かさ比:                      {span/us:.2f}")
        print()
        if span > 3 * us:
            print("  → 変化は不確かさに対して有意。Di_p の Re 依存性を主張できる。")
        else:
            print("  → 変化が不確かさに埋もれる。Di_p の Re 依存性は主張できない。")
        print()
        print(f"  参考: 判定の目安とされた 0.005 との比較")
        for Re in (100, 500):
            u = table[(Re, "u")]
            v = "十分小さい" if u < 0.005 / 3 else ("同程度" if u < 0.005 else "大きすぎる")
            print(f"    Re = {Re:4d}: ±{u:.5f}  → 0.005 に対して {v}")

    os.makedirs("results/periodic", exist_ok=True)
    json.dump({f"Re{k[0]}_{k[1]}": (v if not isinstance(v, float) else v)
               for k, v in table.items()},
              open("results/periodic/di_uncertainty.json", "w"), indent=1, default=str)
    print("\nsaved results/di_uncertainty.json")


if __name__ == "__main__":
    main()
