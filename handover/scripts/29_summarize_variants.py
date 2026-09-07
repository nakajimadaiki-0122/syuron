#!/usr/bin/env python3
"""
29_summarize_variants.py -- 変種の計算結果をまとめる
====================================================

`28_run_variants.py` が出した JSON（`results/single/variant_*.json`）を集め、
表と図にする。2 次元と 3 次元、こぶの大きさ、形状の違いを並べて見るためのもの。

出力
  results/single/variants_summary.json
  figures/flow/variants_di.png
"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import plotstyle

RES = os.path.join(HERE, "..", "results", "single")
SIZE_ORDER = ["small", "std", "large"]
KIND_LABEL = {"stair4": "段々 4 段", "side4": "直線 + 側面こぶ 4 個"}


def collect():
    """向き別のファイルもまとめて、形状ごとに 1 件へ集約する。"""
    rows = {}
    for p in sorted(glob.glob(os.path.join(RES, "variant_*.json"))):
        d = json.load(open(p))
        key = (d["shape"], d["dim"], int(d["Re"]), d["cpm"])
        r = rows.setdefault(key, dict(shape=d["shape"], dim=d["dim"],
                                      Re=d["Re"], cpm=d["cpm"],
                                      geometry=d.get("geometry", {})))
        for lab in ("fwd", "rev"):
            if lab in d:
                r[lab] = d[lab]
    out = []
    for r in rows.values():
        if "fwd" in r and "rev" in r:
            r["Di"] = r["rev"]["dp_Pa"] / r["fwd"]["dp_Pa"]
            # 誤差伝播（各 Δp の標準偏差から）
            f, rv = r["fwd"], r["rev"]
            if f["dp_Pa"] > 0:
                r["Di_std"] = r["Di"] * np.hypot(
                    rv["dp_Pa_std"] / max(rv["dp_Pa"], 1e-12),
                    f["dp_Pa_std"] / max(f["dp_Pa"], 1e-12))
        out.append(r)
    return sorted(out, key=lambda r: (r["dim"], r["shape"]))


def main():
    plotstyle.use_jp()
    rows = collect()
    if not rows:
        print("結果がまだ無い")
        return
    print(f"{'形状':16s}{'次元':>4}{'格子':>6}{'Δp 順 [Pa]':>14}{'Δp 逆 [Pa]':>14}"
          f"{'Di':>9}{'収束':>16}")
    for r in rows:
        f, rv = r.get("fwd"), r.get("rev")
        di = r.get("Di")
        st = []
        for lab, x in (("順", f), ("逆", rv)):
            if x is None:
                st.append(f"{lab}:未")
            elif not x["converged"]:
                st.append(f"{lab}:未収束")
            elif x["unsteady"]:
                st.append(f"{lab}:非定常")
            else:
                st.append(f"{lab}:OK")
        dpf = f"{f['dp_Pa']:.3f}" if f else "—"
        dpr = f"{rv['dp_Pa']:.3f}" if rv else "—"
        dis = f"{di:.4f}" if di else "—"
        print(f"{r['shape']:16s}{r['dim']:>4}{r['cpm']:>6}{dpf:>14}{dpr:>14}"
              f"{dis:>9}{' '.join(st):>16}")

    j = os.path.join(RES, "variants_summary.json")
    json.dump([{k: v for k, v in r.items() if k != "geometry"} for r in rows],
              open(j, "w"), indent=1, default=float)
    print(f"\nsaved {os.path.normpath(j)}")

    # --- 図: こぶの大きさ 対 Di ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, dim in zip(axes, (2, 3)):
        any_pt = False
        for kind, mark in (("stair4", "o-"), ("side4", "s--")):
            xs, ys, es = [], [], []
            for i, size in enumerate(SIZE_ORDER):
                m = [r for r in rows if r["dim"] == dim
                     and r["shape"] == f"{kind}_{size}" and "Di" in r]
                if m:
                    xs.append(i)
                    ys.append(m[0]["Di"])
                    es.append(m[0].get("Di_std", 0.0))
            if xs:
                any_pt = True
                ax.errorbar(xs, ys, yerr=es, fmt=mark, capsize=3,
                            label=KIND_LABEL[kind])
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xticks(range(3))
        ax.set_xticklabels([f"{s}\nR={r}" for s, r in
                            zip(SIZE_ORDER, (1.9, 2.35, 3.0))])
        ax.set_ylabel("Di = Δp 逆 / Δp 順")
        ax.set_title(f"{dim} 次元" + ("" if any_pt else "（結果なし）"))
        ax.grid(alpha=0.3)
        if any_pt:
            ax.legend(fontsize=9)
    fig.suptitle("こぶの大きさと方向依存性（Re = 100）", fontsize=11)
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "flow", "variants_di.png")
    fig.savefig(png, dpi=130)
    print(f"saved {os.path.normpath(png)}")


if __name__ == "__main__":
    main()
