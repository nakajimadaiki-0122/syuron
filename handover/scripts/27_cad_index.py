#!/usr/bin/env python3
"""
27_cad_index.py -- cad/ の中身を一覧にする
==========================================

`cad/*_fluid.stl` を読んで、**ファイル名と形の対応が一目で分かる図**と
一覧表（`cad/INDEX.md`）を作る。ファイル名だけでは形が分からない、という
問題への対処。

命名規則
--------

    {形状}_{こぶの大きさ}_{役割}.{拡張子}

| 形状 | 意味 |
|---|---|
| `valve1` | Gamboa のバルブ 1 個（プレナム込み。CFD の計算領域そのもの） |
| `valve1bare` | 同、プレナムを外したバルブ本体だけ |
| `stair{n}` | 段々（n 段。1 段おきに反転して平均方向を水平にした構成） |
| `side{n}fwd` | 直線流路の側面にこぶ n 個（こぶが順方向を向く） |
| `side{n}rev` | 同、こぶが逆方向を向く（x 反転） |

| 大きさ | ループ外半径 R |
|---|---|
| `small` | 1.9 w_v |
| `std` | 2.35 w_v（Gamboa Table 2 の最適値） |
| `large` | 3.0 w_v |

| 役割 | 中身 |
|---|---|
| `_outline.dxf` | 2D 輪郭（外形 + 島）。CAD で押し出す用 |
| `_fluid.stl` | 流体体積。固体から引く用 |
| `_plate.stl` | 溝を彫った板（裏面にポート穴）。実物モデル |
| `_lid.stl` | 蓋（平板） |
"""
import os, struct, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import plotstyle
CAD = os.path.join(HERE, "..", "cad")

FAMILY = {
    "valve1bare": "バルブ 1 個（プレナム無し）",
    "valve1": "バルブ 1 個（プレナム込み、CFD 領域）",
    "stair": "段々",
    "sidefwd": "直線 + 側面こぶ（順方向）",
    "siderev": "直線 + 側面こぶ（逆方向）",
}
SIZE_R = {"small": 1.9, "std": 2.35, "large": 3.0}


def read_stl(path):
    b = open(path, "rb").read()
    n = struct.unpack("<I", b[80:84])[0]
    a = np.frombuffer(b[84:], dtype=np.uint8).reshape(n, 50)
    v = a[:, 12:48].copy().view("<f4").reshape(n, 3, 3)
    return v


def describe(name):
    """ファイル名から（形状の説明, 大きさ）を作る。"""
    stem = name.rsplit("_", 1)[0]
    parts = stem.split("_")
    fam, size = parts[0], (parts[1] if len(parts) > 1 else "")
    if fam.startswith("stair"):
        return f"段々 {fam[5:]} 段", size
    if fam.startswith("side"):
        n = "".join(c for c in fam[4:] if c.isdigit())
        d = "順方向" if fam.endswith("fwd") else "逆方向"
        return f"直線 + 側面こぶ {n} 個（{d}）", size
    return FAMILY.get(fam, fam), size


def main():
    fnt = plotstyle.use_jp()
    print(f"図のフォント: {fnt}")
    files = sorted(f for f in os.listdir(CAD) if f.endswith("_fluid.stl"))
    if not files:
        print("cad/ に *_fluid.stl がない")
        return
    print(f"{'ファイル':32s}{'形状':30s}{'大きさ':8s}{'外形 [mm]':>18}{'体積 [mm^3]':>13}")
    rows = []
    shapes = []
    for f in files:
        v = read_stl(os.path.join(CAD, f))
        P = v.reshape(-1, 3)
        lo, hi = P.min(0), P.max(0)
        vol = float(np.einsum("ij,ij->i", v[:, 0],
                              np.cross(v[:, 1], v[:, 2])).sum() / 6.0)
        desc, size = describe(f)
        R = SIZE_R.get(size)
        rows.append((f, desc, size, R, hi - lo, vol))
        bot = v[np.abs(v[:, :, 2]).max(1) == 0][:, :, :2]
        shapes.append((f, desc, size, bot, lo, hi))
        print(f"{f:32s}{desc:30s}{size:8s}"
              f"{f'{hi[0]-lo[0]:.1f} x {hi[1]-lo[1]:.1f} x {hi[2]-lo[2]:.1f}':>18}"
              f"{vol:13.2f}")

    # --- 一覧図 ---
    n = len(shapes)
    ncol = 2
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 1.55 * nrow))
    for ax, (f, desc, size, bot, lo, hi) in zip(np.ravel(axes), shapes):
        ax.add_collection(PolyCollection(bot, facecolors="#9fc9cb",
                                         edgecolors="none"))
        ax.set_xlim(lo[0] - 2, lo[0] + 92)
        ax.set_ylim(lo[1] - 1, hi[1] + 1)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#d5dee1")
        ax.set_title(f"{f}   —   {desc}"
                     + (f" / {size}" if size else "")
                     + f"   {hi[0]-lo[0]:.1f} x {hi[1]-lo[1]:.1f} mm",
                     fontsize=8, loc="left", pad=2)
    for ax in np.ravel(axes)[n:]:
        ax.axis("off")
    fig.suptitle("cad/ の流体形状（_fluid.stl）— 名前と形の対応", fontsize=11)
    fig.tight_layout()
    png = os.path.join(HERE, "..", "figures", "cad_index.png")
    fig.savefig(png, dpi=140)
    print(f"\nsaved {os.path.normpath(png)}")

    # --- 一覧表 ---
    md = [
        "# cad/ の一覧",
        "",
        "`scripts/27_cad_index.py` が生成。図は `figures/cad_index.png`。",
        "",
        "## 命名規則",
        "",
        "    {形状}_{こぶの大きさ}_{役割}.{拡張子}",
        "",
        "| 形状 | 意味 |",
        "|---|---|",
        "| `valve1` | Gamboa のバルブ 1 個（プレナム込み。CFD の計算領域そのもの） |",
        "| `valve1bare` | 同、プレナムを外したバルブ本体だけ |",
        "| `stair{n}` | 段々（n 段。1 段おきに反転して平均方向を水平にした構成） |",
        "| `side{n}fwd` | 直線流路の側面にこぶ n 個（こぶが順方向を向く） |",
        "| `side{n}rev` | 同、こぶが逆方向を向く（x 反転） |",
        "",
        "| 大きさ | ループ外半径 R | こぶ頂部 |",
        "|---|---:|---:|",
        "| `small` | 1.9 w_v | 4.3 w_v |",
        "| `std` | 2.35 w_v（Gamboa Table 2） | 5.2 w_v |",
        "| `large` | 3.0 w_v | 6.5 w_v |",
        "",
        "| 役割 | 中身 |",
        "|---|---|",
        "| `_outline.dxf` | 2D 輪郭（外形 + 島）。CAD で押し出す用 |",
        "| `_fluid.stl` | 流体体積。固体から引く用 |",
        "| `_plate.stl` | 溝を彫った板（裏面にポート穴）。実物モデル |",
        "| `_lid.stl` | 蓋（平板） |",
        "",
        "## 形状の一覧",
        "",
        "| ファイル | 形状 | 大きさ | R [w_v] | 外形 [mm] | 流体体積 [mm³] |",
        "|---|---|---|---:|---|---:|",
    ]
    for f, desc, size, R, ext, vol in rows:
        md.append(f"| `{f}` | {desc} | {size} | "
                  f"{R if R else '—'} | "
                  f"{ext[0]:.1f} × {ext[1]:.1f} × {ext[2]:.1f} | {vol:.2f} |")
    md += ["",
           "各形状には同じ名前で `_outline.dxf` と `_plate.stl` がある",
           "（`valve1_std` には `_lid.stl` も）。",
           "",
           "寸法は流路幅 w_v = 1.0 mm、深さ 1.0 mm。変えるときは生成し直す。",
           "",
           "```",
           "scripts/15_gamboa_full_geom.py --wv-mm 1.0      # valve1 の DXF",
           "scripts/20_make_3d_model.py    --wv-mm 1.0      # valve1 の STL",
           "scripts/22_make_chain.py -n 8 --size std        # 段々（任意の段数）",
           "scripts/26_make_variants.py                     # 変種一式（大きさ 3 通り）",
           "scripts/27_cad_index.py                         # この一覧を作り直す",
           "```",
           ""]
    idx = os.path.join(CAD, "INDEX.md")
    open(idx, "w", encoding="utf-8").write("\n".join(md))
    print(f"saved {os.path.normpath(idx)}")


if __name__ == "__main__":
    main()
