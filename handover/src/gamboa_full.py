"""
gamboa_full.py -- Gamboa 2005 の**単発バルブ全体**（プレナム込み）の形状生成
=============================================================================

`gamboa.py` との違い（重要）
----------------------------
`gamboa.py` は周期配列にするために出口区間を「閉じ壁」に置き換え、主流路を
接合部の先へまっすぐ延長していた。2026-09-06 に Fig.2 を実測した結果、
**この置き換えが Tesla バルブの中核機構を壊していた**ことが判明した。

  実形状（Fig.2 実測、`scripts/13_fig2_topology.py`）
    - 主流路は接合部（x ≈ 7.9 w_v）で**終わる**。まっすぐ先へ続かない
    - ループ下流枝と出口区間は**同一直線**
      （右壁の角度差 0.026 度、法線オフセット 0.003 w_v、残差 0.18 px）
    - したがって逆流時、出口区間を上ってきた流れは**折れ角ゼロでループへ入る**

  周期版（gamboa.py）
    - 逆流はまっすぐ主流路を素通りでき、ループへ入るには 132 度曲がる必要がある
    - 実測 φ_loop は順流 4.0 % > 逆流 2.9 %（Re=500）で、設計原理と逆
    - Di_valve = 0.9925（Re=100）、0.9678（Re=500）。Gamboa は 1.02、1.37

本モジュールは Fig.1 / Fig.2 どおりの単発構成
（左プレナム → 入口区間 → ループ → 出口区間 → 右プレナム）を組む。

作図規則の来歴
--------------
[論文記載] Table 2: X2, n, Y3, R, alpha, LENOUT。Fig.1: R_p = 5、原点の定義。
[確認済み事実] 以下は Fig.2 の画素実測（`13_fig2_topology.py`, `14_fig2_plenum.py`）。

  - beta = theta - 90（主流路法線から測る）。Fig.2 実測 161.245 度 対 Table 2 導出 161.607 度
  - alpha も**法線から**測る。Fig.2 実測: 出口区間は主流路軸から 48.03 度
    = 法線から 41.97 度（Table 2: 41.9）
  - 外側円弧は y = +w/2 と戻り流路外壁の両方に接する。中心 (1.98, 2.85)、
    Fig.2 実測 (2.08, 2.86)、円フィット R = 2.32〜2.36（Table 2: 2.35）
  - 左プレナム: 中心 (-5.464, -0.029) w_v、半径 5.060 w_v。
    幾何から予想される中心 (-5.475, 0)・R_p = 5 と一致（**スケールの独立検証**）
  - 右プレナム: 中心 (12.849, -6.275) w_v、半径 5.166 w_v。
    出口区間の中心線上にある（法線距離 0.071 w_v）

[未確定] LENOUT の起終点。出口中心線が主流路中心線 y = 0 と交わる点 J から
右プレナムの円弧までは **3.27 w_v** で、Table 2 の 2.94 と 11 % 違う。
本モジュールは Fig.2 実測を優先し、`lenout_axis` 既定値を 3.274 とする。
差は直線区間 0.33 w_v ぶんで、順逆どちらにも同じだけ効く。

CFD 用の逸脱（`plenum_out="vertical"`、既定）
--------------------------------------------
Gamboa の右プレナムは直線壁が**出口区間に垂直**（水平から 42 度傾く）。
Cartesian 格子の LBM で傾いた境界面に入口・出口 BC を課すのは階段誤差が大きい。
そこで右プレナムを「半径 R_p、中心は同じ、直線壁は**鉛直**」の半円に置き換える。

  - 流入口（円弧上の点）から直線壁までの流路長は元と同じ R_p
  - 面積・スケールも同じ。変わるのは半円をどの向きに切るかだけ
  - 順流・逆流とも BC 面が鉛直になるので、x 方向反転だけで逆流が解ける

`plenum_out="faithful"` にすると Fig.2 どおり（傾いた直線壁）を返す。CAD 用。
"""
from __future__ import annotations
import numpy as np
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union

# Table 2（論文記載値、w_v で無次元化）
OPTIMIZED = dict(X2=1.60, n=0.797, Y3=0.608, R=2.35, alpha=41.9, LENOUT=2.94)
REFERENCE = dict(X2=1.50, n=0.990, Y3=0.600, R=2.50, alpha=45.0, LENOUT=2.00)

RP = 5.0                # [論文記載] プレナム半径 R_p = 5 w_v
LENOUT_AXIS_FIG2 = 3.274   # [確認済み事実] Fig.2 実測（J から円弧まで、軸に沿う）
ARC_PTS = 720


def theta_from_params(X2, n, Y3, w=1.0):
    """A=(X2, w/2) から B=(n*X2, Y3) への向き [deg]。beta = theta - 90。"""
    return float(np.degrees(np.arctan2(Y3 - 0.5 * w, n * X2 - X2)))


def _skeleton(p, w=1.0, arc_pts=ARC_PTS):
    """バルブの骨格点を返す。gamboa.py の `_build` と同じ規則。"""
    h = 0.5 * w
    th = np.radians(theta_from_params(p["X2"], p["n"], p["Y3"], w))
    d = np.array([np.cos(th), np.sin(th)])          # 主流路 -> ループ
    nrm = np.array([-np.sin(th), np.cos(th)])
    R = p["R"] * w
    A = np.array([p["X2"] * w, h])

    Cx = p["X2"] * w + R * (1.0 + np.cos(th)) / np.sin(th)
    C = np.array([Cx, h + R])
    T1 = A + float(np.dot(C - A, d)) * d

    al = np.radians(p["alpha"])
    e = np.array([np.sin(al), -np.cos(al)])         # 出口区間の向き（下流・下向き）
    m = np.array([np.cos(al), np.sin(al)])          # e を +90 度回した法線
    T2 = C + R * m

    # 外側円弧 T1 -> T2。y=h との接点（角度 -90 度）を通らない側を選ぶ。
    a1 = np.arctan2(*(T1 - C)[::-1])
    a2 = np.arctan2(*(T2 - C)[::-1])
    ccw = a2 + 2 * np.pi if a2 <= a1 else a2
    cw = a2 - 2 * np.pi if a2 >= a1 else a2

    def passes_tangent(a_from, a_to, target=-0.5 * np.pi):
        for j in range(-2, 3):
            t = target + 2 * np.pi * j
            if min(a_from, a_to) < t < max(a_from, a_to):
                return True
        return False

    ang = (np.linspace(a1, cw, arc_pts) if passes_tangent(a1, ccw)
           else np.linspace(a1, ccw, arc_pts))
    arc = C + R * np.c_[np.cos(ang), np.sin(ang)]

    # 外壁 L_out: T2 を通り方向 e。y=+h と y=-h を横切る点
    E = T2 + (h - T2[1]) / e[1] * e            # 主流路上壁の高さ
    Q = T2 + (-h - T2[1]) / e[1] * e           # 主流路下壁の高さ（= 接合部の角）

    # 島（内側輪郭）: 外側を w だけ内へオフセット
    ri = R - w
    T1i = C + ri * (T1 - C) / R
    T2i = C + ri * (T2 - C) / R
    arci = C + ri * np.c_[np.cos(ang), np.sin(ang)]
    P1 = T1i - (h - T1i[1]) / (-d[1]) * d       # 戻り流路内壁が y=+h に達する点
    P2 = T2i + (h - T2i[1]) / e[1] * e          # 島の下流端（尖点）

    # 出口区間の中心線が y=0 と交わる点 J（LENOUT の起点として使う）
    cmid = P2 + 0.5 * w * m
    J = cmid + (0.0 - cmid[1]) / e[1] * e
    return dict(h=h, th=th, d=d, nrm=nrm, e=e, m=m, R=R, C=C, A=A,
                T1=T1, T2=T2, E=E, Q=Q, P1=P1, P2=P2, J=J,
                arc=arc, arci=arci, island=np.vstack([P1[None], arci, P2[None]]))


def build(which="optimized", w=1.0, lenout_axis=LENOUT_AXIS_FIG2, Rp=RP,
          plenum_out="vertical", arc_pts=ARC_PTS, include_plenums=True, **over):
    """単発バルブ全体の流体領域を返す。

    Parameters
    ----------
    which       : "optimized" | "reference"
    w           : 流路幅（出力座標の単位。w=1 なら無次元、w=1e-3 なら m）
    lenout_axis : 出口中心線が y=0 と交わる点 J から右プレナム円弧までの長さ [w_v]
    Rp          : プレナム半径 [w_v]
    plenum_out  : "vertical"（CFD 用、直線壁が鉛直）| "faithful"（Fig.2 どおり）
    include_plenums : False にするとプレナムを外し、入口流路の左端 x = -0.5 w_v から
                  出口区間がプレナム円弧に入る点までのバルブ本体だけを返す（CAD 用）

    Returns
    -------
    poly : shapely Polygon（島は内側リング）
    info : dict（骨格点、BC 面の位置など）
    """
    p = dict(OPTIMIZED if which.lower().startswith("opt") else REFERENCE)
    p.update({k: v for k, v in over.items() if k in p})
    s = _skeleton(p, w=w, arc_pts=arc_pts)
    h, e, m = s["h"], s["e"], s["m"]
    Rp = Rp * w
    lo = lenout_axis * w

    # --- 左プレナム: 中心は入口流路中心線上、円弧が (x=-0.5w, y=±h) を通る ---
    Lx = -0.5 * w - np.sqrt(Rp ** 2 - h ** 2)
    Lc = np.array([Lx, 0.0])
    plenL = (Point(*Lc).buffer(Rp, resolution=256)
             .intersection(Polygon([(Lx, -Rp - w), (Lx + Rp + w, -Rp - w),
                                    (Lx + Rp + w, Rp + w), (Lx, Rp + w)])))

    # --- 右プレナム: 中心は出口中心線上、円弧が J から lenout_axis の点を通る ---
    entry = s["J"] + lo * e
    Rc = entry + Rp * e
    disc = Point(*Rc).buffer(Rp, resolution=256)
    if plenum_out == "vertical":
        # 直線壁は鉛直（x = Rc_x）。流入口のある左半分を残す。
        half = Polygon([(Rc[0] - Rp - w, Rc[1] - Rp - w), (Rc[0], Rc[1] - Rp - w),
                        (Rc[0], Rc[1] + Rp + w), (Rc[0] - Rp - w, Rc[1] + Rp + w)])
        out_plane = ("x", float(Rc[0]), float(Rc[1] - Rp), float(Rc[1] + Rp))
    elif plenum_out == "faithful":
        # 直線壁は出口区間に垂直で中心を通る。**流入口のある上流側**の半分を残す
        # （直線壁が出口 BC 面なので、流入口から壁までが流路になる）。
        q = Rc + (Rp + w) * m
        r = Rc - (Rp + w) * m
        half = Polygon([q, r, r - (Rp + w) * e, q - (Rp + w) * e])
        out_plane = ("perp", float(Rc[0]), float(Rc[1]), float(p["alpha"]))
    else:
        raise ValueError(plenum_out)
    plenR = disc.intersection(half)

    # --- 主流路（右端は外壁 L_out で斜めに切られる）---
    x_start = Lx if include_plenums else -0.5 * w
    chan = Polygon([(x_start, -h), s["Q"], s["E"], (x_start, h)])

    # --- ループ（y=+h より上）---
    loop = Polygon(np.vstack([s["A"][None], s["arc"], s["E"][None],
                              [[s["E"][0], h]], [[s["A"][0], h]]]))

    # --- 出口区間（P2 から e 方向へ、幅 w の帯）---
    # 帯はプレナムの内部で止める。長すぎるとプレナムを突き抜けて外へ伸びる。
    P2 = s["P2"]
    t_J = float(np.dot(s["J"] - P2, e))       # P2 から J までの軸方向距離
    Lband = t_J + lo + 0.5 * Rp
    band = Polygon([P2, P2 + w * m, P2 + w * m + Lband * e, P2 + Lband * e])

    parts = ([plenL, chan, loop, band, plenR] if include_plenums
             else [chan, loop, band_valve_only(P2, m, e, w, t_J, lo)])
    poly = unary_union(parts)
    poly = poly.difference(Polygon(s["island"]))
    if poly.geom_type != "Polygon":
        raise ValueError(f"流体領域が単一多角形にならない: {poly.geom_type}")

    info = dict(params=p, w=float(w), which=which,
                theta_deg=float(np.degrees(s["th"])),
                beta_deg=float(np.degrees(s["th"]) - 90.0),
                alpha_deg=float(p["alpha"]),
                centre=s["C"].tolist(), splitter_tip=s["A"].tolist(),
                arc_tangent=s["T1"].tolist(), outlet_tangent=s["T2"].tolist(),
                E=s["E"].tolist(), Q=s["Q"].tolist(),
                island_tip=s["P2"].tolist(), island_heel=s["P1"].tolist(),
                J=s["J"].tolist(), outlet_entry=entry.tolist(),
                plenum_left_centre=Lc.tolist(), plenum_right_centre=Rc.tolist(),
                Rp=float(Rp), lenout_axis=float(lenout_axis),
                inlet_plane=("x", float(Lx), float(-Rp), float(Rp)),
                outlet_plane=out_plane, plenum_out=plenum_out,
                bounds=[float(v) for v in poly.bounds])
    return poly, info


def band_valve_only(P2, m, e, w, t_J, lo):
    """プレナムを外すとき用: 出口区間を円弧の入口までで切った帯。"""
    L = t_J + lo
    return Polygon([P2, P2 + w * m, P2 + w * m + L * e, P2 + L * e])


def rasterize(poly, cells_per_w, w=1.0, pad=2):
    """流体マスクを作る。

    x 方向の両端は入口・出口の BC 面なので流体でよい。
    y 方向の両端は必ず固体にする（np.roll による周期化を防ぐ。README §8 B-1）。

    Returns
    -------
    mask : (nx, ny) bool
    xs, ys : セル中心座標
    """
    from matplotlib.path import Path
    xmin, ymin, xmax, ymax = poly.bounds
    res = cells_per_w / w
    nx = int(np.ceil((xmax - xmin) * res))
    ny = int(np.ceil((ymax - ymin) * res))
    xs = xmin + (np.arange(nx) + 0.5) / res
    ys = ymin + (np.arange(ny) + 0.5) / res
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel()])
    inside = Path(np.asarray(poly.exterior.coords)).contains_points(pts)
    for ring in poly.interiors:
        inside &= ~Path(np.asarray(ring.coords)).contains_points(pts)
    core = inside.reshape(nx, ny)
    # 端の空列・空行を落とす。こうすると i=0 と i=nx-1 が必ず BC 面
    # （プレナムの直線壁に最も近いセル中心）になる。
    ci = np.flatnonzero(core.any(1))
    cj = np.flatnonzero(core.any(0))
    core = core[ci.min():ci.max() + 1, cj.min():cj.max() + 1]
    xs = xs[ci.min():ci.max() + 1]
    ys = ys[cj.min():cj.max() + 1]
    nx, ny = core.shape
    mask = np.zeros((nx, ny + 2 * pad), bool)
    mask[:, pad:-pad] = core
    ys = np.concatenate([ys[0] - np.arange(pad, 0, -1) / res, ys,
                         ys[-1] + np.arange(1, pad + 1) / res])
    return mask, xs, ys


def straight_duct(length, w=1.0):
    """検証用の直線流路（同じ BC で f*Re = 96 を確認するため）。"""
    return Polygon([(0, -0.5 * w), (length, -0.5 * w),
                    (length, 0.5 * w), (0, 0.5 * w)])
