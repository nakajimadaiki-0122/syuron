"""
gamboa.py -- Gamboa 2005 型 Tesla バルブの形状生成
==================================================

作図規則の来歴
--------------
すべて `docs/gamboa2005.pdf` の Fig.1（記号定義）と Fig.2（実形状）に基づく。
Fig.1 は模式図で Table 2 の最適値どおりに描かれていないため、寸法は Fig.1 から
逆算していない。未記載の作図規則は Fig.2 を画素計測して決めた（docs/gamboa2005_geometry.md §5）。

確定した規則

1. 主流路は y = -w/2 .. +w/2。順流は +x。
2. 分流板は**平行平板ではなく楔**。下面が主流路上壁 y = +w/2、上面が戻り流路の外壁で、
   先端 A = (X2, +w/2) で交わる。
3. 戻り流路の向き theta は A から B = (n*X2, Y3) へ向かう向き。
   Gamboa の beta は主流路法線から測るので beta = theta - 90。
   （Table 2 の 2 行を再現し、Fig.2 実測 161.245 度とも 0.36 度で一致）
4. ループ外側円弧（半径 R）は **y = +w/2 と戻り流路外壁の両方に接する**。
   これで中心が一意に決まる。中心の y は w/2 + R。
   （理論 Cx = 1.981、Fig.2 実測 2.075。図の線幅の範囲内）
5. 島（内側輪郭）は外側境界を w だけ内側へオフセットしたもの。
   左端は中心を共有する半径 R - w の円弧。
   （Fig.2 実測: 島左端の中心からの距離 1.386 w_v、理論 1.35）
6. alpha も**主流路法線から測る**。基準は出口区間の上流側の壁。
   （Fig.2 実測 41.91 度、Table 2 の 41.9 度と一致）

周期構成について
----------------
Gamboa のバルブは単発構成で、出口区間が主流路軸から 48.1 度傾いているため
そのままでは x 方向に並べられない。本モジュールは **出口区間を「閉じ壁」に置き換える**。
閉じ壁は Gamboa の出口区間上壁と同じ角度 alpha を保ったまま、y = +w/2 で主流路に戻す。
ループ下流側開口の局所形状は保存される。

この置き換えにより LENOUT は形状に現れない。Fig.7 と不一致だった場合の疑い先として
`docs/gamboa2005_geometry.md` §5.4 に記録してある。
"""
from __future__ import annotations
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

# Gamboa 2005 Table 2（論文記載値。線寸法は w_v で無次元化）
OPTIMIZED = dict(X2=1.60, n=0.797, Y3=0.608, R=2.35, alpha=41.9, LENOUT=2.94)
REFERENCE = dict(X2=1.50, n=0.990, Y3=0.600, R=2.50, alpha=45.0, LENOUT=2.00)

ARC_PTS = 720          # 円弧の分割数。全形状で共通にして再現性を担保する


def theta_from_params(X2: float, n: float, Y3: float, w: float = 1.0) -> float:
    """A=(X2, w/2) から B=(n*X2, Y3) への向き [deg]。beta = theta - 90。"""
    return float(np.degrees(np.arctan2(Y3 - 0.5 * w, n * X2 - X2)))


def _build(theta_deg: float, R: float, alpha_deg: float, X2: float, w: float):
    """作図の核。外側境界と島の頂点列を返す。

    戻り値
    ------
    outer : (m,2) 主流路上壁 y=+w/2 より上のループ外側境界。
            先頭が分流板先端 A、末尾が閉じ壁が y=+w/2 に達する点 E。
    isl   : (k,2) 島の閉多角形
    info  : dict
    """
    h = 0.5 * w
    th = np.radians(theta_deg)
    d = np.array([np.cos(th), np.sin(th)])          # 主流路 -> ループ（上流向き）
    nrm = np.array([-np.sin(th), np.cos(th)])       # d を +90 度回した法線
    A = np.array([X2, h])

    # --- 規則 4: 円弧中心は y=h+R 上にあり、戻り流路外壁に接する ---
    # (C - A)·nrm = -R  を Cx について解く
    Cx = X2 + R * (1.0 + np.cos(th)) / np.sin(th)
    C = np.array([Cx, h + R])
    T1 = A + float(np.dot(C - A, d)) * d             # 外壁と円弧の接点

    # --- 規則 6: 閉じ壁（Gamboa の出口区間上壁と同じ角度、法線から alpha） ---
    al = np.radians(alpha_deg)
    e = np.array([np.sin(al), -np.cos(al)])          # 下流・下向き
    m = np.array([np.cos(al), np.sin(al)])           # e を +90 度回した法線
    # 円弧の外側に接する直線: (P - C)·m = R
    T2 = C + R * m                                   # 接点
    # 閉じ壁が y = h に達する点
    t = (h - T2[1]) / e[1]
    E = T2 + t * e

    # --- 円弧（T1 -> T2）---
    # 2 通りの回り方のうち、y = h との接点（角度 -90 度）を**通らない**方を選ぶ。
    # 接点を通る側を選ぶとループが主流路壁と 1 点で接して幅がゼロになる。
    a1 = np.arctan2(*(T1 - C)[::-1])
    a2 = np.arctan2(*(T2 - C)[::-1])

    def _contains(a_from, a_to, target=-0.5 * np.pi):
        s = np.sign(a_to - a_from)
        k = (target - a_from) / (a_to - a_from) if a_to != a_from else 0.0
        # target を 2pi 周期でずらして区間に入るものがあるか
        for j in range(-2, 3):
            t = target + 2 * np.pi * j
            if min(a_from, a_to) < t < max(a_from, a_to):
                return True
        return False

    ccw = a2 + 2 * np.pi if a2 <= a1 else a2         # 反時計回り
    cw = a2 - 2 * np.pi if a2 >= a1 else a2          # 時計回り
    ang = (np.linspace(a1, cw, ARC_PTS) if _contains(a1, ccw)
           else np.linspace(a1, ccw, ARC_PTS))
    arc = C + R * np.c_[np.cos(ang), np.sin(ang)]

    outer = np.vstack([A[None, :], arc, E[None, :]])

    # --- 島: 外側境界を w だけ内側へオフセット ---
    Ai = A + w * nrm                                  # 戻り流路の内壁上の点
    ri = R - w
    T1i = C + ri * (T1 - C) / R
    T2i = C + ri * (T2 - C) / R
    arci = C + ri * np.c_[np.cos(ang), np.sin(ang)]
    # 内壁 2 本が y = h に達する点
    t1 = (h - T1i[1]) / (-d[1]);  P1 = T1i - t1 * d   # 戻り流路内壁を下流へ延長
    t2 = (h - T2i[1]) / e[1];     P2 = T2i + t2 * e   # 閉じ壁内壁を下流へ延長
    isl = np.vstack([P1[None, :], arci, P2[None, :]])

    info = dict(theta_deg=float(theta_deg), beta_deg=float(theta_deg - 90.0),
                R=float(R), alpha_deg=float(alpha_deg), X2=float(X2), w=float(w),
                center=C.tolist(), splitter_tip=A.tolist(),
                arc_tangent_main_wall=[float(Cx), float(h)],
                outer_close_x=float(E[0]),
                island_bottom_x=[float(P1[0]), float(P2[0])],
                loop_x=[float(Cx - R), float(E[0])])
    return outer, isl, info


def valve_polygon(theta_deg: float, R: float = OPTIMIZED["R"],
                  alpha_deg: float = OPTIMIZED["alpha"],
                  X2: float = OPTIMIZED["X2"], w: float = 1.0,
                  x0: float = None, x1: float = None):
    """1 段分の流体領域（主流路 + ループ）を Polygon で返す。

    x0, x1 は主流路の左右端。None ならループ範囲に 0.5*w の余白を付ける。
    """
    outer, isl, info = _build(theta_deg, R, alpha_deg, X2, w)
    h = 0.5 * w
    lx, rx = info["loop_x"]
    if x0 is None:
        x0 = lx - 0.5 * w
    if x1 is None:
        x1 = rx + 0.5 * w
    if x0 > lx or x1 < rx:
        raise ValueError(f"x0={x0:.3f}, x1={x1:.3f} がループ範囲 "
                         f"[{lx:.3f},{rx:.3f}] を覆っていない")

    main = Polygon([(x0, -h), (x1, -h), (x1, h), (x0, h)])
    loop = Polygon(np.vstack([outer, [[outer[-1, 0], h]], [[outer[0, 0], h]]]))
    poly = unary_union([main, loop]).difference(Polygon(isl))
    if poly.geom_type != "Polygon":
        raise ValueError(f"流体領域が単一多角形にならない: {poly.geom_type}")
    info.update(x0=float(x0), x1=float(x1), cell_length=float(x1 - x0))
    return poly, Polygon(isl), info


def valve_contours(theta_deg: float, R: float = OPTIMIZED["R"],
                   alpha_deg: float = OPTIMIZED["alpha"],
                   X2: float = OPTIMIZED["X2"], w: float = 1.0):
    """外側境界（ループ側）と島の頂点列を返す。接合角の検証に使う。"""
    return _build(theta_deg, R, alpha_deg, X2, w)


def gamboa_valve(which: str = "optimized", w: float = 1.0, **kw):
    """Gamboa Table 2 の無次元パラメータから組む。"""
    p = dict(OPTIMIZED if which.lower().startswith("opt") else REFERENCE)
    p.update({k: v for k, v in kw.items() if k in p})
    th = theta_from_params(p["X2"], p["n"], p["Y3"], w)
    poly, isl, info = valve_polygon(th, R=p["R"] * w, alpha_deg=p["alpha"],
                                    X2=p["X2"] * w, w=w,
                                    x0=kw.get("x0"), x1=kw.get("x1"))
    info.update(source=f"Gamboa 2005 Table 2 ({which})", params=p)
    return poly, isl, info


def tesla_valve_theta(theta_deg: float, w: float = 1.0,
                      R: float = None, alpha_deg: float = None,
                      X2: float = None, **kw):
    """theta を陽に指定して組む（beta 族用）。

    R, alpha, X2 を省略すると Gamboa optimized の値を使う。
    theta を振れば「接合角の符号と大きさが方向依存性を支配する」という仮説を
    単一パラメータで検証できる。
    """
    R = (OPTIMIZED["R"] if R is None else R) * w
    alpha_deg = OPTIMIZED["alpha"] if alpha_deg is None else alpha_deg
    X2 = (OPTIMIZED["X2"] if X2 is None else X2) * w
    poly, isl, info = valve_polygon(theta_deg, R=R, alpha_deg=alpha_deg,
                                    X2=X2, w=w, x0=kw.get("x0"), x1=kw.get("x1"))
    info.update(source=f"parametric theta = {theta_deg:.3f} deg")
    return poly, isl, info


def theta_family(theta_deg: float, w: float = 1.0, R: float = None,
                 alpha_deg: float = None, Cx: float = None, **kw):
    """beta 族: **円弧中心 C を固定**して theta だけを振る。

    接線条件 Cx = X2 + R(1+cos theta)/sin theta を X2 について解くので、
    R・alpha・外側円弧・閉じ壁・島の下流側が**絶対位置ごと同一**になり、
    変わるのは戻り流路（と、それに伴う分流板先端 X2）だけになる。

    避けられない共変量
    ------------------
    接線条件により **弧の掃引 = theta + 90 - alpha** が theta に縛られる。
    したがって theta を下げるとループが短くなる。
    Gamboa optimized の R, alpha で theta = 161.6 -> 90 度のとき

        弧の掃引  209.7 -> 138.1 度
        ループ長  14.25 -> 13.28（-6.8 %）
        ループ面積 13.47 -> 11.63（-13.7 %）

    **Di_p(theta) の変化を theta 由来と断定する前に、theta を固定して
    ループ長だけを変えた対照計算が要る。** どの構成を選んでもこの共変量は残る
    （戻り流路を円弧に接させない構成にすれば消せるが、壁に折れ角が入る）。
    """
    R = (OPTIMIZED["R"] if R is None else R) * w
    alpha_deg = OPTIMIZED["alpha"] if alpha_deg is None else alpha_deg
    if Cx is None:      # Gamboa optimized の円弧中心
        th0 = np.radians(theta_from_params(OPTIMIZED["X2"], OPTIMIZED["n"],
                                           OPTIMIZED["Y3"], w))
        Cx = (OPTIMIZED["X2"] * w
              + R * (1.0 + np.cos(th0)) / np.sin(th0))
    th = np.radians(theta_deg)
    X2 = Cx - R * (1.0 + np.cos(th)) / np.sin(th)
    poly, isl, info = valve_polygon(theta_deg, R=R, alpha_deg=alpha_deg,
                                    X2=X2, w=w, x0=kw.get("x0"), x1=kw.get("x1"))
    info.update(source=f"theta family (C fixed) theta = {theta_deg:.3f} deg",
                arc_sweep_deg=float(theta_deg + 90.0 - alpha_deg),
                Cx_fixed=float(Cx))
    return poly, isl, info


def periodic_cell(poly: Polygon, info: dict, gap: float, w: float = 1.0):
    """バルブの前後に直線区間を足して周期単位セルにする。

    gap : 追加する直線区間の合計長さ [同じ単位]。左右に gap/2 ずつ振り分ける。
    """
    h = 0.5 * w
    lx, rx = info["loop_x"]
    x0 = lx - 0.5 * gap
    x1 = rx + 0.5 * gap
    p2, isl, i2 = valve_polygon(info["theta_deg"], info["R"], info["alpha_deg"],
                                info["X2"], w, x0=x0, x1=x1)
    i2.update(gap=float(gap), source=info.get("source", ""))
    return p2, isl, i2


def straight_cell(length: float, w: float = 1.0):
    """接続部のみ（直線流路）の周期セル。案 2 の補正の対称性検証に使う。"""
    h = 0.5 * w
    return Polygon([(0, -h), (length, -h), (length, h), (0, h)])


def rasterize_cell(poly: Polygon, cells_per_w: int, w: float = 1.0,
                   pad: int = 2, y_margin: float = 0.0):
    """周期セルを、x 方向に厳密に周期となるようラスタ化する。

    x はセル長を cells_per_w/w 刻みで割り切る前提。左右端の列は主流路のみなので
    周期境界で厳密に一致する。y 方向には pad 個の固体セルを付ける
    （geometry.py 冒頭の np.roll の注意を参照）。
    """
    from matplotlib.path import Path
    xmin, ymin, xmax, ymax = poly.bounds
    res = cells_per_w / w
    nx = int(round((xmax - xmin) * res))
    ymin -= y_margin; ymax += y_margin
    ny = int(round((ymax - ymin) * res))
    xs = xmin + (np.arange(nx) + 0.5) / res
    ys = ymin + (np.arange(ny) + 0.5) / res
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel()])
    inside = Path(np.asarray(poly.exterior.coords)).contains_points(pts)
    for ring in poly.interiors:
        inside &= ~Path(np.asarray(ring.coords)).contains_points(pts)
    core = inside.reshape(nx, ny)
    mask = np.zeros((nx, ny + 2 * pad), bool)
    mask[:, pad:-pad] = core
    ys = np.concatenate([ys[0] - np.arange(pad, 0, -1) / res, ys,
                         ys[-1] + np.arange(1, pad + 1) / res])
    return mask, xs, ys
