"""
xia2025.py -- Xia et al. (2025) の TVM / SYMTVM 流路の形状生成
=================================================================

出典: Xia, Wu, Deng, Yuan, Zhang, "Numerical study on heat transfer and flow
characteristics of symmetric Tesla-type microchannel heat sinks",
Applied Thermal Engineering 258 (2025) 124611.  Fig. 2–4, Table 2。

**[論文記載]**（Fig. 4 の寸法、単位 mm）
  流路幅 w = 0.30（trunk・helix とも）、島の円弧 R_in = 0.58、外側円弧 R_out = 0.88、
  島の先端角 30°、入口から最初の島の円弧中心まで 2.12、島のピッチ 4.72、
  出口直線の傾き 165°（= 水平から 15°）、SYMTVM の入口幅 0.62 / TVM 0.31、
  SYMTVM の菱形の辺 1.84、入口壁が外側円弧に当たる位置 1.29、外壁長 3.75 / 4.36。
  Table 2: Ly = 50、W = 3.54、H = 2、Hch = 1.5。

**[逆算]**（Fig. 4 の画素計測と幾何関係から決めたもの。`docs/xia2025_geometry.md`）
  - 主流路は ±15° のジグザグ帯。帯の壁は島の円に接する。
  - 島の中心は y = ±cy、cy = (4.72/2)·tan15° = 0.632。
    これは「帯の壁が次の島の縁と一直線になる」条件から決まり、菱形の辺が
    (2 cy cos15° − 0.30)/(2 sin15° cos15°) = 1.84 になることで裏付けられる。
  - SYMTVM の島は、軸側の角が −15° の弦で切られている（Fig. 4 の 0.87 の辺）。
    弦長 0.87 から位置を決める（中心から 0.384）。前段の −15° 帯の外壁（中心から
    0.342）とほぼ同じ線で、帯はこの島の下を幅 0.26 で通る。
  - 出口は x = 50 で終わる（Table 2 の Ly）。Fig. 4 の出口長 1.00 は 0.8 になる。

Forward / Reverse は同じ流路の 180° 回転（入口を左に戻したもの）。
"""
from __future__ import annotations
import numpy as np
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from shapely import affinity

# ---- 論文記載値 ----
W_CH = 0.30          # trunk / helix の幅
R_IN = 0.58          # 島の円弧
R_OUT = 0.88         # 外側円弧
TIP_DEG = 30.0       # 島の先端角
PITCH = 4.72         # 島のピッチ（x）
X_C1 = 2.12          # 最初の島の円弧中心 x
ANG = 15.0           # 帯の傾き [deg]
LY = 50.0            # 流路長
W_UNIT = 3.54        # 計算単位の幅
H_CH = 1.5           # 流路深さ
H_SOLID = 2.0        # 板の厚さ
N_STAGE = 10
CHORD_LEN = 0.87     # SYMTVM の島の軸側を切る弦の長さ（Fig. 4）

# ---- 逆算値 ----
CY = 0.5 * PITCH * np.tan(np.radians(ANG))      # 0.632

_s, _c = np.sin(np.radians(ANG)), np.cos(np.radians(ANG))
D_UP = np.array([_c, _s]);   N_UP = np.array([-_s, _c])     # +15° とその左法線
D_DN = np.array([_c, -_s]);  N_DN = np.array([_s, _c])      # −15° とその左法線
S_ENTRY = float(np.sqrt(R_OUT ** 2 - R_IN ** 2))            # 接線が外側円に入るまでの距離


def _line_pt(C, n, off, x=None, y=None):
    """(P−C)·n = off の直線上で x または y を指定した点。"""
    C = np.asarray(C, float)
    if x is not None:
        return np.array([x, C[1] + (off - (x - C[0]) * n[0]) / n[1]])
    return np.array([C[0] + (off - (y - C[1]) * n[1]) / n[0], y])


def _isect(C1, n1, o1, C2, n2, o2):
    """2 直線 (P−C1)·n1=o1, (P−C2)·n2=o2 の交点。"""
    A = np.array([n1, n2]); b = np.array([o1 + np.dot(C1, n1), o2 + np.dot(C2, n2)])
    return np.linalg.solve(A, b)


def _arc(C, r, a0, a1, step=2.0):
    """角度 a0 → a1 [deg] の円弧点列（符号で向きを決める）。両端を含む。"""
    n = max(2, int(abs(a1 - a0) / step) + 1)
    a = np.radians(np.linspace(a0, a1, n))
    return np.c_[C[0] + r * np.cos(a), C[1] + r * np.sin(a)]


def _ang(v):
    return float(np.degrees(np.arctan2(v[1], v[0])))


def teardrop(C, upper=True, chord_off=None):
    """島: 半径 R_IN の円 + 先端角 30° の接線。tip は水平方向。

    chord_off : None 以外なら、軸側の角を (P−C)·n = chord_off の直線で切る
                （上側の島は n = N_DN、下側は n = N_UP）。SYMTVM 用。
    """
    C = np.asarray(C, float)
    half = np.radians(0.5 * TIP_DEG)
    tip = C + np.array([R_IN / np.sin(half), 0.0])
    # 接点: 上の縁は +75°、下の縁は −75°
    t_up = C + R_IN * np.array([np.sin(half), np.cos(half)])
    t_lo = C + R_IN * np.array([np.sin(half), -np.cos(half)])
    a_up, a_lo = 90 - 0.5 * TIP_DEG, -(90 - 0.5 * TIP_DEG)
    pts = np.vstack([tip[None], t_up[None], _arc(C, R_IN, a_up, 360 + a_lo), t_lo[None]])
    poly = Polygon(pts)
    if chord_off is not None:
        n = N_DN if upper else N_UP
        # 弦の半平面（中心側を残す）: (P−C)·n >= chord_off (upper) / <= (lower)
        big = 10.0
        d = np.array([n[1], -n[0]])
        P0 = C + chord_off * n
        if upper:
            hp = Polygon([P0 - big * d, P0 + big * d, P0 + big * d + big * n, P0 - big * d + big * n])
        else:
            hp = Polygon([P0 - big * d, P0 + big * d, P0 + big * d - big * n, P0 - big * d - big * n])
        poly = poly.intersection(hp)
    return poly


# ---------------------------------------------------------------------------
# TVM（古典的な多段 Tesla バルブ）
# ---------------------------------------------------------------------------
def tvm():
    """TVM-Forward の流体領域（Polygon）と情報。"""
    Cs = [np.array([X_C1 + PITCH * k, CY if k % 2 == 0 else -CY]) for k in range(N_STAGE)]
    # 入口: 帯 A（+15°、上壁は島 1 の下縁）が x=0.66 で水平に折れる。縦幅 0.31
    x_bend = 0.66
    y_in_top = _line_pt(Cs[0], N_UP, -R_IN, x=x_bend)[1]
    y_in_bot = _line_pt(Cs[0], N_UP, -R_OUT, x=x_bend)[1]
    # 出口: 最後の帯（+15°、上壁は島 10 の下縁）が y=0 に達する x
    x_out = _line_pt(Cs[-1], N_UP, -R_IN, y=0.0)[0]
    y_out_bot = _line_pt(Cs[-1], N_UP, -R_OUT, x=x_out)[1]

    top = [np.array([0.0, y_in_top]), np.array([x_bend, y_in_top])]
    # 上側経路: 上の島 (k 偶数) ごとに 「+15° 線 → 外側円弧 → −15° 壁 → 角」
    for k in range(0, N_STAGE, 2):
        C = Cs[k]
        foot = C - R_IN * N_UP
        entry = foot - S_ENTRY * D_UP
        top.append(entry)
        # 入口点（左下）から時計回りに左・上を通って +75° の接点へ
        top.extend(_arc(C, R_OUT, _ang(entry - C), 90 - 0.5 * TIP_DEG - 360)[1:])
        if k + 1 < N_STAGE:
            # −15° 壁（offset +R_OUT）と、下の島 k+1 の下縁（+15°、offset −R_IN）の角
            Q = _isect(C, N_DN, R_OUT, Cs[k + 1], N_UP, -R_IN)
            top.append(Q)
    top.append(np.array([x_out, 0.0]))
    top.append(np.array([LY, 0.0]))

    bot = [np.array([LY, y_out_bot]), np.array([x_out, y_out_bot])]
    # 下側経路（右→左）: 下の島 (k 奇数) ごとに「+15° 壁（offset −R_OUT）→ 外側円弧 → −15° 線 → 角」
    for k in range(N_STAGE - 1, 0, -2):
        C = Cs[k]
        # 円弧: 下の接点 (−75°) から左回りに、上縁 (−15°, offset +R_IN) の入口まで
        foot = C + R_IN * N_DN
        entry = foot - S_ENTRY * D_DN
        a0 = -(90 - 0.5 * TIP_DEG)
        a1 = _ang(entry - C)
        bot.extend(_arc(C, R_OUT, a0, a1 - 360 if a1 > a0 else a1)[0:])
        # −15° 線を左へ、上の島 k−1 の下壁 (+15°, offset −R_OUT) との角
        Q = _isect(Cs[k - 1], N_UP, -R_OUT, C, N_DN, R_IN)
        bot.append(Q)
    bot.append(np.array([x_bend, y_in_bot]))
    bot.append(np.array([0.0, y_in_bot]))

    outer = Polygon(np.vstack([np.array(top), np.array(bot)]))
    islands = [teardrop(C, upper=(k % 2 == 0)) for k, C in enumerate(Cs)]
    poly = outer.difference(unary_union(islands))
    info = dict(name="TVM", centres=[c.tolist() for c in Cs], x_bend=x_bend,
                x_out=float(x_out), inlet_y=[float(y_in_bot), float(y_in_top)],
                outlet_y=[float(y_out_bot), 0.0], w_in=float(y_in_top - y_in_bot))
    return _check(poly), info


# ---------------------------------------------------------------------------
# SYMTVM（対称型）
# ---------------------------------------------------------------------------
def symtvm():
    """SYMTVM-Forward の流体領域（Polygon）と情報。"""
    Cu = [np.array([X_C1 + PITCH * k, CY]) for k in range(N_STAGE)]
    Cl = [np.array([c[0], -c[1]]) for c in Cu]
    h_in = 0.31                                   # 入口半幅 (Wch = 0.62)
    # 入口上壁 y=h_in が島 1 の外側円弧に当たる x
    x_in = X_C1 - np.sqrt(R_OUT ** 2 - (h_in - CY) ** 2)
    # 最終段の −15° 壁が y=h_in に達する x（= 出口の始まり）
    x_out = _line_pt(Cu[-1], N_DN, R_OUT, y=h_in)[0]

    top = [np.array([0.0, h_in]), np.array([x_in, h_in])]
    for k, C in enumerate(Cu):
        # 外側円弧: 入口点から時計回りに +75° の接点まで
        if k == 0:
            entry = np.array([x_in, h_in])
        else:
            # 前段の −15° 外壁 (offset +R_OUT from Cu[k-1]) が円 (C, R_OUT) に入る点
            T = Cu[k - 1] + R_OUT * N_DN
            rel = T - C
            b = 2 * np.dot(rel, D_DN); c0 = np.dot(rel, rel) - R_OUT ** 2
            t = (-b - np.sqrt(b * b - 4 * c0)) / 2
            entry = T + t * D_DN
        a0 = _ang(entry - C)
        a1 = 90 - 0.5 * TIP_DEG - 360          # 時計回りに左・上を通る
        top.extend(_arc(C, R_OUT, a0, a1)[1:] if k == 0 else _arc(C, R_OUT, a0, a1))
    top.append(np.array([x_out, h_in]))
    top.append(np.array([LY, h_in]))
    bot = [np.array([p[0], -p[1]]) for p in top[::-1]]
    outer = Polygon(np.vstack([np.array(top), np.array(bot)]))

    # 島: 涙滴の軸側の角を −15° の弦で切る（Fig. 4 の 0.87 の辺）。島 1 も同じ。
    # 弦の位置は中心からの距離で与える（弦長 0.87 → 0.384）。前段の −15° 外壁
    # （中心から 0.342）とほぼ平行・同位置だが、完全に同一線にすると境界の 3 点が
    # 共線になり、三角形分割で面積ゼロの三角形が出て非多様体になる。
    chord_off = -float(np.sqrt(R_IN ** 2 - (0.5 * CHORD_LEN) ** 2))   # −0.384
    islands = [teardrop(C, True, chord_off) for C in Cu]
    islands += [teardrop(C, False, -chord_off) for C in Cl]
    # 菱形: 上左辺 = 上の島の下縁から 0.30 外、上右辺 = 上の島の上縁の延長
    diamonds = []
    for C in Cu:
        Dl = _line_pt(C, N_UP, -R_OUT, y=0.0)
        Dr = _line_pt(C, N_DN, R_IN, y=0.0)
        Vt = _isect(C, N_UP, -R_OUT, C, N_DN, R_IN)
        diamonds.append(Polygon([Dl, Vt, Dr, [Vt[0], -Vt[1]]]))
    poly = outer.difference(unary_union(islands + diamonds))
    side = float(np.hypot(*(Vt - Dl)))
    info = dict(name="SYMTVM", centres_upper=[c.tolist() for c in Cu], x_in=float(x_in),
                x_out=float(x_out), inlet_y=[-h_in, h_in], outlet_y=[-h_in, h_in],
                diamond_side=side, chord_off=float(chord_off),
                wall_len_first=float(np.linalg.norm((Cu[0] + R_OUT * N_DN) - _entry_pt(Cu, 1))))
    return _check(poly), info


def _entry_pt(Cu, k):
    T = Cu[k - 1] + R_OUT * N_DN; C = Cu[k]
    rel = T - C; b = 2 * np.dot(rel, D_DN); c0 = np.dot(rel, rel) - R_OUT ** 2
    return T + ((-b - np.sqrt(b * b - 4 * c0)) / 2) * D_DN


def _check(poly):
    if poly.geom_type != "Polygon":
        # 小さなかけらが出たら最大のものを採る
        poly = max(poly.geoms, key=lambda g: g.area)
    if not poly.is_valid:
        poly = poly.buffer(0)
    return poly


def reverse(poly):
    """Reverse 構成: 同じ流路を逆向きに流す = 流路を 180° 回転して入口を左に戻す。

    Fig. 3 の TVM-Reverse は最初の島が「上側・丸い端が右」で、出口が中心線より上。
    これは Forward の 180° 回転（x 反転 + y 反転）で、x 方向の鏡像ではない
    （鏡像だと最初の島が下側になる）。SYMTVM は上下対称なので両者は同じ。
    """
    return affinity.rotate(poly, 180.0, origin=(0.5 * LY, 0.0))


def all_four():
    """Fig. 3 の 4 構成。名前 → (Polygon, info)。"""
    t, ti = tvm(); s, si = symtvm()
    return {
        "SYMTVM-Forward": (s, si),
        "SYMTVM-Reverse": (reverse(s), dict(si, name="SYMTVM-Reverse")),
        "TVM-Forward": (t, ti),
        "TVM-Reverse": (reverse(t), dict(ti, name="TVM-Reverse")),
    }
