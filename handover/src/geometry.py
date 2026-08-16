"""
geometry.py -- DXF から流体マスク（2次元ブール配列）を生成する
==============================================================

重要な注意
----------
`rasterize()` は必ず `pad` 引数（既定 2）で固体セルの縁を付ける。
LBM のストリーミングに numpy.roll を使うため、流体セルが配列端に接すると
その方向が周期境界になり **壁が消滅する**。2026-08-09 のセッションで
直線流路が 36,000 反復収束しなかった原因はこれだった。
"""
from __future__ import annotations
import numpy as np
import ezdxf
from shapely.geometry import Polygon
from shapely.ops import unary_union
from matplotlib.path import Path


def load_dxf(path: str) -> Polygon:
    """DXF の閉じた LWPOLYLINE 群を読み、最大面積を外形・残りを島として差し引く。

    戻り値の単位は DXF の単位（本プロジェクトでは mm, $INSUNITS=4）。
    """
    doc = ezdxf.readfile(path)
    rings = []
    for e in doc.modelspace():
        pts = [(x, y) for x, y, *_ in e.get_points()]
        p = Polygon(pts)
        rings.append((p, p.area))
    if not rings:
        raise ValueError(f"no closed polyline found in {path}")
    rings.sort(key=lambda t: -t[1])
    outer = rings[0][0]
    holes = [r[0] for r in rings[1:]]
    return outer.difference(unary_union(holes)) if holes else outer


def rasterize(poly: Polygon, cells_per_mm: int, pad: int = 2):
    """ポリゴンをセル中心判定でラスタ化し、固体の縁を付けて返す。

    Parameters
    ----------
    poly          : load_dxf() の戻り値
    cells_per_mm  : 解像度。7 月の計算および 2026-08-09 の再計算は 12
    pad           : 上下左右に付ける固体セル数。0 にしてはいけない（冒頭の注意）

    Returns
    -------
    mask : (nx, ny) bool   True = 流体
    xs   : (nx,) float     セル中心の x 座標 [mm]
    ys   : (ny,) float     セル中心の y 座標 [mm]
    """
    xmin, ymin, xmax, ymax = poly.bounds
    nx = int(round((xmax - xmin) * cells_per_mm))
    ny = int(round((ymax - ymin) * cells_per_mm))
    xs = xmin + (np.arange(nx) + 0.5) / cells_per_mm
    ys = ymin + (np.arange(ny) + 0.5) / cells_per_mm
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel()])

    inside = Path(np.asarray(poly.exterior.coords)).contains_points(pts)
    for ring in poly.interiors:
        inside &= ~Path(np.asarray(ring.coords)).contains_points(pts)
    core = inside.reshape(nx, ny)

    if pad <= 0:
        return core, xs, ys
    mask = np.zeros((nx, ny + 2 * pad), bool)
    mask[:, pad:-pad] = core
    ys = np.concatenate([ys[0] - np.arange(pad, 0, -1) / cells_per_mm,
                         ys,
                         ys[-1] + np.arange(1, pad + 1) / cells_per_mm])
    return mask, xs, ys


def mirror_x(mask: np.ndarray) -> np.ndarray:
    """逆流ケース用に幾何を x 方向に鏡像反転する。

    入口・出口の数値処理を順流と完全に同一に保つため、境界条件を
    入れ替えるのではなく幾何を反転する。ソルバ由来の系統誤差が
    順流・逆流の比で相殺される。
    """
    return mask[::-1, :].copy()


def periodic_windows(mask: np.ndarray, period_cells: int) -> np.ndarray:
    """流れ方向に厳密周期となる窓の開始インデックスを列挙する。

    mask[i0 : i0+period_cells] を切り出したとき、その左右が
    連続的に繋がる（= 列パターンが period_cells ごとに一致する）i0 を返す。
    7 月形状（ピッチ 6 mm）では 12 cells/mm で 118 通りが該当する。
    切り出し位置の任意性が結果に効かないことの確認に使う。
    """
    nx = mask.shape[0]
    s = int(period_cells)
    if nx < 2 * s:
        return np.array([], int)
    same = np.array([np.array_equal(mask[i], mask[i + s]) for i in range(nx - s)])
    ok = [i0 for i0 in range(nx - 2 * s + 1) if same[i0:i0 + s].all()]
    return np.array(ok, int)


def seam_windows(mask: np.ndarray, period_cells: int):
    """継ぎ目が厳密に一致する窓を、内部の不一致列数とともに列挙する。

    周期境界で必要なのは「窓の左端の列と、その 1 周期先の列が一致すること」
    だけである。窓の内部が完全な周期でなくても、継ぎ目さえ合っていれば
    壁は連続する（単位セルが隣の段と厳密に同一である必要はない）。

    7 月形状の DXF は 4 つのループが同じ曲線の**異なる折れ線近似**になっており
    （ループ 1 の外壁だけ他と 13.7 um 違い、島 3 の下流端だけ 4.5 um 浮いている）、
    48 cells/mm では厳密周期の窓が存在しない。そこで直線区間で切る。

    Returns
    -------
    list of (i0, n_mismatch)  n_mismatch は窓内部で 1 周期先と異なる列の数
    """
    nx = mask.shape[0]
    s = int(period_cells)
    if nx < 2 * s:
        return []
    same = np.array([np.array_equal(mask[i], mask[i + s]) for i in range(nx - s)])
    out = []
    for i0 in range(nx - 2 * s + 1):
        if same[i0]:
            out.append((i0, int((~same[i0:i0 + s]).sum())))
    return out


def plain_columns(mask: np.ndarray, ys: np.ndarray, main_half_width_mm: float = 0.5):
    """主流路帯の外に流体が無い列（= 直線区間）を True で返す。"""
    main = np.abs(ys) <= main_half_width_mm + 1e-9
    return ~(mask & ~main).any(axis=1)


def choose_unit_cell(mask: np.ndarray, ys: np.ndarray, period_cells: int):
    """単位セルの切り出し位置を選ぶ。

    優先順位
      1. 継ぎ目が一致し、かつ両端が直線区間の列（ループを跨がない）
      2. 継ぎ目が一致する窓のうち、内部の不一致が最小のもの
    戻り値 (i0, n_mismatch, note)
    """
    s = int(period_cells)
    cand = seam_windows(mask, s)
    if not cand:
        raise ValueError(f"継ぎ目が一致する窓がない。period={s} を確認すること")
    plain = plain_columns(mask, ys)
    good = [(i0, n) for i0, n in cand if plain[i0] and plain[i0 + s]]
    if good:
        good.sort(key=lambda t: (t[1], t[0]))
        i0, n = good[len(good) // 2] if good[0][1] == good[-1][1] else good[0]
        return i0, n, "直線区間で切断"
    cand.sort(key=lambda t: (t[1], t[0]))
    return cand[0][0], cand[0][1], "ループ内で切断（直線区間に候補なし）"


def cut_unit_cell(mask: np.ndarray, xs: np.ndarray, i0: int, period_cells: int):
    """周期単位セルを切り出す。x 方向は周期境界で使うので pad を付けない。

    y 方向の固体縁は元の mask が既に持っている（rasterize の pad）。
    """
    s = int(period_cells)
    sub = mask[i0:i0 + s].copy()
    if sub[:, 0].any() or sub[:, -1].any():
        raise ValueError("y 端に流体セルがある。rasterize の pad を確認すること。")
    return sub, xs[i0:i0 + s].copy()


def check_rasterization(poly: Polygon, mask: np.ndarray, cells_per_mm: int) -> dict:
    """ラスタ化の妥当性チェック。面積誤差が 0.1 % を超えたら解像度を上げる。"""
    area_raster = mask.sum() / cells_per_mm ** 2
    return dict(area_polygon_mm2=poly.area,
                area_raster_mm2=area_raster,
                error_pct=100 * (area_raster - poly.area) / poly.area,
                fluid_cells=int(mask.sum()),
                solid_cells=int((~mask).sum()),
                shape=tuple(mask.shape))
