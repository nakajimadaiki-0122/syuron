"""fluid.stl / solid_copper.stl を境界ゾーンごとの STL に分割する。

Fluent Meshing は STL をファイル名のゾーンとして読むので、読み込むだけで
inlet / outlet / fluid_top ... の名前付けが終わる。面の分離を Fluent 側でやらない。

使い方:  python split_zones.py            （4 構成すべて）
         python split_zones.py tvm_fwd    （1 構成だけ）
出力:    <構成>/zones/*.stl（10 ファイル。固体側の界面は流体側の面が兼ねるので出さない）
"""
import struct
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import Polygon
from scipy.spatial import Delaunay

HERE = Path(__file__).resolve().parent
CONFIGS = ["symtvm_fwd", "symtvm_rev", "tvm_fwd", "tvm_rev"]
TOL = 1e-4  # mm
SKIP = {"solid_interface"}  # 流体側と同じ面。二重になるので出さない
# 平らな蓋は元の耳刈り分割に高さ 1e-4 mm の細長い三角形が多く、Fluent が読み込み時に
# 削って穴になる。輪郭の点はそのままに、制約付き Delaunay で分割し直す
CAPS = {"fluid_top", "fluid_bottom"}
SPACING = 0.1   # 面の内側に足す点の間隔 [mm]
METHODS = {}
# 銅ブロックの外面。外側の辺は点が 5 mm おきで、流路の輪郭（0.05 mm おき）と結ぶと
# 細長い三角形になる。外側の辺に STEP 間隔で点を足し、隣り合う面どうしで同じ点を使う
OUTER = {"bottom_heat": 1, "solid_top": 1, "sym_left": 2, "sym_right": 2, "solid_ends": 0}
STEP = 0.25  # mm

# 幾何の約束（Y が鉛直）。README と同じ
X0, X1 = 0.0, 50.0
Y_BOT, Y_CH, Y_TOP = 0.0, 0.5, 2.0
Z_HALF = 1.77


def read_stl(path):
    data = Path(path).read_bytes()
    n = struct.unpack_from("<I", data, 80)[0]
    rec = np.frombuffer(data, dtype=np.dtype([
        ("n", "<3f4"), ("v", "<3,3f4"), ("attr", "<u2")]), count=n, offset=84)
    return rec["v"].astype(float)  # (n, 3, 3)


def write_stl(path, tris):
    tris = np.asarray(tris, float)
    with open(path, "wb") as fh:
        fh.write(b"\0" * 80)
        fh.write(struct.pack("<I", len(tris)))
        for v in tris:
            n = np.cross(v[1] - v[0], v[2] - v[0])
            ln = np.linalg.norm(n)
            n = n / ln if ln > 0 else np.zeros(3)
            fh.write(struct.pack("<3f", *n))
            for p in v:
                fh.write(struct.pack("<3f", *p))
            fh.write(struct.pack("<H", 0))


def on_plane(tris, axis, value):
    """3 頂点すべてが axis = value の平面上にある三角形"""
    return np.all(np.abs(tris[:, :, axis] - value) < TOL, axis=1)


def boundary_loops(tris):
    """蓋の三角形群から外周・穴の閉じた輪郭（頂点列）を取り出す。1 回しか使われない辺を辿る"""
    v = tris.reshape(-1, 3)
    uniq, inv = np.unique(v, axis=0, return_inverse=True)
    idx = inv.reshape(-1, 3)
    use = {}
    for a, b, c in idx:
        for e in ((a, b), (b, c), (c, a)):
            use[e] = use.get(e, 0) + 1
    # 向き付き辺 e が 1 回で、逆向きが無いものが境界。向きは三角形の向きに従う
    bnd = {e for e, n in use.items() if n == 1 and (e[1], e[0]) not in use}
    nxt = {}
    for a, b in bnd:
        assert a not in nxt, "境界が分岐している"
        nxt[a] = b
    loops = []
    seen = set()
    for start in list(nxt):
        if start in seen:
            continue
        loop, cur = [], start
        while cur not in seen:
            seen.add(cur)
            loop.append(cur)
            cur = nxt[cur]
        loops.append(uniq[loop])
    return loops


def key(p):
    return tuple(np.round(np.asarray(p, float), 6))


def on_box_plane(a, b):
    """2 点が同じ外面（x = 0 / 50、z = ±Z_HALF、y = 0）に乗っているか"""
    for ax, vals in ((0, (X0, X1)), (2, (-Z_HALF, Z_HALF)), (1, (Y_BOT,))):
        for v in vals:
            if abs(a[ax] - v) < TOL and abs(b[ax] - v) < TOL:
                return True
    return False


def densify_loop(L, step, fixed):
    """外面に乗る長い辺に等間隔の点を足す。点は端点の並びから決めるので、
    同じ辺を共有する面どうしで同じ点になる。fixed（入口・出口の角）どうしの辺は触らない"""
    out = []
    n = len(L)
    for i in range(n):
        a, b = L[i], L[(i + 1) % n]
        out.append(a)
        if not on_box_plane(a, b) or (key(a) in fixed and key(b) in fixed):
            continue
        length = float(np.linalg.norm(b - a))
        k = int(np.ceil(length / step))
        if k < 2:
            continue
        lo, hi = (a, b) if key(a) < key(b) else (b, a)
        pts = [lo + (hi - lo) * (j / k) for j in range(1, k)]
        if key(a) > key(b):
            pts = pts[::-1]
        out.extend(pts)
    return np.asarray(out)


def triangulate_polygon(poly, spacing):
    """輪郭の点はそのまま、内側に spacing 間隔の点を足して Delaunay 分割する。
    輪郭の辺がすべて残り、面積が一致したときだけ採用。だめなら制約付き Delaunay、
    それもだめなら耳刈り（mapbox_earcut）に落とす。戻り値は (N, 3, 2) と手法名"""
    rings = [np.asarray(poly.exterior.coords)[:-1]] + [np.asarray(r.coords)[:-1] for r in poly.interiors]
    bpts = np.vstack(rings)
    bedges = set()
    for R in rings:
        n = len(R)
        for i in range(n):
            bedges.add(tuple(sorted((key(R[i]), key(R[(i + 1) % n])))))
    if spacing:
        x0, y0, x1, y1 = poly.bounds
        gx = np.arange(x0 + spacing / 2, x1, spacing)
        gy = np.arange(y0 + spacing / 2, y1, spacing)
        G = np.array([(x, y) for x in gx for y in gy]) if len(gx) and len(gy) else np.zeros((0, 2))
        if len(G):
            inside = shapely.contains_xy(poly, G[:, 0], G[:, 1])
            G = G[inside]
        if len(G):
            far = shapely.dwithin(poly.boundary, shapely.points(G), 0.6 * spacing)
            G = G[~far]
        if len(G):
            # 輪郭の辺 (a, b) の直径円の中に点があると Delaunay がその辺を壊す。
            # 長い辺（TVM の直線壁）の直径円に入る格子点を捨てる
            drop = np.zeros(len(G), bool)
            for R in rings:
                A, B = R, np.roll(R, -1, axis=0)
                for a, b in zip(A, B):
                    drop |= ((a - G) * (b - G)).sum(axis=1) < 0
            G = G[~drop]
        pts = np.vstack([bpts, G]) if len(G) else bpts
        if len(pts) >= 3:
            try:
                dt = Delaunay(pts)
                T = pts[dt.simplices]
                c = T.mean(axis=1)
                keep = shapely.contains_xy(poly, c[:, 0], c[:, 1])
                T = T[keep]
                got = set()
                for tri in T:
                    k = [key(v) for v in tri]
                    for e in ((k[0], k[1]), (k[1], k[2]), (k[2], k[0])):
                        got.add(tuple(sorted(e)))
                d1, d2 = T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]
                area = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]).sum()
                if bedges <= got and abs(area - poly.area) < 1e-6 * max(poly.area, 1):
                    return T, "delaunay+grid"
                print(f"    delaunay+grid 不採用: 欠けた辺 {len(bedges - got)}、面積 {area:.4f} / {poly.area:.4f}")
            except Exception as e:  # noqa: BLE001
                print("    delaunay+grid 失敗:", type(e).__name__, e)
    try:
        tri = shapely.constrained_delaunay_triangles(poly)
        return np.array([np.asarray(g.exterior.coords)[:3] for g in tri.geoms]), "cdt"
    except Exception as e:  # noqa: BLE001
        print("    cdt 失敗:", e, "| 頂点数", len(bpts), "穴", len(rings) - 1, "面積", round(poly.area, 4))
        import mapbox_earcut
        idx = mapbox_earcut.triangulate_float64(bpts, np.cumsum([len(R) for R in rings]))
        return bpts[idx.reshape(-1, 3)], "earcut"


def retriangulate_flat(tris, axis=1, step=None, fixed=frozenset(), spacing=None):
    """axis = 一定 の平らな面を、同じ輪郭点で制約付き Delaunay 分割し直す。向きは元に合わせる。
    step を与えると外面の辺に点を足す（densify_loop）"""
    if len(tris) == 0:
        return tris
    n0 = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]).sum(axis=0)
    sign = np.sign(n0[axis])
    uv = [i for i in range(3) if i != axis]          # 平面内の 2 軸
    loops = boundary_loops(tris)
    if step:
        loops = [densify_loop(L, step, fixed) for L in loops]
    rings = [Polygon(L[:, uv]) for L in loops]
    rings = [r if r.is_valid else r.buffer(0) for r in rings]
    # 入れ子の深さ。偶数が外周、その 1 つ内側が穴
    # 判定点は輪郭上の頂点。representative_point だと穴（島）の中に落ちることがある。
    # 入れ子は同じ高さ（axis の値）の輪郭どうしでだけ判定する。端面は x = 0 と x = 50 の
    # 2 枚が投影で重なるので、混ぜると片方が穴に見える
    pts = [shapely.Point(L[0, uv]) for L in loops]
    lev = [round(float(L[0, axis]), 4) for L in loops]
    depth = [sum(1 for j, q in enumerate(rings) if j != i and lev[j] == lev[i] and q.contains(pts[i]))
             for i in range(len(rings))]
    out = []
    for i, r in enumerate(rings):
        if depth[i] % 2:
            continue
        holes = [loops[j][:, uv] for j, q in enumerate(rings)
                 if depth[j] == depth[i] + 1 and lev[j] == lev[i] and r.contains(pts[j])]
        poly = Polygon(loops[i][:, uv], holes)
        level = float(loops[i][0, axis])          # 輪郭ごと（端面は x = 0 と x = 50 の 2 枚）
        T, how = triangulate_polygon(poly, spacing)
        METHODS[how] = METHODS.get(how, 0) + 1
        for c in T:
            p = np.zeros((3, 3))
            p[:, uv] = c
            p[:, axis] = level
            n = np.cross(p[1] - p[0], p[2] - p[0])[axis]
            if n * sign < 0:
                p = p[::-1]
            out.append(p)
    out = np.asarray(out)
    a_old = 0.5 * np.linalg.norm(np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1).sum()
    a_new = 0.5 * np.linalg.norm(np.cross(out[:, 1] - out[:, 0], out[:, 2] - out[:, 0]), axis=1).sum()
    assert abs(a_old - a_new) < 1e-6 * max(a_old, 1), ("面積が変わった", a_old, a_new)
    if not step:
        v_old = {tuple(x) for x in np.unique(tris.reshape(-1, 3), axis=0)}
        v_new = {tuple(x) for x in np.unique(out.reshape(-1, 3), axis=0)}
        assert v_new <= v_old, "新しい頂点が増えた（壁と合わなくなる）"
    return out


def split(tris, rules, rest_name):
    """rules = [(name, mask_fn), ...] を順に適用。先に当たった名前が優先"""
    out = {}
    taken = np.zeros(len(tris), bool)
    for name, fn in rules:
        m = fn(tris) & ~taken
        out[name] = tris[m]
        taken |= m
    out[rest_name] = tris[~taken]
    return out


def split_fluid(tris):
    return split(tris, [
        ("inlet",        lambda t: on_plane(t, 0, X0)),
        ("outlet",       lambda t: on_plane(t, 0, X1)),
        ("fluid_top",    lambda t: on_plane(t, 1, Y_TOP)),
        ("fluid_bottom", lambda t: on_plane(t, 1, Y_CH)),
    ], "fluid_walls")


def split_solid(tris):
    return split(tris, [
        ("bottom_heat", lambda t: on_plane(t, 1, Y_BOT)),
        ("solid_top",   lambda t: on_plane(t, 1, Y_TOP)),
        ("sym_left",    lambda t: on_plane(t, 2, -Z_HALF)),
        ("sym_right",   lambda t: on_plane(t, 2, +Z_HALF)),
        ("solid_ends",  lambda t: on_plane(t, 0, X0) | on_plane(t, 0, X1)),
    ], "solid_interface")


def main(configs):
    for cfg in configs:
        d = HERE / cfg
        out = d / "zones"
        out.mkdir(exist_ok=True)
        fixed = frozenset()
        for src, fn in (("fluid.stl", split_fluid), ("solid_copper.stl", split_solid)):
            tris = read_stl(d / src)
            parts = fn(tris)
            if src == "fluid.stl":
                fixed = frozenset(key(p) for n in ("inlet", "outlet") for p in parts[n].reshape(-1, 3))
            total = sum(len(v) for v in parts.values())
            assert total == len(tris), (cfg, src, total, len(tris))
            print(f"{cfg}/{src}: {len(tris)} 三角形")
            for name, t in parts.items():
                if name in SKIP:
                    # 流体側の fluid_walls / fluid_bottom が共有面を兼ねる。
                    # 固体側にも同じ面を入れると二重になり Fluent で free faces が出る
                    print(f"  {name:16s} {len(t):6d}  (書き出さない。流体側の面と重複)")
                    continue
                if name in CAPS:
                    t = retriangulate_flat(t, axis=1, spacing=SPACING)
                elif name in OUTER:
                    t = retriangulate_flat(t, axis=OUTER[name], step=STEP, fixed=fixed, spacing=SPACING)
                write_stl(out / f"{name}.stl", t)
                lo = t.min(axis=(0, 1)) if len(t) else [np.nan] * 3
                hi = t.max(axis=(0, 1)) if len(t) else [np.nan] * 3
                print(f"  {name:16s} {len(t):6d}  x {lo[0]:6.2f}..{hi[0]:6.2f}"
                      f"  y {lo[1]:5.2f}..{hi[1]:5.2f}  z {lo[2]:6.2f}..{hi[2]:6.2f}")


if __name__ == "__main__":
    main(sys.argv[1:] or CONFIGS)
    print("分割手法の内訳:", METHODS)
