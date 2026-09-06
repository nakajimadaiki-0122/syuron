#!/usr/bin/env python3
"""
21_make_report.py -- 3D モデルと結果をまとめた閲覧用 HTML を作る
================================================================

`report/gamboa_report.html` を出力する。外部依存なしの 1 ファイル
（STL は base64 の float32、図は data URI）。ブラウザで開けばそのまま見られる。

3D 表示は canvas の自前ソフトウェアレンダラ（画家のアルゴリズム + フラットシェーディング）。
ドラッグで回転、ホイールで拡大縮小。
"""
import base64, json, os, struct, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "report", "gamboa_report.html")


def stl_vertices(path):
    b = open(path, "rb").read()
    n = struct.unpack("<I", b[80:84])[0]
    a = np.frombuffer(b[84:], dtype=np.uint8).reshape(n, 50)
    v = a[:, 12:48].copy().view("<f4").reshape(n * 9)
    return n, v.astype(np.float32)


def b64(arr):
    return base64.b64encode(arr.tobytes()).decode("ascii")


def img_data_uri(path):
    with open(path, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode("ascii")


def load_json(p, default=None):
    p = os.path.join(ROOT, p)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else default


def fmt(x, n=4):
    return "—" if x is None else f"{x:.{n}f}"


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    models = []
    for key, fn, label, note in (
        ("fluid", "gamboa_fluid_3d.stl", "流体体積（プレナム込み）",
         "CFD が解いている領域そのもの"),
        ("valve", "gamboa_valve_only_fluid_3d.stl", "流体体積（バルブ本体）",
         "プレナムを外した版"),
        ("plate", "gamboa_plate_3d.stl", "溝を彫った板",
         "実物モデル。裏面にポート穴"),
    ):
        p = os.path.join(ROOT, "cad", fn)
        n, v = stl_vertices(p)
        models.append(dict(key=key, file=fn, label=label, note=note,
                           tris=n, data=b64(v)))
        print(f"  {fn}: {n} 三角形, base64 {len(models[-1]['data'])/1024:.0f} KB")

    figs = {}
    for key, fn in (("overlay", "gamboa_full.png"),
                    ("fields", "gamboa_full_fields_Re100_cpm16.png"),
                    ("topology", "fig2_topology.png")):
        p = os.path.join(ROOT, "figures", fn)
        if os.path.exists(p):
            figs[key] = img_data_uri(p)
            print(f"  {fn}: {os.path.getsize(p)/1024:.0f} KB")

    cad = load_json("results/cad_models.json", {})
    geom = load_json("results/gamboa_full_geom.json", {})
    phi = load_json("results/gamboa_full_phi_Re100_cpm16.json", {})
    top = load_json("results/fig2_topology.json", {})

    runs = []
    for Re, cpm in ((100, 16), (300, 16), (500, 32)):
        d = {}
        for lab in ("fwd", "rev"):
            j = load_json(f"results/gamboa_full_Re{Re}_cpm{cpm}_{lab}.json")
            if j and lab in j:
                d[lab] = j[lab]
            elif j is None:
                j2 = load_json(f"results/gamboa_full_Re{Re}_cpm{cpm}.json")
                if j2 and lab in j2:
                    d[lab] = j2[lab]
        if d:
            runs.append((Re, cpm, d))

    ref = {100: 1.02, 200: 1.07, 300: 1.12, 500: 1.37, 1000: 1.73, 2000: 1.92}
    rows = []
    for Re, cpm, d in runs:
        f = d.get("fwd", {}).get("dp_Pa")
        r = d.get("rev", {}).get("dp_Pa")
        di = (r / f) if (f and r) else None
        rows.append(dict(Re=Re, cpm=cpm, dpf=f, dpr=r, di=di, ref=ref.get(Re),
                         conv=[d.get(k, {}).get("converged") for k in ("fwd", "rev")]))

    html = TEMPLATE.format(
        models_json=json.dumps([{k: m[k] for k in ("key", "label", "note",
                                                   "tris", "data")}
                                for m in models]),
        rows_html=results_rows(rows),
        fig_overlay=figs.get("overlay", ""),
        fig_fields=figs.get("fields", ""),
        fig_topology=figs.get("topology", ""),
        iou=fmt(geom.get("overlay", {}).get("iou"), 4),
        phi_f=fmt(100 * phi.get("fwd", {}).get("phi_mean", float("nan")), 1),
        phi_r=fmt(100 * phi.get("rev", {}).get("phi_mean", float("nan")), 1),
        phi_ratio=fmt((phi.get("rev", {}).get("phi_mean", 1)
                       / phi.get("fwd", {}).get("phi_mean", 1)), 2),
        coll_ang=fmt(top.get("collinearity", {}).get("right", {})
                     .get("dangle_deg", float("nan")), 3),
        coll_off=fmt(top.get("collinearity", {}).get("right", {})
                     .get("offset_wv", float("nan")), 3),
        di100=fmt(next((r["di"] for r in rows if r["Re"] == 100), None), 4),
        plate=plate_spec(cad),
    )
    open(OUT, "w", encoding="utf-8").write(html)
    print(f"\nsaved {os.path.normpath(OUT)}  "
          f"({os.path.getsize(OUT)/1024/1024:.2f} MB)")


def plate_spec(cad):
    ext = cad.get("plate_extent", [0, 0, 0, 0])
    return (f"{ext[1]-ext[0]:.1f} × {ext[3]-ext[2]:.1f} × "
            f"{cad.get('wall_mm',0)+cad.get('depth_mm',0):.1f} mm、"
            f"ポート径 {cad.get('port_mm','—')} mm")


def results_rows(rows):
    out = []
    for r in rows:
        di = r["di"]
        ref = r["ref"]
        if di and ref:
            d = f"{100*(di/ref-1):+.1f} %"
            cls = "ok" if abs(di / ref - 1) < 0.1 else "warn"
        else:
            d, cls = "計算中", "pending"
        out.append(
            f'<tr><td class="num">{r["Re"]}</td>'
            f'<td class="num">{r["cpm"]}</td>'
            f'<td class="num">{fmt(r["dpf"], 3)}</td>'
            f'<td class="num">{fmt(r["dpr"], 3)}</td>'
            f'<td class="num strong">{fmt(r["di"], 4)}</td>'
            f'<td class="num">{fmt(r["ref"], 2)}</td>'
            f'<td class="num"><span class="tag {cls}">{d}</span></td></tr>')
    return "\n".join(out)


TEMPLATE = r"""<title>Tesla バルブ再現ノート</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=Noto+Sans+JP:wght@400;500;700&family=Zen+Kaku+Gothic+New:wght@700;900&display=swap">
<style>
:root {{
  --ground:#f2f5f6; --surface:#ffffff; --surface-2:#e9eef0;
  --ink:#111a1e; --muted:#5b6b72; --rule:#d5dee1;
  --fwd:#0e7c7b; --rev:#c0521c; --ok:#1d6f4f; --warn:#a4611a;
  --shadow:0 1px 2px rgba(17,26,30,.06), 0 8px 24px -16px rgba(17,26,30,.30);
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ground:#0e1417; --surface:#151e22; --surface-2:#1b262b;
    --ink:#e4edf0; --muted:#93a5ac; --rule:#26333a;
    --fwd:#3fbfb4; --rev:#e8874e; --ok:#4cc08d; --warn:#d99a4e;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -18px rgba(0,0,0,.8);
  }}
}}
:root[data-theme="dark"] {{
  --ground:#0e1417; --surface:#151e22; --surface-2:#1b262b;
  --ink:#e4edf0; --muted:#93a5ac; --rule:#26333a;
  --fwd:#3fbfb4; --rev:#e8874e; --ok:#4cc08d; --warn:#d99a4e;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -18px rgba(0,0,0,.8);
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--ground); color:var(--ink);
  font-family:"Noto Sans JP", system-ui, sans-serif;
  font-size:15px; line-height:1.85; -webkit-font-smoothing:antialiased;
}}
.wrap {{ max-width:960px; margin:0 auto; padding:0 24px 96px; }}
h1, h2, h3 {{ font-family:"Zen Kaku Gothic New","Noto Sans JP",sans-serif;
  text-wrap:balance; line-height:1.3; }}
h1 {{ font-size:clamp(28px,4.4vw,42px); font-weight:900; margin:0 0 8px;
  letter-spacing:-.01em; }}
h2 {{ font-size:22px; font-weight:700; margin:0 0 4px; }}
h3 {{ font-size:16px; font-weight:700; margin:28px 0 6px; }}
p {{ margin:0 0 14px; max-width:70ch; }}
a {{ color:var(--fwd); }}
code, .mono, table, .num {{ font-family:"IBM Plex Mono", ui-monospace, monospace;
  font-variant-numeric:tabular-nums; }}
code {{ background:var(--surface-2); padding:1px 5px; border-radius:3px;
  font-size:.88em; }}

header {{ padding:56px 0 28px; border-bottom:1px solid var(--rule); }}
.eyebrow {{ font-family:"IBM Plex Mono",monospace; font-size:11px;
  letter-spacing:.16em; text-transform:uppercase; color:var(--muted);
  margin:0 0 14px; }}
.lede {{ font-size:17px; color:var(--muted); max-width:62ch; margin:10px 0 0; }}

.readout {{ display:grid; gap:1px; background:var(--rule);
  grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  border:1px solid var(--rule); border-radius:6px; overflow:hidden;
  margin:32px 0 0; }}
.readout div {{ background:var(--surface); padding:14px 16px; }}
.readout dt {{ font-family:"IBM Plex Mono",monospace; font-size:11px;
  letter-spacing:.06em; color:var(--muted); margin:0 0 6px; }}
.readout dd {{ margin:0; font-family:"IBM Plex Mono",monospace;
  font-size:24px; font-weight:600; font-variant-numeric:tabular-nums; }}
.readout small {{ display:block; font-size:11.5px; color:var(--muted);
  font-family:"Noto Sans JP",sans-serif; margin-top:2px; line-height:1.5; }}

section {{ padding:44px 0 0; }}
.sec-head {{ display:flex; align-items:baseline; gap:12px;
  border-bottom:1px solid var(--rule); padding-bottom:10px; margin-bottom:20px; }}
.sec-head .n {{ font-family:"IBM Plex Mono",monospace; font-size:11px;
  letter-spacing:.14em; color:var(--muted); padding-bottom:2px; }}

.viewer {{ background:var(--surface); border:1px solid var(--rule);
  border-radius:8px; box-shadow:var(--shadow); overflow:hidden; }}
.viewer canvas {{ display:block; width:100%; height:460px; touch-action:none;
  cursor:grab; background:
   radial-gradient(120% 90% at 50% 0%, var(--surface-2), var(--surface) 70%); }}
.viewer canvas:active {{ cursor:grabbing; }}
.bar {{ display:flex; flex-wrap:wrap; gap:8px; align-items:center;
  padding:12px 14px; border-top:1px solid var(--rule); background:var(--surface); }}
.chip {{ font:500 12.5px/1 "Noto Sans JP",sans-serif; padding:8px 13px;
  border:1px solid var(--rule); border-radius:999px; background:transparent;
  color:var(--muted); cursor:pointer; }}
.chip[aria-pressed="true"] {{ background:var(--ink); color:var(--ground);
  border-color:var(--ink); }}
.chip:focus-visible {{ outline:2px solid var(--fwd); outline-offset:2px; }}
.bar .spacer {{ flex:1 1 auto; }}
.stat {{ font-family:"IBM Plex Mono",monospace; font-size:11.5px;
  color:var(--muted); }}
.hint {{ font-size:12px; color:var(--muted); margin:10px 2px 0; }}

table {{ border-collapse:collapse; width:100%; font-size:13px; }}
.scroll {{ overflow-x:auto; border:1px solid var(--rule); border-radius:6px;
  background:var(--surface); }}
th, td {{ padding:9px 12px; text-align:left; border-bottom:1px solid var(--rule);
  white-space:nowrap; }}
th {{ font-size:11px; letter-spacing:.06em; text-transform:uppercase;
  color:var(--muted); font-weight:600; background:var(--surface-2); }}
tr:last-child td {{ border-bottom:none; }}
td.num {{ text-align:right; }}
td.strong {{ font-weight:600; }}
.tag {{ font-size:11px; padding:2px 7px; border-radius:3px;
  border:1px solid currentColor; }}
.tag.ok {{ color:var(--ok); }}
.tag.warn {{ color:var(--warn); }}
.tag.pending {{ color:var(--muted); }}

figure {{ margin:22px 0; }}
figure img {{ width:100%; max-width:100%; display:block; border-radius:6px;
  border:1px solid var(--rule); background:#fff; }}
figcaption {{ font-size:12.5px; color:var(--muted); margin-top:9px;
  max-width:70ch; }}

.cols {{ display:grid; gap:18px; grid-template-columns:1fr 1fr; }}
@media (max-width:720px) {{ .cols {{ grid-template-columns:1fr; }} }}
.card {{ background:var(--surface); border:1px solid var(--rule);
  border-radius:6px; padding:16px 18px; }}
.card h3 {{ margin:0 0 6px; font-size:14px; }}
.card p {{ margin:0; font-size:13px; color:var(--muted); }}
.bad {{ border-left:3px solid var(--rev); }}
.good {{ border-left:3px solid var(--fwd); }}

ul {{ padding-left:1.15em; max-width:70ch; }}
li {{ margin-bottom:7px; }}
.foot {{ margin-top:56px; padding-top:20px; border-top:1px solid var(--rule);
  font-size:12.5px; color:var(--muted); }}
@media (prefers-reduced-motion:reduce) {{ * {{ animation:none !important; }} }}
</style>

<div class="wrap">
<header>
  <p class="eyebrow">Gamboa 2005 の再現 — 2026-09-06</p>
  <h1>Tesla バルブ再現ノート</h1>
  <p class="lede">周期セルで再現できなかった原因を論文図の実測で突き止め、
  プレナム込みの単発モデルを組み直しました。3D モデルと現時点の数値をまとめています。</p>
  <dl class="readout">
    <div><dt>Di（Re = 100）</dt><dd>{di100}</dd>
      <small>Gamboa Fig.7 は 1.02 ± 0.02。周期版は 0.9925 で 1 を割っていました</small></div>
    <div><dt>形状の一致（IoU）</dt><dd>{iou}</dd>
      <small>Fig.2 の実形状との重ね合わせ</small></div>
    <div><dt>φ_loop 逆 / 順</dt><dd>{phi_ratio}</dd>
      <small>逆流のほうがループへ流れる。設計原理どおり</small></div>
  </dl>
</header>

<section>
  <div class="sec-head"><span class="n">モデル</span><h2>3D モデル</h2></div>
  <p>ドラッグで回転、ホイールで拡大縮小、ダブルクリックで初期姿勢に戻ります。
  寸法は流路幅 w_v = 1.0 mm、深さ 1.0 mm（正方形断面）。</p>
  <div class="viewer">
    <canvas id="cv"></canvas>
    <div class="bar" id="bar">
      <span class="spacer"></span>
      <span class="stat" id="stat"></span>
    </div>
  </div>
  <p class="hint">STL は <code>handover/cad/</code> にあります。すべて閉じたメッシュ
  （辺を共有する三角形がちょうど 2 枚）であることを確認済みです。
  板モデルは {plate}。</p>
</section>

<section>
  <div class="sec-head"><span class="n">診断</span><h2>何が間違っていたか</h2></div>
  <p>それまでは、バルブを周期配列にするために出口区間を「閉じ壁」に置き換え、
  主流路を接合部の先へまっすぐ延長していました。局所形状は保っているつもりでしたが、
  <strong>位相が変わっていました</strong>。</p>
  <div class="cols">
    <div class="card good"><h3>実形状（Fig.2 実測）</h3>
      <p>主流路は接合部で終わる。ループ下流枝と出口区間は同一直線
      （角度差 {coll_ang} 度、法線オフセット {coll_off} w_v）。
      逆流は<strong>折れ角ゼロでループへ入る</strong>。</p></div>
    <div class="card bad"><h3>周期版（閉じ壁）</h3>
      <p>主流路がまっすぐ続き、ループ下流枝は 48 度の側枝。
      逆流はループへ入るのに 132 度曲がる必要があり、
      <strong>主流路を素通りできる</strong>。</p></div>
  </div>
  <p style="margin-top:16px">結果、周期版では φ_loop が順流 4.0 % ＞ 逆流 2.9 %（Re = 500）と
  設計原理と逆になり、Di は Re とともに 1 を割って下がりました
  （0.9925 → 0.9678）。単発モデルではこれが逆転します。</p>
  <figure><img src="{fig_topology}" alt="Fig.2 の流体マスクと壁の直線当てはめ">
    <figcaption>Fig.2（997 × 714 px、再標本化なし）から取り出した流体領域。
    赤がループ下流枝、青が出口区間の壁の当てはめ。外壁は残差 0.18 px で 1 本の直線に載ります。</figcaption></figure>
</section>

<section>
  <div class="sec-head"><span class="n">形状</span><h2>組み直した形状</h2></div>
  <p>左プレナム → 入口区間 → ループ → 出口区間 → 右プレナムという Fig.1 / Fig.2 どおりの
  単発構成です。プレナム半径は実測 5.06 w_v（論文 R_p = 5）、左プレナムの中心は
  実測 (−5.464, −0.029) w_v で、幾何から予想される (−5.475, 0) と一致します。
  これはスケールの独立検証にもなっています。</p>
  <figure><img src="{fig_overlay}" alt="生成形状と Fig.2 の重ね合わせ">
    <figcaption>上: 生成形状（青は Fig.2 どおり、橙は CFD 用に出口プレナムの直線壁を鉛直にした版）。
    下: Fig.2 との重ね（青 = 一致、水色 = Fig.2 のみ、赤紫 = 生成のみ）。
    はみ出しは Fig.2 側 5.4 % / 生成側 2.0 % で、塗り潰しが線幅ぶん太る向きと一致します。</figcaption></figure>
</section>

<section>
  <div class="sec-head"><span class="n">結果</span><h2>いまの結果</h2></div>
  <div class="scroll"><table>
    <thead><tr><th>Re</th><th>格子</th><th>Δp 順 [Pa]</th><th>Δp 逆 [Pa]</th>
      <th>Di</th><th>Gamboa Fig.7</th><th>差</th></tr></thead>
    <tbody>{rows_html}</tbody>
  </table></div>
  <p class="hint">Δp はプレナムの直線壁どうしの流量重み平均圧力差。
  Di = Δp<sub>逆</sub> / Δp<sub>順</sub>（同一流量）。格子は流路幅あたりのセル数で、
  Gamboa 自身も 16 要素を使っています。</p>

  <h3>機構が働いているか</h3>
  <p>ループを通る質量流量の割合 φ_loop は、順流 {phi_f} %、逆流 {phi_r} %。
  比は {phi_ratio} で、<strong>逆流のほうがループへ流れて</strong>います。
  Gamboa 3.1 節が高 Di の前提とする条件です。周期版では 0.72 と逆でした。</p>
  <figure><img src="{fig_fields}" alt="順流・逆流の速度場と phi_loop">
    <figcaption>Re = 100、格子 16 cells/w_v。上が順流、中が逆流、下が断面ごとの φ_loop。
    逆流では出口区間から折れ角ゼロでループへ入り、戻り流路から主流に対向して吐き出されています。</figcaption></figure>
</section>

<section>
  <div class="sec-head"><span class="n">検証</span><h2>ソルバの検証</h2></div>
  <p>単発モデルには入口・出口境界条件が要ります（周期化できないため）。
  Zou-He で組み直し、平行平板の解析解で確認しました。</p>
  <div class="scroll"><table>
    <thead><tr><th>Re</th><th>格子</th><th>τ_p</th><th>u_max/u_avg</th>
      <th>f·Re（補正後）</th><th>判定</th></tr></thead>
    <tbody>
      <tr><td class="num">100</td><td class="num">16</td><td class="num">0.5480</td>
        <td class="num">1.48916</td><td class="num strong">95.690</td>
        <td><span class="tag ok">合格</span></td></tr>
      <tr><td class="num">300</td><td class="num">16</td><td class="num">0.5160</td>
        <td class="num">1.48907</td><td class="num strong">95.694</td>
        <td><span class="tag ok">合格</span></td></tr>
      <tr><td class="num">400</td><td class="num">16</td><td class="num">0.5120</td>
        <td class="num">1.48908</td><td class="num strong">95.691</td>
        <td><span class="tag ok">合格</span></td></tr>
      <tr><td class="num">500</td><td class="num">16</td><td class="num">0.5096</td>
        <td class="num">—</td><td class="num">—</td>
        <td><span class="tag warn">発散</span></td></tr>
      <tr><td class="num">500</td><td class="num">20</td><td class="num">0.5120</td>
        <td class="num">1.49257</td><td class="num strong">95.824</td>
        <td><span class="tag ok">合格</span></td></tr>
    </tbody>
  </table></div>
  <p class="hint">離散理論値は 16 セルで約 95.6、20 セルで 95.76。誤差 0.1 %。
  発散は τ が小さいこと自体ではなく、<strong>セル Reynolds 数</strong>
  （Re<sub>cell</sub> = 0.75 Re / 格子）が 20 を超えると起きます。
  同じ τ でも周期版（体積力駆動）は安定でした。</p>
</section>

<section>
  <div class="sec-head"><span class="n">この先</span><h2>次にやること</h2></div>
  <ul>
    <li>Re = 300 と Re = 500 を出して Di(Re) 曲線を Fig.7 と比べる。ここが再現の判定</li>
    <li>格子依存（16 → 24）と収束判定の確認。Gamboa 自身は「要素数倍増で 4 % 未満」と報告</li>
    <li>7 月形状も同じ単発枠組みで解き、Gamboa と同条件で比較する</li>
    <li>3 次元は当面やらない。Gamboa のベンチマーク自体が 2D なので不要で、
      同じ格子解像度なら計算量が約 20 倍になる。3D が要るのは Porwal（正方形断面）と
      Han（共役伝熱）の段階。そのときは今回の STL を外部ソルバへ渡すのが現実的</li>
    <li>熱はまだ入れない</li>
  </ul>
</section>

<p class="foot">Gamboa, A. R., Morris, C. J., Forster, F. K. (2005),
“Improvements in Fixed-Valve Micropump Performance Through Shape Optimization of Valves,”
ASME J. Fluids Eng. 127(2), 339–346.<br>
数値はすべて自前の計算・実測です。論文記載値と自分の計算値は本文で区別しています。</p>
</div>

<script>
const MODELS = {models_json};
const cv = document.getElementById('cv');
const ctx = cv.getContext('2d', {{ alpha:false }});
const bar = document.getElementById('bar');
const stat = document.getElementById('stat');
let cur = 0, yaw = 0.38, pitch = -0.42, zoom = 1, spin = true, drag = null;

function decode(b64) {{
  const bin = atob(b64), u = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
  return new Float32Array(u.buffer);
}}
MODELS.forEach(m => {{
  const v = decode(m.data);
  let lo = [1e9, 1e9, 1e9], hi = [-1e9, -1e9, -1e9];
  for (let i = 0; i < v.length; i += 3)
    for (let k = 0; k < 3; k++) {{
      if (v[i+k] < lo[k]) lo[k] = v[i+k];
      if (v[i+k] > hi[k]) hi[k] = v[i+k];
    }}
  const c = [0,1,2].map(k => (lo[k]+hi[k])/2);
  const s = Math.max(hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]);
  for (let i = 0; i < v.length; i += 3)
    for (let k = 0; k < 3; k++) v[i+k] = (v[i+k]-c[k]) / s;
  m.v = v; m.n = v.length / 9;
  m.tmp = new Float32Array(v.length);
  m.order = new Int32Array(m.n);
  m.zc = new Float32Array(m.n);
}});

MODELS.forEach((m, i) => {{
  const b = document.createElement('button');
  b.className = 'chip'; b.textContent = m.label;
  b.setAttribute('aria-pressed', i === 0 ? 'true' : 'false');
  b.onclick = () => {{ cur = i; sync(); }};
  bar.insertBefore(b, bar.querySelector('.spacer'));
}});
function sync() {{
  [...bar.querySelectorAll('.chip')].forEach((b, i) =>
    b.setAttribute('aria-pressed', i === cur ? 'true' : 'false'));
  const m = MODELS[cur];
  stat.textContent = m.n.toLocaleString() + ' 三角形 · ' + m.note;
  draw();
}}

function resize() {{
  const r = cv.getBoundingClientRect(), d = Math.min(devicePixelRatio || 1, 2);
  cv.width = Math.round(r.width * d); cv.height = Math.round(r.height * d);
  draw();
}}
addEventListener('resize', resize);

function css(v) {{ return getComputedStyle(document.body).getPropertyValue(v).trim(); }}

function draw() {{
  const m = MODELS[cur], W = cv.width, H = cv.height;
  ctx.fillStyle = css('--surface'); ctx.fillRect(0, 0, W, H);
  const cy = Math.cos(yaw), sy = Math.sin(yaw);
  const cp = Math.cos(pitch), sp = Math.sin(pitch);
  const S = Math.min(W, H) * 0.78 * zoom, ox = W/2, oy = H/2;
  const v = m.v, t = m.tmp;
  for (let i = 0; i < v.length; i += 3) {{
    const x = v[i], y = v[i+1], z = v[i+2];
    const X = cy*x + sy*z, Z = -sy*x + cy*z;
    const Y = cp*y - sp*Z, Z2 = sp*y + cp*Z;
    t[i] = ox + X*S; t[i+1] = oy - Y*S; t[i+2] = Z2;
  }}
  const n = m.n, zc = m.zc, order = m.order;
  let cnt = 0;
  for (let f = 0; f < n; f++) {{
    const a = f*9;
    const ax = t[a], ay = t[a+1], bx = t[a+3], by = t[a+4], cx2 = t[a+6], cy2 = t[a+7];
    if ((bx-ax)*(cy2-ay) - (by-ay)*(cx2-ax) >= 0) continue;   // 裏面を捨てる
    zc[cnt] = t[a+2] + t[a+5] + t[a+8]; order[cnt] = f; cnt++;
  }}
  const pairs = [];
  for (let i = 0; i < cnt; i++) pairs.push([zc[i], order[i]]);
  pairs.sort((p, q) => p[0] - q[0]);
  const base = css('--fwd');
  const rgb = hexToRgb(base);
  for (let i = 0; i < pairs.length; i++) {{
    const f = pairs[i][1], a = f*9;
    const ux = t[a+3]-t[a], uy = t[a+4]-t[a+1], uz = t[a+5]-t[a+2];
    const wx = t[a+6]-t[a], wy = t[a+7]-t[a+1], wz = t[a+8]-t[a+2];
    let nx = uy*wz - uz*wy, ny = uz*wx - ux*wz, nz = ux*wy - uy*wx;
    const L = Math.hypot(nx, ny, nz) || 1; nx/=L; ny/=L; nz/=L;
    // canvas は y が下向きで左手系になるので、法線は視点側（nz > 0）へ揃える
    if (nz < 0) {{ nx = -nx; ny = -ny; nz = -nz; }}
    const lam = Math.max(0, 0.42*nx - 0.55*ny + 0.72*nz);
    const k = Math.min(1, 0.24 + 1.02*lam);
    ctx.fillStyle = 'rgb(' + Math.round(rgb[0]*k) + ',' + Math.round(rgb[1]*k)
      + ',' + Math.round(rgb[2]*k) + ')';
    ctx.beginPath();
    ctx.moveTo(t[a], t[a+1]); ctx.lineTo(t[a+3], t[a+4]); ctx.lineTo(t[a+6], t[a+7]);
    ctx.closePath(); ctx.fill();
  }}
}}
function hexToRgb(h) {{
  h = h.replace('#','');
  if (h.length === 3) h = h.split('').map(c => c+c).join('');
  const v = parseInt(h, 16);
  return [(v>>16)&255, (v>>8)&255, v&255];
}}

cv.addEventListener('pointerdown', e => {{
  drag = [e.clientX, e.clientY]; spin = false; cv.setPointerCapture(e.pointerId);
}});
cv.addEventListener('pointermove', e => {{
  if (!drag) return;
  yaw += (e.clientX - drag[0]) * 0.01;
  pitch = Math.max(-1.55, Math.min(1.55, pitch + (e.clientY - drag[1]) * 0.01));
  drag = [e.clientX, e.clientY]; draw();
}});
addEventListener('pointerup', () => {{ drag = null; }});
cv.addEventListener('wheel', e => {{
  e.preventDefault();
  zoom = Math.max(0.4, Math.min(6, zoom * (e.deltaY < 0 ? 1.12 : 0.89))); draw();
}}, {{ passive:false }});
cv.addEventListener('dblclick', () => {{
  yaw = 0.38; pitch = -0.42; zoom = 1; spin = true; draw();
}});

const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
function tick() {{
  if (spin && !reduce) {{ yaw += 0.0045; draw(); }}
  requestAnimationFrame(tick);
}}
resize(); sync(); tick();
</script>
"""


if __name__ == "__main__":
    main()
