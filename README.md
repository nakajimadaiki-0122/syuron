# 修論 — Tesla 型マイクロ流路

半導体冷却向けマイクロ流路について、Tesla 型構造を出発点に、
先行研究の再現 → ギャップの定量化 → ボトルネックの特定 → 新規構造の提案、
という順序で進めています。

## どこから読むか

1. [目標.md](目標.md) — 研究全体の 5 段階
2. [handover/HANDOFF.md](handover/HANDOFF.md) — **現在地・地図・次の一手。まずこれ**
3. [handover/progress/](handover/progress/) — セッションごとの作業記録（最新が現状）
4. [handover/README.md](handover/README.md) — 技術詳細（定義・踏んだバグ・全結果）

作業はすべて [handover/](handover/) の中で行います。

## いまの段階

目標.md の第 1 段階（既存の数値を再現できるか）と第 2 段階（その数値でモデル作成）の境目です。

- Gamboa 2005 の幾何は復元済み（Fig.2 と IoU 0.9285）
- 単発モデル（プレナム込み）で Re = 100 の Di = 1.0038（論文 1.02 ± 0.02）
- Re = 300 / 500 を計算中。ここが再現の判定
- CAD 用の DXF / STL は [handover/cad/](handover/cad/)、閲覧用ページは
  [handover/report/gamboa_report.html](handover/report/gamboa_report.html)

熱解析・新構造へはベンチマークが閉じるまで進みません（理由は handover/README.md §3）。
