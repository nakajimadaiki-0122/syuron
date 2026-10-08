# 文献調査: 半導体冷却 × 流路形状 × Tesla バルブ

作成日: 2026-09-26

対象を「半導体（チップ・パワーデバイス）の冷却 → 流路形状による強化 → Tesla バルブ」に
絞った Web 検索の結果。**すべて [Web検索・未読] であり、原典は 1 本も確認していない。**
引用する前に必ず本文を読み、`literature_findings.md` に `[原典確認]` で移すこと。

既存の中核文献（Gamboa 2005 / Porwal 2016 / Han 2024）は `literature_findings.md` 参照。

---

## 0. 全体像（[推測]）

Tesla バルブを半導体冷却に使う研究は 2021 年頃から急増し、大きく 3 群に分かれる。

| 群 | 何をしているか | 本研究との関係 |
|---|---|---|
| **A. 単相・形状研究** | Tesla 型マイクロチャネルの Nu・f・PEC を通常流路や他形状と比較。形状パラメータを振る | **本研究の直接の土俵。** Han 2024 はここ |
| **B. チップ実装** | マニホールド・埋め込み・チップレットなど実装系に Tesla パターンを組み込み、R_th を実測 | 「あるべき姿」の設定と、ボトルネック（R_cond / R_cal）議論の材料 |
| **C. 二相（沸騰）** | 蒸気の逆流抑制・圧力振動抑制・CHF 向上 | Tesla の実績が主に二相にあることの裏付け。**本研究の対象外**だが、区別して書く必要がある |

**[推測]** A 群で「接続角・ループ幅・こぶ半径」を系統的に振った論文は少なく、
Han 2024 と対称型（2025）が中心。Gamboa の β を伝熱の文脈で検証した論文は
検索では見つからなかった。ここが本研究の位置になりうる。

---

## A. 単相・形状研究（本研究の直接の比較対象）

| 文献 | 内容（検索結果の要約） | 見るべき点 |
|---|---|---|
| Han et al. (2024) *A comparative study of enhanced thermal performance in Tesla-type microchannels*, ATE 239, 122157 | 既存中核文献。Nu 10–25（Re 100–450）、Di_p ≈ 1.9、Di_t ≈ 1.3 | `literature_findings.md` §5 |
| *Numerical study on heat transfer and flow characteristics of symmetric Tesla-type microchannel heat sinks*, ATE (2025) [S1359431124022798](https://www.sciencedirect.com/science/article/abs/pii/S1359431124022798) | Han 系の続編と思われる。対称型 Tesla（SYMTVM）。**接続角 30° の逆流配置**で Nu が矩形流路の 5 倍、総合性能 2.59 倍 | **接続角を振っている**可能性。V1（β）の直接の比較対象になりうる |
| *A comprehensive analysis of flow and heat transfer performance in a novel Tesla valve microchannel* [ResearchGate](https://www.researchgate.net/publication/399653824_A_comprehensive_analysis_of_flow_and_heat_transfer_performance_in_a_novel_Tesla_valve_microchannel) | 新型 Tesla マイクロチャネルの流動・伝熱 | 形状定義を確認 |
| *Study on thermal-hydraulic characteristics of tesla valve-type microchannels in liquid-cooled plates: Experiment and numerical simulation* (2025) [S0735193325011480](https://www.sciencedirect.com/science/article/abs/pii/S0735193325011480) | 実験 + 数値。Di・熱的 Di（TDi）・PEC で評価 | 実測に進むときの手順の参照。Di と TDi の定義 |
| *Multi-objective optimization and ANN models for … microchannel heat sink with fins inspired Tesla valve profile* (2024) [S2214157X24010049](https://www.sciencedirect.com/science/article/pii/S2214157X24010049) | Tesla 形状のフィンを多目的最適化 | 設計変数の取り方 |
| *Design and optimization of a novel tesla valve air-cooled heat sink* ATE (2025) [S1359431125025700](https://www.sciencedirect.com/science/article/abs/pii/S1359431125025700) | 空冷ヒートシンクへの応用 | 対象外だが形状最適化の手法は参考 |
| *Performance analysis on the liquid cooling plate with the new Tesla valve capillary channel based on the fluid solid coupling simulation* ATE (2023) [S1359431123010062](https://www.sciencedirect.com/science/article/abs/pii/S1359431123010062) | 流体–固体連成（CHT）で液冷板を評価 | R_cond を含めた評価の先行例 |
| Porwal & Thompson (2018?) *Heat transfer and fluid flow characteristics in multistaged Tesla valves* [ResearchGate](https://www.researchgate.net/publication/324032178_Heat_transfer_and_fluid_flow_characteristics_in_multistaged_Tesla_valves) | Porwal 修論の論文版と思われる | **修論より査読誌を引用する。** 書誌を確定させること |

---

## B. チップ実装（マニホールド・埋め込み・パワーデバイス）

| 文献 | 内容（検索結果の要約） | 本研究との関係 |
|---|---|---|
| *Experiment on enhanced heat transfer of embedded manifold Tesla-patterned microchannel heat sink*, ICHMT (2025) [S0735193325004762](https://www.sciencedirect.com/science/article/abs/pii/S0735193325004762) | **マニホールド × Tesla パターン**を製作・実測。ピンフィン・矩形と比較し、h が +55〜88 %、R_th が −22〜32 %、チップ表面 −20 °C。単相・二相の両方 | `literature_findings.md` §1「分野の到達解はマニホールド化」と Tesla を結ぶ。**ボトルネック議論の中心文献候補** |
| *Thermal and hydraulic optimization of embedded manifold microchannel heat sink for chip cooling based on Tesla-pattern*, ATE (2025) [S1359431125029461](https://www.sciencedirect.com/science/article/abs/pii/S1359431125029461) | 同グループの最適化版。R_th 0.21–0.44 K/W、ポンプ動力 0.032–0.57 W。AI チップ向け。楕円基底 NN で設計変数→応答を予測。Nu +13.8 / 56.8 / 161 % | 設計変数の選び方と応答量の定義 |
| *High-performance near-substrate heat sink with Tesla-like rotor-wing microchannel for chiplet cooling application*, ATE (2024) [S1359431124024438](https://www.sciencedirect.com/science/article/abs/pii/S1359431124024438) | 3D-IC 埋め込み。R_th 0.02 K/W。非一様アスペクト比でホットスポット −16.9 % | 埋め込み型での R_cond の扱い |
| *Numerical investigation of tesla valve-based liquid cooling system with diamond substrate and TIM for high-power GaN device* (2026) [S2590123026032822](https://www.sciencedirect.com/science/article/pii/S2590123026032822) | GaN + ダイヤモンド基板 + Tesla 流路。逆流配置でチップ最高温度 −8〜17 °C | **R_cond（基板・TIM）と流路を同時に扱った例。** 支配抵抗の議論に直接使える |
| *The 3D Tesla Valve Manifold Two-Phase Cold Plate for the SiC Power Device Module*, IEEE (2022) [ieeexplore 9899661](https://ieeexplore.ieee.org/abstract/document/9899661/) | SiC パワーモジュール向け二重層マニホールド二相コールドプレート。順逆多段 Tesla を交互配置。375 W/cm² | パワーデバイス側の実装例（二相） |
| *Effect of the Tesla Valve on the heat transfer performance and the suppression of pressure drop oscillation in a liquid cooling loop*, IJTS (2024) [S1290072924004782](https://www.sciencedirect.com/science/article/abs/pii/S1290072924004782) | 冷却ループ全体での圧力振動抑制 | システム側の効果 |

---

## C. 二相（沸騰）— 対象外だが区別のために把握

| 文献 | 要点 |
|---|---|
| *The role of Tesla valves in microchannel flow boiling*, IJHMT (2024) [S0017931024009785](https://www.sciencedirect.com/science/article/abs/pii/S0017931024009785) | 総説的。逆流抑制の機構 |
| *Enhance flow boiling in Tesla-type microchannels by inhibiting two-phase backflow*, IJHMT (2023) [S0017931023006166](https://www.sciencedirect.com/science/article/abs/pii/S0017931023006166) | Han 系。CHF +26.2 %、HTC +120 %、壁温の標準偏差 1/10 |
| *Flow boiling in a relatively large copper heat sink comprised of Tesla microchannels*, IJHMT (2024) [S0017931024011955](https://www.sciencedirect.com/science/article/abs/pii/S0017931024011955) | IGBT 向け大型銅ヒートシンク。1000 W/cm² |
| *Thermo-hydrodynamic performances of flow boiling in Tesla-microchannels*, IJHMT (2025) [S0017931025010245](https://www.sciencedirect.com/science/article/abs/pii/S0017931025010245) | |
| *Flow boiling in parallel copper microchannels with asymmetric Tesla valves*, ATE (2026) [S1359431126012809](https://www.sciencedirect.com/science/article/abs/pii/S1359431126012809) | 非対称配置 |
| *Enhanced flow boiling by manipulating two-phase flow in Tesla channel heat sink using HFE-7100*, IJTS (2025) [S1290072924006938](https://www.sciencedirect.com/science/article/abs/pii/S1290072924006938) | 誘電性冷媒 |
| *Research on single-phase flow and two-phase flow boiling cooling performance of microchannel TMS with novel Tesla Valve design*, IJHMT (2024) [S001793102400591X](https://www.sciencedirect.com/science/article/abs/pii/S001793102400591X) | **単相と二相を同じ形状で比較**。単相での優位性の有無を見るのに有用 |

---

## D. 隣接領域（電池・機構論）

半導体ではないが、形状パラメータの振り方の参考になる。

| 文献 | 要点 |
|---|---|
| *A numerical analysis on multi-stage Tesla valve based cold plate for cooling of pouch type Li-ion batteries*, IJHMT (2021) [S0017931021006633](https://www.sciencedirect.com/science/article/abs/pii/S0017931021006633) | **段数・外側曲率半径・ユニットの向き**を振っている。V2（R）の参考 |
| *Augmentation of multi-stage Tesla valve design cold plate with reverse flow …*, IJHMT (2023) [S0017931023005847](https://www.sciencedirect.com/science/article/abs/pii/S0017931023005847) | 逆流配置 |
| *Performance optimisation of Tesla valve-type channel for cooling lithium-ion batteries*, ATE (2022) [S1359431122005312](https://www.sciencedirect.com/science/article/abs/pii/S1359431122005312) | |
| *Numerical Analysis and Multi-Objective Optimization of Tesla Valve Cold Plates for Li-Ion BTMS* (2026) [S1555256X26000512](https://www.sciencedirect.com/org/science/article/pii/S1555256X26000512) | |
| *Numerical Study of Diodicity Mechanism in Different Tesla-Type Microvalves* [S1665642313715943](https://www.sciencedirect.com/science/article/pii/S1665642313715943) | 機構論。側枝角・内側曲線接線角・幅を設計変数に。**本研究の 3 変数と同じ軸** |
| *Numerical calculation of forward and reverse flow in Tesla valves with different longitudinal width-to-narrow ratios*, Sci. Rep. (2023) [s41598-023-39758-3](https://www.nature.com/articles/s41598-023-39758-3) | 幅比を振った層流 CFD。V3 の参考 |
| *Research Progress and Application Review of Tesla Valve* [drpress](https://drpress.org/ojs/index.php/ajst/article/download/29618/29059/43138) | 総説（出版元の質は要確認） |

---

## 次にやること

1. **A 群の対称型（2025）と B 群のマニホールド Tesla（2025）を最優先で読む。**
   前者は接続角を振っている可能性があり V1 の比較対象、後者はボトルネック議論の中心
2. Porwal & Thompson の論文版の書誌を確定し、修論の代わりに引用する
3. IEEE Xplore（ITherm / IEEE TCPMT）で `Tesla` を直接検索する。
   Google 検索では半導体パッケージ系の会議論文がほとんど拾えなかった
4. 読んだものは `literature_findings.md` に `[原典確認]` で移す
