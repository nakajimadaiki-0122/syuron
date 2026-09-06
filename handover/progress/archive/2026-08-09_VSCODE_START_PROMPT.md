# VS Code Chatに最初に渡すプロンプト

このリポジトリはTesla型マイクロチャネル研究の解析コードです。まず `HANDOFF.md` と `configs/benchmark_specs.yaml` と `results/current_state.json` を読んでください。

重要な制約:

1. 現時点では熱解析・新構造設計へ進まないでください。
2. 最優先は2D flow solverのbenchmarkです。
3. straight validation → grid check → Gamboa 2005 の順で進めてください。
4. Gamboaの論文値を勝手に推定しないでください。Re=100のDi targetが必要ならPDF Fig.7を確認・digitizeしてください。
5. mass flowは `rho*u_n` の断面積分で評価してください。
6. `np.roll` による意図しない周期境界が再発していないか確認してください。
7. 変更前に既存コードを読んで、最小変更で進めてください。
8. 計算ごとに `templates/case_summary.example.json` 形式で結果を保存してください。
9. Gamboa benchmarkが通るまでHan/Porwalの熱解析へ進まないでください。

最初のタスク:
- straight caseで `u_max/u_avg`, `f_D*Re`, mass imbalance を再計算し、合否を報告してください。
- 次に7月Tesla形状を12/24 cells/mmでRe=100のみ比較してください。
- その後Gamboa optimized geometryの実装に入ってください。
