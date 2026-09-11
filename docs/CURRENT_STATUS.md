# CURRENT STATUS

## 当前阶段

Q1 已冻结；Q2 FINALIZATION ROUND 已完成并通过 final review，结论为 **Q2 可以封板**。尚未进入 Q3/Q4、频道调度或总路径。当前等待用户审核 Q2 结果与 Git 提交范围。

## 统一模型口径

- `S1=(0,0)`、第一次示向读数 `0°`、误差 `±1°`。
- 目标 `g=(r cos alpha,r sin alpha)`，`5<r<=1500 m`、`alpha in [-1,+1]°`；alpha 已包含第一次误差，不重复扫描。
- 第二次误差 `epsilon2 in [-1,+1]°`；全部两测定位调用 Q1 exact geometry。
- clear-ready：区域非空有界、diameter `<=40 m`、至少一个最大直径点对圆覆盖整个区域。
- `R_min=max(1000,r)`；`||S2-g||<=R_min` 才是 guaranteed receive。5 m 内为 near/direct-clear，不构造第二方向角。
- primary coverage 按面积元 `r dr d alpha` 加权；uniform `(r,alpha)` 仅作敏感性，不是物理面积或成功概率。

## 已确认结论

- Q1 核心行为保持冻结，未发现新 bug。
- 35 组 `d1×d2` exact 实验中，robust sampled minimax beta 均为 91.00°；nominal optimum 为 90.03°–90.13°。相对 beta=90° 的最坏直径改善仅 0.243%–0.892%。偏移主要来自 epsilon 端点 minimax，不是 0.05°网格误差或 unbounded 状态切换。论文应写“90°解析指导 + 约1°数值修正”，不得写成解析恒等式。
- 完整未知目标下不存在严格万能固定 S2；既有 near/far 回归 2564 点共同可行数仍为 0。
- 最终推荐代表点：`S2=(691,475) m`，移动 838.514 m。
- 稳定候选区：cross-grid sampled coverage `>=99.5%` 最大值时，`x=686..691 m, y=430..490 m`；99.95%核心为 `x=688..691 m, y=466..479 m`。向上补扫到 y=520 后峰值仍在内部。
- 推荐点 ultra offset validation 面积覆盖 25.5213%；另一高密度网格为 25.8785%。最终可信表述：约 25.5%–25.9%，不是连续域精确上界。
- epsilon step 0.1/0.05/0.02/0.01°下，所有重点候选的 robust coverage 分类不变。
- 推荐点 target-space 成功集主要为两个径向带：`r=275..662.5 m` 与 `700..825 m`；`675..687.5 m` 是窄的圆覆盖失败带。面积上约 72.62% 因 diameter>40 失败，1.50% 因 diameter-circle failure，25.88% clear-ready；guaranteed receive failure、unbounded、empty 均为 0。
- Pareto 备选：中间点 `(648,386)`，移动 754.255 m、面积覆盖 25.0964%；90%阈值最短点 `(620,320)`，移动 697.711 m、面积覆盖 23.5690%。
- 正交 baseline `(806.98,399.94)` 面积覆盖 22.6433%，直径统计较好但移动更长；正交只能解释大致走廊，不能替代全目标集优化。
- 权重敏感：uniform 最佳约 `(550,400)`，uniform coverage 47.19%，但面积覆盖仅 22.88%；两种权重的高值区明显不同，主结论必须使用面积权重。

## 最近实验与入口

- 主封板入口：`experiments/q2/run_q2_finalization.py`。
- 前置关键验证入口：`experiments/q2/run_q2_key_validation.py`。
- 最终报告：`results/q2/Q2_EXPERIMENT_REPORT.md`。
- 封板复核：`results/q2/Q2_FINAL_REVIEW.md`。
- 主要数据：orthogonal correction、2/1/0.5 m refinement、cross-grid/upper-extension、ultra validation、epsilon convergence、target-space map/summary、weighting sensitivity、top/Pareto CSV，均在 `results/q2/`。
- 论文候选图：`fig_q2_orthogonal_correction.png`、`fig_q2_coverage_heatmap.png`、`fig_q2_candidate_region_final.png`、`fig_q2_pareto_final.png`、`fig_q2_target_space_failure.png`。
- 长阶段使用签名 meta cache；脚本完全确定性，无随机 seed 依赖。

## 已知局限 / 风险

- 覆盖率是目标集合几何面积测度，不是概率。
- 推荐点随目标网格有小幅漂移，因此封板对象是稳定区域与代表点，不是唯一连续域最优点。
- 91°修正和目标径向带是高精度有限网格结论，尚非解析定理。
- `q2_target_space_map.csv` 为逐点解释数据，提交前需要检查体积并决定是否纳入版本控制；不要擅自删除。

## 测试与 Git

- 本轮新增 `tests/q2/test_finalization.py`，覆盖 near 优先级、guaranteed receive failure 优先级和正常几何分类入口。
- 最新完整 pytest：25 passed、0 failed（既有 22 项 + 本轮新增 3 项）。
- 当前基线 commit 为 `ece84e4`；Q2 工作尚未 commit、未 push。

## 下一步

1. 向用户提交“Q2 可以封板”的简要结论并建议一次 checkpoint commit。
2. 用户确认后再确定 checkpoint 文件范围；未经确认，不 commit、不 push。
3. 未经用户确认，不开始 Q3。
