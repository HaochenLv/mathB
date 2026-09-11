# Q2 FINAL REVIEW

## 封板结论

**Q2 可以封板。**

理由：模型口径未变；Q1 exact geometry 始终作为最终裁判；正交修正、目标失效结构、局部坐标与 epsilon 收敛、目标权重敏感性和 Pareto baseline 均已完成；结论在多套错位网格上保持一致，没有发现 Q1 bug 或会推翻 Q2 主结论的证据。

## 最终交付

- 推荐代表点：`S2=(691,475) m`，移动 838.514 m。
- 推荐稳定区：cross-grid sampled coverage 不低于最大值 99.5% 的 `x=686..691 m, y=430..490 m`。
- 更严格核心：99.95% sampled maximum 的 `x=688..691 m, y=466..479 m`。
- robust clear-ready 面积覆盖：ultra offset grid 25.5213%，另一高密度 grid 25.8785%；最终表述为约 25.5%–25.9%。
- 备选中间点：`(648,386) m`，移动 754.255 m，高密度面积覆盖 25.0964%。
- 备选短距离点：`(620,320) m`，移动 697.711 m，高密度面积覆盖 23.5690%。
- 正交 baseline：`(806.98,399.94) m`，移动 900.650 m，高密度面积覆盖 22.6433%。

## A–D 完成检查

- A：35 个 `(d1,d2)` 组合、0.01° beta/epsilon 精细验证完成。robust sampled optimum 全为 91°；nominal optimum 90.03°–90.13°；91°收益仅 0.243%–0.892%，偏移主要由 epsilon 端点 minimax 导致，而非粗网格或状态切换。
- B：推荐点完整 target-space robust/fail、worst diameter、guaranteed receive 和 failure reason 图/CSV 已生成。成功集为两个主要径向带；72.62% 目标面积因直径超限失败，1.50% 因圆覆盖失败。
- C：2→1→0.5 m 局部细化、目标网格错位验证、上边界补扫和 ultra validation 完成。epsilon step 0.1→0.01°不改变候选覆盖分类。
- D：最大面积、90%阈值最短、中间 Pareto、正交和 global uniform-weight baseline 已在同一高密度网格比较。面积权重和 `(r,alpha)` 等权会产生明显不同最优区。

## 一致性审核

- 数学口径：与任务书一致；没有重复扫描第一次误差，没有修改 clear-ready。
- 接收口径：使用 `R_min=max(1000,r)`；near 不进入 Q1。
- 数据口径：primary coverage 使用权重 r；uniform 指标明确仅为敏感性。
- 优化/验证分离：包含原优化网格、错开 validation、cross-grid、边界扩展和 ultra offset validation。
- 结论等级：90°指导与已有空集冲突推理单独表述；坐标、覆盖率、条带边界均标为有限网格数值结论。
- Q1：未修改核心源码；没有发现需报告的反例。
- 工程：长扫描提供签名 cache/resume；所有配置集中；无随机过程。
- 测试：完整 pytest 为 25 passed、0 failed。

## 已知局限

- 25.5%–25.9% 是目标集合的面积测度，不是干扰源出现概率。
- 推荐坐标在不同目标网格间有数十米级 y 方向漂移，因此只推荐稳定区域和代表点，不宣称连续域唯一最优。
- 径向成功带和 91°修正尚未被证明为连续域解析定理。
- 本轮只解决固定第二检测点的 Q2 几何，不包含后续调度或路径。

## Git 与下一步

这是一组值得保留的稳定结果，建议作为 Q2 封板 checkpoint 提交。提交前应保留 finalization 代码、测试、报告、关键 CSV 和论文候选图，并检查是否需要排除体积过大的逐点 target-map CSV。未经用户明确授权不得 push。

用户确认 Q2 封板及 Git 范围之前，不开始 Q3。
