# Q2 几何实验最终报告

## 1. 范围、模型与证据等级

本报告只研究第二检测点 `S2` 的几何选择，不涉及频道调度、总路径或 Q3/Q4。所有两测定位结论均直接调用冻结的 Q1 exact geometry：角域半平面求交、区域状态、polygon diameter、全部最大直径点对及其圆覆盖判定。实验中没有以条带近似替代 Q1 真值。

统一模型为：

- `S1=(0,0)`，第一次读数 `theta1=0 deg`，误差 `delta=1 deg`；
- 真实目标 `g=(r cos(alpha),r sin(alpha))`，其中 `5<r<=1500 m`、`alpha in [-1,+1] deg`；
- `alpha` 已包含第一次测向误差，不再重复扫描；
- 第二次读数中心为 `atan2(g_y-s_y,g_x-s_x)+epsilon2`，其中 `epsilon2 in [-1,+1] deg`；
- `R_min(g)=max(1000,r)`，`||s-g||<=R_min(g)` 才是 guaranteed receive；
- `||s-g||<=5 m` 记为 near/direct-clear，不构造第二方向角；
- clear-ready 必须同时满足：区域非空有界、直径不超过 40 m、至少一个以最大直径点对为直径的圆覆盖整个区域。

本文严格区分：Q1 对单个已采样情景的结果是 exact geometry 判定；候选坐标、覆盖率及稳定区域是有限网格数值结论；面积覆盖率是目标位置集合的几何测度，不是真实成功概率。

## 2. A：91° 偏移的最终解释

选择 `d1={200,400,600,800,1000,1200,1500} m` 和 `d2={200,400,600,800,1000} m`，共 35 组实际有意义的视距组合。beta 先以 0.05 deg 扫描正交邻域，再以 0.01 deg 局部细化；robust epsilon 最终步长为 0.01 deg，共 201 个误差值。

结果高度一致：

- 35/35 组的 sampled robust minimax 点均为 `beta=91.00 deg`；
- 不考虑第二次误差、只取 `epsilon2=0` 时，nominal optimum 位于 `90.03..90.13 deg`；
- beta=90 与 robust optimum 的最坏直径差仅为 0.243%–0.892%；
- robust 曲线在 91° 附近由 `epsilon2=-1 deg` 或 `+1 deg` 端点控制，35 组中最坏端点分别出现 11 次和 24 次；
- 所有组合在最优点均保持 bounded，没有几何状态切换；
- 20/35 组对全部 epsilon 满足直径不超过 40 m，但 0/35 组对全部 epsilon 满足最大直径圆覆盖，因此 0/35 组完整 clear-ready。

由此可以排除“0.05° beta 网格误差”和“unbounded 状态切换”这两个解释。nominal optimum 仍非常接近 90°，而 robust minimax 在对完整 `[-1,+1]°` 取最大值后稳定移到 91°附近，说明主要修正来自第二次误差端点的 minimax 平衡，并叠加有限角宽多边形造成的很小非对称。

论文建议表述为：90°是解析几何指导；在当前有限角宽与 Q1 clear-ready 口径下，exact numerical minimax 出现约 +1° 修正，但收益不足 0.9%。不能把“恒等于 91°”写成未经证明的解析定理。

## 3. B/C：目标空间结构与候选区收敛

### 3.1 搜索与交叉验证

在此前完整大矩形粗扫基础上，对 `x=630..730 m, y=360..480 m` 完成：

- 2 m 全矩形网格 3111 点；
- 面积权重与均匀权重共同选种子的 1 m 网格 569 点；
- 0.5 m 网格 569 点；
- 目标网格错位的 cross-grid：143 个粗点、向 `y=520 m` 补扫的 88 个边界点和 168 个 1 m 精细点；
- 7 个候选的 ultra validation：38240 个目标样本，每目标 101 个 epsilon 值。

0.5 m 坐标细化没有产生可信的单像素尖峰。原优化网格上多个点量化为相同覆盖率；换用错开的目标网格后，最优点会在 `(686,430)`、`(690,468)`、`(690,471)`、`(691,475)` 等邻近位置间移动。因此最终结论必须是“稳定区域 + 代表点”，而不是声称一个连续域精确最优坐标。

cross-grid 上以 `coverage >= 99.5%` sampled maximum 定义的稳定候选区为：

`C_0.5%: x=686..691 m, y=430..490 m`。

更严格的 `99.95%` sampled core 为：

`x=688..691 m, y=466..479 m`。

向上补扫至 `y=520 m` 后，最高值仍在内部的 `(690,471)` 附近，而不是新边界，排除了原局部窗口截断最优带的风险。

最终推荐代表点取：

`S2=(691,475) m`，移动距离 838.514 m。

ultra offset validation 给出：

- guaranteed receive area coverage：100%；
- robust clear-ready area coverage：25.5213%；
- uniform `(r,alpha)` coverage：35.7008%；
- worst-diameter p50/p90/p95：68.798/149.563/161.039 m；
- unbounded fraction：0。

另一套高密度对齐网格在同一点给出面积覆盖 25.8785%。因此最终稳妥表述为：推荐区面积覆盖约 25.5%–25.9%，而不是把任一网格值当作连续域精确最优。

### 3.2 epsilon 收敛

对最大面积、旧验证点、短距离点、中间 Pareto 点和均匀权重点分别采用 `epsilon step={0.1,0.05,0.02,0.01} deg`。所有候选的 robust clear-ready coverage 在四档 epsilon 网格上保持不变；直径均值只有远低于决策精度的变化。

因此当前坐标波动来自目标 `(r,alpha)` 离散边界，而不是遗漏 epsilon 内部失效点。0.01° epsilon 已足以支持本轮封板，不需要继续追逐更密误差小数位。

### 3.3 为什么只有约 25%

在推荐点 `(691,475)` 的高密度 target map 上：

- robust clear-ready 面积比例：25.8785%；
- diameter `>40 m`：72.6201%；
- diameter `<=40 m` 但最大直径圆覆盖失败：1.5014%；
- guaranteed receive failure、unbounded、empty、near：均为 0。

成功集不是零散噪声，而是近乎横跨全部 `alpha` 的两个连续径向带：

- 全 alpha 成功主带：`r=275..662.5 m`；
- 全 alpha 成功次带：`r=700..825 m`；
- `r=262.5 m` 仅部分 alpha 成功；
- `r=675..687.5 m` 是窄的 diameter-circle failure 分隔带；
- 更近和更远目标主要因 worst diameter 超过 40 m 失败。

这解释了面积覆盖率为何只有约四分之一：接收范围和有界性在推荐点并未造成损失，主要限制是随目标距离增长/变化的定位直径门槛；圆覆盖条件只切掉一条较窄的内部带。上述条带是高密度数值结构，尚不是解析边界定理。

## 4. D：最终 baseline、敏感性与 Pareto

下表全部来自同一高密度 target map（`r step=12.5 m`、`alpha step=0.025 deg`、`epsilon step=0.02 deg`）。所有点的 guaranteed receive 为 100%、unbounded fraction 为 0。

| 候选 | S2 (m) | 移动 (m) | 面积覆盖 | uniform 覆盖 | p50 (m) | p90 (m) | p95 (m) | 圆覆盖失败面积 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 推荐代表点 | (691,475) | 838.51 | 25.8785% | 36.1420% | 69.43 | 151.25 | 162.89 | 1.5014% |
| 旧优化点 | (690,468) | 833.74 | 25.8749% | 36.1317% | 69.88 | 153.00 | 164.92 | 1.5014% |
| 90%阈值最短点 | (620,320) | 697.71 | 23.5690% | 35.9671% | 96.80 | 234.50 | 256.38 | 1.3636% |
| 中间 Pareto 点 | (648,386) | 754.25 | 25.0964% | 37.4074% | 81.91 | 189.76 | 205.57 | 1.4187% |
| uniform 权重最佳 | (550,400) | 680.07 | 22.8808% | 47.1914% | 97.26 | 213.65 | 230.73 | 1.1983% |
| 正交 baseline | (806.98,399.94) | 900.65 | 22.6433% | 24.5782% | 60.99 | 142.19 | 154.58 | 2.6446% |

移动–质量关系明确：推荐点比 90% 阈值最短点多移动 140.80 m，换来约 2.31 个百分点的面积覆盖提升，同时 p95 从 256.38 m 降到 162.89 m；中间点提供连续折中。正交 baseline 的直径统计最好，但移动更长且覆盖率较低，说明单目标正交几何不能替代完整不确定集优化。

权重敏感性同样明确：若错误地对 `(r,alpha)` 等权，最优区移向 `(550,400)`，uniform 指标达到约 47.2%，但物理面积覆盖降至约 22.9%。这是因为等权会高估小 r 环带，而面积元必须乘 r。两种权重下高质量区域明显变化，所以论文主结果必须坚持面积权重，uniform 只作为敏感性分析。

最终四指标 Pareto 文件同时考虑面积覆盖、guaranteed receive、p95 和移动距离；不人为指定 lambda。前沿点较多是因为指标离散和四维非支配关系，不意味着每个点都值得作为策略候选。

## 5. 严格全域回归与 Q2 逻辑闭环

此前对 `g_near=(6,0)`、`g_far=(1500,0)` 和 2564 个候选的严格回归结果仍为：near feasible 353、far feasible 52、both feasible 0。它与已有解析冲突结论一致，但本身只是有限候选一致性检查。

Q2 现形成闭环：

1. 两次角域交会全部由 Q1 exact geometry 判定；
2. 解析正交直觉说明接近 90°有利；
3. 35 组 exact 数值实验确认有限角宽下 minimax 位于正交邻域，当前采样修正约 +1°且收益很小；
4. 严格全域固定 S2 不存在；
5. 因而改为最大化完整不确定位置集上的 robust clear-ready 面积覆盖；
6. 多级、错位和 ultra 网格得到稳定候选区 `x=686..691, y=430..490 m`；
7. 推荐点 `(691,475)`、中间点 `(648,386)` 和短距离点 `(620,320)` 给出清晰 Pareto 选择；
8. epsilon 收敛、失效结构和目标权重敏感性均已解释。

## 6. 输出与复现

主入口：

```text
python experiments/q2/run_q2_finalization.py --workers 8
```

脚本完全确定性。正交、候选细化、交叉验证、epsilon 收敛和 target map 均有 CSV；长阶段带签名元数据并支持默认 resume/cache。主要最终文件包括：

- `q2_orthogonal_correction.csv`
- `q2_candidate_refinement.csv`
- `q2_candidate_cross_grid.csv`
- `q2_candidate_cross_grid_upper_extension.csv`
- `q2_candidate_cross_grid_refined.csv`
- `q2_candidate_ultra_validation.csv`
- `q2_epsilon_convergence.csv`
- `q2_target_space_map.csv`
- `q2_target_space_summary.csv`
- `q2_weighting_sensitivity.csv`
- `q2_top_candidates.csv`
- `q2_pareto_candidates.csv`
- `q2_finalization_summary.json`

论文候选图包括正交修正、完整 S2 热图、最终候选区/权重敏感性、Pareto 和推荐点 target-space failure structure。

## 7. 最终判断

**Q2 可以封板。**

剩余不确定性主要是连续域边界与有限目标网格之间的亚百分点差异，不影响模型选择、推荐区域、失效机理或 Pareto 排序的核心叙事。继续为 0.05% 做全域搜索不具备足够收益。后续只应在论文写作时准确保留“有限网格数值结论”的限定，不再继续优化 Q2，也不在未获确认前进入 Q3。
