# CURRENT STATUS

## 当前阶段

Q1 数值鲁棒性修复完成，建议冻结。

## 当前实现

示向度观测仍按正东 0°、逆时针为正和 `theta ± 1°` 构造两个归一化半平面。定位区域继续使用边界线求交、可行性筛选和真实退化锥判定来区分 `EMPTY`、`UNBOUNDED`、`BOUNDED`，不使用人工边界框；有界后继续增量裁剪凸多边形。

## 数值鲁棒性修复

- 审查复现了 near-tangent bearing wedges 的状态误判：远处边界交点的大数相消可令理论为零的残差略超绝对 `EPS`，导致真实 `UNBOUNDED` 被判成 `EMPTY`。
- 快速浮点候选和公共 API 保持不变；仅在没有浮点可行候选且即将判空时，对已有二进制浮点半平面系数执行精确有理数可行性复核。
- 精确复核沿用原 `eps`，不扩大几何容差；近临界但方向区间确实分离的 `EMPTY` 回归案例仍为 `EMPTY`。

## 回归验证

- `python -m pytest -q`：17 passed，0 failed。
- near-tangent 严格正重叠 0.1°、0.01°、0.001°：全部为 `UNBOUNDED`。
- Monte Carlo，seed 20260910，500 例：passed 500，failed 0，unbounded/skipped 0。
- Monte Carlo，seed 20260911，5000 例：passed 5000，failed 0，unbounded/skipped 0。
- Monte Carlo 现显式检查：已知可行目标不得分类为 `EMPTY`；所有有界顶点满足全部原始观测；增量 polygon 包含关系、面积和直径均保持单调。
- Q1 命令行演示正常：状态 `BOUNDED`，直径和直径圆覆盖正常，图像生成成功。

## 冻结结论

YES
