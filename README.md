# 2026 高教社杯 B 题实现

当前仓库只实现 Question 1：将多个带 ±1° 确定性误差的示向度表示为前向角域，计算其交集状态与有界凸多边形，并分析区域直径及对应直径圆的覆盖性。

## 运行

```powershell
python -m pytest
python -m src.q1.solve_q1
python experiments/q1/run_q1_validation.py
```

核心算法位于 `src/geometry/`，Q1 的可复现实例和绘图入口位于 `src/q1/`。算法不使用斜率、不使用人工大边界框，也不依赖 GIS 或计算几何库。

