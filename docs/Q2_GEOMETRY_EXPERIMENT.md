# Q2 Geometry-Only Candidate-Region Experiment

## Scope

This is the first geometry-only exploration of the second detector position. It does not include reception range, travel cost, timing, channels, or an operational strategy. The first observation is normalized to `S1=(0,0)`, `theta1=0°`, with the unchanged Q1 error model `±1°`.

## Method

- Target uncertainty set: `r=50,100,...,1500 m` and `alpha ∈ {-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1}°`, for 270 deterministic target samples across the full first wedge.
- Second-measurement errors: `e2 ∈ {-1,-0.5,0,0.5,1}°`.
- Coarse S2 grid: `x ∈ [-500,1800] m`, `y ∈ [-1800,1800] m`, step 100 m, for 888 candidates.
- Local refinement: `x ∈ [350,750] m`, `|y| ∈ [250,850] m`, step 50 m, for 234 candidates.
- Total sampled scenarios: 1,514,700. There were 75 grid-coincidence scenarios where `S2=G` and the bearing is undefined; these were explicitly marked invalid rather than assigning a synthetic bearing. Q1 localization was evaluated 1,514,625 times.
- For every valid `(S2,G,e2)`, the experiment directly calls Q1's `construct_bounded_polygon`, `polygon_diameter`, and `diameter_circle_coverage`. A target is robust clear-ready only if all five sampled errors produce a bounded nonempty polygon, diameter at most 40 m, and at least one covering diameter circle.
- `valid_polygon_ratio`, `diameter_pass_ratio`, and `circle_coverage_pass_ratio` use all individual `(G,e2)` outcomes as their denominator. `coverage_ratio` is the fraction of target samples for which all five errors are clear-ready; it is a discrete uncertainty-set coverage ratio, not a probability.

Run with:

```powershell
python experiments/q2_geometry_experiment.py --grid-step 100
```

## Findings

- No overall strict candidate exists: the maximum robust uncertainty-set coverage ratio is `117/270 = 0.433333`.
- The coarse best is `(600,-500) m`; refinement finds a tied nearest best at `(550,-400) m`, with a mirror at `(550,400) m`. Several points near `x=550..600 m`, `|y|=400..650 m` attain the same maximum.
- As a descriptive, non-prescriptive band, sampled points with ratio at least 0.425 lie at `x=500..600 m`, `|y|=400..650 m`. No point reaches 0.90 or 0.95.
- All four ratio metrics are exactly symmetric on the sampled `y` grid: maximum and mean mirror errors are both zero.
- Every sampled centerline point has zero robust coverage. Every coarse candidate at distance at least 2000 m also has zero robust coverage, so moving arbitrarily far does not improve this criterion.
- Across the coarse grid, 98.7602% of individual outcomes form bounded polygons, but only 7.1431% pass `D<=40`; diameter is the main filter. Circle coverage rejects a further 5,392 diameter-passing outcomes, or 6.2967% of all diameter passes, so it is a real but secondary restriction.

## Distance slices

| Target r | Maximum ratio | Strict candidates | Approximate `ratio>=0.90` bounds |
|---:|---:|---:|---|
| 300 m | 1.0 | 102 | `x=0..700`, `y=-1000..1000` m |
| 600 m | 1.0 | 52 | `x=400..800`, `y=-800..800` m |
| 900 m | 1.0 | 10 | `x=800..1000`, `y=-500..500` m |
| 1200 m | 0.0 | 0 | empty |
| 1500 m | 0.0 | 0 | empty |

The useful S2 band moves outward and contracts as target distance increases. The abrupt loss at 1200 m is consistent with the first-wedge transverse width crossing the 40 m clear diameter: `2r tan(1°)` is about 41.89 m at 1200 m and 52.37 m at 1500 m. This is an empirical explanation under the sampled second-error set, not a continuous worst-case proof.

A non-monotone detail also appears for individual S2 points: the refined best `(550,-400)` passes every sampled `alpha` at several shorter radii, fails at 550 m, passes again at 600–700 m, and then fails from 750 m onward. Therefore target distance changes not only the scale but also the topology of the favorable S2 region.

## Outputs

Data:

- `results/q2_geometry_experiment/q2_s2_grid_metrics.csv`
- `results/q2_geometry_experiment/q2_s2_refined_metrics.csv`
- `results/q2_geometry_experiment/q2_top_candidates.csv`
- `results/q2_geometry_experiment/q2_distance_slice_summary.csv`
- `results/q2_geometry_experiment/q2_distance_slice_grid.csv`
- `results/q2_geometry_experiment/q2_experiment_summary.json`

Figures:

- `figures/q2_geometry_experiment/figure1_first_bearing_uncertainty.png`
- `figures/q2_geometry_experiment/figure2_s2_coverage_heatmap.png`
- `figures/q2_geometry_experiment/figure3_candidate_regions.png`
- `figures/q2_geometry_experiment/figure4_distance_slices.png`

## Interpretation and next step

The experiment supports the hypothesis that useful second-detector positions form off-axis, mirror-symmetric regions and that unknown target distance materially changes their location. It rejects the stronger hypothesis that one fixed geometry-only S2 can robustly clear the entire sampled first-wedge uncertainty set.

The results are meaningful and reproducible enough to be worth a dedicated Git commit after review. A new conversation is recommended before adding any new modeling layer, so this experiment remains a clean checkpoint. No next modeling round should begin without explicit confirmation.

本轮实验只验证 Q2 的 geometry-only hypothesis，未形成最终 Q2 策略。
