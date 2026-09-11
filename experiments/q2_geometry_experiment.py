"""Deterministic geometry-only exploration for the second Q2 detector.

The experiment deliberately excludes reception range, travel cost, timing,
channels, and operational strategy.  It reuses the Q1 half-plane intersection,
diameter, and diameter-circle coverage implementations without changing their
mathematical conventions.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import statistics
import sys
import tempfile
import time
from typing import Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.geometry.coverage import diameter_circle_coverage
from src.geometry.diameter import polygon_diameter
from src.geometry.localization import (
    EPS,
    Measurement,
    Point,
    RegionStatus,
    construct_bounded_polygon,
)


FIRST_MEASUREMENT = Measurement(0.0, 0.0, 0.0)
TARGET_RADII = tuple(float(value) for value in range(50, 1501, 50))
TARGET_ANGLES_DEG = (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
SECOND_ERRORS_DEG = (-1.0, -0.5, 0.0, 0.5, 1.0)
DISTANCE_SLICES = (300.0, 600.0, 900.0, 1200.0, 1500.0)
DIAMETER_LIMIT = 40.0
COINCIDENCE_EPS = 1e-12
RATIO_FIELDS = (
    "coverage_ratio",
    "valid_polygon_ratio",
    "diameter_pass_ratio",
    "circle_coverage_pass_ratio",
)


@dataclass(frozen=True, slots=True)
class Outcome:
    """Q1-derived metrics for one target and one sampled second-angle error."""

    valid_polygon: bool
    diameter_pass: bool
    circle_coverage_pass: bool
    clear_ready: bool
    diameter: float | None


@dataclass(slots=True)
class Accumulator:
    """Counts accumulated over individual (target, error) outcomes."""

    evaluations: int = 0
    valid_polygons: int = 0
    diameter_passes: int = 0
    circle_coverage_passes: int = 0
    clear_ready_passes: int = 0
    diameters: list[float] | None = None

    def __post_init__(self) -> None:
        if self.diameters is None:
            self.diameters = []

    def add(self, outcome: Outcome) -> None:
        self.evaluations += 1
        self.valid_polygons += int(outcome.valid_polygon)
        self.diameter_passes += int(outcome.diameter_pass)
        self.circle_coverage_passes += int(outcome.circle_coverage_pass)
        self.clear_ready_passes += int(outcome.clear_ready)
        if outcome.diameter is not None:
            self.diameters.append(outcome.diameter)

    def add_undefined(self, count: int) -> None:
        self.evaluations += count


def _inclusive_grid(start: float, stop: float, step: float) -> tuple[float, ...]:
    count = int(round((stop - start) / step))
    return tuple(start + index * step for index in range(count + 1))


def _target_point(radius: float, angle_deg: float) -> Point:
    angle = math.radians(angle_deg)
    return Point(radius * math.cos(angle), radius * math.sin(angle))


def _bearing(origin: Point, target: Point) -> float:
    return math.degrees(math.atan2(target.y - origin.y, target.x - origin.x)) % 360.0


def _percentile(values: Sequence[float], probability: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _evaluate_outcome(sensor: Point, target: Point, error_deg: float) -> Outcome:
    true_bearing = _bearing(sensor, target)
    second_measurement = Measurement(
        sensor.x,
        sensor.y,
        (true_bearing + error_deg) % 360.0,
    )
    construction = construct_bounded_polygon(
        (FIRST_MEASUREMENT, second_measurement),
        delta_deg=1.0,
        eps=EPS,
    )
    if construction.status is not RegionStatus.BOUNDED or not construction.vertices:
        return Outcome(False, False, False, False, None)

    diameter = polygon_diameter(construction.vertices, EPS)
    coverage = diameter_circle_coverage(construction.vertices, diameter, EPS)
    diameter_pass = diameter.length <= DIAMETER_LIMIT + EPS
    circle_pass = coverage.exists_covering_diameter_circle
    return Outcome(
        True,
        diameter_pass,
        circle_pass,
        diameter_pass and circle_pass,
        diameter.length,
    )


def _metric_row(
    sensor: Point,
    accumulator: Accumulator,
    robust_passes: int,
    target_count: int,
    undefined_target_count: int,
) -> dict[str, float | int]:
    assert accumulator.diameters is not None
    denominator = accumulator.evaluations
    return {
        "x": sensor.x,
        "y": sensor.y,
        "distance_from_s1": math.hypot(sensor.x, sensor.y),
        "angle_from_first_bearing": math.degrees(math.atan2(sensor.y, sensor.x)),
        "coverage_ratio": robust_passes / target_count,
        "valid_polygon_ratio": accumulator.valid_polygons / denominator,
        "diameter_pass_ratio": accumulator.diameter_passes / denominator,
        "circle_coverage_pass_ratio": (
            accumulator.circle_coverage_passes / denominator
        ),
        "clear_ready_outcome_ratio": accumulator.clear_ready_passes / denominator,
        "median_diameter": (
            statistics.median(accumulator.diameters)
            if accumulator.diameters
            else math.nan
        ),
        "p90_diameter": _percentile(accumulator.diameters, 0.90),
        "max_valid_diameter": (
            max(accumulator.diameters) if accumulator.diameters else math.nan
        ),
        "robust_pass_count": robust_passes,
        "target_count": target_count,
        "outcome_count": denominator,
        "valid_polygon_count": accumulator.valid_polygons,
        "diameter_pass_count": accumulator.diameter_passes,
        "circle_coverage_pass_count": accumulator.circle_coverage_passes,
        "clear_ready_outcome_count": accumulator.clear_ready_passes,
        "circle_blocked_diameter_pass_count": (
            accumulator.diameter_passes - accumulator.clear_ready_passes
        ),
        "undefined_target_count": undefined_target_count,
    }


def evaluate_candidate(
    sensor: Point,
    radii: Sequence[float] = TARGET_RADII,
    target_angles_deg: Sequence[float] = TARGET_ANGLES_DEG,
    errors_deg: Sequence[float] = SECOND_ERRORS_DEG,
    distance_slices: Sequence[float] = DISTANCE_SLICES,
) -> tuple[dict[str, float | int], list[dict[str, float | int]]]:
    """Evaluate one S2 point over the structured target and error grids."""

    total = Accumulator()
    slice_accumulators = {radius: Accumulator() for radius in distance_slices}
    slice_robust_passes = {radius: 0 for radius in distance_slices}
    slice_undefined_targets = {radius: 0 for radius in distance_slices}
    robust_passes = 0
    undefined_target_count = 0

    for radius in radii:
        for angle_deg in target_angles_deg:
            target = _target_point(radius, angle_deg)
            slice_accumulator = slice_accumulators.get(radius)
            if math.hypot(target.x - sensor.x, target.y - sensor.y) <= COINCIDENCE_EPS:
                total.add_undefined(len(errors_deg))
                undefined_target_count += 1
                if slice_accumulator is not None:
                    slice_accumulator.add_undefined(len(errors_deg))
                    slice_undefined_targets[radius] += 1
                continue

            robust_clear_ready = True
            for error_deg in errors_deg:
                outcome = _evaluate_outcome(sensor, target, error_deg)
                total.add(outcome)
                if slice_accumulator is not None:
                    slice_accumulator.add(outcome)
                robust_clear_ready = robust_clear_ready and outcome.clear_ready
            robust_passes += int(robust_clear_ready)
            if slice_accumulator is not None:
                slice_robust_passes[radius] += int(robust_clear_ready)

    target_count = len(radii) * len(target_angles_deg)
    overall_row = _metric_row(
        sensor,
        total,
        robust_passes,
        target_count,
        undefined_target_count,
    )
    slice_rows: list[dict[str, float | int]] = []
    for radius in distance_slices:
        row = _metric_row(
            sensor,
            slice_accumulators[radius],
            slice_robust_passes[radius],
            len(target_angles_deg),
            slice_undefined_targets[radius],
        )
        row["target_radius"] = radius
        slice_rows.append(row)
    return overall_row, slice_rows


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _extent(rows: Sequence[dict[str, object]], threshold: float) -> dict[str, object]:
    selected = [row for row in rows if float(row["coverage_ratio"]) >= threshold]
    if not selected:
        return {"count": 0, "bounds": None}
    return {
        "count": len(selected),
        "bounds": {
            "x_min": min(float(row["x"]) for row in selected),
            "x_max": max(float(row["x"]) for row in selected),
            "y_min": min(float(row["y"]) for row in selected),
            "y_max": max(float(row["y"]) for row in selected),
        },
    }


def _symmetry_summary(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    lookup = {
        (float(row["x"]), float(row["y"])): row
        for row in rows
    }
    differences = {field: [] for field in RATIO_FIELDS}
    for (x, y), row in lookup.items():
        if y <= 0.0:
            continue
        mirror = lookup[(x, -y)]
        for field in RATIO_FIELDS:
            differences[field].append(abs(float(row[field]) - float(mirror[field])))
    return {
        field: {
            "max_absolute_error": max(values, default=0.0),
            "mean_absolute_error": statistics.fmean(values) if values else 0.0,
        }
        for field, values in differences.items()
    }


def _distance_slice_summary(
    slice_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    for radius in DISTANCE_SLICES:
        rows = [
            row for row in slice_rows if float(row["target_radius"]) == radius
        ]
        ordered = sorted(
            rows,
            key=lambda row: (
                -float(row["coverage_ratio"]),
                -float(row["valid_polygon_ratio"]),
                float(row["distance_from_s1"]),
            ),
        )
        best = ordered[0]
        high_90 = _extent(rows, 0.90)
        bounds = high_90["bounds"] or {}
        summaries.append(
            {
                "target_radius": radius,
                "max_coverage_ratio": best["coverage_ratio"],
                "best_x": best["x"],
                "best_y": best["y"],
                "best_distance_from_s1": best["distance_from_s1"],
                "candidate_count_ge_0_90": high_90["count"],
                "candidate_count_ge_0_95": _extent(rows, 0.95)["count"],
                "strict_candidate_count": sum(
                    float(row["coverage_ratio"]) == 1.0 for row in rows
                ),
                "ge_0_90_x_min": bounds.get("x_min", ""),
                "ge_0_90_x_max": bounds.get("x_max", ""),
                "ge_0_90_y_min": bounds.get("y_min", ""),
                "ge_0_90_y_max": bounds.get("y_max", ""),
            }
        )
    return summaries


def _matrix(
    rows: Sequence[dict[str, object]],
    x_values: Sequence[float],
    y_values: Sequence[float],
    field: str,
) -> list[list[float]]:
    lookup = {
        (float(row["x"]), float(row["y"])): float(row[field])
        for row in rows
    }
    return [[lookup[(x, y)] for x in x_values] for y in y_values]


def _add_first_wedge(axes, x_limit: float) -> None:
    boundary = x_limit * math.tan(math.radians(1.0))
    axes.plot((0.0, x_limit), (0.0, boundary), color="#d95f02", linewidth=1.0)
    axes.plot((0.0, x_limit), (0.0, -boundary), color="#d95f02", linewidth=1.0)
    axes.plot((0.0, x_limit), (0.0, 0.0), color="#444444", linewidth=0.8)
    axes.scatter((0.0,), (0.0,), marker="*", s=70, color="#d73027", zorder=5)


def _contours(axes, x_values, y_values, values, add_labels: bool = True) -> None:
    minimum = min(min(row) for row in values)
    maximum = max(max(row) for row in values)
    levels = [
        value
        for value in (0.25, 0.40, 0.50, 0.75, 0.90, 0.95)
        if minimum < value < maximum
    ]
    if levels:
        contours = axes.contour(
            x_values,
            y_values,
            values,
            levels=levels,
            colors="white",
            linewidths=0.7,
        )
        if add_labels:
            axes.clabel(contours, inline=True, fontsize=7, fmt="%.2f")


def create_figures(
    grid_rows: Sequence[dict[str, object]],
    slice_rows: Sequence[dict[str, object]],
    x_values: Sequence[float],
    y_values: Sequence[float],
    figure_dir: Path,
    grid_step: float,
) -> list[Path]:
    """Create the four requested static scientific figures."""

    matplotlib_cache = Path(tempfile.gettempdir()) / "mathb-matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Polygon

    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    # Figure 1: the first-measurement uncertainty wedge.
    figure, axes = plt.subplots(figsize=(10, 4.8))
    radius = 1500.0
    angles = [math.radians(-1.0 + index * 2.0 / 200) for index in range(201)]
    sector = [(0.0, 0.0)] + [
        (radius * math.cos(angle), radius * math.sin(angle)) for angle in angles
    ]
    axes.add_patch(Polygon(sector, closed=True, facecolor="#9ecae1", alpha=0.55))
    _add_first_wedge(axes, radius)
    width = 2.0 * radius * math.tan(math.radians(1.0))
    half_width = width / 2.0
    axes.annotate(
        "",
        xy=(radius, half_width),
        xytext=(radius, -half_width),
        arrowprops={"arrowstyle": "<->", "color": "#2166ac", "linewidth": 1.5},
    )
    axes.text(radius - 20.0, -43.0, f"w(1500) = {width:.2f} m", ha="right", va="center")
    clear_x = 1125.0
    axes.annotate(
        "",
        xy=(clear_x, 20.0),
        xytext=(clear_x, -20.0),
        arrowprops={"arrowstyle": "<->", "color": "#7b3294", "linewidth": 2.0},
    )
    axes.text(clear_x - 15.0, 34.0, "clear diameter = 40 m", ha="right", va="center")
    axes.text(15.0, 3.0, "S1", ha="left", va="bottom")
    axes.set_xlim(-50.0, 1580.0)
    axes.set_ylim(-70.0, 70.0)
    axes.set_aspect("equal", adjustable="box")
    axes.set_xlabel("x (m)")
    axes.set_ylabel("y (m)")
    axes.set_title("Figure 1. First-bearing uncertainty set F1")
    axes.grid(alpha=0.2)
    figure.tight_layout()
    path = figure_dir / "figure1_first_bearing_uncertainty.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(path)

    # Figure 2: robust coverage heatmap.
    coverage = _matrix(grid_rows, x_values, y_values, "coverage_ratio")
    extent = (
        x_values[0] - grid_step / 2.0,
        x_values[-1] + grid_step / 2.0,
        y_values[0] - grid_step / 2.0,
        y_values[-1] + grid_step / 2.0,
    )
    figure, axes = plt.subplots(figsize=(9, 9))
    image = axes.imshow(
        coverage,
        origin="lower",
        extent=extent,
        vmin=0.0,
        vmax=1.0,
        cmap="viridis",
        interpolation="nearest",
        aspect="equal",
    )
    _contours(axes, x_values, y_values, coverage)
    _add_first_wedge(axes, max(x_values))
    axes.set_xlabel("S2 x (m)")
    axes.set_ylabel("S2 y (m)")
    axes.set_title("Figure 2. Robust uncertainty-set coverage ratio")
    figure.colorbar(image, ax=axes, label="coverage ratio")
    figure.tight_layout()
    path = figure_dir / "figure2_s2_coverage_heatmap.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(path)

    # Figure 3: exploratory threshold regions.
    categories = [
        [
            3 if value == 1.0 else 2 if value >= 0.95 else 1 if value >= 0.90 else 0
            for value in row
        ]
        for row in coverage
    ]
    colors = ListedColormap(("#f0f0f0", "#9ecae1", "#3182bd", "#54278f"))
    norm = BoundaryNorm((-0.5, 0.5, 1.5, 2.5, 3.5), colors.N)
    figure, axes = plt.subplots(figsize=(9, 9))
    image = axes.imshow(
        categories,
        origin="lower",
        extent=extent,
        cmap=colors,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
    )
    _add_first_wedge(axes, max(x_values))
    axes.set_xlabel("S2 x (m)")
    axes.set_ylabel("S2 y (m)")
    axes.set_title("Figure 3. Geometry-only candidate thresholds (exploratory)")
    colorbar = figure.colorbar(image, ax=axes, ticks=(0, 1, 2, 3))
    colorbar.ax.set_yticklabels(("< 0.90", ">= 0.90", ">= 0.95", "= 1.00"))
    if not any(float(row["coverage_ratio"]) == 1.0 for row in grid_rows):
        axes.text(
            0.02,
            0.98,
            "strict candidate region is empty under the current discretization\n"
            f"maximum coverage ratio = {max(max(row) for row in coverage):.3f}",
            transform=axes.transAxes,
            ha="left",
            va="top",
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
    figure.tight_layout()
    path = figure_dir / "figure3_candidate_regions.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(path)

    # Figure 4: robust coverage by target-distance slice.
    figure = plt.figure(figsize=(16, 10), constrained_layout=True)
    layout = figure.add_gridspec(2, 4, width_ratios=(1.0, 1.0, 1.0, 0.06))
    map_axes = [
        figure.add_subplot(layout[0, 0]),
        figure.add_subplot(layout[0, 1]),
        figure.add_subplot(layout[0, 2]),
        figure.add_subplot(layout[1, 0]),
        figure.add_subplot(layout[1, 1]),
    ]
    summary_axes = figure.add_subplot(layout[1, 2])
    colorbar_axes = figure.add_subplot(layout[:, 3])
    last_image = None
    maxima: list[float] = []
    for axes, radius in zip(map_axes, DISTANCE_SLICES):
        rows = [
            row for row in slice_rows if float(row["target_radius"]) == radius
        ]
        values = _matrix(rows, x_values, y_values, "coverage_ratio")
        maxima.append(max(max(row) for row in values))
        last_image = axes.imshow(
            values,
            origin="lower",
            extent=extent,
            vmin=0.0,
            vmax=1.0,
            cmap="viridis",
            interpolation="nearest",
            aspect="equal",
        )
        _contours(axes, x_values, y_values, values, add_labels=False)
        _add_first_wedge(axes, max(x_values))
        axes.set_title(f"Target distance r = {radius:.0f} m")
        axes.set_xlabel("S2 x (m)")
        axes.set_ylabel("S2 y (m)")
    summary_axes.plot(DISTANCE_SLICES, maxima, marker="o", color="#2166ac")
    summary_axes.set_ylim(0.0, 1.05)
    summary_axes.set_xlabel("Target distance r (m)")
    summary_axes.set_ylabel("Maximum slice coverage ratio")
    summary_axes.set_title("Best sampled coverage by distance")
    summary_axes.grid(alpha=0.25)
    assert last_image is not None
    figure.colorbar(last_image, cax=colorbar_axes, label="slice coverage ratio")
    path = figure_dir / "figure4_distance_slices.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(path)
    return paths


def _scan_candidates(
    x_values: Sequence[float],
    y_values: Sequence[float],
    label: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grid_rows: list[dict[str, object]] = []
    slice_rows: list[dict[str, object]] = []
    candidate_count = len(x_values) * len(y_values)
    completed = 0
    for y in y_values:
        for x in x_values:
            row, candidate_slice_rows = evaluate_candidate(Point(x, y))
            grid_rows.append(row)
            slice_rows.extend(candidate_slice_rows)
            completed += 1
        if completed % (10 * len(x_values)) == 0 or completed == candidate_count:
            print(f"{label} candidates: {completed}/{candidate_count}", flush=True)
    return grid_rows, slice_rows


def run_experiment(
    grid_step: float,
    output_dir: Path,
    figure_dir: Path,
    top_count: int,
) -> dict[str, object]:
    """Run the deterministic coarse grid and write all data products."""

    started = time.perf_counter()
    x_values = _inclusive_grid(-500.0, 1800.0, grid_step)
    y_values = _inclusive_grid(-1800.0, 1800.0, grid_step)
    grid_rows, slice_rows = _scan_candidates(x_values, y_values, "coarse")
    candidate_count = len(x_values) * len(y_values)

    refine_step = 50.0
    refine_x_values = _inclusive_grid(350.0, 750.0, refine_step)
    refine_y_values = (
        *_inclusive_grid(-850.0, -250.0, refine_step),
        *_inclusive_grid(250.0, 850.0, refine_step),
    )
    refined_rows, _ = _scan_candidates(
        refine_x_values,
        refine_y_values,
        "refined",
    )
    refined_candidate_count = len(refine_x_values) * len(refine_y_values)

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "q2_s2_grid_metrics.csv", grid_rows)
    _write_csv(output_dir / "q2_s2_refined_metrics.csv", refined_rows)
    ordered = sorted(
        refined_rows,
        key=lambda row: (
            -float(row["coverage_ratio"]),
            -float(row["valid_polygon_ratio"]),
            float(row["distance_from_s1"]),
        ),
    )
    top_rows = [dict(rank=index + 1, **row) for index, row in enumerate(ordered[:top_count])]
    _write_csv(output_dir / "q2_top_candidates.csv", top_rows)
    _write_csv(output_dir / "q2_distance_slice_grid.csv", slice_rows)
    distance_summaries = _distance_slice_summary(slice_rows)
    _write_csv(output_dir / "q2_distance_slice_summary.csv", distance_summaries)

    figures = create_figures(
        grid_rows,
        slice_rows,
        x_values,
        y_values,
        figure_dir,
        grid_step,
    )
    all_candidate_rows = [*grid_rows, *refined_rows]
    strict = [
        row for row in all_candidate_rows if float(row["coverage_ratio"]) == 1.0
    ]
    centerline = [row for row in grid_rows if float(row["y"]) == 0.0]
    far = [row for row in grid_rows if float(row["distance_from_s1"]) >= 2000.0]
    best = ordered[0]
    total_outcomes = sum(int(row["outcome_count"]) for row in grid_rows)
    total_valid = sum(int(row["valid_polygon_count"]) for row in grid_rows)
    total_diameter_passes = sum(int(row["diameter_pass_count"]) for row in grid_rows)
    total_circle_passes = sum(
        int(row["circle_coverage_pass_count"]) for row in grid_rows
    )
    total_clear_ready = sum(
        int(row["clear_ready_outcome_count"]) for row in grid_rows
    )
    undefined_bearing_scenarios = sum(
        int(row["undefined_target_count"]) * len(SECOND_ERRORS_DEG)
        for row in grid_rows
    )
    coarse_scenario_count = (
        candidate_count
        * len(TARGET_RADII)
        * len(TARGET_ANGLES_DEG)
        * len(SECOND_ERRORS_DEG)
    )
    refined_scenario_count = (
        refined_candidate_count
        * len(TARGET_RADII)
        * len(TARGET_ANGLES_DEG)
        * len(SECOND_ERRORS_DEG)
    )
    summary: dict[str, object] = {
        "configuration": {
            "grid_step": grid_step,
            "x_range": [x_values[0], x_values[-1]],
            "y_range": [y_values[0], y_values[-1]],
            "candidate_count": candidate_count,
            "refined_x_range": [refine_x_values[0], refine_x_values[-1]],
            "refined_abs_y_range": [250.0, 850.0],
            "refined_grid_step": refine_step,
            "refined_candidate_count": refined_candidate_count,
            "target_radii": list(TARGET_RADII),
            "target_angles_deg": list(TARGET_ANGLES_DEG),
            "target_count": len(TARGET_RADII) * len(TARGET_ANGLES_DEG),
            "second_errors_deg": list(SECOND_ERRORS_DEG),
            "coarse_scenario_count": coarse_scenario_count,
            "refined_scenario_count": refined_scenario_count,
            "total_scenario_count": coarse_scenario_count + refined_scenario_count,
            "undefined_bearing_scenario_count": undefined_bearing_scenarios,
            "q1_localization_evaluation_count": (
                coarse_scenario_count
                + refined_scenario_count
                - undefined_bearing_scenarios
            ),
        },
        "coarse_best_candidate": max(
            grid_rows,
            key=lambda row: float(row["coverage_ratio"]),
        ),
        "best_candidate": best,
        "strict_candidate_count": len(strict),
        "strict_candidate_extent": _extent(all_candidate_rows, 1.0),
        "candidate_extent_ge_0_90": _extent(all_candidate_rows, 0.90),
        "candidate_extent_ge_0_95": _extent(all_candidate_rows, 0.95),
        "candidate_extent_ge_0_40": _extent(all_candidate_rows, 0.40),
        "candidate_extent_ge_0_425": _extent(all_candidate_rows, 0.425),
        "coarse_outcome_totals": {
            "outcome_count": total_outcomes,
            "valid_polygon_ratio": total_valid / total_outcomes,
            "diameter_pass_ratio": total_diameter_passes / total_outcomes,
            "circle_coverage_pass_ratio": total_circle_passes / total_outcomes,
            "clear_ready_outcome_ratio": total_clear_ready / total_outcomes,
            "diameter_passes_blocked_by_circle": (
                total_diameter_passes - total_clear_ready
            ),
            "circle_block_fraction_among_diameter_passes": (
                (total_diameter_passes - total_clear_ready) / total_diameter_passes
                if total_diameter_passes
                else math.nan
            ),
        },
        "symmetry": _symmetry_summary(grid_rows),
        "centerline": {
            "max_coverage_ratio": max(float(row["coverage_ratio"]) for row in centerline),
            "mean_coverage_ratio": statistics.fmean(
                float(row["coverage_ratio"]) for row in centerline
            ),
        },
        "far_distance_ge_2000": {
            "candidate_count": len(far),
            "max_coverage_ratio": max(
                (float(row["coverage_ratio"]) for row in far),
                default=math.nan,
            ),
            "mean_coverage_ratio": (
                statistics.fmean(float(row["coverage_ratio"]) for row in far)
                if far
                else math.nan
            ),
        },
        "distance_slices": distance_summaries,
        "outputs": {
            "figures": [str(path.resolve()) for path in figures],
            "csv": [
                str((output_dir / name).resolve())
                for name in (
                    "q2_s2_grid_metrics.csv",
                    "q2_s2_refined_metrics.csv",
                    "q2_top_candidates.csv",
                    "q2_distance_slice_summary.csv",
                    "q2_distance_slice_grid.csv",
                )
            ],
        },
        "elapsed_seconds": time.perf_counter() - started,
    }
    summary_path = output_dir / "q2_experiment_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, allow_nan=True)
        handle.write("\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid-step", type=float, default=100.0)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/q2_geometry_experiment"),
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("figures/q2_geometry_experiment"),
    )
    parser.add_argument("--top-count", type=int, default=100)
    arguments = parser.parse_args()
    if arguments.grid_step <= 0.0:
        parser.error("--grid-step must be positive")
    run_experiment(
        arguments.grid_step,
        arguments.output_dir,
        arguments.figure_dir,
        arguments.top_count,
    )


if __name__ == "__main__":
    main()
