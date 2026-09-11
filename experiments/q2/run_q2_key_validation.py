"""Q2 key geometry validation using the exact Q1 implementation.

This deterministic experiment studies only two-bearing geometry, guaranteed
reception, and movement distance.  It deliberately excludes scheduling,
multi-source routing, and any final operational strategy.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import sys
import tempfile
import time
from typing import Iterable, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.geometry.coverage import diameter_circle_coverage
from src.geometry.diameter import polygon_diameter
from src.geometry.localization import (
    EPS,
    Measurement,
    Point,
    RegionStatus,
    construct_bounded_polygon,
)


CACHE_VERSION = 1
FIRST_MEASUREMENT = Measurement(0.0, 0.0, 0.0)


@dataclass(frozen=True, slots=True)
class ExperimentConfig:
    delta_deg: float = 1.0
    diameter_limit: float = 40.0
    near_radius: float = 5.0
    min_effective_radius: float = 1000.0
    max_effective_radius: float = 1500.0
    target_r_start: float = 25.0
    target_r_stop: float = 1500.0
    target_r_step: float = 25.0
    target_alpha_start: float = -1.0
    target_alpha_stop: float = 1.0
    target_alpha_step: float = 0.1
    epsilon_start: float = -1.0
    epsilon_stop: float = 1.0
    epsilon_step: float = 0.25
    coarse_x_start: float = -1500.0
    coarse_x_stop: float = 3000.0
    coarse_y_start: float = 0.0
    coarse_y_stop: float = 3000.0
    coarse_step: float = 100.0
    orthogonal_target_x: float = 800.0
    orthogonal_d2: tuple[float, ...] = (200.0, 400.0, 600.0, 800.0, 1000.0)


CONFIG = ExperimentConfig()


@dataclass(frozen=True, slots=True)
class ReceiverState:
    distance_to_target: float
    guaranteed_rx: bool
    possible_rx: bool
    near: bool
    direction_available: bool


@dataclass(frozen=True, slots=True)
class GeometryOutcome:
    bounded: bool
    unbounded: bool
    empty: bool
    diameter: float | None
    diameter_pass: bool
    circle_pass: bool
    clear_ready: bool


def inclusive_grid(start: float, stop: float, step: float) -> tuple[float, ...]:
    count = int(round((stop - start) / step))
    return tuple(round(start + index * step, 10) for index in range(count + 1))


def target_point(radius: float, alpha_deg: float) -> Point:
    alpha = math.radians(alpha_deg)
    return Point(radius * math.cos(alpha), radius * math.sin(alpha))


def bearing(origin: Point, target: Point) -> float:
    return math.degrees(math.atan2(target.y - origin.y, target.x - origin.x)) % 360.0


def receiver_state(sensor: Point, target: Point, radius: float) -> ReceiverState:
    distance = math.hypot(sensor.x - target.x, sensor.y - target.y)
    near = distance <= CONFIG.near_radius
    guaranteed = distance <= max(CONFIG.min_effective_radius, radius)
    possible = distance <= CONFIG.max_effective_radius
    return ReceiverState(
        distance,
        guaranteed,
        possible,
        near,
        guaranteed and not near,
    )


def evaluate_q1(sensor: Point, target: Point, epsilon_deg: float) -> GeometryOutcome:
    theta2 = (bearing(sensor, target) + epsilon_deg) % 360.0
    construction = construct_bounded_polygon(
        (FIRST_MEASUREMENT, Measurement(sensor.x, sensor.y, theta2)),
        delta_deg=CONFIG.delta_deg,
        eps=EPS,
    )
    if construction.status is not RegionStatus.BOUNDED or not construction.vertices:
        return GeometryOutcome(
            False,
            construction.status is RegionStatus.UNBOUNDED,
            construction.status is RegionStatus.EMPTY,
            None,
            False,
            False,
            False,
        )
    diameter = polygon_diameter(construction.vertices, EPS)
    coverage = diameter_circle_coverage(construction.vertices, diameter, EPS)
    diameter_pass = diameter.length <= CONFIG.diameter_limit + EPS
    circle_pass = coverage.exists_covering_diameter_circle
    return GeometryOutcome(
        True,
        False,
        False,
        diameter.length,
        diameter_pass,
        circle_pass,
        diameter_pass and circle_pass,
    )


def intersection_angle_deg(sensor: Point, target: Point, alpha_deg: float) -> float:
    difference = abs((bearing(sensor, target) - alpha_deg + 180.0) % 360.0 - 180.0)
    return min(difference, 180.0 - difference)


def weighted_mean(values: Sequence[tuple[float, float]]) -> float:
    if not values:
        return math.nan
    total_weight = sum(weight for _, weight in values)
    return sum(value * weight for value, weight in values) / total_weight


def weighted_percentile(
    values: Sequence[tuple[float, float]], probability: float
) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    threshold = probability * sum(weight for _, weight in ordered)
    cumulative = 0.0
    for value, weight in ordered:
        cumulative += weight
        if cumulative >= threshold:
            return value
    return ordered[-1][0]


def evaluate_orthogonal_case(
    d2: float,
    beta_deg: float,
    epsilons: Sequence[float],
    stage: str,
) -> dict[str, object]:
    target = Point(CONFIG.orthogonal_target_x, 0.0)
    beta = math.radians(beta_deg)
    sensor = Point(
        target.x - d2 * math.cos(beta),
        target.y - d2 * math.sin(beta),
    )
    outcomes = [evaluate_q1(sensor, target, epsilon) for epsilon in epsilons]
    finite = [
        (outcome.diameter, epsilon)
        for outcome, epsilon in zip(outcomes, epsilons)
        if outcome.diameter is not None
    ]
    unbounded_epsilons = [
        epsilon
        for outcome, epsilon in zip(outcomes, epsilons)
        if outcome.unbounded or outcome.empty
    ]
    if unbounded_epsilons:
        worst_diameter = math.inf
        worst_epsilon = unbounded_epsilons[0]
    else:
        worst_diameter, worst_epsilon = max(finite)
    return {
        "stage": stage,
        "d2": d2,
        "beta_deg": beta_deg,
        "sensor_x": sensor.x,
        "sensor_y": sensor.y,
        "epsilon_step": epsilons[1] - epsilons[0],
        "epsilon_count": len(epsilons),
        "worst_case_diameter": worst_diameter,
        "mean_diameter": (
            statistics.fmean(value for value, _ in finite) if finite else math.inf
        ),
        "worst_epsilon": worst_epsilon,
        "bounded_count": sum(outcome.bounded for outcome in outcomes),
        "unbounded_or_empty_count": len(outcomes) - sum(
            outcome.bounded for outcome in outcomes
        ),
        "diameter_pass_all": all(outcome.diameter_pass for outcome in outcomes),
        "circle_pass_all": all(outcome.circle_pass for outcome in outcomes),
        "clear_ready_all": all(outcome.clear_ready for outcome in outcomes),
    }


def run_orthogonal_experiment(output_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    coarse_betas = inclusive_grid(3.5, 176.5, 0.5)
    coarse_epsilons = inclusive_grid(-1.0, 1.0, 0.1)
    refined_epsilons = inclusive_grid(-1.0, 1.0, 0.02)
    rows: list[dict[str, object]] = []
    optima: list[dict[str, object]] = []
    for d2 in CONFIG.orthogonal_d2:
        coarse = [
            evaluate_orthogonal_case(d2, beta, coarse_epsilons, "coarse")
            for beta in coarse_betas
        ]
        rows.extend(coarse)
        finite_coarse = [
            row for row in coarse if math.isfinite(float(row["worst_case_diameter"]))
        ]
        coarse_best = min(
            finite_coarse,
            key=lambda row: (float(row["worst_case_diameter"]), abs(float(row["beta_deg"]) - 90.0)),
        )
        center = float(coarse_best["beta_deg"])
        refined_betas = inclusive_grid(max(3.5, center - 1.0), min(176.5, center + 1.0), 0.05)
        refined = [
            evaluate_orthogonal_case(d2, beta, refined_epsilons, "refined")
            for beta in refined_betas
        ]
        rows.extend(refined)
        best = min(
            refined,
            key=lambda row: (float(row["worst_case_diameter"]), abs(float(row["beta_deg"]) - 90.0)),
        )
        optimum = dict(best)
        optimum["deviation_from_90_deg"] = float(best["beta_deg"]) - 90.0
        optima.append(optimum)
    write_csv(output_dir / "q2_orthogonal_validation.csv", rows)
    write_csv(output_dir / "q2_orthogonal_optima.csv", optima)
    return rows, optima


def _evaluate_sensor_on_grid(
    sensor: Point,
    radii: Sequence[float],
    alphas: Sequence[float],
    epsilons: Sequence[float],
    resolution_m: float,
    stage: str,
) -> dict[str, object]:
    total_weight = 0.0
    guaranteed_weight = 0.0
    possible_weight = 0.0
    near_weight = 0.0
    direction_weight = 0.0
    clear_weight = 0.0
    operational_weight = 0.0
    unbounded_weight = 0.0
    orthogonal_band_weight = 0.0
    guaranteed_count = 0
    clear_count = 0
    operational_count = 0
    target_count = 0
    worst_diameters: list[tuple[float, float]] = []
    orthogonality_errors: list[tuple[float, float]] = []

    for radius in radii:
        weight = radius
        for alpha_deg in alphas:
            target_count += 1
            total_weight += weight
            target = target_point(radius, alpha_deg)
            state = receiver_state(sensor, target, radius)
            possible_weight += weight * int(state.possible_rx)
            if not state.guaranteed_rx:
                continue
            guaranteed_count += 1
            guaranteed_weight += weight
            if state.near:
                near_weight += weight
                operational_weight += weight
                operational_count += 1
                continue

            direction_weight += weight
            crossing_angle = intersection_angle_deg(sensor, target, alpha_deg)
            orthogonality_errors.append((abs(crossing_angle - 90.0), weight))
            if crossing_angle >= 80.0:
                orthogonal_band_weight += weight

            outcomes = [evaluate_q1(sensor, target, epsilon) for epsilon in epsilons]
            all_clear = all(outcome.clear_ready for outcome in outcomes)
            clear_count += int(all_clear)
            operational_count += int(all_clear)
            if all_clear:
                clear_weight += weight
                operational_weight += weight
            if any(not outcome.bounded for outcome in outcomes):
                unbounded_weight += weight
            else:
                worst_diameters.append(
                    (max(float(outcome.diameter) for outcome in outcomes), weight)
                )

    diameter_mean = weighted_mean(worst_diameters)
    return {
        "stage": stage,
        "resolution_m": resolution_m,
        "x": sensor.x,
        "y": sensor.y,
        "movement_distance": math.hypot(sensor.x, sensor.y),
        "guaranteed_receive_coverage": guaranteed_weight / total_weight,
        "possible_receive_coverage": possible_weight / total_weight,
        "near_area_coverage": near_weight / total_weight,
        "direction_available_coverage": direction_weight / total_weight,
        "clear_ready_coverage": clear_weight / total_weight,
        "operational_good_coverage": operational_weight / total_weight,
        "uniform_guaranteed_receive_coverage": guaranteed_count / target_count,
        "uniform_clear_ready_coverage": clear_count / target_count,
        "uniform_operational_good_coverage": operational_count / target_count,
        "worst_diameter_median": weighted_percentile(worst_diameters, 0.50),
        "worst_diameter_mean": diameter_mean,
        "worst_diameter_p90": weighted_percentile(worst_diameters, 0.90),
        "worst_diameter_p95": weighted_percentile(worst_diameters, 0.95),
        "worst_diameter_max": (
            max((value for value, _ in worst_diameters), default=math.nan)
        ),
        "unbounded_fraction_of_direction_area": (
            unbounded_weight / direction_weight if direction_weight else math.nan
        ),
        "mean_abs_orthogonality_error_deg": weighted_mean(orthogonality_errors),
        "orthogonal_80_100_area_coverage": orthogonal_band_weight / total_weight,
        "target_sample_count": target_count,
        "epsilon_sample_count": len(epsilons),
    }


def optimization_grids() -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    return (
        inclusive_grid(CONFIG.target_r_start, CONFIG.target_r_stop, CONFIG.target_r_step),
        inclusive_grid(
            CONFIG.target_alpha_start,
            CONFIG.target_alpha_stop,
            CONFIG.target_alpha_step,
        ),
        inclusive_grid(CONFIG.epsilon_start, CONFIG.epsilon_stop, CONFIG.epsilon_step),
    )


def _worker(payload: tuple[float, float, float, str]) -> dict[str, object]:
    x, y, resolution_m, stage = payload
    radii, alphas, epsilons = optimization_grids()
    return _evaluate_sensor_on_grid(
        Point(x, y), radii, alphas, epsilons, resolution_m, stage
    )


def config_signature(extra: object) -> str:
    payload = {
        "cache_version": CACHE_VERSION,
        "config": asdict(CONFIG),
        "extra": extra,
    }
    encoded = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def scan_candidates(
    candidates: Sequence[Point],
    resolution_m: float,
    stage: str,
    output_path: Path,
    workers: int,
    resume: bool,
) -> list[dict[str, object]]:
    signature = config_signature(
        {
            "stage": stage,
            "resolution_m": resolution_m,
            "candidates": [(point.x, point.y) for point in candidates],
        }
    )
    metadata_path = output_path.with_suffix(".meta.json")
    if resume and output_path.exists() and metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("signature") == signature:
            print(f"cache hit: {output_path}", flush=True)
            return [dict(row) for row in read_csv(output_path)]

    payloads = [(point.x, point.y, resolution_m, stage) for point in candidates]
    started = time.perf_counter()
    if workers == 1:
        iterator: Iterable[dict[str, object]] = map(_worker, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        iterator = executor.map(_worker, payloads, chunksize=1)
    rows: list[dict[str, object]] = []
    try:
        for index, row in enumerate(iterator, 1):
            rows.append(row)
            if index % 100 == 0 or index == len(candidates):
                print(f"{stage}: {index}/{len(candidates)}", flush=True)
    finally:
        if executor is not None:
            executor.shutdown()
    write_csv(output_path, rows)
    metadata_path.write_text(
        json.dumps(
            {
                "signature": signature,
                "candidate_count": len(candidates),
                "elapsed_seconds": time.perf_counter() - started,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return rows


def _number(row: dict[str, object], field: str) -> float:
    try:
        return float(row[field])
    except (TypeError, ValueError):
        return math.nan


def pareto_front(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    objectives = [
        (
            -_number(row, "clear_ready_coverage"),
            -_number(row, "guaranteed_receive_coverage"),
            _number(row, "worst_diameter_p95"),
            _number(row, "movement_distance"),
        )
        for row in rows
    ]
    front: list[dict[str, object]] = []
    for index, values in enumerate(objectives):
        comparable = tuple(math.inf if math.isnan(value) else value for value in values)
        dominated = False
        for other_index, other_values in enumerate(objectives):
            if other_index == index:
                continue
            other = tuple(
                math.inf if math.isnan(value) else value for value in other_values
            )
            if all(left <= right for left, right in zip(other, comparable)) and any(
                left < right for left, right in zip(other, comparable)
            ):
                dominated = True
                break
        if not dominated:
            front.append(rows[index])
    return sorted(front, key=lambda row: _number(row, "movement_distance"))


def top_rows(rows: Sequence[dict[str, object]], count: int) -> list[dict[str, object]]:
    return sorted(
        rows,
        key=lambda row: (
            -_number(row, "clear_ready_coverage"),
            -_number(row, "operational_good_coverage"),
            _number(row, "worst_diameter_p95"),
            _number(row, "movement_distance"),
        ),
    )[:count]


def representative_rows(rows: Sequence[dict[str, object]], count: int) -> list[dict[str, object]]:
    if len(rows) <= count:
        return list(rows)
    return [rows[round(index * (len(rows) - 1) / (count - 1))] for index in range(count)]


def rapid_change_rows(
    rows: Sequence[dict[str, object]], spacing: float, count: int
) -> list[dict[str, object]]:
    lookup = {(_number(row, "x"), _number(row, "y")): row for row in rows}
    scored: list[tuple[float, dict[str, object]]] = []
    for row in rows:
        x, y = _number(row, "x"), _number(row, "y")
        differences = [
            abs(
                _number(row, "clear_ready_coverage")
                - _number(neighbor, "clear_ready_coverage")
            )
            for point in ((x - spacing, y), (x + spacing, y), (x, y - spacing), (x, y + spacing))
            if (neighbor := lookup.get(point)) is not None
        ]
        scored.append((max(differences, default=0.0), row))
    return [row for _, row in sorted(scored, key=lambda item: -item[0])[:count]]


def neighborhood_candidates(
    seeds: Sequence[dict[str, object]], step: float, radius: float
) -> list[Point]:
    offsets = inclusive_grid(-radius, radius, step)
    points = {
        (round(_number(seed, "x") + dx, 10), round(_number(seed, "y") + dy, 10))
        for seed in seeds
        for dx in offsets
        for dy in offsets
        if CONFIG.coarse_x_start <= _number(seed, "x") + dx <= CONFIG.coarse_x_stop
        and CONFIG.coarse_y_start <= _number(seed, "y") + dy <= CONFIG.coarse_y_stop
    }
    return [Point(x, y) for x, y in sorted(points)]


def merge_finest(rows_by_stage: Sequence[Sequence[dict[str, object]]]) -> list[dict[str, object]]:
    selected: dict[tuple[float, float], dict[str, object]] = {}
    for rows in rows_by_stage:
        for row in rows:
            key = (_number(row, "x"), _number(row, "y"))
            existing = selected.get(key)
            if existing is None or _number(row, "resolution_m") < _number(
                existing, "resolution_m"
            ):
                selected[key] = row
    return list(selected.values())


def run_strict_regression(
    candidates: Sequence[Point], output_dir: Path
) -> list[dict[str, object]]:
    epsilons = inclusive_grid(-1.0, 1.0, 0.05)
    rows: list[dict[str, object]] = []
    for sensor in candidates:
        values: dict[str, object] = {"x": sensor.x, "y": sensor.y}
        both_operational = True
        for label, radius in (("near", 6.0), ("far", 1500.0)):
            target = Point(radius, 0.0)
            state = receiver_state(sensor, target, radius)
            if state.near:
                robust_clear = False
                operational = True
                worst_diameter = math.nan
            elif state.direction_available:
                outcomes = [evaluate_q1(sensor, target, epsilon) for epsilon in epsilons]
                robust_clear = all(outcome.clear_ready for outcome in outcomes)
                operational = robust_clear
                worst_diameter = (
                    max(float(outcome.diameter) for outcome in outcomes)
                    if all(outcome.diameter is not None for outcome in outcomes)
                    else math.inf
                )
            else:
                robust_clear = False
                operational = False
                worst_diameter = math.nan
            values.update(
                {
                    f"{label}_guaranteed_rx": state.guaranteed_rx,
                    f"{label}_near": state.near,
                    f"{label}_direction_available": state.direction_available,
                    f"{label}_robust_clear_ready": robust_clear,
                    f"{label}_operational_good": operational,
                    f"{label}_worst_diameter": worst_diameter,
                }
            )
            both_operational = both_operational and operational
        values["both_operational_good"] = both_operational
        rows.append(values)
    write_csv(output_dir / "q2_strict_regression.csv", rows)
    return rows


def dense_validation(
    candidates: Sequence[dict[str, object]], output_dir: Path
) -> list[dict[str, object]]:
    radii = inclusive_grid(18.75, 1493.75, 12.5)
    alphas = inclusive_grid(-0.975, 0.975, 0.05)
    epsilons = inclusive_grid(-0.95, 0.95, 0.1)
    rows = [
        _evaluate_sensor_on_grid(
            Point(_number(candidate, "x"), _number(candidate, "y")),
            radii,
            alphas,
            epsilons,
            0.0,
            "offset_dense_validation",
        )
        for candidate in candidates
    ]
    write_csv(output_dir / "q2_final_validation.csv", rows)
    return rows


def create_figures(
    output_dir: Path,
    orthogonal_rows: Sequence[dict[str, object]],
    optima: Sequence[dict[str, object]],
    coarse_rows: Sequence[dict[str, object]],
    finest_rows: Sequence[dict[str, object]],
    pareto: Sequence[dict[str, object]],
) -> None:
    matplotlib_cache = Path(tempfile.gettempdir()) / "mathb-matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    output_dir.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(figsize=(10, 6))
    for d2 in CONFIG.orthogonal_d2:
        rows = [
            row
            for row in orthogonal_rows
            if str(row["stage"]) == "coarse" and _number(row, "d2") == d2
        ]
        axes.semilogy(
            [_number(row, "beta_deg") for row in rows],
            [_number(row, "worst_case_diameter") for row in rows],
            label=f"d2={d2:.0f} m",
        )
    axes.axvline(90.0, color="#444444", linestyle="--", linewidth=1.0, label="90 deg")
    axes.axhline(CONFIG.diameter_limit, color="#d73027", linestyle=":", linewidth=1.2, label="D=40 m")
    axes.set_xlim(3.5, 176.5)
    axes.set_xlabel("True line-intersection angle beta (deg)")
    axes.set_ylabel("Worst sampled diameter (m, log scale)")
    axes.set_title("Q1 exact geometry: orthogonal-intersection validation")
    axes.grid(alpha=0.25)
    axes.legend(ncol=2)
    figure.tight_layout()
    figure.savefig(output_dir / "fig_q2_beta_validation.png", dpi=180)
    plt.close(figure)

    x_values = inclusive_grid(CONFIG.coarse_x_start, CONFIG.coarse_x_stop, CONFIG.coarse_step)
    y_values = inclusive_grid(CONFIG.coarse_y_start, CONFIG.coarse_y_stop, CONFIG.coarse_step)
    lookup = {
        (_number(row, "x"), _number(row, "y")): _number(row, "clear_ready_coverage")
        for row in coarse_rows
    }
    matrix = [[lookup[(x, y)] for x in x_values] for y in y_values]
    extent = (
        x_values[0] - CONFIG.coarse_step / 2.0,
        x_values[-1] + CONFIG.coarse_step / 2.0,
        y_values[0] - CONFIG.coarse_step / 2.0,
        y_values[-1] + CONFIG.coarse_step / 2.0,
    )
    figure, axes = plt.subplots(figsize=(11, 7))
    image = axes.imshow(
        matrix,
        origin="lower",
        extent=extent,
        aspect="equal",
        cmap="viridis",
        vmin=0.0,
        vmax=max(max(row) for row in matrix),
        interpolation="nearest",
    )
    axes.scatter((0.0,), (0.0,), marker="*", color="#d73027", s=70, label="S1")
    axes.set_xlabel("S2 x (m)")
    axes.set_ylabel("S2 y (m)")
    axes.set_title("Area-weighted robust clear-ready coverage")
    figure.colorbar(image, ax=axes, label="target-area coverage")
    axes.legend()
    figure.tight_layout()
    figure.savefig(output_dir / "fig_q2_coverage_heatmap.png", dpi=180)
    plt.close(figure)

    best = top_rows(finest_rows, 1)[0]
    top_threshold = _number(best, "clear_ready_coverage") * 0.95
    high = [
        row for row in finest_rows if _number(row, "clear_ready_coverage") >= top_threshold
    ]
    figure, axes = plt.subplots(figsize=(10, 7))
    axes.scatter(
        [_number(row, "x") for row in coarse_rows],
        [_number(row, "y") for row in coarse_rows],
        color="#9e9e9e",
        s=14,
        alpha=0.30,
        label="coarse scan",
    )
    axes.scatter(
        [_number(row, "x") for row in high],
        [_number(row, "y") for row in high],
        facecolors="none",
        edgecolors="#d73027",
        s=45,
        label=">=95% of sampled maximum",
    )
    displayed_pareto = representative_rows(pareto, 80)
    axes.scatter(
        [_number(row, "x") for row in displayed_pareto],
        [_number(row, "y") for row in displayed_pareto],
        marker="x",
        color="#111111",
        s=24,
        label="representative Pareto candidates",
    )
    axes.scatter(
        (_number(best, "x"),),
        (_number(best, "y"),),
        marker="*",
        color="#fdae61",
        edgecolors="#111111",
        s=150,
        label="numerical maximum",
    )
    axes.axvline(
        CONFIG.orthogonal_target_x,
        color="#2166ac",
        linestyle="--",
        linewidth=1.0,
        label="orthogonal line for g=(800,0)",
    )
    axes.set_xlim(-200.0, 1300.0)
    axes.set_ylim(0.0, 1200.0)
    axes.set_aspect("equal", adjustable="box")
    axes.set_xlabel("S2 x (m)")
    axes.set_ylabel("S2 y (m)")
    axes.set_title("Coarse-to-fine candidate region and Pareto set")
    axes.legend(loc="upper right")
    axes.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(output_dir / "fig_q2_candidate_region.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(figsize=(10, 6))
    finite_p95 = sorted(
        _number(row, "worst_diameter_p95")
        for row in finest_rows
        if math.isfinite(_number(row, "worst_diameter_p95"))
        and _number(row, "worst_diameter_p95") > 0.0
    )
    color_cap = finite_p95[round(0.95 * (len(finite_p95) - 1))]
    color_floor = max(1.0, finite_p95[0])
    color_values = [
        min(_number(row, "worst_diameter_p95"), color_cap)
        if math.isfinite(_number(row, "worst_diameter_p95"))
        else math.nan
        for row in finest_rows
    ]
    scatter = axes.scatter(
        [_number(row, "movement_distance") for row in finest_rows],
        [_number(row, "clear_ready_coverage") for row in finest_rows],
        c=color_values,
        cmap="plasma_r",
        norm=LogNorm(vmin=color_floor, vmax=color_cap),
        s=14,
        alpha=0.35,
    )
    axes.scatter(
        [_number(row, "movement_distance") for row in pareto],
        [_number(row, "clear_ready_coverage") for row in pareto],
        facecolors="none",
        edgecolors="#111111",
        linewidths=0.7,
        s=24,
        label="four-metric Pareto set",
    )
    axes.set_xlabel("Movement distance from S1 (m)")
    axes.set_ylabel("Area-weighted clear-ready coverage")
    axes.set_title("Coverage-distance trade-off")
    axes.grid(alpha=0.25)
    axes.legend()
    figure.colorbar(
        scatter,
        ax=axes,
        label="p95 worst-case diameter (m, log scale; upper 5% clipped)",
    )
    figure.tight_layout()
    figure.savefig(output_dir / "fig_q2_pareto.png", dpi=180)
    plt.close(figure)


def benchmark(evaluations: int = 5000) -> dict[str, float]:
    sensor = Point(550.0, 500.0)
    target = Point(800.0, 0.0)
    epsilons = inclusive_grid(-1.0, 1.0, 0.1)
    started = time.perf_counter()
    for index in range(evaluations):
        evaluate_q1(sensor, target, epsilons[index % len(epsilons)])
    elapsed = time.perf_counter() - started
    return {
        "evaluations": evaluations,
        "elapsed_seconds": elapsed,
        "evaluations_per_second": evaluations / elapsed,
    }


def run_all(output_dir: Path, workers: int, resume: bool) -> dict[str, object]:
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_result = benchmark()
    print(json.dumps({"benchmark": benchmark_result}, indent=2), flush=True)

    orthogonal_rows, optima = run_orthogonal_experiment(output_dir)
    coarse_candidates = [
        Point(x, y)
        for y in inclusive_grid(
            CONFIG.coarse_y_start, CONFIG.coarse_y_stop, CONFIG.coarse_step
        )
        for x in inclusive_grid(
            CONFIG.coarse_x_start, CONFIG.coarse_x_stop, CONFIG.coarse_step
        )
    ]
    coarse = scan_candidates(
        coarse_candidates,
        CONFIG.coarse_step,
        "coarse_100m",
        output_dir / "q2_coarse_scan.csv",
        workers,
        resume,
    )

    top_five_percent = top_rows(coarse, max(1, math.ceil(0.05 * len(coarse))))
    coarse_pareto = pareto_front(coarse)
    stage_50_seeds = [
        *top_five_percent,
        *representative_rows(coarse_pareto, 8),
        *rapid_change_rows(coarse, CONFIG.coarse_step, 10),
    ]
    candidates_50 = neighborhood_candidates(stage_50_seeds, 50.0, 50.0)
    stage_50 = scan_candidates(
        candidates_50,
        50.0,
        "refine_50m",
        output_dir / "q2_refine_50m.csv",
        workers,
        resume,
    )

    stage_50_pareto = pareto_front(stage_50)
    stage_10_seeds = [
        *top_rows(stage_50, 5),
        *representative_rows(stage_50_pareto, 3),
        *rapid_change_rows(stage_50, 50.0, 3),
    ]
    candidates_10 = neighborhood_candidates(stage_10_seeds, 10.0, 30.0)
    stage_10 = scan_candidates(
        candidates_10,
        10.0,
        "refine_10m",
        output_dir / "q2_refine_10m.csv",
        workers,
        resume,
    )

    stage_10_pareto = pareto_front(stage_10)
    stage_2_seeds = [
        *top_rows(stage_10, 3),
        *representative_rows(stage_10_pareto, 2),
    ]
    candidates_2 = neighborhood_candidates(stage_2_seeds, 2.0, 8.0)
    stage_2 = scan_candidates(
        candidates_2,
        2.0,
        "refine_2m",
        output_dir / "q2_refine_2m.csv",
        workers,
        resume,
    )

    finest = merge_finest((coarse, stage_50, stage_10, stage_2))
    top = top_rows(finest, 100)
    top_ranked = [dict(rank=index + 1, **row) for index, row in enumerate(top)]
    write_csv(output_dir / "q2_top_candidates.csv", top_ranked)
    pareto = pareto_front(finest)
    pareto_ranked = [dict(pareto_index=index + 1, **row) for index, row in enumerate(pareto)]
    write_csv(output_dir / "q2_pareto_candidates.csv", pareto_ranked)

    dense = dense_validation(top[:5], output_dir)
    strict_candidate_points = {
        (_number(row, "x"), _number(row, "y")) for row in finest
    }
    strict_candidate_points.update(
        (float(x), float(y))
        for center in (6.0, 1500.0)
        for x in range(round(center) - 5, round(center) + 6)
        for y in range(0, 6)
    )
    strict_rows = run_strict_regression(
        [Point(x, y) for x, y in sorted(strict_candidate_points)],
        output_dir,
    )
    create_figures(output_dir, orthogonal_rows, optima, coarse, finest, pareto)

    best = top[0]
    max_coverage = _number(best, "clear_ready_coverage")
    short_candidates = [
        row for row in finest if _number(row, "clear_ready_coverage") >= 0.90 * max_coverage
    ]
    shortest_high = min(short_candidates, key=lambda row: _number(row, "movement_distance"))
    radii, alphas, epsilons = optimization_grids()
    baseline_rows: list[dict[str, object]] = []
    for optimum in optima:
        baseline = _evaluate_sensor_on_grid(
            Point(_number(optimum, "sensor_x"), _number(optimum, "sensor_y")),
            radii,
            alphas,
            epsilons,
            0.0,
            "orthogonal_reference",
        )
        baseline_rows.append(
            dict(baseline=f"orthogonal_d2_{_number(optimum, 'd2'):.0f}m", **baseline)
        )
    baseline_rows.extend(
        (
            dict(baseline="global_numerical_best", **best),
            dict(baseline="shortest_at_least_90pct_max", **shortest_high),
            dict(
                baseline="highest_coverage_pareto",
                **max(pareto, key=lambda row: _number(row, "clear_ready_coverage")),
            ),
        )
    )
    write_csv(output_dir / "q2_baseline_comparison.csv", baseline_rows)
    summary: dict[str, object] = {
        "configuration": asdict(CONFIG),
        "benchmark": benchmark_result,
        "workers": workers,
        "candidate_counts": {
            "coarse_100m": len(coarse),
            "refine_50m": len(stage_50),
            "refine_10m": len(stage_10),
            "refine_2m": len(stage_2),
            "unique_candidates": len(finest),
            "pareto_candidates": len(pareto),
        },
        "orthogonal_optima": optima,
        "numerical_best": best,
        "shortest_candidate_at_least_90pct_max_coverage": shortest_high,
        "baseline_comparison": baseline_rows,
        "dense_offset_validation": dense,
        "strict_regression_both_operational_count": sum(
            str(row["both_operational_good"]).lower() == "true" for row in strict_rows
        ),
        "strict_regression_candidate_count": len(strict_rows),
        "elapsed_seconds": time.perf_counter() - started,
    }
    (output_dir / "q2_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, allow_nan=True), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("results/q2"))
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
    )
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--benchmark-only", action="store_true")
    arguments = parser.parse_args()
    if arguments.workers <= 0:
        parser.error("--workers must be positive")
    if arguments.benchmark_only:
        print(json.dumps(benchmark(), indent=2))
        return
    run_all(arguments.output_dir, arguments.workers, not arguments.no_resume)


if __name__ == "__main__":
    main()
