"""Final Q2 geometry experiments built on the frozen Q1 exact geometry.

This module closes the remaining Q2 numerical questions: the near-orthogonal
minimax correction, target-space failure structure, local candidate
convergence, and weighting/Pareto sensitivity.  It does not model scheduling,
routing, or Q3/Q4.
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
import sys
import tempfile
import time
from typing import Iterable, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.q2.run_q2_key_validation import (
    CONFIG,
    Point,
    _evaluate_sensor_on_grid,
    _number,
    bearing,
    evaluate_q1,
    inclusive_grid,
    merge_finest,
    pareto_front,
    read_csv,
    receiver_state,
    scan_candidates,
    target_point,
    top_rows,
    weighted_mean,
    weighted_percentile,
    write_csv,
)


FINALIZATION_CACHE_VERSION = 1


@dataclass(frozen=True, slots=True)
class FinalizationConfig:
    orthogonal_d1: tuple[float, ...] = (
        200.0,
        400.0,
        600.0,
        800.0,
        1000.0,
        1200.0,
        1500.0,
    )
    orthogonal_d2: tuple[float, ...] = (200.0, 400.0, 600.0, 800.0, 1000.0)
    orthogonal_beta_start: float = 88.0
    orthogonal_beta_stop: float = 93.0
    orthogonal_beta_step: float = 0.05
    orthogonal_epsilon_step: float = 0.02
    orthogonal_refine_half_width: float = 0.10
    orthogonal_refine_beta_step: float = 0.01
    orthogonal_refine_epsilon_step: float = 0.01
    local_x_start: float = 630.0
    local_x_stop: float = 730.0
    local_y_start: float = 360.0
    local_y_stop: float = 480.0
    local_step_2m: float = 2.0
    local_step_1m: float = 1.0
    local_step_05m: float = 0.5
    local_seed_count: int = 12
    local_radius_1m: float = 5.0
    local_radius_05m: float = 2.5
    stable_fraction: float = 0.99
    convergence_r_start: float = 18.75
    convergence_r_stop: float = 1493.75
    convergence_r_step: float = 12.5
    convergence_alpha_start: float = -0.975
    convergence_alpha_stop: float = 0.975
    convergence_alpha_step: float = 0.05
    convergence_epsilon_steps: tuple[float, ...] = (0.1, 0.05, 0.02, 0.01)
    target_map_r_start: float = 12.5
    target_map_r_stop: float = 1500.0
    target_map_r_step: float = 12.5
    target_map_alpha_start: float = -1.0
    target_map_alpha_stop: float = 1.0
    target_map_alpha_step: float = 0.025
    target_map_epsilon_step: float = 0.02


FINAL = FinalizationConfig()


@dataclass(frozen=True, slots=True)
class CrossValidationConfig:
    x_start: float = 680.0
    x_stop: float = 700.0
    x_step: float = 2.0
    y_start: float = 420.0
    y_stop: float = 480.0
    y_step: float = 5.0
    refine_step: float = 1.0
    refine_radius: float = 3.0
    seed_count: int = 5
    r_start: float = 12.5
    r_stop: float = 1500.0
    r_step: float = 12.5
    alpha_start: float = -1.0
    alpha_stop: float = 1.0
    alpha_step: float = 0.05
    epsilon_step: float = 0.05
    ultra_r_start: float = 9.375
    ultra_r_stop: float = 1496.875
    ultra_r_step: float = 6.25
    ultra_alpha_start: float = -0.99375
    ultra_alpha_stop: float = 0.99375
    ultra_alpha_step: float = 0.0125
    ultra_epsilon_step: float = 0.02
    ultra_candidate_count: int = 5


CROSS = CrossValidationConfig()


FAILURE_REASONS = (
    "near_direct_clear",
    "guaranteed_rx_failure",
    "localization_unbounded",
    "localization_empty",
    "diameter_over_40",
    "diameter_circle_failure",
    "robust_clear_ready",
)


def _signature(kind: str, extra: object) -> str:
    payload = {
        "cache_version": FINALIZATION_CACHE_VERSION,
        "base_config": asdict(CONFIG),
        "final_config": asdict(FINAL),
        "kind": kind,
        "extra": extra,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _cache_hit(path: Path, signature: str, resume: bool) -> bool:
    metadata = path.with_suffix(".meta.json")
    if not resume or not path.exists() or not metadata.exists():
        return False
    return json.loads(metadata.read_text(encoding="utf-8")).get("signature") == signature


def _write_metadata(path: Path, signature: str, started: float, count: int) -> None:
    path.with_suffix(".meta.json").write_text(
        json.dumps(
            {
                "signature": signature,
                "row_count": count,
                "elapsed_seconds": time.perf_counter() - started,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _orthogonal_metrics(
    d1: float, d2: float, beta_deg: float, epsilons: Sequence[float]
) -> dict[str, object]:
    target = Point(d1, 0.0)
    beta = math.radians(beta_deg)
    sensor = Point(target.x - d2 * math.cos(beta), -d2 * math.sin(beta))
    outcomes = [evaluate_q1(sensor, target, epsilon) for epsilon in epsilons]
    if any(not outcome.bounded for outcome in outcomes):
        worst_diameter = math.inf
        worst_epsilon = next(
            epsilon
            for outcome, epsilon in zip(outcomes, epsilons)
            if not outcome.bounded
        )
    else:
        worst_diameter, worst_epsilon = max(
            (float(outcome.diameter), epsilon)
            for outcome, epsilon in zip(outcomes, epsilons)
        )
    return {
        "d1": d1,
        "d2": d2,
        "beta_deg": beta_deg,
        "sensor_x": sensor.x,
        "sensor_y": sensor.y,
        "worst_diameter": worst_diameter,
        "mean_diameter": (
            sum(float(outcome.diameter) for outcome in outcomes) / len(outcomes)
            if all(outcome.diameter is not None for outcome in outcomes)
            else math.inf
        ),
        "worst_epsilon": worst_epsilon,
        "bounded_all": all(outcome.bounded for outcome in outcomes),
        "diameter_pass_all": all(outcome.diameter_pass for outcome in outcomes),
        "circle_pass_all": all(outcome.circle_pass for outcome in outcomes),
        "clear_ready_all": all(outcome.clear_ready for outcome in outcomes),
        "epsilon_count": len(epsilons),
    }


def _orthogonal_worker(payload: tuple[float, float]) -> tuple[list[dict[str, object]], dict[str, object]]:
    d1, d2 = payload
    coarse_epsilons = inclusive_grid(-1.0, 1.0, FINAL.orthogonal_epsilon_step)
    coarse_rows = [
        _orthogonal_metrics(d1, d2, beta, coarse_epsilons)
        for beta in inclusive_grid(
            FINAL.orthogonal_beta_start,
            FINAL.orthogonal_beta_stop,
            FINAL.orthogonal_beta_step,
        )
    ]
    coarse_best = min(coarse_rows, key=lambda row: _number(row, "worst_diameter"))
    center = _number(coarse_best, "beta_deg")
    refined_epsilons = inclusive_grid(-1.0, 1.0, FINAL.orthogonal_refine_epsilon_step)
    refined_rows = [
        _orthogonal_metrics(d1, d2, beta, refined_epsilons)
        for beta in inclusive_grid(
            center - FINAL.orthogonal_refine_half_width,
            center + FINAL.orthogonal_refine_half_width,
            FINAL.orthogonal_refine_beta_step,
        )
    ]
    robust_best = min(refined_rows, key=lambda row: _number(row, "worst_diameter"))
    at_90 = _orthogonal_metrics(d1, d2, 90.0, refined_epsilons)
    nominal_rows = [
        _orthogonal_metrics(d1, d2, beta, (0.0,))
        for beta in inclusive_grid(
            FINAL.orthogonal_beta_start,
            FINAL.orthogonal_beta_stop,
            FINAL.orthogonal_refine_beta_step,
        )
    ]
    nominal_best = min(nominal_rows, key=lambda row: _number(row, "worst_diameter"))
    summary = {
        "d1": d1,
        "d2": d2,
        "nominal_optimum_beta": _number(nominal_best, "beta_deg"),
        "robust_optimum_beta": _number(robust_best, "beta_deg"),
        "deviation_from_90_deg": _number(robust_best, "beta_deg") - 90.0,
        "beta_90_worst_diameter": _number(at_90, "worst_diameter"),
        "optimum_worst_diameter": _number(robust_best, "worst_diameter"),
        "diameter_improvement_from_90": _number(at_90, "worst_diameter")
        - _number(robust_best, "worst_diameter"),
        "relative_improvement_from_90": (
            (_number(at_90, "worst_diameter") - _number(robust_best, "worst_diameter"))
            / _number(at_90, "worst_diameter")
        ),
        "worst_epsilon": robust_best["worst_epsilon"],
        "bounded_all_at_optimum": robust_best["bounded_all"],
        "diameter_pass_all_at_optimum": robust_best["diameter_pass_all"],
        "circle_pass_all_at_optimum": robust_best["circle_pass_all"],
        "clear_ready_all_at_optimum": robust_best["clear_ready_all"],
        "refined_beta_step": FINAL.orthogonal_refine_beta_step,
        "refined_epsilon_step": FINAL.orthogonal_refine_epsilon_step,
    }
    curves = [dict(stage="coarse", **row) for row in coarse_rows]
    curves.extend(dict(stage="refined", **row) for row in refined_rows)
    return curves, summary


def run_orthogonal_correction(
    output_dir: Path, workers: int, resume: bool
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    summary_path = output_dir / "q2_orthogonal_correction.csv"
    curves_path = output_dir / "q2_orthogonal_correction_curves.csv"
    pairs = [(d1, d2) for d1 in FINAL.orthogonal_d1 for d2 in FINAL.orthogonal_d2]
    signature = _signature("orthogonal_correction", pairs)
    if _cache_hit(summary_path, signature, resume) and curves_path.exists():
        print("cache hit: orthogonal correction", flush=True)
        return [dict(row) for row in read_csv(curves_path)], [dict(row) for row in read_csv(summary_path)]
    started = time.perf_counter()
    if workers == 1:
        results: Iterable[tuple[list[dict[str, object]], dict[str, object]]] = map(
            _orthogonal_worker, pairs
        )
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        results = executor.map(_orthogonal_worker, pairs, chunksize=1)
    curves: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    try:
        for index, (pair_curves, summary) in enumerate(results, 1):
            curves.extend(pair_curves)
            summaries.append(summary)
            print(f"orthogonal correction: {index}/{len(pairs)}", flush=True)
    finally:
        if executor is not None:
            executor.shutdown()
    write_csv(curves_path, curves)
    write_csv(summary_path, summaries)
    _write_metadata(summary_path, signature, started, len(summaries))
    return curves, summaries


def _candidate_grid(
    x_start: float,
    x_stop: float,
    y_start: float,
    y_stop: float,
    step: float,
) -> list[Point]:
    return [
        Point(x, y)
        for y in inclusive_grid(y_start, y_stop, step)
        for x in inclusive_grid(x_start, x_stop, step)
    ]


def _seeded_grid(
    seeds: Sequence[dict[str, object]], step: float, radius: float
) -> list[Point]:
    offsets = inclusive_grid(-radius, radius, step)
    return [
        Point(x, y)
        for x, y in sorted(
            {
                (
                    round(_number(seed, "x") + dx, 10),
                    round(_number(seed, "y") + dy, 10),
                )
                for seed in seeds
                for dx in offsets
                for dy in offsets
                if FINAL.local_x_start <= _number(seed, "x") + dx <= FINAL.local_x_stop
                and FINAL.local_y_start <= _number(seed, "y") + dy <= FINAL.local_y_stop
            }
        )
    ]


def _seeded_grid_with_bounds(
    seeds: Sequence[dict[str, object]],
    step: float,
    radius: float,
    x_start: float,
    x_stop: float,
    y_start: float,
    y_stop: float,
) -> list[Point]:
    offsets = inclusive_grid(-radius, radius, step)
    return [
        Point(x, y)
        for x, y in sorted(
            {
                (
                    round(_number(seed, "x") + dx, 10),
                    round(_number(seed, "y") + dy, 10),
                )
                for seed in seeds
                for dx in offsets
                for dy in offsets
                if x_start <= _number(seed, "x") + dx <= x_stop
                and y_start <= _number(seed, "y") + dy <= y_stop
            }
        )
    ]


def _dual_metric_seeds(rows: Sequence[dict[str, object]], count: int) -> list[dict[str, object]]:
    area = sorted(rows, key=lambda row: -_number(row, "clear_ready_coverage"))[:count]
    uniform = sorted(rows, key=lambda row: -_number(row, "uniform_clear_ready_coverage"))[:count]
    selected: dict[tuple[float, float], dict[str, object]] = {}
    for row in (*area, *uniform):
        selected[(_number(row, "x"), _number(row, "y"))] = row
    return list(selected.values())


def run_local_refinement(
    output_dir: Path, workers: int, resume: bool
) -> tuple[list[dict[str, object]], dict[str, object]]:
    stage_2 = scan_candidates(
        _candidate_grid(
            FINAL.local_x_start,
            FINAL.local_x_stop,
            FINAL.local_y_start,
            FINAL.local_y_stop,
            FINAL.local_step_2m,
        ),
        FINAL.local_step_2m,
        "final_refine_2m",
        output_dir / "q2_candidate_refine_2m.csv",
        workers,
        resume,
    )
    seeds_1 = _dual_metric_seeds(stage_2, FINAL.local_seed_count)
    stage_1 = scan_candidates(
        _seeded_grid(seeds_1, FINAL.local_step_1m, FINAL.local_radius_1m),
        FINAL.local_step_1m,
        "final_refine_1m",
        output_dir / "q2_candidate_refine_1m.csv",
        workers,
        resume,
    )
    seeds_05 = _dual_metric_seeds(stage_1, FINAL.local_seed_count)
    stage_05 = scan_candidates(
        _seeded_grid(seeds_05, FINAL.local_step_05m, FINAL.local_radius_05m),
        FINAL.local_step_05m,
        "final_refine_0p5m",
        output_dir / "q2_candidate_refine_0p5m.csv",
        workers,
        resume,
    )
    merged = merge_finest((stage_2, stage_1, stage_05))
    best_area = max(merged, key=lambda row: _number(row, "clear_ready_coverage"))
    best_uniform = max(
        merged, key=lambda row: _number(row, "uniform_clear_ready_coverage")
    )
    area_threshold = FINAL.stable_fraction * _number(best_area, "clear_ready_coverage")
    uniform_threshold = FINAL.stable_fraction * _number(
        best_uniform, "uniform_clear_ready_coverage"
    )
    enriched = []
    for row in merged:
        item = dict(row)
        item["area_rank"] = 0
        item["uniform_rank"] = 0
        item["in_area_99pct_region"] = _number(row, "clear_ready_coverage") >= area_threshold
        item["in_uniform_99pct_region"] = (
            _number(row, "uniform_clear_ready_coverage") >= uniform_threshold
        )
        enriched.append(item)
    for rank, row in enumerate(
        sorted(enriched, key=lambda item: -_number(item, "clear_ready_coverage")), 1
    ):
        row["area_rank"] = rank
    for rank, row in enumerate(
        sorted(
            enriched,
            key=lambda item: -_number(item, "uniform_clear_ready_coverage"),
        ),
        1,
    ):
        row["uniform_rank"] = rank
    write_csv(output_dir / "q2_candidate_refinement.csv", enriched)
    stage_best = []
    for stage, rows in (("2m", stage_2), ("1m", stage_1), ("0.5m", stage_05)):
        row = max(rows, key=lambda item: _number(item, "clear_ready_coverage"))
        stage_best.append(
            {
                "stage": stage,
                "x": _number(row, "x"),
                "y": _number(row, "y"),
                "clear_ready_coverage": _number(row, "clear_ready_coverage"),
                "uniform_clear_ready_coverage": _number(
                    row, "uniform_clear_ready_coverage"
                ),
            }
        )
    write_csv(output_dir / "q2_candidate_refinement_best_by_stage.csv", stage_best)
    area_stable = [row for row in enriched if row["in_area_99pct_region"]]
    uniform_stable = [row for row in enriched if row["in_uniform_99pct_region"]]

    def bounds(rows: Sequence[dict[str, object]]) -> dict[str, object]:
        return {
            "count": len(rows),
            "x_min": min(_number(row, "x") for row in rows),
            "x_max": max(_number(row, "x") for row in rows),
            "y_min": min(_number(row, "y") for row in rows),
            "y_max": max(_number(row, "y") for row in rows),
        }

    summary = {
        "best_area": best_area,
        "best_uniform": best_uniform,
        "area_99pct_threshold": area_threshold,
        "uniform_99pct_threshold": uniform_threshold,
        "area_99pct_bounds": bounds(area_stable),
        "uniform_99pct_bounds": bounds(uniform_stable),
        "best_by_stage": stage_best,
        "candidate_counts": {
            "2m": len(stage_2),
            "1m": len(stage_1),
            "0.5m": len(stage_05),
            "merged": len(enriched),
        },
    }
    return enriched, summary


def _convergence_worker(
    payload: tuple[str, float, float, float]
) -> dict[str, object]:
    label, x, y, epsilon_step = payload
    radii = inclusive_grid(
        FINAL.convergence_r_start,
        FINAL.convergence_r_stop,
        FINAL.convergence_r_step,
    )
    alphas = inclusive_grid(
        FINAL.convergence_alpha_start,
        FINAL.convergence_alpha_stop,
        FINAL.convergence_alpha_step,
    )
    epsilons = inclusive_grid(-1.0, 1.0, epsilon_step)
    return dict(
        candidate=label,
        epsilon_step=epsilon_step,
        **_evaluate_sensor_on_grid(
            Point(x, y), radii, alphas, epsilons, 0.0, "epsilon_convergence"
        ),
    )


def run_epsilon_convergence(
    candidates: Sequence[tuple[str, Point]],
    output_dir: Path,
    workers: int,
    resume: bool,
) -> list[dict[str, object]]:
    path = output_dir / "q2_epsilon_convergence.csv"
    payloads = [
        (label, point.x, point.y, step)
        for label, point in candidates
        for step in FINAL.convergence_epsilon_steps
    ]
    signature = _signature("epsilon_convergence", payloads)
    if _cache_hit(path, signature, resume):
        print("cache hit: epsilon convergence", flush=True)
        return [dict(row) for row in read_csv(path)]
    started = time.perf_counter()
    if workers == 1:
        iterator: Iterable[dict[str, object]] = map(_convergence_worker, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        iterator = executor.map(_convergence_worker, payloads, chunksize=1)
    rows: list[dict[str, object]] = []
    try:
        for index, row in enumerate(iterator, 1):
            rows.append(row)
            print(f"epsilon convergence: {index}/{len(payloads)}", flush=True)
    finally:
        if executor is not None:
            executor.shutdown()
    write_csv(path, rows)
    _write_metadata(path, signature, started, len(rows))
    return rows


def _cross_grid_worker(payload: tuple[float, float, float, str]) -> dict[str, object]:
    x, y, resolution_m, stage = payload
    return _evaluate_sensor_on_grid(
        Point(x, y),
        inclusive_grid(CROSS.r_start, CROSS.r_stop, CROSS.r_step),
        inclusive_grid(CROSS.alpha_start, CROSS.alpha_stop, CROSS.alpha_step),
        inclusive_grid(-1.0, 1.0, CROSS.epsilon_step),
        resolution_m,
        stage,
    )


def _scan_cross_grid(
    candidates: Sequence[Point],
    resolution_m: float,
    stage: str,
    path: Path,
    workers: int,
    resume: bool,
) -> list[dict[str, object]]:
    signature = _signature(
        "cross_grid",
        {
            "config": asdict(CROSS),
            "stage": stage,
            "resolution_m": resolution_m,
            "candidates": [(point.x, point.y) for point in candidates],
        },
    )
    if _cache_hit(path, signature, resume):
        print(f"cache hit: {path}", flush=True)
        return [dict(row) for row in read_csv(path)]
    started = time.perf_counter()
    payloads = [(point.x, point.y, resolution_m, stage) for point in candidates]
    if workers == 1:
        iterator: Iterable[dict[str, object]] = map(_cross_grid_worker, payloads)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        iterator = executor.map(_cross_grid_worker, payloads, chunksize=1)
    rows: list[dict[str, object]] = []
    try:
        for index, row in enumerate(iterator, 1):
            rows.append(row)
            if index % 25 == 0 or index == len(payloads):
                print(f"{stage}: {index}/{len(payloads)}", flush=True)
    finally:
        if executor is not None:
            executor.shutdown()
    write_csv(path, rows)
    _write_metadata(path, signature, started, len(rows))
    return rows


def _ultra_worker(payload: tuple[float, float]) -> dict[str, object]:
    x, y = payload
    return _evaluate_sensor_on_grid(
        Point(x, y),
        inclusive_grid(CROSS.ultra_r_start, CROSS.ultra_r_stop, CROSS.ultra_r_step),
        inclusive_grid(
            CROSS.ultra_alpha_start, CROSS.ultra_alpha_stop, CROSS.ultra_alpha_step
        ),
        inclusive_grid(-1.0, 1.0, CROSS.ultra_epsilon_step),
        0.0,
        "ultra_cross_validation",
    )


def run_cross_grid_validation(
    output_dir: Path, workers: int, resume: bool
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    coarse_candidates = [
        Point(x, y)
        for y in inclusive_grid(CROSS.y_start, CROSS.y_stop, CROSS.y_step)
        for x in inclusive_grid(CROSS.x_start, CROSS.x_stop, CROSS.x_step)
    ]
    coarse = _scan_cross_grid(
        coarse_candidates,
        CROSS.x_step,
        "cross_grid_coarse",
        output_dir / "q2_candidate_cross_grid.csv",
        workers,
        resume,
    )
    extension_y_stop = 520.0
    extension = _scan_cross_grid(
        [
            Point(x, y)
            for y in inclusive_grid(CROSS.y_stop + CROSS.y_step, extension_y_stop, CROSS.y_step)
            for x in inclusive_grid(CROSS.x_start, CROSS.x_stop, CROSS.x_step)
        ],
        CROSS.x_step,
        "cross_grid_upper_extension",
        output_dir / "q2_candidate_cross_grid_upper_extension.csv",
        workers,
        resume,
    )
    pre_refinement = merge_finest((coarse, extension))
    seeds = top_rows(pre_refinement, CROSS.seed_count)
    refined_candidates = _seeded_grid_with_bounds(
        seeds,
        CROSS.refine_step,
        CROSS.refine_radius,
        CROSS.x_start,
        CROSS.x_stop,
        CROSS.y_start,
        extension_y_stop,
    )
    refined = _scan_cross_grid(
        refined_candidates,
        CROSS.refine_step,
        "cross_grid_refined",
        output_dir / "q2_candidate_cross_grid_refined.csv",
        workers,
        resume,
    )
    combined = merge_finest((coarse, extension, refined))
    ultra_seed_rows = top_rows(combined, CROSS.ultra_candidate_count)
    ultra_points = {
        (_number(row, "x"), _number(row, "y")) for row in ultra_seed_rows
    }
    ultra_points.update({(686.0, 430.0), (690.0, 462.0), (690.0, 468.0)})
    payloads = sorted(ultra_points)
    ultra_path = output_dir / "q2_candidate_ultra_validation.csv"
    signature = _signature(
        "ultra_cross_validation", {"config": asdict(CROSS), "candidates": payloads}
    )
    if _cache_hit(ultra_path, signature, resume):
        print("cache hit: ultra cross validation", flush=True)
        ultra = [dict(row) for row in read_csv(ultra_path)]
    else:
        started = time.perf_counter()
        if workers == 1:
            iterator: Iterable[dict[str, object]] = map(_ultra_worker, payloads)
            executor = None
        else:
            executor = ProcessPoolExecutor(max_workers=min(workers, len(payloads)))
            iterator = executor.map(_ultra_worker, payloads, chunksize=1)
        ultra = []
        try:
            for index, row in enumerate(iterator, 1):
                ultra.append(row)
                print(f"ultra cross validation: {index}/{len(payloads)}", flush=True)
        finally:
            if executor is not None:
                executor.shutdown()
        write_csv(ultra_path, ultra)
        _write_metadata(ultra_path, signature, started, len(ultra))
    best = top_rows(ultra, 1)[0]
    summary = {
        "config": asdict(CROSS),
        "coarse_candidate_count": len(coarse),
        "upper_extension_candidate_count": len(extension),
        "refined_candidate_count": len(refined),
        "ultra_candidate_count": len(ultra),
        "best": best,
        "coverage_range_ultra": {
            "min": min(_number(row, "clear_ready_coverage") for row in ultra),
            "max": max(_number(row, "clear_ready_coverage") for row in ultra),
        },
    }
    return combined, ultra, summary


def classify_target(
    sensor: Point,
    radius: float,
    alpha_deg: float,
    epsilons: Sequence[float],
) -> dict[str, object]:
    target = target_point(radius, alpha_deg)
    state = receiver_state(sensor, target, radius)
    base: dict[str, object] = {
        "r": radius,
        "alpha_deg": alpha_deg,
        "target_x": target.x,
        "target_y": target.y,
        "distance_to_target": state.distance_to_target,
        "guaranteed_rx": state.guaranteed_rx,
        "possible_rx": state.possible_rx,
        "near": state.near,
        "direction_available": state.direction_available,
        "operational_good": False,
        "robust_clear_ready": False,
        "worst_diameter": math.nan,
        "first_failure_epsilon": math.nan,
    }
    if state.near:
        return {**base, "failure_reason": "near_direct_clear", "operational_good": True}
    if not state.guaranteed_rx:
        return dict(failure_reason="guaranteed_rx_failure", **base)
    outcomes = [evaluate_q1(sensor, target, epsilon) for epsilon in epsilons]
    finite = [float(outcome.diameter) for outcome in outcomes if outcome.diameter is not None]
    base["worst_diameter"] = max(finite) if finite else math.nan
    for outcome, epsilon in zip(outcomes, epsilons):
        if outcome.unbounded:
            base["first_failure_epsilon"] = epsilon
            return dict(failure_reason="localization_unbounded", **base)
    for outcome, epsilon in zip(outcomes, epsilons):
        if outcome.empty:
            base["first_failure_epsilon"] = epsilon
            return dict(failure_reason="localization_empty", **base)
    for outcome, epsilon in zip(outcomes, epsilons):
        if not outcome.diameter_pass:
            base["first_failure_epsilon"] = epsilon
            return dict(failure_reason="diameter_over_40", **base)
    for outcome, epsilon in zip(outcomes, epsilons):
        if not outcome.circle_pass:
            base["first_failure_epsilon"] = epsilon
            return dict(failure_reason="diameter_circle_failure", **base)
    base["operational_good"] = True
    base["robust_clear_ready"] = True
    return dict(failure_reason="robust_clear_ready", **base)


def _target_map_worker(
    payload: tuple[str, float, float]
) -> tuple[list[dict[str, object]], dict[str, object]]:
    label, x, y = payload
    sensor = Point(x, y)
    radii = inclusive_grid(
        FINAL.target_map_r_start, FINAL.target_map_r_stop, FINAL.target_map_r_step
    )
    alphas = inclusive_grid(
        FINAL.target_map_alpha_start,
        FINAL.target_map_alpha_stop,
        FINAL.target_map_alpha_step,
    )
    epsilons = inclusive_grid(-1.0, 1.0, FINAL.target_map_epsilon_step)
    rows = [
        dict(candidate=label, sensor_x=x, sensor_y=y, **classify_target(sensor, r, a, epsilons))
        for r in radii
        for a in alphas
    ]
    total_area_weight = sum(_number(row, "r") for row in rows)
    total_uniform = len(rows)
    summary: dict[str, object] = {
        "candidate": label,
        "x": x,
        "y": y,
        "movement_distance": math.hypot(x, y),
        "target_count": len(rows),
        "epsilon_count": len(epsilons),
    }
    for reason in FAILURE_REASONS:
        matching = [row for row in rows if row["failure_reason"] == reason]
        summary[f"area_fraction_{reason}"] = (
            sum(_number(row, "r") for row in matching) / total_area_weight
        )
        summary[f"uniform_fraction_{reason}"] = len(matching) / total_uniform
    diameter_values = [
        (_number(row, "worst_diameter"), _number(row, "r"))
        for row in rows
        if math.isfinite(_number(row, "worst_diameter"))
    ]
    summary.update(
        {
            "guaranteed_receive_area_coverage": sum(
                _number(row, "r") for row in rows if row["guaranteed_rx"]
            )
            / total_area_weight,
            "robust_clear_ready_area_coverage": summary[
                "area_fraction_robust_clear_ready"
            ],
            "operational_good_area_coverage": sum(
                _number(row, "r") for row in rows if row["operational_good"]
            )
            / total_area_weight,
            "robust_clear_ready_uniform_coverage": summary[
                "uniform_fraction_robust_clear_ready"
            ],
            "worst_diameter_p50": weighted_percentile(diameter_values, 0.50),
            "worst_diameter_mean": weighted_mean(diameter_values),
            "worst_diameter_p90": weighted_percentile(diameter_values, 0.90),
            "worst_diameter_p95": weighted_percentile(diameter_values, 0.95),
            "worst_diameter_max": max((value for value, _ in diameter_values), default=math.nan),
            "unbounded_fraction_of_target_area": summary[
                "area_fraction_localization_unbounded"
            ],
            "circle_failure_fraction_of_target_area": summary[
                "area_fraction_diameter_circle_failure"
            ],
        }
    )
    return rows, summary


def run_target_maps(
    candidates: Sequence[tuple[str, Point]],
    output_dir: Path,
    workers: int,
    resume: bool,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    map_path = output_dir / "q2_target_space_map.csv"
    summary_path = output_dir / "q2_target_space_summary.csv"
    payloads = [(label, point.x, point.y) for label, point in candidates]
    signature = _signature("target_maps", payloads)
    if _cache_hit(summary_path, signature, resume) and map_path.exists():
        print("cache hit: target maps", flush=True)
        return [dict(row) for row in read_csv(map_path)], [dict(row) for row in read_csv(summary_path)]
    started = time.perf_counter()
    if workers == 1:
        iterator: Iterable[tuple[list[dict[str, object]], dict[str, object]]] = map(
            _target_map_worker, payloads
        )
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=min(workers, len(payloads)))
        iterator = executor.map(_target_map_worker, payloads, chunksize=1)
    map_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    try:
        for index, (rows, summary) in enumerate(iterator, 1):
            map_rows.extend(rows)
            summaries.append(summary)
            print(f"target maps: {index}/{len(payloads)}", flush=True)
    finally:
        if executor is not None:
            executor.shutdown()
    write_csv(map_path, map_rows)
    write_csv(summary_path, summaries)
    _write_metadata(summary_path, signature, started, len(summaries))
    return map_rows, summaries


def _point_from_row(row: dict[str, object]) -> Point:
    return Point(_number(row, "x"), _number(row, "y"))


def _dedupe_labeled(points: Sequence[tuple[str, Point]]) -> list[tuple[str, Point]]:
    result: list[tuple[str, Point]] = []
    seen: set[tuple[float, float]] = set()
    for label, point in points:
        key = (point.x, point.y)
        if key not in seen:
            result.append((label, point))
            seen.add(key)
    return result


def select_candidates(
    local_rows: Sequence[dict[str, object]],
) -> tuple[list[tuple[str, Point]], dict[str, dict[str, object]]]:
    best = top_rows(local_rows, 1)[0]
    best_uniform = max(
        local_rows, key=lambda row: _number(row, "uniform_clear_ready_coverage")
    )
    max_coverage = _number(best, "clear_ready_coverage")
    short = min(
        (
            row
            for row in local_rows
            if _number(row, "clear_ready_coverage") >= 0.90 * max_coverage
        ),
        key=lambda row: _number(row, "movement_distance"),
    )
    middle = min(
        (
            row
            for row in local_rows
            if _number(row, "clear_ready_coverage") >= 0.97 * max_coverage
        ),
        key=lambda row: _number(row, "movement_distance"),
    )
    rows = {
        "local_area_best": best,
        "local_uniform_best": best_uniform,
        "short_90pct": short,
        "middle_97pct": middle,
    }
    candidates = _dedupe_labeled(
        [
            ("local_area_best", _point_from_row(best)),
            ("local_uniform_best", _point_from_row(best_uniform)),
            ("previous_validation_best", Point(690.0, 462.0)),
            ("previous_optimization_best", Point(690.0, 468.0)),
            ("short_90pct", _point_from_row(short)),
            ("middle_97pct", _point_from_row(middle)),
        ]
    )
    return candidates, rows


def create_sensitivity_data(
    local_rows: Sequence[dict[str, object]], output_dir: Path
) -> list[dict[str, object]]:
    area_sorted = sorted(
        local_rows, key=lambda row: -_number(row, "clear_ready_coverage")
    )
    uniform_sorted = sorted(
        local_rows, key=lambda row: -_number(row, "uniform_clear_ready_coverage")
    )
    union = {
        (_number(row, "x"), _number(row, "y")): row
        for row in (*area_sorted[:100], *uniform_sorted[:100])
    }
    area_rank = {
        (_number(row, "x"), _number(row, "y")): rank
        for rank, row in enumerate(area_sorted, 1)
    }
    uniform_rank = {
        (_number(row, "x"), _number(row, "y")): rank
        for rank, row in enumerate(uniform_sorted, 1)
    }
    rows = []
    for key, row in union.items():
        rows.append(
            {
                "x": key[0],
                "y": key[1],
                "movement_distance": _number(row, "movement_distance"),
                "area_weighted_coverage": _number(row, "clear_ready_coverage"),
                "uniform_r_alpha_coverage": _number(
                    row, "uniform_clear_ready_coverage"
                ),
                "area_rank": area_rank[key],
                "uniform_rank": uniform_rank[key],
            }
        )
    rows.sort(key=lambda row: min(int(row["area_rank"]), int(row["uniform_rank"])))
    write_csv(output_dir / "q2_weighting_sensitivity.csv", rows)
    return rows


def create_final_pareto(
    local_rows: Sequence[dict[str, object]], output_dir: Path
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    existing_paths = (
        "q2_coarse_scan.csv",
        "q2_refine_50m.csv",
        "q2_refine_10m.csv",
        "q2_refine_2m.csv",
    )
    old_stages = [
        [dict(row) for row in read_csv(output_dir / name)]
        for name in existing_paths
        if (output_dir / name).exists()
    ]
    local_metric_rows = [
        {
            key: value
            for key, value in row.items()
            if key
            not in {
                "area_rank",
                "uniform_rank",
                "in_area_99pct_region",
                "in_uniform_99pct_region",
            }
        }
        for row in local_rows
    ]
    combined = merge_finest((*old_stages, local_metric_rows))
    top = top_rows(combined, 100)
    write_csv(
        output_dir / "q2_top_candidates.csv",
        [dict(rank=index, **row) for index, row in enumerate(top, 1)],
    )
    pareto = pareto_front(combined)
    write_csv(
        output_dir / "q2_pareto_candidates.csv",
        [dict(pareto_index=index, **row) for index, row in enumerate(pareto, 1)],
    )
    return combined, pareto


def create_figures(
    output_dir: Path,
    orthogonal_summary: Sequence[dict[str, object]],
    local_rows: Sequence[dict[str, object]],
    local_summary: dict[str, object],
    cross_rows: Sequence[dict[str, object]],
    ultra_rows: Sequence[dict[str, object]],
    combined: Sequence[dict[str, object]],
    pareto: Sequence[dict[str, object]],
    target_rows: Sequence[dict[str, object]],
    recommended_label: str,
    recommended_point: Point,
) -> None:
    matplotlib_cache = Path(tempfile.gettempdir()) / "mathb-matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.lines import Line2D

    d1s = list(FINAL.orthogonal_d1)
    d2s = list(FINAL.orthogonal_d2)
    lookup = {
        (_number(row, "d1"), _number(row, "d2")): row for row in orthogonal_summary
    }
    deviation = [
        [_number(lookup[(d1, d2)], "deviation_from_90_deg") for d2 in d2s]
        for d1 in d1s
    ]
    improvement = [
        [100.0 * _number(lookup[(d1, d2)], "relative_improvement_from_90") for d2 in d2s]
        for d1 in d1s
    ]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for axis, matrix, title, color_label, cmap in (
        (axes[0], deviation, "Exact minimax correction", "optimal beta - 90 deg", "coolwarm"),
        (axes[1], improvement, "Gain versus beta=90 deg", "worst-diameter reduction (%)", "viridis"),
    ):
        image = axis.imshow(matrix, origin="lower", aspect="auto", cmap=cmap)
        axis.set_xticks(range(len(d2s)), [f"{value:.0f}" for value in d2s])
        axis.set_yticks(range(len(d1s)), [f"{value:.0f}" for value in d1s])
        axis.set_xlabel("Second range d2 (m)")
        axis.set_ylabel("First range d1 (m)")
        axis.set_title(title)
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                axis.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)
        figure.colorbar(image, ax=axis, label=color_label, shrink=0.85)
    figure.savefig(output_dir / "fig_q2_orthogonal_correction.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    cross_best = max(_number(row, "clear_ready_coverage") for row in cross_rows)
    stable = [
        row
        for row in cross_rows
        if _number(row, "clear_ready_coverage") >= 0.995 * cross_best
    ]
    scatter = axes[0].scatter(
        [_number(row, "x") for row in cross_rows],
        [_number(row, "y") for row in cross_rows],
        c=[_number(row, "clear_ready_coverage") for row in cross_rows],
        cmap="viridis",
        s=24,
    )
    axes[0].scatter(
        [_number(row, "x") for row in stable],
        [_number(row, "y") for row in stable],
        color="#d73027",
        s=10,
        alpha=0.8,
        label=">=99.5% on cross grid",
    )
    axes[0].scatter(
        [_number(row, "x") for row in ultra_rows],
        [_number(row, "y") for row in ultra_rows],
        facecolors="none",
        edgecolors="#111111",
        marker="s",
        s=45,
        label="ultra-validated points",
    )
    axes[0].scatter(
        [recommended_point.x],
        [recommended_point.y],
        marker="*",
        color="#fdae61",
        edgecolors="#111111",
        s=150,
        label="recommended point",
    )
    axes[0].set_title("Area-weighted cross-grid refinement")
    axes[0].legend(loc="lower right", fontsize=7)
    figure.colorbar(scatter, ax=axes[0], label="area-weighted coverage")

    sensitivity_region = [
        row
        for row in combined
        if 400.0 <= _number(row, "x") <= 800.0
        and 250.0 <= _number(row, "y") <= 650.0
    ]
    uniform_best = max(
        combined, key=lambda row: _number(row, "uniform_clear_ready_coverage")
    )
    scatter = axes[1].scatter(
        [_number(row, "x") for row in sensitivity_region],
        [_number(row, "y") for row in sensitivity_region],
        c=[_number(row, "uniform_clear_ready_coverage") for row in sensitivity_region],
        cmap="viridis",
        s=8,
    )
    axes[1].scatter(
        [_number(uniform_best, "x")],
        [_number(uniform_best, "y")],
        marker="D",
        color="#d73027",
        edgecolors="#111111",
        s=70,
        label="uniform-grid best",
    )
    axes[1].scatter(
        [recommended_point.x],
        [recommended_point.y],
        marker="*",
        color="#fdae61",
        edgecolors="#111111",
        s=150,
        label="area-weighted recommendation",
    )
    axes[1].set_title("Uniform (r, alpha) sensitivity")
    axes[1].legend(loc="lower right", fontsize=8)
    figure.colorbar(scatter, ax=axes[1], label="uniform-grid coverage")
    for axis in axes:
        axis.set_xlabel("S2 x (m)")
        axis.set_ylabel("S2 y (m)")
        axis.set_aspect("equal", adjustable="box")
    figure.savefig(output_dir / "fig_q2_candidate_region_final.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    axis.scatter(
        [_number(row, "movement_distance") for row in combined],
        [_number(row, "clear_ready_coverage") for row in combined],
        color="#bdbdbd",
        s=8,
        alpha=0.35,
        label="evaluated candidates",
    )
    axis.scatter(
        [_number(row, "movement_distance") for row in pareto],
        [_number(row, "clear_ready_coverage") for row in pareto],
        facecolors="none",
        edgecolors="#2166ac",
        s=22,
        linewidths=0.7,
        label="four-metric Pareto set",
    )
    axis.set_xlabel("Movement distance from S1 (m)")
    axis.set_ylabel("Area-weighted robust clear-ready coverage")
    axis.set_title("Final coverage-distance Pareto view")
    axis.grid(alpha=0.2)
    axis.legend()
    figure.savefig(output_dir / "fig_q2_pareto_final.png", dpi=180)
    plt.close(figure)

    selected = [row for row in target_rows if row["candidate"] == recommended_label]
    radii = sorted({_number(row, "r") for row in selected})
    alphas = sorted({_number(row, "alpha_deg") for row in selected})
    by_key = {
        (_number(row, "r"), _number(row, "alpha_deg")): row for row in selected
    }
    reason_index = {reason: index for index, reason in enumerate(FAILURE_REASONS)}
    robust = [
        [int(str(by_key[(r, a)]["robust_clear_ready"]).lower() == "true") for r in radii]
        for a in alphas
    ]
    diameter = [
        [
            min(_number(by_key[(r, a)], "worst_diameter"), 200.0)
            if math.isfinite(_number(by_key[(r, a)], "worst_diameter"))
            else math.nan
            for r in radii
        ]
        for a in alphas
    ]
    guaranteed = [
        [int(str(by_key[(r, a)]["guaranteed_rx"]).lower() == "true") for r in radii]
        for a in alphas
    ]
    reasons = [
        [reason_index[str(by_key[(r, a)]["failure_reason"])] for r in radii]
        for a in alphas
    ]
    extent = (
        radii[0] - FINAL.target_map_r_step / 2,
        radii[-1] + FINAL.target_map_r_step / 2,
        alphas[0] - FINAL.target_map_alpha_step / 2,
        alphas[-1] + FINAL.target_map_alpha_step / 2,
    )
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    binary_cmap = ListedColormap(("#f0f0f0", "#1b9e77"))
    image = axes[0, 0].imshow(
        robust, origin="lower", extent=extent, aspect="auto", cmap=binary_cmap, vmin=0, vmax=1
    )
    axes[0, 0].set_title("Robust clear-ready")
    figure.colorbar(image, ax=axes[0, 0], ticks=(0, 1), label="0 fail / 1 pass")
    image = axes[0, 1].imshow(
        diameter, origin="lower", extent=extent, aspect="auto", cmap="magma", vmin=0, vmax=200
    )
    axes[0, 1].set_title("Worst diameter over epsilon")
    figure.colorbar(image, ax=axes[0, 1], label="diameter (m; clipped at 200)")
    image = axes[1, 0].imshow(
        guaranteed, origin="lower", extent=extent, aspect="auto", cmap=binary_cmap, vmin=0, vmax=1
    )
    axes[1, 0].set_title("Guaranteed reception")
    figure.colorbar(image, ax=axes[1, 0], ticks=(0, 1), label="0 no / 1 yes")
    reason_colors = (
        "#66c2a5",
        "#8da0cb",
        "#e78ac3",
        "#a6d854",
        "#fc8d62",
        "#ffd92f",
        "#1b9e77",
    )
    categorical = ListedColormap(reason_colors)
    norm = BoundaryNorm(range(len(FAILURE_REASONS) + 1), categorical.N)
    axes[1, 1].imshow(
        reasons, origin="lower", extent=extent, aspect="auto", cmap=categorical, norm=norm
    )
    axes[1, 1].set_title("Dominant robust outcome")
    handles = [
        Line2D((0,), (0,), marker="s", linestyle="", color=color, label=reason.replace("_", " "))
        for reason, color in zip(FAILURE_REASONS, reason_colors)
    ]
    axes[1, 1].legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8)
    for axis in axes.flat:
        axis.set_xlabel("Target radius r (m)")
        axis.set_ylabel("Target angle alpha (deg)")
    figure.suptitle(f"Target-space structure for {recommended_label}")
    figure.savefig(output_dir / "fig_q2_target_space_failure.png", dpi=180)
    plt.close(figure)


def run_all(output_dir: Path, workers: int, resume: bool) -> dict[str, object]:
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    orthogonal_curves, orthogonal_summary = run_orthogonal_correction(
        output_dir, workers, resume
    )
    local_rows, local_summary = run_local_refinement(output_dir, workers, resume)
    combined, pareto = create_final_pareto(local_rows, output_dir)
    selected, selected_rows = select_candidates(combined)
    convergence = run_epsilon_convergence(selected, output_dir, workers, resume)
    cross_rows, ultra_rows, cross_summary = run_cross_grid_validation(
        output_dir, workers, resume
    )
    recommended_convergence = top_rows(ultra_rows, 1)[0]
    recommended_label = "ultra_cross_validation_best"
    recommended_point = Point(
        _number(recommended_convergence, "x"), _number(recommended_convergence, "y")
    )
    global_uniform_best = max(
        combined, key=lambda row: _number(row, "uniform_clear_ready_coverage")
    )
    selected_for_maps = _dedupe_labeled(
        [
            ("recommended", recommended_point),
            ("previous_optimization_best", Point(690.0, 468.0)),
            ("short_90pct", _point_from_row(selected_rows["short_90pct"])),
            ("middle_97pct", _point_from_row(selected_rows["middle_97pct"])),
            ("uniform_weight_best", _point_from_row(global_uniform_best)),
            ("orthogonal_baseline", Point(806.9809625749134, 399.9390780625565)),
        ]
    )
    target_rows, target_summaries = run_target_maps(
        selected_for_maps, output_dir, workers, resume
    )
    write_csv(output_dir / "q2_baseline_comparison.csv", target_summaries)
    sensitivity = create_sensitivity_data(combined, output_dir)
    create_figures(
        output_dir,
        orthogonal_summary,
        local_rows,
        local_summary,
        cross_rows,
        ultra_rows,
        combined,
        pareto,
        target_rows,
        "recommended",
        recommended_point,
    )
    summary: dict[str, object] = {
        "configuration": asdict(FINAL),
        "orthogonal_pair_count": len(orthogonal_summary),
        "orthogonal_summary": orthogonal_summary,
        "local_refinement": local_summary,
        "selected_candidates": {
            label: {"x": point.x, "y": point.y} for label, point in selected
        },
        "global_uniform_weight_best": global_uniform_best,
        "final_recommended": {
            "source_label": recommended_label,
            "x": recommended_point.x,
            "y": recommended_point.y,
            "metrics": recommended_convergence,
        },
        "epsilon_convergence": convergence,
        "cross_grid_validation": cross_summary,
        "cross_grid_candidate_count": len(cross_rows),
        "ultra_cross_validation": ultra_rows,
        "target_space_summary": target_summaries,
        "sensitivity_row_count": len(sensitivity),
        "combined_candidate_count": len(combined),
        "pareto_candidate_count": len(pareto),
        "elapsed_seconds": time.perf_counter() - started,
    }
    (output_dir / "q2_finalization_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, allow_nan=True), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("results/q2"))
    parser.add_argument(
        "--workers", type=int, default=max(1, min(8, os.cpu_count() or 1))
    )
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("--workers must be positive")
    run_all(args.output_dir, args.workers, not args.no_resume)


if __name__ == "__main__":
    main()
