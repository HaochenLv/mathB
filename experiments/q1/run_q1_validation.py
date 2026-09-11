"""Lightweight reproducible Monte Carlo validation for Q1."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import random
import sys
from dataclasses import dataclass

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.geometry.localization import (
    EPS,
    LocalizationRegion,
    Point,
    RegionStatus,
    point_in_convex_polygon,
    point_satisfies_measurement,
)


@dataclass
class ValidationSummary:
    total: int = 0
    passed: int = 0
    failed: int = 0
    unbounded_or_skipped: int = 0


def _bearing(origin: Point, target: Point) -> float:
    return math.degrees(math.atan2(target.y - origin.y, target.x - origin.x)) % 360.0


def run_validation(cases: int = 500, seed: int = 20260910) -> ValidationSummary:
    """Validate containment and monotonicity on reproducible random cases."""

    randomizer = random.Random(seed)
    summary = ValidationSummary(total=cases)
    for case_index in range(cases):
        target = Point(
            randomizer.uniform(-800.0, 800.0),
            randomizer.uniform(-800.0, 800.0),
        )
        count = randomizer.randint(3, 7)
        base_angle = randomizer.uniform(0.0, 360.0)
        sensors: list[Point] = []
        measurements: list[tuple[Point, float]] = []
        for index in range(count):
            angle = math.radians(base_angle + index * 360.0 / count + randomizer.uniform(-18.0, 18.0))
            distance = randomizer.uniform(700.0, 1500.0)
            sensor = Point(
                target.x + distance * math.cos(angle),
                target.y + distance * math.sin(angle),
            )
            theta = (_bearing(sensor, target) + randomizer.uniform(-1.0, 1.0)) % 360.0
            sensors.append(sensor)
            measurements.append((sensor, theta))

        region = LocalizationRegion(delta_deg=1.0)
        previous_vertices: tuple[Point, ...] | None = None
        previous_area: float | None = None
        previous_diameter: float | None = None
        failure = ""
        for step, (sensor, theta) in enumerate(measurements, 1):
            region.add_measurement(sensor.x, sensor.y, theta)
            target_is_feasible = all(
                point_satisfies_measurement(
                    target,
                    measurement,
                    delta_deg=region.delta_deg,
                    eps=region.eps,
                )
                for measurement in region.measurements
            )
            if not target_is_feasible:
                failure = f"true target rejected at step {step}"
                break
            if region.status is RegionStatus.EMPTY:
                failure = f"known feasible target classified EMPTY at step {step}"
                break
            if region.status is RegionStatus.BOUNDED:
                area = region.area()
                diameter = region.diameter().length
                if not all(
                    point_satisfies_measurement(
                        vertex,
                        measurement,
                        delta_deg=region.delta_deg,
                        eps=100.0 * EPS,
                    )
                    for vertex in region.vertices
                    for measurement in region.measurements
                ):
                    failure = f"bounded vertex violates a measurement at step {step}"
                    break
                if previous_vertices is not None and not all(
                    point_in_convex_polygon(vertex, previous_vertices, 100.0 * EPS)
                    for vertex in region.vertices
                ):
                    failure = f"set containment failed at step {step}"
                    break
                if previous_area is not None and area > previous_area + 1000.0 * EPS:
                    failure = f"area increased at step {step}"
                    break
                if previous_diameter is not None and diameter > previous_diameter + 1000.0 * EPS:
                    failure = f"diameter increased at step {step}"
                    break
                previous_vertices = region.vertices
                previous_area = area
                previous_diameter = diameter

        if failure:
            summary.failed += 1
            print(f"FAIL seed={seed} case={case_index} reason={failure}")
            print(f"  target={target}")
            print(f"  sensors={sensors}")
            print(f"  measurements={measurements}")
        elif previous_vertices is None:
            summary.unbounded_or_skipped += 1
        else:
            summary.passed += 1
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260910)
    arguments = parser.parse_args()
    summary = run_validation(arguments.cases, arguments.seed)
    print(f"seed: {arguments.seed}")
    print(f"total: {summary.total}")
    print(f"passed: {summary.passed}")
    print(f"failed: {summary.failed}")
    print(f"unbounded/skipped: {summary.unbounded_or_skipped}")
    raise SystemExit(1 if summary.failed else 0)


if __name__ == "__main__":
    main()
