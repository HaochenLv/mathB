import math

import pytest

from src.geometry.localization import (
    LocalizationRegion,
    Measurement,
    Point,
    RegionStatus,
    point_in_convex_polygon,
    point_satisfies_measurement,
)


DEMO = (
    Measurement(0.0, 0.0, 45.4),
    Measurement(1000.0, 0.0, 134.6),
    Measurement(0.0, 1000.0, 315.7),
)


def test_angle_wraparound_uses_vector_halfplanes() -> None:
    point = Point(100.0, 0.0)
    assert point_satisfies_measurement(point, Measurement(0, 0, 0.2))
    assert point_satisfies_measurement(point, Measurement(0, 0, 359.8))
    assert not point_satisfies_measurement(Point(-100, 0), Measurement(0, 0, 359.8))


def test_single_and_nearly_parallel_measurements_are_unbounded() -> None:
    region = LocalizationRegion()
    assert region.add_measurement(0, 0, 0.0) is RegionStatus.UNBOUNDED
    assert region.add_measurement(0, 10, 0.001) is RegionStatus.UNBOUNDED
    assert region.vertices == ()


def test_near_tangent_regression_is_unbounded() -> None:
    region = LocalizationRegion()
    region.add_measurement(-954.0845491620161, 513.6663198127981, 135.69)
    region.add_measurement(937.646961779461, 1334.822136238413, 133.70)

    # The angular intervals [134.69, 136.69] and [132.70, 134.70]
    # overlap with strictly positive width, so they share a recession direction.
    assert region.status is RegionStatus.UNBOUNDED
    assert region.vertices == ()


@pytest.mark.parametrize("overlap_deg", [0.1, 0.01, 0.001])
def test_near_tangent_positive_overlaps_remain_unbounded(
    overlap_deg: float,
) -> None:
    first_theta = 135.69
    second_theta = first_theta - (2.0 - overlap_deg)
    region = LocalizationRegion()
    region.add_measurement(-954.0845491620161, 513.6663198127981, first_theta)
    region.add_measurement(937.646961779461, 1334.822136238413, second_theta)

    assert region.status is RegionStatus.UNBOUNDED
    assert region.vertices == ()


def test_near_tangent_disjoint_directions_remain_empty() -> None:
    first_theta = 135.69
    region = LocalizationRegion()
    region.add_measurement(-954.0845491620161, 513.6663198127981, first_theta)
    region.add_measurement(
        937.646961779461,
        1334.822136238413,
        first_theta - 2.001,
    )

    assert region.status is RegionStatus.EMPTY
    assert region.vertices == ()


def test_contradictory_measurements_are_empty_without_crash() -> None:
    region = LocalizationRegion()
    region.add_measurement(0, 0, 0.0)
    assert region.add_measurement(10, 10, 180.0) is RegionStatus.EMPTY
    assert region.vertices == ()


def test_manual_localization_contains_target_and_shrinks() -> None:
    target = Point(500.0, 500.0)
    region = LocalizationRegion()
    previous_vertices = None
    previous_area = None
    previous_diameter = None
    bounded_count = 0

    for measurement in DEMO:
        region.add_measurement(measurement.x, measurement.y, measurement.theta_deg)
        assert region.contains(target)
        if region.status is RegionStatus.BOUNDED:
            bounded_count += 1
            assert all(
                point_satisfies_measurement(
                    vertex,
                    measurement,
                    delta_deg=region.delta_deg,
                    eps=1e-7,
                )
                for vertex in region.vertices
                for measurement in region.measurements
            )
            if previous_vertices is not None:
                assert all(
                    point_in_convex_polygon(vertex, previous_vertices, 1e-7)
                    for vertex in region.vertices
                )
                assert region.area() <= previous_area + 1e-7
                assert region.diameter().length <= previous_diameter + 1e-7
            previous_vertices = region.vertices
            previous_area = region.area()
            previous_diameter = region.diameter().length

    assert bounded_count == 2
    assert len(region.vertices) >= 3
    assert region.area() > 0.0


def test_redundant_measurement_leaves_region_unchanged() -> None:
    region = LocalizationRegion()
    for measurement in DEMO[:2]:
        region.add_measurement(measurement.x, measurement.y, measurement.theta_deg)
    before = region.vertices
    before_area = region.area()
    region.add_measurement(DEMO[0].x, DEMO[0].y, DEMO[0].theta_deg)
    assert region.status is RegionStatus.BOUNDED
    assert region.area() == pytest.approx(before_area, abs=1e-7)
    assert len(region.vertices) == len(before)


def test_small_scale_region_is_finite_and_stable() -> None:
    target = Point(5e-4, 5e-4)
    sensors = (Point(0, 0), Point(1e-3, 0), Point(0, 1e-3))
    region = LocalizationRegion(delta_deg=0.1, eps=1e-12)
    for sensor in sensors:
        bearing = math.degrees(math.atan2(target.y - sensor.y, target.x - sensor.x))
        region.add_measurement(sensor.x, sensor.y, bearing)
    assert region.status is RegionStatus.BOUNDED
    assert region.contains(target)
    assert math.isfinite(region.area())
    assert math.isfinite(region.diameter().length)


def test_point_on_wedge_boundary_is_included() -> None:
    angle = math.radians(1.0)
    boundary_point = Point(100.0 * math.cos(angle), 100.0 * math.sin(angle))
    assert point_satisfies_measurement(boundary_point, Measurement(0, 0, 0.0))
