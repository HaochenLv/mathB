"""Coverage tests for circles whose diameters attain polygon diameter."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

from .diameter import DiameterResult, polygon_diameter
from .localization import EPS, Point


@dataclass(frozen=True, slots=True)
class CoverageResult:
    """Coverage result for one selected diameter pair."""

    covered: bool
    center: Point
    radius: float
    max_vertex_distance: float
    margin: float
    point_a: Point
    point_b: Point


@dataclass(frozen=True, slots=True)
class DiameterCircleCoverage:
    """Coverage results over every maximum-diameter vertex pair."""

    diameter: DiameterResult
    circles: tuple[CoverageResult, ...]
    exists_covering_diameter_circle: bool
    all_diameter_circles_cover: bool


def diameter_circle_coverage(
    vertices: Sequence[Point],
    diameter: DiameterResult | None = None,
    eps: float = EPS,
) -> DiameterCircleCoverage:
    """Check whether each diametral circle contains every polygon vertex."""

    if not vertices:
        raise ValueError("coverage is undefined for an empty polygon")
    diameter = diameter or polygon_diameter(vertices, eps)
    results: list[CoverageResult] = []
    for point_a, point_b in diameter.diametral_pairs:
        center = Point(
            (point_a.x + point_b.x) / 2.0,
            (point_a.y + point_b.y) / 2.0,
        )
        radius = math.hypot(point_a.x - point_b.x, point_a.y - point_b.y) / 2.0
        max_distance = max(
            math.hypot(vertex.x - center.x, vertex.y - center.y)
            for vertex in vertices
        )
        margin = radius - max_distance
        results.append(
            CoverageResult(
                covered=max_distance <= radius + eps,
                center=center,
                radius=radius,
                max_vertex_distance=max_distance,
                margin=margin,
                point_a=point_a,
                point_b=point_b,
            )
        )

    return DiameterCircleCoverage(
        diameter=diameter,
        circles=tuple(results),
        exists_covering_diameter_circle=any(item.covered for item in results),
        all_diameter_circles_cover=all(item.covered for item in results),
    )
