"""Diameter computation for bounded convex polygons."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

from .localization import EPS, Point


@dataclass(frozen=True, slots=True)
class DiameterResult:
    """The polygon diameter and every vertex pair attaining it."""

    length: float
    point_a: Point
    point_b: Point
    diametral_pairs: tuple[tuple[Point, Point], ...]


def polygon_diameter(
    vertices: Sequence[Point], eps: float = EPS
) -> DiameterResult:
    """Return the convex polygon diameter by enumerating all vertex pairs."""

    if not vertices:
        raise ValueError("diameter is undefined for an empty polygon")
    if len(vertices) == 1:
        pair = (vertices[0], vertices[0])
        return DiameterResult(0.0, *pair, (pair,))

    best_squared = -1.0
    best_pairs: list[tuple[Point, Point]] = []
    for first_index, first in enumerate(vertices):
        for second in vertices[first_index + 1 :]:
            dx = first.x - second.x
            dy = first.y - second.y
            squared = dx * dx + dy * dy
            tolerance = eps * max(1.0, best_squared, squared)
            if squared > best_squared + tolerance:
                best_squared = squared
                best_pairs = [(first, second)]
            elif abs(squared - best_squared) <= tolerance:
                best_pairs.append((first, second))
    point_a, point_b = best_pairs[0]
    return DiameterResult(
        math.sqrt(max(0.0, best_squared)),
        point_a,
        point_b,
        tuple(best_pairs),
    )

