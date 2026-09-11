"""Localization by intersecting deterministic bearing-error wedges.

The public angle convention is degrees counter-clockwise from the positive
x-axis.  Each measurement is represented by two normalized half-planes; no
slope form or artificial bounding box is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
import math
from typing import Iterable, Sequence


EPS = 1e-9


@dataclass(frozen=True, slots=True)
class Point:
    """A point or vector in the Euclidean plane."""

    x: float
    y: float

    def __add__(self, other: "Point") -> "Point":
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Point") -> "Point":
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: float) -> "Point":
        return Point(self.x * scalar, self.y * scalar)

    __rmul__ = __mul__


@dataclass(frozen=True, slots=True)
class Measurement:
    """A detector position and its measured bearing, in degrees."""

    x: float
    y: float
    theta_deg: float

    @property
    def sensor(self) -> Point:
        return Point(self.x, self.y)


class RegionStatus(str, Enum):
    """Possible states of a half-plane intersection."""

    UNBOUNDED = "UNBOUNDED"
    BOUNDED = "BOUNDED"
    EMPTY = "EMPTY"


@dataclass(frozen=True, slots=True)
class _HalfPlane:
    """Normalized inequality ``a*x + b*y <= c``."""

    a: float
    b: float
    c: float

    def value(self, point: Point) -> float:
        return self.a * point.x + self.b * point.y - self.c

    def contains(self, point: Point, eps: float = EPS) -> bool:
        return self.value(point) <= eps


@dataclass(frozen=True, slots=True)
class PolygonConstruction:
    """Result returned while constructing an intersection from scratch."""

    status: RegionStatus
    vertices: tuple[Point, ...]


def cross(a: Point, b: Point) -> float:
    """Return the scalar two-dimensional cross product."""

    return a.x * b.y - a.y * b.x


def _direction(angle_deg: float) -> Point:
    angle_rad = math.radians(angle_deg)
    return Point(math.cos(angle_rad), math.sin(angle_rad))


def _measurement_halfplanes(
    measurement: Measurement, delta_deg: float
) -> tuple[_HalfPlane, _HalfPlane]:
    low = _direction(measurement.theta_deg - delta_deg)
    high = _direction(measurement.theta_deg + delta_deg)
    sensor = measurement.sensor

    # cross(low, P-S) >= 0
    first_a, first_b = low.y, -low.x
    # cross(high, P-S) <= 0
    second_a, second_b = -high.y, high.x
    return (
        _HalfPlane(
            first_a,
            first_b,
            first_a * sensor.x + first_b * sensor.y,
        ),
        _HalfPlane(
            second_a,
            second_b,
            second_a * sensor.x + second_b * sensor.y,
        ),
    )


def point_satisfies_measurement(
    point: Point,
    measurement: Measurement,
    delta_deg: float = 1.0,
    eps: float = EPS,
) -> bool:
    """Return whether ``point`` lies in a measurement's forward wedge."""

    return all(
        halfplane.contains(point, eps)
        for halfplane in _measurement_halfplanes(measurement, delta_deg)
    )


def _line_intersection(
    first: _HalfPlane, second: _HalfPlane, eps: float
) -> Point | None:
    determinant = first.a * second.b - second.a * first.b
    if abs(determinant) <= eps:
        return None
    return Point(
        (first.c * second.b - second.c * first.b) / determinant,
        (first.a * second.c - second.a * first.c) / determinant,
    )


def _exact_feasible_intersections(
    halfplanes: Sequence[_HalfPlane], eps: float
) -> list[Point]:
    """Recover feasible vertices when floating-point cancellation hides them.

    This slow path is used only before declaring a system empty.  It treats the
    already-computed binary floating-point coefficients as exact rationals, so
    it removes intersection/evaluation roundoff without widening ``eps``.
    """

    exact_halfplanes = [
        tuple(Fraction.from_float(value) for value in (hp.a, hp.b, hp.c))
        for hp in halfplanes
    ]
    exact_eps = Fraction(eps)
    feasible: list[Point] = []
    for first_index, (first_a, first_b, first_c) in enumerate(exact_halfplanes):
        for second_a, second_b, second_c in exact_halfplanes[first_index + 1 :]:
            determinant = first_a * second_b - second_a * first_b
            if determinant == 0:
                continue
            x = (first_c * second_b - second_c * first_b) / determinant
            y = (first_a * second_c - second_a * first_c) / determinant
            if all(
                a * x + b * y - c <= exact_eps
                for a, b, c in exact_halfplanes
            ):
                point = Point(float(x), float(y))
                if math.isfinite(point.x) and math.isfinite(point.y):
                    feasible.append(point)
    return feasible


def _distance_squared(first: Point, second: Point) -> float:
    dx = first.x - second.x
    dy = first.y - second.y
    return dx * dx + dy * dy


def _deduplicate_points(points: Iterable[Point], eps: float) -> list[Point]:
    unique: list[Point] = []
    for point in points:
        if not any(_distance_squared(point, seen) <= eps * eps for seen in unique):
            unique.append(point)
    return unique


def _signed_area(vertices: Sequence[Point]) -> float:
    return 0.5 * sum(
        cross(vertices[index], vertices[(index + 1) % len(vertices)])
        for index in range(len(vertices))
    )


def _convex_hull(points: Iterable[Point], eps: float) -> list[Point]:
    ordered = sorted(_deduplicate_points(points, eps), key=lambda p: (p.x, p.y))
    if len(ordered) <= 1:
        return ordered

    def turn(origin: Point, first: Point, second: Point) -> float:
        return cross(first - origin, second - origin)

    lower: list[Point] = []
    for point in ordered:
        while len(lower) >= 2 and turn(lower[-2], lower[-1], point) <= eps:
            lower.pop()
        lower.append(point)

    upper: list[Point] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and turn(upper[-2], upper[-1], point) <= eps:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def _all_normals_parallel(halfplanes: Sequence[_HalfPlane], eps: float) -> bool:
    first = halfplanes[0]
    return all(abs(first.a * hp.b - first.b * hp.a) <= eps for hp in halfplanes[1:])


def _parallel_system_is_feasible(
    halfplanes: Sequence[_HalfPlane], eps: float
) -> bool:
    reference = halfplanes[0]
    lower = -math.inf
    upper = math.inf
    for halfplane in halfplanes:
        coefficient = halfplane.a * reference.a + halfplane.b * reference.b
        if coefficient > eps:
            upper = min(upper, halfplane.c / coefficient)
        elif coefficient < -eps:
            lower = max(lower, halfplane.c / coefficient)
        elif halfplane.c < -eps:
            return False
    return lower <= upper + eps


def _has_nonzero_recession_direction(
    halfplanes: Sequence[_HalfPlane], eps: float
) -> bool:
    """Test whether ``a*d.x+b*d.y <= 0`` has a nonzero solution."""

    boundaries: list[float] = []
    full_turn = 2.0 * math.pi
    for halfplane in halfplanes:
        normal_angle = math.atan2(halfplane.b, halfplane.a)
        boundaries.extend(
            ((normal_angle + math.pi / 2.0) % full_turn,
             (normal_angle - math.pi / 2.0) % full_turn)
        )
    boundaries.sort()

    candidates = list(boundaries)
    for index, angle in enumerate(boundaries):
        following = boundaries[(index + 1) % len(boundaries)]
        if index == len(boundaries) - 1:
            following += full_turn
        candidates.append(((angle + following) / 2.0) % full_turn)

    for angle in candidates:
        dx, dy = math.cos(angle), math.sin(angle)
        if all(hp.a * dx + hp.b * dy <= eps for hp in halfplanes):
            return True
    return False


def construct_bounded_polygon(
    observations: Sequence[Measurement],
    delta_deg: float = 1.0,
    eps: float = EPS,
) -> PolygonConstruction:
    """Construct and classify the exact wedge intersection.

    Boundary lines are intersected pairwise and feasible intersections become
    polygon candidates.  A separate homogeneous recession-cone test determines
    boundedness, so a finite artificial box is never introduced.
    """

    if not observations:
        return PolygonConstruction(RegionStatus.UNBOUNDED, ())
    halfplanes = [
        halfplane
        for measurement in observations
        for halfplane in _measurement_halfplanes(measurement, delta_deg)
    ]

    feasible: list[Point] = []
    for first_index, first in enumerate(halfplanes):
        for second in halfplanes[first_index + 1 :]:
            candidate = _line_intersection(first, second, eps)
            if candidate is not None and all(
                halfplane.contains(candidate, eps) for halfplane in halfplanes
            ):
                feasible.append(candidate)

    # These witnesses cover coincident/parallel systems that have no vertices.
    witnesses = [Point(0.0, 0.0), *(item.sensor for item in observations)]
    has_feasible_witness = any(
        all(halfplane.contains(point, eps) for halfplane in halfplanes)
        for point in witnesses
    )

    if not feasible and not has_feasible_witness:
        # A near-tangent pair can intersect very far from the sensors.  Cramer's
        # rule then loses enough absolute precision that a truly active
        # constraint may miss ``eps`` by a few ulps.  Recheck exact signs before
        # making the irreversible EMPTY classification.
        feasible = _exact_feasible_intersections(halfplanes, eps)

    if not feasible and not has_feasible_witness:
        if _all_normals_parallel(halfplanes, eps):
            if _parallel_system_is_feasible(halfplanes, eps):
                return PolygonConstruction(RegionStatus.UNBOUNDED, ())
        return PolygonConstruction(RegionStatus.EMPTY, ())

    if _has_nonzero_recession_direction(halfplanes, eps):
        return PolygonConstruction(RegionStatus.UNBOUNDED, ())

    hull = _convex_hull(feasible, eps)
    if not hull and has_feasible_witness:
        # This is only reachable for a bounded, zero-dimensional coincident
        # system; preserve a valid witness as the degenerate polygon.
        hull = [
            point
            for point in witnesses
            if all(halfplane.contains(point, eps) for halfplane in halfplanes)
        ][:1]
    return PolygonConstruction(RegionStatus.BOUNDED, tuple(hull))


def _segment_boundary_intersection(
    start: Point, end: Point, halfplane: _HalfPlane, eps: float
) -> Point:
    start_value = halfplane.value(start)
    end_value = halfplane.value(end)
    denominator = start_value - end_value
    if abs(denominator) <= eps:
        return start
    ratio = start_value / denominator
    return start + (end - start) * ratio


def _normalize_polygon(vertices: Iterable[Point], eps: float) -> list[Point]:
    points: list[Point] = []
    for point in vertices:
        if not points or _distance_squared(point, points[-1]) > eps * eps:
            points.append(point)
    if len(points) > 1 and _distance_squared(points[0], points[-1]) <= eps * eps:
        points.pop()
    if len(points) <= 2:
        return _deduplicate_points(points, eps)

    changed = True
    while changed and len(points) >= 3:
        changed = False
        cleaned: list[Point] = []
        count = len(points)
        for index, current in enumerate(points):
            previous = points[index - 1]
            following = points[(index + 1) % count]
            if abs(cross(current - previous, following - current)) <= eps:
                if (current.x - previous.x) * (following.x - current.x) + (
                    current.y - previous.y
                ) * (following.y - current.y) >= -eps:
                    changed = True
                    continue
            cleaned.append(current)
        points = cleaned
    if len(points) >= 3 and _signed_area(points) < 0.0:
        points.reverse()
    return points


def _clip_polygon_with_halfplane(
    vertices: Sequence[Point], halfplane: _HalfPlane, eps: float
) -> list[Point]:
    if not vertices:
        return []
    if len(vertices) == 1:
        return list(vertices) if halfplane.contains(vertices[0], eps) else []
    if len(vertices) == 2:
        start, end = vertices
        start_inside = halfplane.contains(start, eps)
        end_inside = halfplane.contains(end, eps)
        if start_inside and end_inside:
            return list(vertices)
        if not start_inside and not end_inside:
            return []
        crossing = _segment_boundary_intersection(start, end, halfplane, eps)
        return [start, crossing] if start_inside else [crossing, end]

    output: list[Point] = []
    previous = vertices[-1]
    previous_inside = halfplane.contains(previous, eps)
    for current in vertices:
        current_inside = halfplane.contains(current, eps)
        if previous_inside and current_inside:
            output.append(current)
        elif previous_inside and not current_inside:
            output.append(
                _segment_boundary_intersection(previous, current, halfplane, eps)
            )
        elif not previous_inside and current_inside:
            output.append(
                _segment_boundary_intersection(previous, current, halfplane, eps)
            )
            output.append(current)
        previous, previous_inside = current, current_inside
    return _normalize_polygon(output, eps)


def point_in_convex_polygon(
    point: Point, vertices: Sequence[Point], eps: float = EPS
) -> bool:
    """Return whether a point lies inside or on a CCW convex polygon."""

    if not vertices:
        return False
    if len(vertices) == 1:
        return _distance_squared(point, vertices[0]) <= eps * eps
    if len(vertices) == 2:
        segment = vertices[1] - vertices[0]
        relative = point - vertices[0]
        return abs(cross(segment, relative)) <= eps and (
            relative.x * (point.x - vertices[1].x)
            + relative.y * (point.y - vertices[1].y)
            <= eps
        )
    return all(
        cross(vertices[(index + 1) % len(vertices)] - vertices[index],
              point - vertices[index])
        >= -eps
        for index in range(len(vertices))
    )


class LocalizationRegion:
    """Incrementally intersect bearing wedges for one interference source."""

    def __init__(self, delta_deg: float = 1.0, eps: float = EPS) -> None:
        if not math.isfinite(delta_deg) or not 0.0 < delta_deg < 90.0:
            raise ValueError("delta_deg must be finite and in (0, 90)")
        if not math.isfinite(eps) or eps <= 0.0:
            raise ValueError("eps must be a positive finite number")
        self.delta_deg = delta_deg
        self.eps = eps
        self.measurements: list[Measurement] = []
        self.status = RegionStatus.UNBOUNDED
        self.vertices: tuple[Point, ...] = ()

    def add_measurement(self, x: float, y: float, theta_deg: float) -> RegionStatus:
        """Add one bearing observation and return the updated region status."""

        if not all(math.isfinite(value) for value in (x, y, theta_deg)):
            raise ValueError("measurement values must be finite")
        measurement = Measurement(float(x), float(y), float(theta_deg) % 360.0)
        self.measurements.append(measurement)

        if self.status is RegionStatus.EMPTY:
            return self.status
        if self.status is not RegionStatus.BOUNDED:
            construction = construct_bounded_polygon(
                self.measurements, self.delta_deg, self.eps
            )
            self.status = construction.status
            self.vertices = construction.vertices
            return self.status

        clipped = list(self.vertices)
        for halfplane in _measurement_halfplanes(measurement, self.delta_deg):
            clipped = _clip_polygon_with_halfplane(clipped, halfplane, self.eps)
            if not clipped:
                break
        self.vertices = tuple(_normalize_polygon(clipped, self.eps))
        self.status = RegionStatus.BOUNDED if self.vertices else RegionStatus.EMPTY
        return self.status

    def contains(self, point: Point) -> bool:
        """Return whether a point satisfies every measurement constraint."""

        return all(
            point_satisfies_measurement(
                point, measurement, self.delta_deg, self.eps
            )
            for measurement in self.measurements
        )

    def area(self) -> float:
        """Return the bounded polygon area."""

        self._require_bounded()
        if len(self.vertices) < 3:
            return 0.0
        return abs(_signed_area(self.vertices))

    def diameter(self):
        """Return all maximum-distance vertex pairs of the bounded region."""

        self._require_bounded()
        from .diameter import polygon_diameter

        return polygon_diameter(self.vertices, self.eps)

    def diameter_circle_coverage(self):
        """Evaluate the diameter circle associated with every diametral pair."""

        self._require_bounded()
        from .coverage import diameter_circle_coverage

        return diameter_circle_coverage(self.vertices, self.diameter(), self.eps)

    def _require_bounded(self) -> None:
        if self.status is not RegionStatus.BOUNDED:
            raise ValueError(f"operation requires BOUNDED region, got {self.status.value}")
