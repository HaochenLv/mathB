"""Two-dimensional geometry primitives used by all problem stages."""

from .coverage import CoverageResult, DiameterCircleCoverage, diameter_circle_coverage
from .diameter import DiameterResult, polygon_diameter
from .localization import (
    EPS,
    LocalizationRegion,
    Measurement,
    Point,
    RegionStatus,
    construct_bounded_polygon,
    point_in_convex_polygon,
    point_satisfies_measurement,
)

__all__ = [
    "EPS",
    "CoverageResult",
    "DiameterCircleCoverage",
    "DiameterResult",
    "LocalizationRegion",
    "Measurement",
    "Point",
    "RegionStatus",
    "construct_bounded_polygon",
    "diameter_circle_coverage",
    "point_in_convex_polygon",
    "point_satisfies_measurement",
    "polygon_diameter",
]

