import math

from src.geometry.coverage import diameter_circle_coverage
from src.geometry.localization import Point


def test_every_square_diagonal_circle_covers_square() -> None:
    vertices = (Point(0, 0), Point(2, 0), Point(2, 2), Point(0, 2))
    result = diameter_circle_coverage(vertices)
    assert result.exists_covering_diameter_circle
    assert result.all_diameter_circles_cover
    assert len(result.circles) == 2


def test_no_equilateral_triangle_side_circle_covers_triangle() -> None:
    vertices = (Point(0, 0), Point(2, 0), Point(1, math.sqrt(3.0)))
    result = diameter_circle_coverage(vertices)
    assert len(result.circles) == 3
    assert not result.exists_covering_diameter_circle
    assert not result.all_diameter_circles_cover
    assert all(not circle.covered for circle in result.circles)

