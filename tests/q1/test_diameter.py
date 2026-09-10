import math

import pytest

from src.geometry.diameter import polygon_diameter
from src.geometry.localization import Point


def test_three_four_five_rectangle_has_two_diameters() -> None:
    vertices = (Point(0, 0), Point(3, 0), Point(3, 4), Point(0, 4))
    result = polygon_diameter(vertices)
    assert result.length == pytest.approx(5.0)
    assert len(result.diametral_pairs) == 2


def test_square_diameter() -> None:
    vertices = (Point(0, 0), Point(2, 0), Point(2, 2), Point(0, 2))
    result = polygon_diameter(vertices)
    assert result.length == pytest.approx(2.0 * math.sqrt(2.0))
    assert len(result.diametral_pairs) == 2

