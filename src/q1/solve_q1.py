"""Run the deterministic Q1 localization demonstration.

Usage:
    python -m src.q1.solve_q1
"""

from __future__ import annotations

import math
import os
from pathlib import Path
import tempfile

from src.geometry.localization import LocalizationRegion, Measurement, Point, RegionStatus


DEMO_TARGET = Point(500.0, 500.0)
DEMO_MEASUREMENTS = (
    Measurement(0.0, 0.0, 45.4),
    Measurement(1000.0, 0.0, 134.6),
    Measurement(0.0, 1000.0, 315.7),
)


def solve_demo() -> tuple[LocalizationRegion, list[tuple[Point, ...]]]:
    """Solve the three-observation hand-checkable localization example."""

    region = LocalizationRegion(delta_deg=1.0)
    bounded_history: list[tuple[Point, ...]] = []
    for measurement in DEMO_MEASUREMENTS:
        region.add_measurement(
            measurement.x, measurement.y, measurement.theta_deg
        )
        if region.status is RegionStatus.BOUNDED:
            bounded_history.append(region.vertices)
    return region, bounded_history


def _ray_endpoint(measurement: Measurement, angle_deg: float, length: float) -> Point:
    angle_rad = math.radians(angle_deg)
    return Point(
        measurement.x + length * math.cos(angle_rad),
        measurement.y + length * math.sin(angle_rad),
    )


def create_demo_figure(
    region: LocalizationRegion,
    bounded_history: list[tuple[Point, ...]],
    output_path: Path,
) -> Path:
    """Render the measurements, localization polygons, and a diameter circle."""

    matplotlib_cache = Path(tempfile.gettempdir()) / "mathb-matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Polygon

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(figsize=(8, 8))
    colors = ("#1f77b4", "#ff7f0e", "#2ca02c")
    ray_length = 900.0

    for index, (measurement, color) in enumerate(zip(DEMO_MEASUREMENTS, colors), 1):
        axes.scatter(measurement.x, measurement.y, color=color, s=48, zorder=5)
        axes.annotate(
            f"S{index}",
            (measurement.x, measurement.y),
            xytext=(7, 7),
            textcoords="offset points",
        )
        for offset, style, alpha in ((0.0, "-", 0.85), (-1.0, "--", 0.55), (1.0, "--", 0.55)):
            endpoint = _ray_endpoint(
                measurement, measurement.theta_deg + offset, ray_length
            )
            axes.plot(
                (measurement.x, endpoint.x),
                (measurement.y, endpoint.y),
                linestyle=style,
                color=color,
                alpha=alpha,
                linewidth=1.2,
            )

    for index, vertices in enumerate(bounded_history):
        if len(vertices) >= 3:
            axes.add_patch(
                Polygon(
                    [(point.x, point.y) for point in vertices],
                    closed=True,
                    facecolor="#9ecae1" if index == 0 else "#e6550d",
                    edgecolor="#3182bd" if index == 0 else "#a63603",
                    alpha=0.20 if index == 0 else 0.32,
                    linewidth=1.8,
                    label="Initial bounded region" if index == 0 else "Final region",
                )
            )

    axes.scatter(DEMO_TARGET.x, DEMO_TARGET.y, marker="*", color="black", s=110, zorder=6)
    axes.annotate("G", (DEMO_TARGET.x, DEMO_TARGET.y), xytext=(7, 7), textcoords="offset points")

    diameter = region.diameter()
    axes.plot(
        (diameter.point_a.x, diameter.point_b.x),
        (diameter.point_a.y, diameter.point_b.y),
        color="#7a0177",
        linewidth=2.4,
        label="Maximum diameter",
    )
    coverage = region.diameter_circle_coverage().circles[0]
    axes.add_patch(
        Circle(
            (coverage.center.x, coverage.center.y),
            coverage.radius,
            fill=False,
            color="#7a0177",
            linestyle=":",
            linewidth=1.8,
            label="Diameter circle",
        )
    )
    axes.scatter(
        [point.x for point in region.vertices],
        [point.y for point in region.vertices],
        color="#a63603",
        s=20,
        zorder=7,
    )

    zoom_axes = axes.inset_axes([0.58, 0.08, 0.38, 0.38])
    zoom_axes.add_patch(
        Polygon(
            [(point.x, point.y) for point in region.vertices],
            closed=True,
            facecolor="#e6550d",
            edgecolor="#a63603",
            alpha=0.32,
            linewidth=1.8,
        )
    )
    zoom_axes.add_patch(
        Circle(
            (coverage.center.x, coverage.center.y),
            coverage.radius,
            fill=False,
            color="#7a0177",
            linestyle=":",
            linewidth=1.8,
        )
    )
    zoom_axes.plot(
        (diameter.point_a.x, diameter.point_b.x),
        (diameter.point_a.y, diameter.point_b.y),
        color="#7a0177",
        linewidth=2.2,
    )
    zoom_axes.scatter(DEMO_TARGET.x, DEMO_TARGET.y, marker="*", color="black", s=65)
    zoom_axes.scatter(
        [point.x for point in region.vertices],
        [point.y for point in region.vertices],
        color="#a63603",
        s=18,
    )
    x_values = [point.x for point in region.vertices]
    y_values = [point.y for point in region.vertices]
    padding = max(max(x_values) - min(x_values), max(y_values) - min(y_values)) * 0.25
    zoom_axes.set_xlim(min(x_values) - padding, max(x_values) + padding)
    zoom_axes.set_ylim(min(y_values) - padding, max(y_values) + padding)
    zoom_axes.set_aspect("equal", adjustable="box")
    zoom_axes.text(
        0.03,
        0.96,
        "Final region detail",
        transform=zoom_axes.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    zoom_axes.grid(alpha=0.20)
    zoom_axes.tick_params(labelsize=7)

    axes.set_aspect("equal", adjustable="box")
    axes.set_xlabel("x (m)")
    axes.set_ylabel("y (m)")
    axes.set_title("Q1 bearing-wedge intersection and diameter circle")
    axes.grid(alpha=0.22)
    axes.legend(loc="best")
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)
    return output_path


def main() -> None:
    """Print the Q1 demo result and generate its validation figure."""

    region, history = solve_demo()
    output_path = Path("figures/q1/q1_localization_demo.png").resolve()
    create_demo_figure(region, history, output_path)
    diameter = region.diameter()
    coverage = region.diameter_circle_coverage()

    print(f"Observation count: {len(region.measurements)}")
    print(f"Region status: {region.status.value}")
    print("Polygon vertices (CCW):")
    for point in region.vertices:
        print(f"  ({point.x:.6f}, {point.y:.6f})")
    print(f"Polygon area: {region.area():.6f} m^2")
    print(f"Maximum diameter: {diameter.length:.6f} m")
    print(
        "Representative diameter endpoints: "
        f"({diameter.point_a.x:.6f}, {diameter.point_a.y:.6f}) -> "
        f"({diameter.point_b.x:.6f}, {diameter.point_b.y:.6f})"
    )
    print(f"Diametral-pair count: {len(diameter.diametral_pairs)}")
    print(
        "Exists covering diameter circle: "
        f"{coverage.exists_covering_diameter_circle}"
    )
    print(f"All diameter circles cover: {coverage.all_diameter_circles_cover}")
    print(f"Figure: {output_path}")


if __name__ == "__main__":
    main()
