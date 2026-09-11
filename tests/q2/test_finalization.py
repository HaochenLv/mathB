from experiments.q2.run_q2_finalization import classify_target
from experiments.q2.run_q2_key_validation import Point


def test_target_classification_keeps_near_out_of_q1() -> None:
    row = classify_target(Point(9.0, 0.0), 6.0, 0.0, (-1.0, 0.0, 1.0))

    assert row["failure_reason"] == "near_direct_clear"
    assert row["operational_good"]
    assert not row["direction_available"]


def test_target_classification_reports_receive_failure_first() -> None:
    row = classify_target(Point(-1000.0, 0.0), 6.0, 0.0, (-1.0, 0.0, 1.0))

    assert row["failure_reason"] == "guaranteed_rx_failure"
    assert not row["robust_clear_ready"]


def test_orthogonal_target_reaches_geometry_classification() -> None:
    row = classify_target(Point(800.0, -400.0), 800.0, 0.0, (-1.0, 0.0, 1.0))

    assert row["guaranteed_rx"]
    assert row["direction_available"]
    assert row["failure_reason"] in {
        "diameter_over_40",
        "diameter_circle_failure",
        "robust_clear_ready",
    }
