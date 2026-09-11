import pytest

from experiments.q2.run_q2_key_validation import (
    CONFIG,
    Point,
    evaluate_q1,
    receiver_state,
    target_point,
)


def test_receiver_state_uses_observation_consistent_minimum_radius() -> None:
    target = target_point(1200.0, 0.0)
    guaranteed = receiver_state(Point(0.0, 0.0), target, 1200.0)
    not_guaranteed = receiver_state(Point(-1.0, 0.0), target, 1200.0)

    assert guaranteed.guaranteed_rx
    assert not not_guaranteed.guaranteed_rx
    assert not_guaranteed.possible_rx


def test_near_is_not_reported_as_direction_available() -> None:
    target = Point(100.0, 0.0)
    state = receiver_state(Point(104.0, 0.0), target, 100.0)

    assert state.near
    assert state.guaranteed_rx
    assert not state.direction_available


@pytest.mark.parametrize("epsilon", [-1.0, 0.0, 1.0])
def test_exact_q1_orthogonal_construction_is_bounded(epsilon: float) -> None:
    target = Point(CONFIG.orthogonal_target_x, 0.0)
    sensor = Point(target.x, -400.0)
    outcome = evaluate_q1(sensor, target, epsilon)

    assert outcome.bounded
    assert outcome.diameter is not None
