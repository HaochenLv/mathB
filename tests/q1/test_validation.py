from experiments.q1.run_q1_validation import run_validation


def test_seeded_monte_carlo_properties() -> None:
    summary = run_validation(cases=100, seed=20260910)
    assert summary.failed == 0
    assert summary.passed + summary.unbounded_or_skipped == summary.total

