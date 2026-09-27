"""Deterministic budget tests: no test waits for wall-clock time to pass."""

import json
from unittest.mock import patch

import numpy as np

from hyperell_regulator.regulator import _legs
from hyperell_regulator.regulator.diagnostics import (
    IntegrationLimitError, IntegrationMonitor, checkpoint, phase,
)


class Clock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now


def test_progress_is_throttled_and_records_are_serializable():
    clock, events = Clock(), []
    with IntegrationMonitor(clock=clock, callback=events.append) as monitor:
        with monitor.edge("finite", degree=14, za=1+2j, zb=3):
            with phase("side", target=3):
                clock.now = 0.5
                checkpoint("walk_attempts", t=0.1)
                clock.now = 1.1
                checkpoint("walk_attempts", t=0.2, chart="x", sheet=2)
                clock.now = 1.2
                checkpoint("newton_iterations", amount=3)
            clock.now = 1.5
    assert [event["event"] for event in events] == ["start", "progress", "end"]
    assert events[1]["state"]["sheet"] == 2
    assert events[1]["phase"] == "side"
    record = monitor.records[0]
    assert record["seconds"] == 1.5
    assert record["steps"] == 2
    assert record["counters"]["newton_iterations"] == 3
    assert record["phases"][0]["counters"]["walk_attempts"] == 2
    json.dumps(monitor.records, allow_nan=False)
    json.dumps(events, allow_nan=False)


def test_work_limit_does_not_double_count_newton_iterations():
    monitor = IntegrationMonitor(max_edge_seconds=None, max_edge_steps=2, clock=Clock())
    try:
        with monitor, monitor.edge("finite", degree=7):
            with phase("side", target=0):
                checkpoint("walk_attempts", amount=2, t=0.5)
                checkpoint("newton_iterations", amount=100)
                checkpoint("ode_evaluations", t=0.6)
    except IntegrationLimitError as exc:
        assert exc.record["steps"] == 3
        assert exc.record["state"]["t"] == 0.6
        assert exc.record["phase"] == "side"
    else:
        raise AssertionError("the third continuation/ODE operation must stop the edge")
    assert monitor.records[0]["status"] == "limit"
    assert monitor.records[0]["failure"]["state"]["t"] == 0.6
    assert not issubclass(IntegrationLimitError, RuntimeError)
    # A stopped monitor is removed from the current context.
    checkpoint("walk_attempts", amount=1000)


def test_time_budget_interrupts_inside_ode_rhs():
    clock, evaluations = Clock(), []
    monitor = IntegrationMonitor(max_edge_seconds=2, clock=clock)

    def solver(rhs, span, y0, **kwargs):
        # Simulate a solver stuck taking evaluations without making progress.
        for _ in range(100):
            clock.now += 1
            rhs(0.0, y0)
        raise AssertionError("a stalled solver must be interrupted internally")

    def rhs(t, y):
        evaluations.append(t)
        return y

    try:
        with patch.object(_legs, "solve_ivp", solver):
            with monitor, monitor.edge("finite", degree=14):
                _legs._solve(rhs, (0.0, 1.0), np.zeros(1))
    except IntegrationLimitError as exc:
        assert exc.record["counters"]["ode_evaluations"] == 3
        assert exc.record["phase"] == "ode_solve"
    else:
        raise AssertionError("the ODE exceeded its time budget")
    assert len(evaluations) == 2
    assert monitor.records[0]["status"] == "limit"
    assert monitor.records[0]["counters"]["ode_attempts"] == 1


def test_edge_entry_points_start_and_end_monitoring():
    class Cover:
        N = 2

        def fibre(self, w):
            return np.ones(2, complex), np.ones(2, complex)

    def side(*args, **kwargs):
        checkpoint("walk_attempts", t=1.0)
        return _legs.Leg(np.zeros((2, 2), complex), np.zeros((2, 2), complex), None)

    with patch.object(_legs, "side", side):
        with IntegrationMonitor(clock=Clock()) as monitor:
            _legs.edge_leg(Cover(), 0.0, 1.0, 0.0, "P", "simple")
    assert len(monitor.records) == 1
    record = monitor.records[0]
    assert record["edge"] == "finite"
    assert record["degree"] == 2
    assert record["state"]["za"] == 0.0
    assert record["steps"] == 2
    assert record["status"] == "ok"


def test_edge_budget_resets_for_each_edge():
    clock = Clock()
    with IntegrationMonitor(max_edge_seconds=1, max_edge_steps=1, clock=clock) as monitor:
        for index in range(2):
            with monitor.edge("finite", degree=2):
                checkpoint("walk_attempts")
                clock.now += 0.75
    assert [record["seconds"] for record in monitor.records] == [0.75, 0.75]
    assert [record["steps"] for record in monitor.records] == [1, 1]


def test_nested_monitor_contexts_restore_the_parent():
    with IntegrationMonitor(clock=Clock()) as outer:
        with outer.edge("finite", degree=7):
            checkpoint("walk_attempts")
            with IntegrationMonitor(clock=Clock()) as inner:
                with inner.edge("finite", degree=14):
                    checkpoint("walk_attempts", amount=2)
            checkpoint("walk_attempts")
    assert outer.records[0]["steps"] == 2
    assert inner.records[0]["steps"] == 2
