"""Deterministic budget tests: no test waits for wall-clock time to pass."""

import json
from unittest.mock import patch

import numpy as np

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
                checkpoint("quadrature_nodes", t=0.1)
                clock.now = 1.1
                checkpoint("quadrature_nodes", t=0.2, chart="x", sheet=2)
                clock.now = 1.2
                checkpoint("extra_checks", amount=3)
            clock.now = 1.5
    assert [event["event"] for event in events] == ["start", "progress", "end"]
    assert events[1]["state"]["sheet"] == 2
    assert events[1]["phase"] == "side"
    record = monitor.records[0]
    assert record["seconds"] == 1.5
    assert record["steps"] == 2
    assert record["counters"]["extra_checks"] == 3
    assert record["phases"][0]["counters"]["quadrature_nodes"] == 2
    json.dumps(monitor.records, allow_nan=False)
    json.dumps(events, allow_nan=False)


def test_work_limit_does_not_double_count_extra_checks():
    monitor = IntegrationMonitor(max_edge_seconds=None, max_edge_steps=2, clock=Clock())
    try:
        with monitor, monitor.edge("finite", degree=7):
            with phase("side", target=0):
                checkpoint("quadrature_nodes", amount=2, t=0.5)
                checkpoint("extra_checks", amount=100)
                checkpoint("quadrature_nodes", t=0.6)
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
    checkpoint("quadrature_nodes", amount=1000)


def test_edge_budget_resets_for_each_edge():
    clock = Clock()
    with IntegrationMonitor(max_edge_seconds=1, max_edge_steps=1, clock=clock) as monitor:
        for index in range(2):
            with monitor.edge("finite", degree=2):
                checkpoint("quadrature_nodes")
                clock.now += 0.75
    assert [record["seconds"] for record in monitor.records] == [0.75, 0.75]
    assert [record["steps"] for record in monitor.records] == [1, 1]


def test_nested_monitor_contexts_restore_the_parent():
    with IntegrationMonitor(clock=Clock()) as outer:
        with outer.edge("finite", degree=7):
            checkpoint("quadrature_nodes")
            with IntegrationMonitor(clock=Clock()) as inner:
                with inner.edge("finite", degree=14):
                    checkpoint("quadrature_nodes", amount=2)
            checkpoint("quadrature_nodes")
    assert outer.records[0]["steps"] == 2
    assert inner.records[0]["steps"] == 2


def test_experiment_failure_retains_edge_statistics():
    from hyperell_regulator.experiments import beilinson as b

    clock = Clock()
    monitor = IntegrationMonitor(clock=clock)
    result = dict(label='test')
    with patch('hyperell_regulator.regulator.IntegrationMonitor', return_value=monitor), \
            patch.object(b.time, 'monotonic', clock):
        with b._monitored_result(result, False, 15, 100000):
            with monitor.edge('finite', degree=2):
                checkpoint('quadrature_nodes')
                clock.now = 1.5
                raise ValueError('preparation failed')
    assert result['status'] == 'preparation_error'
    assert result['error'] == 'ValueError: preparation failed'
    assert result['edges_attempted'] == 1
    assert result['seconds'] == result['slowest_edge_seconds'] == 1.5
    # A failed experiment must release its monitor before the next run.
    checkpoint('quadrature_nodes', amount=100001)
