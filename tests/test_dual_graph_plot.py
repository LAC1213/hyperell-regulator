"""Stable-reduction plots retain the exact labels of all nodes."""

from types import SimpleNamespace

import pytest
from sage.all import QQ
from sage.plot.text import Text

from hyperell_regulator.dual_graph_plot import plot_dual_graph
from hyperell_regulator.graph_boundary import FibreGraph
from hyperell_regulator.stable_model import StableModel


@pytest.mark.parametrize('genera, edges', [
    ([0], [(0, 0, QQ(1)/3), (0, 0, QQ(2)/3)]),
    ([0, 0], [(0, 1, QQ(1)/3), (0, 1, QQ(2)/3), (0, 1, QQ(2)/3)]),
])
def test_plot_labels_loops_and_parallel_edges(genera, edges):
    prime = SimpleNamespace(smallest_integer=lambda: 2)
    model = StableModel(QQ, prime, FibreGraph(genera, edges), (),
                        length_normalization='v(2)=1')
    plot = plot_dual_graph(model)
    labels = [primitive.string for primitive in plot if isinstance(primitive, Text)]
    assert sorted(t for t in labels if '\n' not in t) == sorted(str(e[2]) for e in edges)
    assert sorted(t for t in labels if '\n' in t) == [
        'C%d\ng = %d' % (i, g) for i, g in enumerate(genera)]
