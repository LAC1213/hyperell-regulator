"""Plot the stable reduction graph of a cached genus-two curve."""

import argparse
from pathlib import Path

from sage.all import Graph, ZZ
from sage.plot.text import Text

from .components import CACHE, curve_from_row, load_curves
from .paths import ensure_parent
from .stable_model import stable_model


def plot_dual_graph(model, *, label=None):
    """Return a Sage plot with component genera and exact node thicknesses."""
    fibre = model.stable_graph
    graph = Graph(loops=True, multiedges=True)
    graph.add_vertices(range(len(fibre.genera)))
    graph.add_edges((i, j, str(n)) for i, j, n in fibre.edges)
    labels = {i: 'C%d\ng = %d' % (i, g) for i, g in enumerate(fibre.genera)}
    colors = {color: [i for i, g in enumerate(fibre.genera) if g == genus]
              for genus, color in enumerate(('#e8edf3', '#c0e0df', '#c9daef'))}
    prime = model.prime.smallest_integer()
    caption = '%s at p = %s' % (label or 'Stable reduction', prime)
    if fibre.edges:
        normalization = model.length_normalization
        if normalization == 'uniformizer':
            normalization = 'v(π) = 1, e = %s' % model.prime.ramification_index()
        caption += '\nNode thicknesses with ' + normalization
    layout = {'pos': {0: (-0.95, 0), 1: (0.95, 0)}} if len(labels) == 2 else {'layout': 'circular'}
    plot = graph.plot(**layout, vertex_labels=labels, vertex_colors=colors,
                      vertex_size=1900, label_fontsize=11, edge_labels=True,
                      edge_labels_background='white', edge_color='#37546c',
                      edge_thickness=1.6, loop_size=0.26, dist=0.22,
                      title=caption, title_pos=(0.5, 1.02))
    for primitive in plot:
        if isinstance(primitive, Text) and 'background_color' in primitive.options():
            primitive.set_options(dict(primitive.options(), zorder=9))
    limits = plot.get_minmax_data()
    plot.set_axes_range(limits['xmin']-0.36, limits['xmax']+0.36,
                        limits['ymin']-0.36, limits['ymax']+0.36)
    return plot


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('curve', help='curve label in the local LMFDB cache')
    parser.add_argument('prime', type=int, help='rational prime')
    parser.add_argument('--cache', type=Path, default=CACHE, help='local LMFDB curve cache')
    parser.add_argument('--out', type=Path, help='output image path, such as graph.png or graph.svg')
    args = parser.parse_args(argv)
    if not ZZ(args.prime).is_prime():
        parser.error('the place must be a rational prime')
    try:
        rows = load_curves(path=args.cache, simple_only=False)
    except OSError as error:
        parser.error(str(error))
    row = next((row for row in rows if row['label'] == args.curve), None)
    if row is None:
        parser.error('curve %s is absent from %s' % (args.curve, args.cache))
    try:
        model = stable_model(curve_from_row(row), args.prime)
    except NotImplementedError as error:
        parser.error(str(error))
    out = args.out or Path('%s-p%s-dual-graph.png' % (args.curve, args.prime))
    plot_dual_graph(model, label=args.curve).save(str(ensure_parent(out)), figsize=(6.4, 3.6), dpi=160)
    print(out)


if __name__ == '__main__':
    main()
