"""Export marked semistable fibre graphs for the static curve explorer."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from sage.all import NumberField, PolynomialRing, QQ, lcm, matrix, number_field_elements_from_algebraics, vector

from ..components import curve_from_row
from ..graph_boundary import FibreGraph
from ..number_field_integrality import _LocalSections, _curve_polynomials
from .._projected_reduction import _split_node_field, nodal_model
from ..stable_model import StableModel, _extend, stable_model
from .geometric_torsion import decode_point


SCHEMA_VERSION = 1


def regular_graph(model):
    """Resolve metric edges into unit edges after a common finite base change.

    A self-node is separated by a further quadratic ramification when
    necessary. Existing marked vertices survive with their original indices.
    For wild disc charts, include the chart field, square completion and
    disc-scale ramification in the common valuation scale.
    """
    graph = model.graph
    denominators = lcm([length.denominator() for _, _, length in graph.edges] or [1])
    scale = denominators
    if model.length_normalization == 'v(2)=1':
        extension = model.prime.ramification_index()
        for chart in model.charts:
            if 'centre_prime' in chart:
                extension *= 2*(1 if chart.get('parameters_in_model_field') else
                                chart['centre_prime'].ramification_index())
            if 'rho_power' in chart:
                extension *= int(chart['rho_power'][0])
        scale = lcm(denominators, extension)
    if any(a == b and scale*length == 1 for a, b, length in graph.edges):
        scale *= 2
    if len(graph.genera)+sum(int(scale*n)-1 for _, _, n in graph.edges) > 400:
        raise NotImplementedError('the regular resolution exceeds 400 components')
    genera, edges = list(graph.genera), []
    for a, b, length in graph.edges:
        count = int(scale*length)
        if count < 1 or count != scale*length:
            raise ValueError('the base-change scale must make node thicknesses integral')
        path = [a]
        for _ in range(count-1):
            path.append(len(genera))
            genera.append(0)
        path.append(b)
        edges.extend((i, j, 1) for i, j in zip(path, path[1:]))
    if len(genera) > 400:
        raise NotImplementedError('the regular resolution exceeds 400 components')
    return FibreGraph(genera, edges), int(scale)


def precise_nodal_specializations(model, local, points):
    """Distinguish the two arms of a resolved node using its tangent lines.

    The projected integrality API intentionally identifies reflected arms
    when they have the same boundary signature. The picture needs the
    actual section: a Newton lift of the critical point determines which
    tangent branch it follows. Orient a self-node by the ordered slopes.
    """
    if model.field is not local.K:
        return
    paths = model.charts[0]['paths']
    nodes = model.charts[0]['nodes']
    resolved, critical = list(model.specializations), {}
    for i, P in enumerate(points):
        if resolved[i] is not None and len(set(resolved[i])) == 1:
            continue
        label, _, steps = local.point(P)
        if label is None or len(steps) < 2:
            continue
        key = (steps[0]['patch'], steps[0]['residue_x'], steps[0]['residue_y'])
        if key not in paths:
            continue
        path, depth = paths[key], len(steps)-1
        n = len(path)-1
        if 2*depth == n:
            resolved[i] = (path[depth],)
            continue
        if 2*depth > n:
            continue
        if key not in critical:
            a, b = local.k(key[1]), local.k(key[2])
            f, h = (local.f, local.h) if key[0] == 'finite' else (local.inf_f, local.inf_h)
            x, y = local.K(local.lift(a)), local.K(local.lift(b))
            for _ in range(10):
                gradient = vector(local.K, [h.derivative()(x)*y-f.derivative()(x), 2*y+h(x)])
                if 2*min(map(local.valuation, gradient)) > n:
                    break
                H = matrix(local.K, [[h.derivative(2)(x)*y-f.derivative(2)(x), h.derivative()(x)],
                                     [h.derivative()(x), 2]])
                x, y = vector(local.K, [x, y])-H.solve_right(gradient)
            else:
                continue
            fb, hb = map(local.reduce_poly, (f, h))
            X = fb.parent().gen()
            g, j = fb(X+a), hb(X+a)
            T = PolynomialRing(local.k, 't')
            t = T.gen()
            slopes = sorted((t*t+j[1]*t+b*j[2]-g[2]).roots(multiplicities=False), key=str)
            if len(slopes) != 2:
                continue
            ends = nodes[key][:2]
            if ends[0] != ends[1]:
                factor = (local.main if key[0] == 'finite' else local.main_infinity)[ends[0]]
                tangent = -(factor[0].derivative()(a)+factor[1].derivative()(a)*b)/factor[1](a)
                slopes.sort(key=lambda r: r != tangent)
            critical[key] = x, y, slopes
        x, y, slopes = critical[key]
        if key[0] == 'infinity':
            px, py = ((local.K(0), P[1]) if P[0] == 'infinity' else (1/P[0], P[1]/P[0]**3))
        else:
            px, py = P
        if px == x or local.valuation(py-y) < local.valuation(px-x):
            continue
        slope = local.reduction((py-y)/(px-x))
        if slope in slopes:
            resolved[i] = (path[depth] if slope == slopes[0] else path[-1-depth],)
    model.specializations = tuple(resolved)


def good_square_chart(local, mark_count):
    """Verify the smooth fifth-root chart for the two square 688 models."""
    models = {
        (1, 2, 1, 4, 0, 2): (1, 1),
        (1, 4, 2, -8, -7, 2): (1, 2, 1),
    }
    key = tuple(QQ(a) for a in local.f)
    if local.h or key not in models or local.k.characteristic() != 2:
        raise NotImplementedError('no verified square-completion chart for this equation')
    T = PolynomialRing(local.K, 's')
    K, q, embed, rho = _extend(local.K, local.prime, T.gen()**5-2, 's')
    R = PolynomialRing(K, 'u')
    u = R.gen()
    f = R([embed(a)*rho**i for i, a in enumerate(local.f)])
    a = R([K(b)*rho**i for i, b in enumerate(models[key])])
    g = (f-a*a)/4
    chart = _LocalSections(g, a, q, max_blowups=64)
    model = nodal_model(chart, [], [])
    if model.graph.genera != (2,) or model.graph.edges:
        raise ArithmeticError('the proposed good-reduction chart was not smooth')
    model.specializations = ((0,),)*mark_count
    return model


def _lift_to(field, value):
    """Follow Sage's relative/absolute field structure maps exactly."""
    if value in QQ:
        return field(QQ(value))
    if value.parent() is field:
        return value
    morphism = field.coerce_map_from(value.parent())
    if morphism is not None:
        return morphism(value)
    structure = field.structure()
    if structure:
        to_absolute = structure[1]
        relative = to_absolute.domain()
        return to_absolute(relative(_lift_to(relative.base_field(), value)))
    raise ValueError('the marked-point field is not in the recorded extension tower')


def precise_tame_specializations(model, C, points):
    """Use normalized ordinates to distinguish split tame components."""
    if not model.charts or 'reciprocal_centre' not in model.charts[0]:
        return
    K, q = model.field, model.prime
    red = q.residue_field().reduction_map()
    f, h = _curve_polynomials(C, K)
    r = model.charts[0]['reciprocal_centre']
    F = h*h+4*f
    if r is not None:
        shifted = F(F.parent().gen()+r)
        F = F.parent()([shifted[6-i] for i in range(7)])
    roots = F.roots(multiplicities=False)
    resolved = list(model.specializations)
    for index, (x, y) in enumerate(points):
        if resolved[index] is not None:
            continue
        y = _lift_to(K, y)
        if x == 'infinity':
            if r is None:
                continue
            t, Y = K(0), 2*y+h[3]
        else:
            x = _lift_to(K, x)
            if r is not None and x == r:
                continue
            t = x if r is None else 1/(x-r)
            Y = (2*y+h(x))*(1 if r is None else t**3)
        charts = [chart for chart in model.charts
                  if (t-chart['centre']).valuation(q) >= chart['scale'].valuation(q)]
        if not charts:
            continue
        chart = max(charts, key=lambda c: c['scale'].valuation(q))
        if len(chart['vertices']) == 1:
            resolved[index] = chart['vertices']
            continue
        a = red((t-chart['centre'])/chart['scale'])
        factors = {}
        for root in roots:
            if (root-chart['centre']).valuation(q) >= chart['scale'].valuation(q):
                residue = red((root-chart['centre'])/chart['scale'])
                factors[residue] = factors.get(residue, 0)+1
        square = red(K(1))
        for residue, count in factors.items():
            square *= (a-residue)**(count//2)
        if not square:
            continue
        ordinate = red(Y/chart['y_scale'])/square
        eta = chart['residual_polynomial'].leading_coefficient().sqrt()
        if ordinate == eta or ordinate == -eta:
            resolved[index] = (chart['vertices'][0 if ordinate == eta else 1],)
    model.specializations = tuple(resolved)


def precise_bolza_specializations(model, points):
    """Put the section coordinates and both elliptic discs in one local field."""
    if model.length_normalization != 'v(2)=1' or not model.charts:
        return
    chart = model.charts[0]
    if 'centre_prime' not in chart or chart['radius'] != 1 or chart['reciprocal']:
        return
    K, q, centre = model.field, model.prime, chart['centre']
    A, place = centre.parent(), chart['centre_prime']
    structure = A.structure()
    if not structure:
        return
    relative, to_absolute = structure[1].domain(), structure[1]
    unit = to_absolute(relative.gen())
    if unit*unit != -1:
        return
    for L, from_K, from_A, _ in K.composite_fields(A, both_maps=True):
        if L.degree() > 256:
            continue
        candidates = L.ideal([from_K(a) for a in q.gens()]+[from_A(a) for a in place.gens()]).prime_factors()
        if not candidates:
            continue
        prime = candidates[0]
        c, other = from_A(centre), from_A(centre*unit)
        e = prime.ramification_index()
        join = (c-other).valuation(prime)/e
        if join != QQ(1)/2:
            continue
        weights = []
        for x, _ in points:
            if x == 'infinity':
                weights.append(QQ(1)/2)
                continue
            x = from_K(x)
            a = max(join, min(QQ(1), (x-c).valuation(prime)/e))
            b = max(join, min(QQ(1), (x-other).valuation(prime)/e))
            weights.append(QQ(1)/2+a-b)
        depths = sorted({QQ(0), QQ(1), *weights})
        model.field, model.prime = L, prime
        model.graph = FibreGraph([1]+[0]*(len(depths)-2)+[1],
                                 [(i, i+1, (b-a)/2) for i, (a, b) in enumerate(zip(depths, depths[1:]))])
        model.specializations = tuple((depths.index(w),) for w in weights)
        model.charts = (dict(chart, centre=c, centre_prime=prime, parameters_in_model_field=True),)
        return


def encode_model(model, support):
    graph, scale = regular_graph(model)
    if len(model.specializations) != len(support):
        raise ValueError('the local model must retain every marked section')
    if sum(graph.genera)+graph.cycles.nrows() != 2:
        raise ArithmeticError('the geometric fibre must have arithmetic genus 2')
    marks = []
    for point, positions in zip(support, model.specializations):
        positions = sorted(set(positions)) if positions is not None else []
        marks.append(dict(point_id=point['id'], components=positions,
                          status='known' if len(positions) == 1 else
                                 'ambiguous' if positions else 'unresolved'))
    return dict(status='available', genera=list(graph.genera),
                edges=[[a, b] for a, b, _ in graph.edges],
                cycles=[[str(x) for x in row] for row in graph.cycles],
                laplacian=[[str(x) for x in row] for row in graph.laplacian],
                specializations=marks,
                source_metric_edges=[[a, b, str(n)] for a, b, n in model.graph.edges],
                source_vertex_count=len(model.graph.genera),
                valuation_normalization=model.length_normalization,
                regularization_scale=scale,
                field_polynomial=list(map(str, model.field.polynomial())),
                prime_ideal=str(model.prime),
                ramification_index=int(model.prime.ramification_index()),
                residue_degree=int(model.prime.residue_class_degree()))


def compute_prime(payload):
    entry, p = payload['entry'], int(payload['prime'])
    C = curve_from_row(entry['original_curve'])
    support = entry['support']
    points = [decode_point(P['lmfdb']['point']) for P in support]
    flat = [a for x, y in points for a in ((y,) if x == 'infinity' else (x, y))]
    K, values, _ = number_field_elements_from_algebraics(flat, embedded=True)
    if K is QQ:
        K = NumberField(PolynomialRing(QQ, 'a').gen(), 'a')
        values = list(map(K, values))
    values = iter(values)
    points = [('infinity' if x == 'infinity' else next(values), next(values)) for x, _ in points]
    places = []
    for prime in K.primes_above(p):
        try:
            f, h = _curve_polynomials(C, K)
            exponent = max([0, *(-a.valuation(prime) for a in list(f)+list(h) if a)])
            scale = K.uniformizer(prime)**exponent if exponent else K(1)
            marks = [(x, scale*y) for x, y in points]
            local = _LocalSections(scale**2*f, scale*h, prime, max_blowups=64)
            try:
                local, marks = _split_node_field(local, marks)
                model = nodal_model(local, marks, [(i,) for i in range(len(marks))])
                precise_nodal_specializations(model, local, marks)
            except (ArithmeticError, NotImplementedError, ValueError):
                try:
                    model = good_square_chart(local, len(points))
                except (ArithmeticError, NotImplementedError, ValueError):
                    model = stable_model(C, prime, field=K, points=points,
                                         max_degree=96, max_blowups=64)
            precision_issue = None
            if any(mark is None for mark in model.specializations):
                try:
                    precise_tame_specializations(model, C, points)
                except (ArithmeticError, NotImplementedError, ValueError, TypeError) as error:
                    precision_issue = str(error)
            if tuple(QQ(a) for a in C.f) == (0, -1, 0, 0, 0, 1) and not C.h:
                if payload.get('stream_partial'):
                    partial = dict(encode_model(model, support), input_prime_ideal=str(prime))
                    print(json.dumps(dict(prime=p, status='available', places=[*places, partial],
                                          partial_specializations=True)), flush=True)
                try:
                    precise_bolza_specializations(model, points)
                except (ArithmeticError, NotImplementedError, ValueError, TypeError) as error:
                    precision_issue = str(error)
            data = encode_model(model, support)
            if precision_issue:
                data['specialization_issue'] = precision_issue
        except (ArithmeticError, NotImplementedError, ValueError, TypeError) as error:
            data = dict(status='unavailable', reason=str(error))
        places.append(dict(data, input_prime_ideal=str(prime)))
    return dict(prime=p, status='available' if any(x['status'] == 'available' for x in places)
                else 'unavailable', places=places)


def _job(entry, p, cache, seconds):
    payload = dict(entry=dict(original_curve=entry['original_curve'], support=entry['support']),
                   prime=p, stream_partial=True)
    sources = [Path(__file__), *[Path(__file__).parents[1]/name for name in (
        'stable_model.py', '_projected_reduction.py', '_dyadic_projected.py', 'graph_boundary.py',
        'number_field_integrality.py', 'components.py', 'experiments/geometric_torsion.py')]]
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()+
                                 b''.join(path.read_bytes() for path in sources)).hexdigest()
    path = cache/(fingerprint+'.json')
    if path.exists():
        saved = json.loads(path.read_text())
        if not saved.get('partial_specializations') or saved.get('export_seconds', 0) >= seconds:
            return saved
    try:
        process = subprocess.run([sys.executable, '-m', __spec__.name, '--worker'],
                                 input=json.dumps(payload), text=True, capture_output=True,
                                 timeout=seconds)
        if process.returncode:
            raise ArithmeticError(process.stderr.strip().splitlines()[-1])
        result = json.loads(process.stdout.strip().splitlines()[-1])
    except subprocess.TimeoutExpired as error:
        partial = error.stdout.decode() if isinstance(error.stdout, bytes) else error.stdout
        if partial:
            result = json.loads(partial.strip().splitlines()[-1])
        else:
            result = dict(prime=p, status='unavailable', places=[],
                          reason='The exact local model exceeded the export time budget.')
    except (ArithmeticError, ValueError, IndexError) as error:
        result = dict(prime=p, status='unavailable', places=[], reason=str(error))
    result['export_seconds'] = seconds
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+'\n')
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog', type=Path, default=Path('docs/_extra/explorer/catalog.json'))
    parser.add_argument('--out', type=Path, default=Path('docs/_extra/explorer/reduction'))
    parser.add_argument('--cache', type=Path, default=Path('results/fibre_explorer'))
    parser.add_argument('--seconds', type=float, default=45)
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--label', action='append', default=[])
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.worker:
        print(json.dumps(compute_prime(json.load(sys.stdin))))
        return
    if args.jobs < 1 or args.seconds <= 0:
        parser.error('jobs and seconds must be positive')
    curves = json.loads(args.catalog.read_text())['curves']
    records = {}
    args.out.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {}
        for entry in curves:
            if args.label and entry['label'] not in args.label:
                continue
            bad = list(map(int, curve_from_row(entry['original_curve']).bad_primes()))
            records[entry['label']] = dict(schema_version=SCHEMA_VERSION, label=entry['label'],
                                          bad_primes=bad, fibres=[])
            for p in bad:
                pending[pool.submit(_job, entry, p, args.cache, args.seconds)] = entry['label']
        for future in as_completed(pending):
            label, result = pending[future], future.result()
            records[label]['fibres'].append(result)
            print('%s p=%s: %s' % (label, result['prime'], result['status']), flush=True)
    for label, record in records.items():
        record['fibres'].sort(key=lambda x: x['prime'])
        path = args.out/(label+'.json')
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
        temporary.replace(path)


if __name__ == '__main__':
    main()
