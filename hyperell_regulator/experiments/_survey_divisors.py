"""Shared exact divisor and integration metadata for tables and the explorer."""

from decimal import Decimal, localcontext
from fractions import Fraction

from .._algebraic import complex_coordinates, decode_algebraic, encode_algebraic


def decimal_value(value, digits=150):
    """Display a ball by its midpoint without discarding its saved endpoints."""
    if value is None:
        return None
    if isinstance(value, dict) and value.get('type') == 'real_ball':
        midpoint = (Fraction(value['lower'])+Fraction(value['upper']))/2
        with localcontext() as context:
            context.prec = digits
            return str(Decimal(midpoint.numerator)/Decimal(midpoint.denominator))
    text = str(value)
    if '/' in text:
        q = Fraction(text)
        with localcontext() as context:
            context.prec = digits
            return str(Decimal(q.numerator)/Decimal(q.denominator))
    return text


def comparison_error(row):
    linear = row.get('linear_comparison', {})
    if len(linear.get('terms', [])) > 1:
        value = linear.get('relative_residual')
    else:
        comparison = row.get('comparison', {})
        value = comparison.get('relative_agreement_upper', comparison.get('relative_agreement'))
    if isinstance(value, dict) and value.get('type') == 'real_ball':
        value = value['upper']
    return decimal_value(value)


def selected_divisors(row):
    """Return point records, selected coefficient rows and their saved IDs."""
    terms = row.get('linear_comparison', {}).get('terms', [])
    if len(terms) > 1:
        search = row['pair_searches'][row['pair_search_index']]
        kernel = row['integral_kernels'][search['kernel_index']]
        indices = sorted({i for term in terms for i in term['basis_indices']})
        return kernel['points'], [kernel['divisor_basis'][i] for i in indices], indices
    return row['points'], row['divisors'], row.get('basis_indices', [0, 1])


def packet_order(row):
    """Exact exponent of packet differences, or None when it was not recorded."""
    if row.get('packet_torsion_order') is not None and row.get('packet_torsion_order_exact'):
        return str(row['packet_torsion_order'])
    for kernel in row.get('integral_kernels', []):
        if kernel.get('points') == row.get('points') and kernel.get('packet_torsion_order_exact'):
            return str(kernel['packet_torsion_order'])
    return None


def endomorphism_type(row, curve_row=None):
    return row.get('end_alg') or row.get('input_curve', {}).get('end_alg') or (curve_row or {}).get('end_alg')


def divisor_metadata(row, curve_row=None):
    """Display divisors in LMFDB coordinates, retaining exact integration images.

    Exact point records and coefficients are retained. Approximate complex
    coordinates are supplied separately for display, never for integration.
    Each support entry retains its integration image in ``point`` and its
    original coordinates in ``lmfdb`` for divisor formulas and support tables.
    """
    from sage.all import PolynomialRing, QQ, QQbar, latex
    from .geometric_torsion import decode_point, encode_point

    raw_points, coefficients, indices = selected_divisors(row)
    coefficients = [[QQ(c) for c in D] for D in coefficients]
    if any(len(D) != len(raw_points) or sum(D) for D in coefficients):
        raise ValueError('invalid degree-zero divisor coefficients on '+row['label'])
    algebraic = row['algebraic_data']
    R = PolynomialRing(QQ, 'x')
    h = R(algebraic['original_h'])
    shift = None if algebraic.get('shift') is None else QQ(algebraic['shift'])
    scale = QQ(algebraic.get('scale', 1))

    def image(record):
        x, y = decode_point(record)
        if shift is None:
            return (('infinity', 2*y+h[3]) if x == 'infinity' else (x, 2*y+h(x)))
        if x == 'infinity':
            return QQbar(0), scale*(2*y+h[3])
        if x == shift:
            return 'infinity', scale*(2*y+h(x))
        return 1/(x-shift), scale*(2*y+h(x))/(x-shift)**3

    originals = [decode_point(P) for P in raw_points]
    images = [image(P) for P in raw_points]

    R = PolynomialRing(QQ, 'x')
    G = R(algebraic['model'].replace('u', 'x'))
    for x, y in images:
        if (y*y != G[6] if x == 'infinity' else y*y != G(x)):
            raise ValueError('transformed divisor point is off the integration model: '+row['label'])
    lmfdb_model = None
    if algebraic.get('original_f') is not None:
        f, h = R(algebraic['original_f']), R(algebraic['original_h'])
        for x, y in originals:
            if y*y + (h[3] if x == 'infinity' else h(x))*y != (f[6] if x == 'infinity' else f(x)):
                raise ValueError('divisor point is off the LMFDB model: '+row['label'])
        discriminant = h*h+4*f
        lhs = 'y^2'+(r'+\left(%s\right)y' % latex(h) if h else '')
        lmfdb_model = dict(f=str(f), h=str(h), discriminant=str(discriminant),
                          lhs_tex=lhs, rhs_tex=str(latex(f)),
                          equation_tex=r'C:\quad '+lhs+'='+str(latex(f)),
                          branch_values=[complex_coordinates(x) for x in
                                         discriminant.roots(QQbar, multiplicities=False)])
    identifiers = ['D%d' % (i+1) for i in indices]
    support = []
    for i, (x, y) in enumerate(originals):
        if not any(D[i] for D in coefficients):
            continue
        infinite = x == 'infinity'
        record = encode_point((x, y))
        rational = record['field_polynomial'] is None
        point_tex = ((r'\infty_{Y=%s}' % latex(y) if infinite else
                      r'\left(%s,%s\right)' % (latex(x), latex(y)))
                     if rational else r'P_{%d}' % (i+1))
        definition = r'P_{%d}=%s' % (i+1, point_tex)
        if not rational:
            R = PolynomialRing(QQ, 't')
            theta = r'\theta_{%d}' % (i+1)
            coordinate_tex = [latex(R(c)).replace('t', theta) for c in
                              (record['coordinates'][1:] if infinite else record['coordinates'])]
            coordinate = (r'\infty_{Y=%s}' % coordinate_tex[0] if infinite else
                          r'\left(%s\right)' % ','.join(coordinate_tex))
            definition = (r'P_{%d}=%s,\quad %s=0\;\text{(root index %d)}' %
                          (i+1, coordinate, latex(R(record['field_polynomial'])).replace('t', theta),
                           record['root_index']))
        u, v = images[i]
        integration_infinity = u == 'infinity'
        support.append(dict(id='P%d' % (i+1), index=i, point=encode_point((u, v)),
                            tex=point_tex, definition_tex=definition,
                            x=None if integration_infinity else complex_coordinates(u), y=complex_coordinates(v),
                            infinity=integration_infinity,
                            real=bool(v.imag() == 0 and (integration_infinity or u.imag() == 0)),
                            lmfdb=dict(point=record, x=None if infinite else complex_coordinates(x),
                                       y=complex_coordinates(y), infinity=infinite,
                                       real=bool(y.imag() == 0 and (infinite or x.imag() == 0))),
                            coefficients={name: str(D[i]) for name, D in zip(identifiers, coefficients)}))
    divisors = []
    for name, D in zip(identifiers, coefficients):
        terms = []
        for P in support:
            c = D[P['index']]
            if not c:
                continue
            body = ('' if abs(c) == 1 else str(latex(abs(c)))+r'\,')+P['tex']
            terms.append(('-' if c < 0 else '+' if terms else '')+body)
        divisors.append(dict(id=name, tex=r'D_{%s}=' % name[1:]+''.join(terms),
                             coefficients=list(map(str, D))))

    recorded = row.get('regulator_input', {}).get('centers')
    if recorded is not None:
        centers = [decode_algebraic(c) for c in recorded]
        center_source = 'computation'
    else:
        kernel = next((k for k in row.get('integral_kernels', []) if k['points'] == raw_points), None)
        cut_rows = (kernel['divisor_basis'] if row.get('settings', {}).get('ignore_integrality') else
                    kernel['basis']) if kernel else coefficients
        # Group by original x before applying the recorded coordinate change,
        # exactly as _divisor_evaluator does. Ordinate cancellation matters.
        original_x = [decode_point(P)[0] for P in raw_points]
        centers = []
        for x in dict.fromkeys(original_x):
            group = [i for i, value in enumerate(original_x) if value == x]
            if any(sum(QQ(D[i]) for i in group) for D in cut_rows):
                transformed = images[group[0]][0]
                if transformed != 'infinity':
                    centers.append(transformed)
        center_source = 'saved_basis' if kernel else 'selected_divisors'
    centers = list(dict.fromkeys(centers))
    center_indices = {x: i for i, x in enumerate(centers)}
    for P in support:
        P['center_index'] = None if P['infinity'] else center_indices.get(images[P['index']][0])
    logarithms = [dict(exact=encode_algebraic(x), position=complex_coordinates(x),
                       tex=latex(x) if x in QQ else r'a_{%d}' % (i+1),
                       label=str(x) if x in QQ else 'a%d' % (i+1)) for i, x in enumerate(centers)]
    matrix_ids = ['D%d' % (i+1) for i in row.get('basis_indices', [0, 1])]
    return dict(divisors=divisors, support=support, logarithmic_centers=logarithms,
                coordinate_system='lmfdb', lmfdb_model=lmfdb_model,
                centers_source=center_source, matrix_divisor_ids=matrix_ids,
                packet_torsion_order=packet_order(row), end_alg=endomorphism_type(row, curve_row),
                packet_size=len(raw_points))
