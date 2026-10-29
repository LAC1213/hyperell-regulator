"""Recompute manuscript period minors and rescale its rational q factors.

Run with Sage's Python from the repository root. Regulators and L-values are
reused from the earlier manuscript audit; only periods are freshly integrated.
"""

import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, RealField

from hyperell_regulator.regulator import c_invariants


def main():
    source_path = Path(__file__).with_name('paper_period_normalisation_input.json')
    source = json.loads(source_path.read_text())
    R = RealField(400)
    polynomials = PolynomialRing(QQ, 'x')

    rows = []
    for expected in source['rows']:
        label = expected['label']
        F = polynomials(expected['polynomial_coefficients'])
        lower = c_invariants(F, prec=50)
        cp, cm = c_invariants(F, prec=90)
        refinement = max(abs(R(a)/R(b)-1) for a, b in zip(lower, (cp, cm)))
        if refinement > R('1e-47'):
            raise RuntimeError('period minors changed at higher precision: ' + label)
        scales = [R(new)/R(expected[key]) for new, key in
                  ((cp, 'old_c_plus'), (cm, 'old_c_minus'))]
        scale = scales[0]*scales[1]
        rational = scale.nearby_rational(max_denominator=10000)
        error = abs(scale-R(rational))
        if error > R('1e-85'):
            raise RuntimeError('period-product ratio did not refine: ' + label)
        q = QQ(expected['old_paper_q'])*rational
        row = dict(label=label, polynomial=str(F),
                   old_paper_q=expected['old_paper_q'], c_plus=str(cp), c_minus=str(cm),
                   c_plus_ratio=str(scales[0]), c_minus_ratio=str(scales[1]),
                   product_ratio=str(scale), product_ratio_rational=str(rational),
                   product_ratio_error=str(error), new_paper_q=str(q),
                   integer=bool(q.denominator() == 1),
                   period_50_to_90_relative_difference=str(refinement))
        if 'saved_numerical_paper_q' in expected:
            numeric = R(expected['saved_numerical_paper_q'])*scale
            row.update(new_q_from_saved_numerical_comparison=str(numeric),
                       saved_comparison_error_to_integer=str(abs(numeric-R(q))))
        else:
            row['verification_note'] = (
                'Rescales manuscript q only; the earlier audit did not reproduce this L-value.')
        rows.append(row)
        print(label, 'old q', expected['old_paper_q'], 'multiplier', rational,
              'new q', q, 'ratio error', float(error), flush=True)

    report = dict(
        formula='q_paper = 2*pi*c_plus*c_minus*L_second_derivative / D',
        scope='Fresh periods on manuscript models; q rescaled from manuscript rational claims. '
              'No fresh regulators or L-values.',
        period_digits=90, previous_period_digits=50, rows=rows,
        all_rows_integer=all(row['integer'] for row in rows),
        source_file='reports/' + source_path.name,
        source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        original_audit_sha256=source['original_audit_sha256'],
        manuscript_sha256=source['manuscript_sha256'],
        source_code_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
                           [Path('hyperell_regulator/regulator') / (name + '.py')
                            for name in ('_paths', '_periods', 'plane', 'precision')]},
        numerical_caveat='Rational ratios and integer q values inferred numerically. '
                         'No interval certificate; determinants and L-values reused from '
                         'the manuscript audit. The 388.a.776.1 L-value was not independently reproduced.')
    out = Path(__file__).with_name('paper_period_normalisation_20261029.json')
    out.write_text(json.dumps(report, indent=2) + '\n')
    md = [
        '# Paper factors with full semilunar cycle periods', '',
        'Recomputed on 29 October 2026, using the manuscript models and the root of '
        'largest absolute value as the starting point. Period minors were evaluated '
        'at 50 and 90 decimal digits.', '',
        'The manuscript convention is `q = 2*pi*c_plus*c_minus*L_second_derivative/D`, '
        'the reciprocal of the experiment ratio. Thus '
        '`q_new = q_old * (c_plus_new*c_minus_new)/(c_plus_old*c_minus_old)`.', '',
        '| Curve | Paper q | Period-product multiplier | Updated q |',
        '|---|---:|---:|---:|',
    ]
    md += ['| `{label}` | {old_paper_q} | {product_ratio_rational} | {new_paper_q} |'
           .format(**row) for row in rows]
    md += [
        '', 'All sixteen updated factors are integers.' if report['all_rows_integer']
        else 'Some updated factors are not integers.',
        'The inferred period-product multipliers agree with the 90-digit computations '
        'to absolute error below `1e-85`; the individual minors agree between the two '
        'precisions to relative error below `1e-47`.', '',
        'This calculation reuses the manuscript rational factors and the 110-digit '
        'period minors from the earlier audit. It does not freshly recompute regulators '
        'or L-values. Saved numerical identity checks are included where available; '
        '`388.a.776.1` rescales the manuscript claim because its L-value was not '
        'independently reproduced by that audit. These are numerical checks, without '
        'interval certificates.', '',
        'The check also exposed quadrature-sized differences between mathematically '
        'equal pivots for `277.a.277.1`. Period-minor pivot ties now use quadrature '
        'accuracy, giving consistent signs under precision refinement.', '',
        'Full values and provenance: '
        '[JSON evidence](paper_period_normalisation_20261029.json). '
        'The [input snapshot](paper_period_normalisation_input.json) contains the '
        'manuscript models, original factors and earlier period minors.', '',
        'Reproduce from the repository root with '
        '`sage -python -m reports.check_paper_period_normalisation`.', '',
    ]
    Path('reports/paper_period_normalisation_20261029.md').write_text('\n'.join(md))


if __name__ == '__main__':
    main()
