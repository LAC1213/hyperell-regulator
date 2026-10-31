"""Regenerate the manuscript tables with the package's full cycle periods.

    sage -python -m reports.remake_paper_tables

Only periods are freshly integrated. The portable input snapshot records the
saved determinants for weights (1, Re(z), |z|^2), L-values, original models
and their precision/provenance. The middle row is doubled for the current
weights (1, 2 Re(z), |z|^2), and q = pi*c_plus*c_minus*L''/D.
The only output is one LaTeX file containing the three replacement tables.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(Path(path).read_text())


def significant(value, digits=10):
    with localcontext() as context:
        context.prec = max(180, digits + 20)
        return format(Decimal(str(value)), f'.{digits}g')


def tex_number(value, digits=10):
    text = significant(value, digits)
    if 'e' in text.lower():
        mantissa, exponent = text.lower().split('e')
        return rf'{mantissa}\times 10^{{{int(exponent)}}}'
    return text


def tex_rational(value):
    q = Fraction(value)
    if q.denominator == 1:
        return str(q.numerator)
    sign = '-' if q < 0 else ''
    return sign + rf'\frac{{{abs(q.numerator)}}}{{{q.denominator}}}'


def curve_link(label):
    return (r'\href{https://www.lmfdb.org/Genus2Curve/Q/' +
            '/'.join(label.split('.')) + r'}{\texttt{' + label + '}}')


def error_tex(row):
    error = row['relative_error']
    return r'\text{--}' if error is None else tex_number(error, 2)


def rounded_match(value, printed):
    with localcontext() as context:
        context.prec = 180
        target = Decimal(printed)
        unit = Decimal(1).scaleb(target.as_tuple().exponent)
        return abs(Decimal(str(value)) - target) <= unit / 2


def compute(source, precision, check_precision):
    from sage.all import PolynomialRing, QQ, RealField
    from hyperell_regulator.regulator import c_invariants

    real = RealField(int(3.33 * (precision + 25)) + 32)
    ring = PolynomialRing(QQ, 'x')
    rows = []
    for original in source['rows']:
        label = original['label']
        polynomial = ring(original['polynomial_coefficients'])
        for x, y in original['points']:
            if QQ(y)**2 != polynomial(QQ(x)):
                raise ArithmeticError('point is off the manuscript model: ' + label)
        lower = c_invariants(polynomial, prec=check_precision)
        cp, cm = c_invariants(polynomial, prec=precision)
        refinement = max(abs(real(a)/real(b)-1) for a, b in zip(lower, (cp, cm)))
        if refinement > real(10)**(-check_precision + 3):
            raise ArithmeticError('periods did not refine consistently: ' + label)
        ratios = [real(new)/real(original[key]) for new, key in
                  ((cp, 'old_c_plus'), (cm, 'old_c_minus'))]
        ratio = ratios[0] * ratios[1]
        multiplier = ratio.nearby_rational(max_denominator=4096)
        ratio_error = abs(ratio - real(multiplier))
        ratio_digits = min(precision, original['old_period_digits'])
        if ratio_error > real(10)**(-ratio_digits + 5):
            raise ArithmeticError('period product multiplier did not converge: ' + label)
        # Old q used 2*pi and an undoubled middle determinant row.
        q = QQ(original['old_paper_q']) * multiplier / 4
        if not rounded_match(original['determinant'], original['original_display']['det']):
            raise ArithmeticError('saved determinant disagrees with manuscript: ' + label)
        with localcontext() as context:
            context.prec = max(180, len(original['determinant']) + 5)
            determinant = str(2 * Decimal(original['determinant']))
        lvalue = original['L_second_derivative']
        error = None
        if lvalue is not None:
            if not rounded_match(lvalue, original['original_display']['L_second_derivative']):
                raise ArithmeticError('saved L-value disagrees with manuscript: ' + label)
            numeric_q = real.pi()*real(cp)*real(cm)*real(lvalue)/real(determinant)
            error = abs(numeric_q/real(q)-1)
            verified_digits = min(precision, original['determinant_digits'], original['lvalue_digits'])
            if error > real(10)**(-verified_digits + 3):
                raise ArithmeticError('rescaled identity exceeds input accuracy: ' + label)
        row = dict(original, determinant=determinant,
                   c_plus=cp.str(truncate=False), c_minus=cm.str(truncate=False),
                   new_paper_q=str(q),
                   relative_error=None if error is None else error.str(truncate=False))
        rows.append(row)
        print(label, 'q =', q, 'period multiplier =', multiplier,
              'identity error =', '--' if error is None else error.str(digits=4), flush=True)
    return rows


def render(rows, output, digits):
    preamble = [r'\begin{table}[htbp]', r'\centering\scriptsize',
                r'\setlength{\tabcolsep}{3pt}', r'\renewcommand{\arraystretch}{1.3}']
    torsion = preamble + [r'\begin{tabular}{@{}lcccrr@{}}', r'\toprule',
                         r'Curve & Point $P$ & Order & $\alpha_P$ integral? & $q$ & Relative error \\',
                         r'\midrule']
    models = preamble + [r'\begin{tabular}{@{}ll@{}}', r'\toprule',
                        r'Curve & $f(x)$ in $y^2=f(x)$ \\', r'\midrule']
    numerical = preamble + [r'\begin{tabular}{@{}lrrrrr@{}}', r'\toprule',
                           r'Curve & $D$ & $c^+$ & $c^-$ & $L^{\prime\prime}(M_{\mathrm{std}},1)$ & Relative error \\',
                           r'\midrule']
    for row in rows:
        label = curve_link(row['label'])
        pending = row['relative_error'] is None
        integrality_pending = any(p['same_component_everywhere'] is not True
                                  for p in row['integrality_component_verdicts'])
        marker = r'^{\dagger}' if integrality_pending else ''
        points = r'\begin{gathered}' + r'\\'.join(
            '(' + ','.join(tex_rational(a) for a in point) + ')' for point in row['points']) + r'\end{gathered}'
        orders = r'\begin{gathered}' + r'\\'.join(map(str, row['orders'])) + r'\end{gathered}'
        integral = r'\begin{gathered}\mathrm{Yes}\\\mathrm{Yes}\end{gathered}'
        if integrality_pending:
            integral += marker
        torsion.append(label + ' & $' + points + '$ & $' + orders + '$ & $' + integral +
                       '$ & $' + tex_rational(row['new_paper_q']) + marker + '$ & $' + error_tex(row) + r'$ \\')
        dagger = r'$^{\dagger}$' if row['reciprocal'] else ''
        models.append(label + dagger + ' & $' + row['polynomial_latex'] + r'$ \\')
        lvalue = row['L_second_derivative'] or row['original_display']['L_second_derivative']
        values = [tex_number(row[k], digits) for k in ('determinant', 'c_plus', 'c_minus')]
        values += [tex_number(lvalue, digits) + (r'^{\dagger}' if pending else ''), error_tex(row)]
        numerical.append(label + ' & ' + ' & '.join('$' + v + '$' for v in values) + r' \\')
    scope = (r'Errors use full saved values and freshly computed periods. Thirteen rows use '
             r'100-digit determinant and $L$-value inputs; \texttt{277.a.277.1} uses '
             r'15-digit determinant and 30-digit $L$-value inputs; \texttt{389.a.389.1} uses '
             r'30-digit inputs. \texttt{388.a.776.1} uses a 15-digit determinant and a '
             r'fresh 15-digit $L$-value, whose functional-equation error is below $10^{-23}$.')
    torsion += [r'\bottomrule\end{tabular}',
                r'\caption{Torsion point classes with the full cycle period convention, '
                r'$\pi c^+c^-L^{\prime\prime}=qD$. Orders are relative to $\mathcal O^+$ '
                r'on the displayed model (the unique infinity for a quintic, or '
                r'$y/x^3\to+1$ for a sextic). The fourth column records integrality. '
                r'For $\dagger$, integrality retains the manuscript claim; the package '
                r'component test at 2 remains undecided. ' + scope + '}',
                r'\label{tab:torsion-integrality}', r'\end{table}']
    models += [r'\bottomrule\end{tabular}',
               r'\caption{The sixteen manuscript models. For unmarked curves, '
               r'$x=x_{\mathrm L}$ and $y=2y_{\mathrm L}+h(x_{\mathrm L})$. '
               r'For curves marked $\dagger$, $x=1/x_{\mathrm L}$ and '
               r'$y=(2y_{\mathrm L}+h(x_{\mathrm L}))/x_{\mathrm L}^3$.}',
               r'\label{tab:curve-models}', r'\end{table}']
    numerical += [r'\bottomrule\end{tabular}',
                  r'\caption{Numerical data with full cycle periods. $D$ is the determinant '
                  r'of the real plane-integral matrix with weights '
                  r'$(d_{11},d_{12},d_{22})=(1,2\operatorname{Re}(z),|z|^2)$. The period columns '
                  r'are signed minors of the real and imaginary parts of the full cycle '
                  r'period matrix. Displayed numbers have ' + str(digits) +
                  r' significant digits; the error is '
                  r'$|\pi c^+c^-L^{\prime\prime}/(qD)-1|$. ' + scope + '}',
                  r'\label{tab:data}', r'\end{table}']
    fragments = dict(torsion=torsion, models=models, numerical=numerical)
    note = '% Generated by reports.remake_paper_tables; requires amsmath, booktabs, hyperref.\n'
    combined = note + '\n\n'.join('\n'.join(fragment) for fragment in fragments.values()) + '\n'
    output.write_text(combined)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=HERE/'paper_full_periods_input.json')
    parser.add_argument('--out', type=Path, default=HERE/'paper_full_periods.tex')
    parser.add_argument('--prec', type=int, default=110)
    parser.add_argument('--check-prec', type=int, default=50)
    parser.add_argument('--display-digits', type=int, default=10)
    args = parser.parse_args()
    if args.prec < 30 or not 15 <= args.check_prec < args.prec:
        parser.error('require prec >= 30 and 15 <= check-prec < prec')
    if not 1 <= args.display_digits <= args.prec:
        parser.error('display-digits must be between 1 and prec')
    source = load(args.input)
    rows = compute(source, args.prec, args.check_prec)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    render(rows, args.out, args.display_digits)
    print('Wrote', args.out)


if __name__ == '__main__':
    main()
