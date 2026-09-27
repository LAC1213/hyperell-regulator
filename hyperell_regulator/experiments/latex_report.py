r"""Export successful Beilinson reports without recomputing numerical integrals.

Usage: ``sage -python -m hyperell_regulator.experiments.latex_report input.json
--out reports/beilinson.tex``. A compact, self-contained JSON snapshot beside
the TeX retains the exact functions, models, and full-precision numbers.
"""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path

from sage.all import QQ, PolynomialRing, latex
from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.regulator import function_with_divisor


def ratio_error(ratio, rational):
    """Absolute error, not relative agreement or an integer-relation residual."""
    if ratio is None or rational is None:
        return None
    with localcontext() as ctx:
        ctx.prec = max(160, len(str(ratio)) + len(str(rational)) + 40)
        r = QQ(rational)
        return str(abs(Decimal(str(ratio)) - Decimal(int(r.numerator())) / Decimal(int(r.denominator()))))


def scientific(value, figures=8):
    if value is None:
        return r'\text{--}'
    with localcontext() as ctx:
        ctx.prec = max(160, len(str(value)) + 20)
        v = Decimal(str(value))
        if not v.is_finite():
            raise ValueError('non-finite report value')
        if not v:
            return '0'
        mantissa, exponent = format(v, '.%dE' % (figures-1)).split('E')
        exponent = int(exponent)
        if -3 <= exponent <= 3:
            return format(Decimal(mantissa).scaleb(exponent), 'f')
        return mantissa + r'\times 10^{' + str(exponent) + '}'


def point_tex(point):
    if point[0] == 'infinity':
        return r'\infty_{%s}' % latex(QQ(point[1]))
    return r'(%s,%s)' % (latex(QQ(point[0])), latex(QQ(point[1])))


def _rational_function_tex(A, B):
    left = latex(A.factor()) if A else ''
    right = (r'\left(%s\right)v' % latex(B.factor())) if B else ''
    return '+'.join(part for part in (left, right) if part) or '0'


def prepare_rows(records, cache=b.CACHE):
    """Recover exact algebraic data; numerical values are copied unchanged."""
    rows = [dict(r) for r in records if r.get('status') == 'ok']
    rows.sort(key=lambda r: (int(r['cond']), r['label']))
    need_cache = any('algebraic_data' not in r for r in rows)
    cached = {}
    if need_cache:
        cached = {r['label']: r for r in json.loads(Path(cache).read_text())}
    for row in rows:
        cmp = row.get('comparison') or {}
        row['absolute_ratio_error'] = ratio_error(cmp.get('ratio'), cmp.get('rational'))
        if 'algebraic_data' in row:
            continue
        print('Exact functions:', row['label'], flush=True)
        C = b.curve_from_row(cached[row['label']])
        triple = dict(row['triple'])
        for name in ('P', 'Q', 'R'):
            point = triple[name]
            triple[name] = (point[0] if point[0] == 'infinity' else QQ(point[0]), QQ(point[1]))
        G, images, shift, pole = b._model_for_triple(C, triple)
        S = PolynomialRing(QQ, 'u')
        if S(str(row['model']).replace('x', 'u')) != S(list(G)):
            raise ValueError('saved and reconstructed models differ for ' + row['label'])
        functions = {}
        for name, image in zip(('P', 'Q'), images):
            A, B = function_with_divisor(G, image, pole, int(triple['order_' + name]))
            # Put every model in u,v for unambiguous typesetting.
            A, B = S.fraction_field()(A), S.fraction_field()(B)
            functions[name] = dict(A=str(A), B=str(B), latex=_rational_function_tex(A, B))
        scale = QQ(1)
        if shift is not None:
            original = S(list(C.F))
            shifted = original(S.gen()+QQ(shift))
            inverse = S(list(shifted)[::-1]) * S.gen()**(2*C.genus()+2-shifted.degree())
            scale = QQ(S(list(G)).leading_coefficient()/inverse.leading_coefficient()).sqrt()
        row['algebraic_data'] = dict(model=str(S(list(G))), model_latex=latex(S(list(G))),
            original_f=str(C.f), original_h=str(C.h),
            image_points=[p if p == 'infinity' else [str(v) for v in p] for p in images],
            shift=None if shift is None else str(shift),
            scale=str(scale), pole=str(pole), functions=functions)
    return rows


def render(rows):
    lines = [r'\documentclass[10pt]{article}',
        r'\usepackage[a4paper,landscape,margin=15mm]{geometry}',
        r'\usepackage{amsmath,amssymb,booktabs,longtable}',
        r'\setlength{\tabcolsep}{5pt}', r'\begin{document}',
        r'\noindent\textbf{Successful Beilinson regulator computations.}',
        r'$q=\det R/(2\pi c^+c^-L^{\prime\prime}(\mathrm{std},1))$, '
        r'$r$ is the saved rational guess, and $\delta=|q-r|$. '
        r'Values are rounded for display; errors use the full saved decimals. '
        r'A successful regulator does not by itself certify a 100-digit comparison.',
        r'\begin{longtable}{lrrrrr}', r'\toprule',
        r'Curve & $\det R$ & $L^{\prime\prime}(\mathrm{std},1)$ & $q$ & $r$ & $\delta$ \\',
        r'\midrule\endhead']
    for row in rows:
        cmp = row.get('comparison') or {}
        rational = latex(QQ(cmp['rational'])) if cmp.get('rational') is not None else r'\text{--}'
        vals = [scientific(row.get('det')), scientific(row.get('L_second_derivative')),
                scientific(cmp.get('ratio'), 10), rational,
                scientific(row.get('absolute_ratio_error'), 3)]
        lines.append(row['label'] + ' & ' + ' & '.join('$'+v+'$' for v in vals) + r' \\')
    lines += [r'\bottomrule\end{longtable}']
    missing = [r for r in rows if not (r.get('comparison') or {}).get('ratio')]
    if missing:
        lines.append(r'\noindent Missing comparisons: ' + '; '.join(
            r['label'] + (' (unsupported Frobenius action)' if 'order > 2' in r.get('L_status', '')
                         else ' (L-value computation failed)') for r in missing) + '.')
    lines += [r'\smallskip',
        r'\noindent $P,Q,R$ below are in the original cached model $y^2+h(x)y=f(x)$; '
        r'$\infty_\epsilon$ retains the saved infinity label. '
        r'The functions are in the integration model $v^2=G(u)$ shown on each row, '
        r'and satisfy $(\phi_P)=n_P(P-R)$ and $(\phi_Q)=n_Q(Q-R)$ after the coordinate change. '
        r'For $s=\text{--}$, $(u,v)=(x,2y+h(x))$; otherwise '
        r'$(u,v)=((x-s)^{-1},d(2y+h(x))/(x-s)^3)$.',
        r'\scriptsize', r'\begin{longtable}{llllrr}', r'\toprule',
        r'Curve & $P$ & $Q$ & $R$ & $n_P$ & $n_Q$ \\', r'\midrule\endhead']
    for row in rows:
        triple, algebra = row['triple'], row['algebraic_data']
        lines.append(row['label'] + ' & ' + ' & '.join('$'+point_tex(triple[n])+'$' for n in ('P','Q','R'))
                     + ' & %s & %s' % (triple['order_P'],triple['order_Q']) + r' \\*')
        shift = r'\text{--}' if algebra['shift'] is None else latex(QQ(algebra['shift']))
        scale = latex(QQ(algebra.get('scale', 1)))
        lines.append(r'\multicolumn{6}{l}{$G(u)=' + algebra['model_latex']
                     + r';\quad s=' + shift + r',\ d=' + scale + r'$} \\*')
        for name in ('P','Q'):
            suffix = r' \\*' if name == 'P' else r' \\ \addlinespace'
            lines.append(r'\multicolumn{6}{l}{$\phi_' + name + '=' + algebra['functions'][name]['latex'] + r'$}' + suffix)
    lines += [r'\bottomrule\end{longtable}', r'\end{document}', '']
    return '\n'.join(lines)


def export(source, out, cache=b.CACHE):
    records = json.loads(Path(source).read_text())
    rows = prepare_rows(records, cache)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render(rows))
    # Keep only the data needed to reproduce this table, not progress traces.
    keys = ('label','cond','status','triple','model','regulator_digits','L_digits',
            'det','c_plus','c_minus','L_second_derivative','L_status','comparison',
            'absolute_ratio_error','algebraic_data')
    snapshot = [{k:r[k] for k in keys if k in r} for r in rows]
    out.with_suffix('.json').write_text(json.dumps(snapshot, indent=2, default=str)+'\n')
    print('%d successful rows: %s' % (len(rows), out), flush=True)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('--out', default='reports/beilinson.tex')
    parser.add_argument('--cache', default=b.CACHE)
    args = parser.parse_args(argv)
    export(args.source, args.out, args.cache)


if __name__ == '__main__':
    main()
