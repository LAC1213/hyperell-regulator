"""Export Beilinson survey divisors, packet orders and numerical comparisons."""

import argparse
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

from ._survey_divisors import (
    comparison_error, decimal_value, divisor_metadata, endomorphism_type,
)


BOUND = 1000
DIGITS = 100
EXCLUDED_DISPLAY_CONDUCTORS = frozenset({1116})
NORMALISATION = "pi c+ c- L''(std,1) / D"


def select_results(results, *, bound=BOUND, digits=DIGITS):
    """Select survey records at the requested precision, including failed attempts."""
    from .beilinson import REPORT_VERSION

    selected = []
    for original in results:
        settings = original.get("settings", {})
        if not (0 < int(original["cond"]) < bound
                and original.get("regulator_method") != "isogeny_transport"
                and settings.get("report_version") == REPORT_VERSION
                and (original.get("status") != "ok" or original.get("divisors") is not None)
                and settings.get("regulator_digits") == digits
                and settings.get("lvalue_digits") == digits
                and settings.get("with_lvalue") is True):
            continue
        row = deepcopy(original)
        if row.get("status") == "ok":
            if row.get("regulator_digits") != digits:
                raise ValueError("Successful regulator has the wrong precision: " + row["label"])
            comparison = row.get("comparison", {})
            if comparison and comparison.get("normalisation") != NORMALISATION:
                raise ValueError("Unexpected comparison convention: " + row["label"])
            if comparison and (int(row.get("L_digits", 0)) < digits
                               or row.get("L_second_derivative") is None):
                    raise ValueError("Comparison lacks a %d-digit L-value: %s" % (digits, row["label"]))
            ratio, rational = comparison.get("ratio"), comparison.get("rational")
            if ratio is not None and rational is not None:
                with localcontext() as ctx:
                    ctx.prec = 180
                    q = Fraction(rational)
                    difference = abs(Decimal(decimal_value(ratio)) - Decimal(q.numerator)/Decimal(q.denominator))
                row["absolute_ratio_error"] = str(difference)
        selected.append(row)
    labels = [row["label"] for row in selected]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate curve labels in the saved results")
    return sorted(selected, key=lambda row: (int(row["cond"]), row["label"]))


def algebraic_data(row, curve_row):
    """Reconstruct the recorded divisor integration chart from the LMFDB model."""
    from sage.all import PolynomialRing, QQ, latex
    from hyperell_regulator.components import curve_from_row
    from ..motivic import integration_chart

    curve = curve_from_row(curve_row)
    F, _, shift = integration_chart(curve.F.change_ring(QQ), (), row.get('shift'))
    ring = PolynomialRing(QQ, 'u')
    model = ring(F.list())
    if model != ring(re.sub(r'\bx\b', 'u', row['model'])):
        raise ValueError('Reconstructed integration model differs: ' + row['label'])
    return dict(model=str(model), model_latex=latex(model), original_f=str(curve.f),
                original_h=str(curve.h), shift=None if shift is None else str(shift),
                scale='1')


def number_tex(value, digits=10):
    if value is None:
        return r"\text{--}"
    with localcontext() as ctx:
        ctx.prec = 180
        text = format(Decimal(decimal_value(value)), f".{digits}g")
    if "e" in text.lower():
        mantissa, exponent = text.lower().split("e")
        return mantissa + r"\times 10^{" + str(int(exponent)) + "}"
    return text


def rational_tex(value):
    if value is None:
        return r"\text{--}"
    q = Fraction(str(value))
    if q.denominator == 1:
        return str(q.numerator)
    return ("-" if q < 0 else "") + rf"\frac{{{abs(q.numerator)}}}{{{q.denominator}}}"


def escape_tex(text):
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
                    "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
                    "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(replacements.get(c, c) for c in str(text))


def curve_link(label, changed_model=False):
    """LMFDB link, with a dagger for a reciprocal coordinate change."""
    url = "https://www.lmfdb.org/Genus2Curve/Q/" + "/".join(label.split("."))
    link = r"\href{" + url + r"}{\texttt{" + escape_tex(label) + "}}"
    return link + (r"$^{\dagger}$" if changed_model else "")


def stacked(values):
    return r"$\begin{gathered}" + r"\\".join(values) + r"\end{gathered}$"


def relation_tex(row):
    terms = row.get('linear_comparison', {}).get('terms', [])
    if len(terms) > 1:
        return '+'.join(rational_tex(term['coefficient'])+r'\det R(D_{%d},D_{%d})' %
                        tuple(i+1 for i in term['basis_indices']) for term in terms).replace('+-', '-')
    return rational_tex(row.get('comparison', {}).get('rational'))


def render(rows, source_hash, *, bound=BOUND, digits=DIGITS):
    """Three multipage tables, with exact selected divisors in LMFDB coordinates."""
    rows = [row for row in rows if int(row.get('cond', 0)) < bound
            and int(row.get('cond', 0)) not in EXCLUDED_DISPLAY_CONDUCTORS]
    ok = [row for row in rows if row.get('status') == 'ok']
    metadata = {row['label']: row.get('divisor_data') or divisor_metadata(row) for row in ok}
    lines = [r'\documentclass[10pt]{article}', r'\usepackage[a4paper,margin=12mm]{geometry}',
             r'\usepackage{amsmath,amssymb,booktabs,longtable,array,pdflscape}',
             r'\usepackage[hidelinks]{hyperref}',
             r'\newcommand{\entryseparator}{\noalign{\vskip3pt\hbox to\linewidth{\leaders\hbox{\kern2pt\rule{3pt}{0.3pt}\kern2pt}\hfill}\nobreak\vskip3pt\nobreak}}',
             r'\begin{document}', r'\begin{landscape}',
             rf'\noindent\textbf{{Beilinson survey: conductor $<{bound}$, {digits}-digit targets.}}',
             r'We use $\pi c^+(C)c^-(C)L^{\prime\prime}(M_{\mathrm{std}}(C),1)=q\det R(D_i,D_j)$,',
             r'where $R(D_i,D_j)$ is the plane-integral matrix with column factors $(L_{D_i},L_{D_j},1)$',
             r'and row weights $(1,2\operatorname{Re}(z),|z|^2)$.',
             r'For $D=\sum_P a_P P$, $L_D(z)=\sum_{P\text{ finite}}a_P\log|z-u(P)|$ in the chart $v^2=F(u)$.',
             r'The periods $c^+(C),c^-(C)$ are the signed maximal minors of $\operatorname{Re}H,\operatorname{Im}H$',
             r'for the full-cycle period matrix $H$ in this chart.',
             r'The divisors below are the actual coefficient rows of the regulator columns,',
             r'in the original LMFDB coordinates of Table~\ref{tab:curve-models}.',
             r'The associated classes $\alpha(D)$ are integral.',
             r'Packet order means $N=\mathop{\mathrm{lcm}}_P\operatorname{ord}(P-P_0)$ in the Jacobian,',
             r'not the number of packet points. A dash means the exact order was not recorded.',
             r'The endomorphism column is the LMFDB rational endomorphism algebra label;',
             r'\texttt{CM} does not by itself assert that the Jacobian is of CM type.',
             r'The column $N_{\mathrm{std}}$ gives the conductor of $L(M_{\mathrm{std}}(C),s)$.',
             r'\begin{longtable}{@{}p{.4\linewidth}>{\centering\arraybackslash}p{.2\linewidth}>{\raggedleft\arraybackslash}p{.16\linewidth}>{\raggedleft\arraybackslash}p{\dimexpr.24\linewidth-6\tabcolsep\relax}@{}}',
             r'\caption{Selected torsion-packet divisors and standard-motive conductors.}\label{tab:torsion-integrality}\\',
             r'\toprule Curve & $\mathrm{End}_{\mathbf Q}^0(J)$ & Packet $N$ & $N_{\mathrm{std}}$ \\',
             r'\midrule\endfirsthead',
             r'\toprule Curve & $\mathrm{End}_{\mathbf Q}^0(J)$ & Packet $N$ & $N_{\mathrm{std}}$ \\',
             r'\midrule\endhead', r'\bottomrule\endfoot']
    for index, row in enumerate(ok):
        if index:
            lines.append(r'\entryseparator')
        data = metadata[row['label']]
        lines.append(' & '.join([curve_link(row['label']),
            r'\texttt{'+escape_tex(data.get('end_alg') or '--')+'}',
            '$'+rational_tex(data.get('packet_torsion_order'))+'$',
            '$'+rational_tex(row.get('standard_conductor'))+'$'])+r' \\*')
        for divisor_index, divisor in enumerate(data['divisors']):
            ending = r' \\*' if divisor_index < len(data['divisors'])-1 else r' \\'
            lines.append(r'\multicolumn{4}{@{}l@{}}{$'+divisor['tex']+'$}'+ending)
    lines += [r'\end{longtable}', r'\clearpage',
              r'\begin{longtable}{@{}p{3.4cm}p{\dimexpr\linewidth-3.4cm-2\tabcolsep\relax}@{}}',
              r'\caption{LMFDB models, integration charts and exact algebraic support definitions.}\label{tab:curve-models}\\',
              r'\toprule Curve & $y^2+h(x)y=f(x)$; integration chart and LMFDB support definitions \\ \midrule\endfirsthead',
              r'\toprule Curve & $y^2+h(x)y=f(x)$; integration chart and LMFDB support definitions \\ \midrule\endhead',
              r'\bottomrule\endfoot']
    for index, row in enumerate(ok):
        if index:
            lines.append(r'\entryseparator')
        a, d = row['algebraic_data'], metadata[row['label']]
        model = re.sub(r'\bu\b', 'x', a['model_latex'])
        display_model = d.get('lmfdb_model')
        descriptions = ([display_model['equation_tex'], r'v^2='+a['model_latex']]
                        if display_model else [model])
        if a.get('shift') is not None:
            descriptions.append(r'(s,d)=('+rational_tex(a['shift'])+','+rational_tex(a.get('scale', '1'))+')')
        descriptions += [P['definition_tex'] for P in d['support'] if P['point']['field_polynomial'] is not None]
        lines.append(curve_link(row['label'], a.get('shift') is not None)+' & '+stacked(descriptions)+r' \\')
    lines += [r'\end{longtable}',
              r'\noindent All displayed divisor points use the original LMFDB equation $y^2+h(x)y=f(x)$.',
              r'The cut model uses $(u,v)=(x,2y+h(x))$ in unmarked rows.',
              r'Rows marked $\dagger$ use $u=(x-s)^{-1}$ and $v=d(2y+h(x))/(x-s)^3$.',
              r'The notation $\infty_{Y=\eta}$ means $Y=\lim_{x\to\infty}y/x^3=\eta$ on the LMFDB model. Algebraic root indices are zero-based,',
              r'ordered by real then imaginary part, and fix the exact embeddings retained in JSON.',
              r'\clearpage', r'\begin{longtable}{@{}lrrrrrr@{}}',
              r'\caption{Numerical data. The factor is a rational guess; a multi-determinant comparison is shown as a sum.}\label{tab:data}\\',
              r'\toprule Curve & $\det R$ & $c^+(C)$ & $c^-(C)$ & $L^{\prime\prime}(M_{\mathrm{std}}(C),1)$ & $q$ or comparison & Relative error \\ \midrule\endfirsthead',
              r'\toprule Curve & $\det R$ & $c^+(C)$ & $c^-(C)$ & $L^{\prime\prime}(M_{\mathrm{std}}(C),1)$ & $q$ or comparison & Relative error \\ \midrule\endhead',
              r'\bottomrule\endfoot']
    for row in ok:
        values = [number_tex(row.get(k), 5) for k in ('det', 'c_plus', 'c_minus', 'L_second_derivative')]
        values += [relation_tex(row), number_tex(comparison_error(row), 5)]
        lines.append(curve_link(row['label'])+' & '+' & '.join('$'+v+'$' for v in values)+r' \\')
    lines += [r'\end{longtable}', r'\noindent Decimal values are displayed to 5 significant digits.',
              r'Certified interval endpoints and analytic assumptions are retained in the JSON;',
              r'the error column uses an interval upper bound when available. Numerical targets alone are not certified bounds.',
              r'Numerical agreement does not prove exact rational equality.']
    issues = [row for row in rows if row.get('status') != 'ok' or not row.get('comparison')]
    if issues:
        lines.append(r'\par\noindent\textbf{Incomplete computations (full records in JSON).}')
        for row in issues:
            reason = row.get('error') or row.get('L_status') or row.get('periods_error') or row['status']
            lines.append(r'\par\noindent\texttt{'+escape_tex(row['label'])+'}: '+escape_tex(str(reason).partition('; state=')[0]))
    lines += [r'\par\noindent\texttt{Source SHA-256: '+source_hash+'}', r'\end{landscape}', r'\end{document}']
    return '\n'.join(lines)+'\n'


def export(source, output, cache, *, bound=BOUND, digits=DIGITS):
    raw = Path(source).read_bytes()
    rows = select_results(json.loads(raw), bound=bound, digits=digits)
    if not rows:
        raise ValueError('No current %d-digit results below conductor %d' % (digits, bound))
    curves = {r['label']: r for r in json.loads(Path(cache).read_text())} if Path(cache).exists() else {}
    for row in rows:
        curve_row = row.get('input_curve') or curves.get(row['label'])
        row['end_alg'] = endomorphism_type(row, curve_row)
        if row.get('status') == 'ok':
            if not row.get('algebraic_data'):
                if curve_row is None:
                    raise ValueError('Missing original curve equation for '+row['label'])
                row['algebraic_data'] = algebraic_data(row, curve_row)
            row['divisor_data'] = divisor_metadata(row, curve_row)
    checksum = hashlib.sha256(raw).hexdigest()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    readme = f"""# Beilinson survey: conductor < {bound}, {digits} digits

- [Saved data, exact divisor coefficients and interval endpoints](beilinson_100_digits.json).
- [Tables of divisors, packet torsion orders, endomorphism types and numerical comparisons](beilinson_100_digits.tex).

The snapshot contains {len(rows)} records and {sum(r.get('status') == 'ok' for r in rows)} completed regulators.
Failed attempts remain in the JSON. Q and CM refer to rational endomorphism algebra labels.
Divisors use the original LMFDB coordinates; the integration charts are recorded separately. Packet order is the exponent of point differences,
not packet cardinality. Unknown exact orders are left blank rather than inferred from Jacobian torsion.

Regulator determinants use weights `(1, 2 Re(z), |z|^2)` and full cycle periods.
Certified interval endpoints and their analytic assumptions are retained; numerical agreement does not prove equality.

Regenerate from a completed experiment without repeating integrations:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report --input {source} --bound {bound} --digits {digits}
```

Source SHA-256: `{checksum}`.
"""
    for name, contents in [('beilinson_100_digits.json', json.dumps(rows, indent=2, allow_nan=False)+'\n'),
                           ('beilinson_100_digits.tex', render(rows, checksum, bound=bound, digits=digits)),
                           ('README.md', readme)]:
        path = output/name
        temporary = path.with_suffix(path.suffix+'.tmp')
        temporary.write_text(contents)
        temporary.replace(path)
    print('Exported %d records to %s' % (len(rows), output))


def explorer_inputs(report, output, *, bound=BOUND, digits=DIGITS, extra_results=()):
    """Export exact divisor support and every center used on the shared cut."""
    raw = Path(report).read_bytes()
    curves = []
    records = [(row, False) for row in json.loads(raw)]
    additional = []
    for path in extra_results:
        contents = Path(path).read_bytes()
        saved = json.loads(contents)
        extra = saved['curves'] if isinstance(saved, dict) else saved
        records.extend((row, True) for row in extra)
        additional.append(dict(snapshot=str(path), sha256=hashlib.sha256(contents).hexdigest(),
                               labels=[row['label'] for row in extra]))
    for row, explicit in records:
        if (row.get('status') != 'ok' or not (0 < int(row['cond']) and (explicit or int(row['cond']) < bound))
                or int(row['cond']) in EXCLUDED_DISPLAY_CONDUCTORS):
            continue
        if any(row.get(k) is None for k in ('matrix', 'det', 'c_plus', 'c_minus', 'L_second_derivative')):
            continue
        data = dict(row.get('integral_basis_display') or row.get('divisor_data') or divisor_metadata(row))
        data['matrix_divisor_ids'] = (row.get('divisor_data') or divisor_metadata(row))['matrix_divisor_ids']
        comparison = row.get('comparison', {})
        error = comparison_error(row)
        verified = bool(comparison.get('rational')) and error is not None and Decimal(error) < Decimal(10)**(-digits+5)
        curves.append(dict(label=row['label'], cond=int(row['cond']),
            F=re.sub(r'\bu\b', 'x', row['algebraic_data']['model']), **data,
            q=comparison.get('rational') if verified else None,
            quotient=decimal_value(comparison.get('ratio')), rational_guess=comparison.get('rational'),
            numerical={k: decimal_value(row.get(k)) for k in ('det', 'c_plus', 'c_minus', 'L_second_derivative')},
            reference_matrix=[[decimal_value(v) for v in line] for line in row.get('matrix', [])],
            reference_intervals={k: row.get(k) for k in ('matrix', 'det', 'c_plus', 'c_minus', 'L_second_derivative')},
            L_num_coeffs=row.get('L_num_coeffs'), feq_check=row.get('feq_check'),
            L_certified=row.get('L_certified', False), L_analytic_assumptions=row.get('L_analytic_assumptions'),
            relative_agreement=error, regulator_digits=row['regulator_digits'], L_digits=row.get('L_digits'),
            comparison_verified=verified, certified=row.get('certified', False),
            integral_rank=row.get('integral_rank'), integral_rank_field=row.get('integral_rank_field', 'Q'),
            integral_rank_space=row.get('integral_rank_space', 'galois_traces'),
            trace_rank_certificate=row.get('trace_rank_certificate'),
            integrality_verified=row.get('integral'),
            linear_comparison=row.get('linear_comparison'), original_curve=row.get('input_curve')))
    labels = [curve['label'] for curve in curves]
    if len(labels) != len(set(labels)):
        raise ValueError('Duplicate curve labels in explorer sources')
    snapshot = dict(schema_version=3, source=dict(title=f'Beilinson survey: conductor < {bound}, {digits}-digit targets',
        input_snapshot=str(output), report_snapshot=str(report), report_sha256=hashlib.sha256(raw).hexdigest(),
        conductor_bound=bound, regulator_digits=digits, lvalue_digits=digits,
        excluded_conductors=sorted(EXCLUDED_DISPLAY_CONDUCTORS)), curves=curves)
    if additional:
        snapshot['source']['additional_results'] = additional
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(snapshot, indent=2, allow_nan=False)+'\n')
    return snapshot


def main(argv=None):
    from .beilinson import CACHE, RESULTS_OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=RESULTS_OUT)
    parser.add_argument("--cache", type=Path, default=CACHE)
    parser.add_argument("--out", type=Path, default=Path("reports"))
    parser.add_argument("--explorer-inputs", type=Path,
                        help="also write an explorer input catalog from the new report")
    parser.add_argument('--extra-results', type=Path, action='append', default=[],
                        help='include a saved single-curve result collection in the explorer')
    parser.add_argument('--bound', type=int, default=BOUND)
    parser.add_argument('--digits', type=int, default=DIGITS)
    args = parser.parse_args(argv)
    if args.extra_results and not args.explorer_inputs:
        parser.error('--extra-results requires --explorer-inputs')
    export(args.input, args.out, args.cache, bound=args.bound, digits=args.digits)
    if args.explorer_inputs:
        explorer_inputs(args.out/"beilinson_100_digits.json", args.explorer_inputs, bound=args.bound, digits=args.digits,
                        extra_results=args.extra_results)


if __name__ == "__main__":
    main()
