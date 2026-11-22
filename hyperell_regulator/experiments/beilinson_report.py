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
                and settings.get("report_version") in (9, 10, 11, 12, REPORT_VERSION)
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


def _parse_point(v):
    r"""Read an affine or labelled infinite point from JSON."""
    if v is None:
        return None
    from sage.all import QQ
    a, b = v
    if a == "infinity":
        return ("infinity", QQ(b))
    return (QQ(a), QQ(b))

def _saved_model_with_pole(C, R, points):
    r"""A model with common pole ``R`` at infinity.

    Return ``(G, images, xR)``, using `u=1/(x-x(R))` if needed.
    ``images`` contains affine pairs or ``"infinity"``.
    """
    from sage.all import PolynomialRing, QQ, ZZ, lcm
    g = C.genus()
    RQ = PolynomialRing(QQ, 'x')
    hQ, FQ = RQ(C.h), RQ(C.F)

    def Ycoord(P):
        return 2 * P[1] + hQ(P[0])

    if R[0] == "infinity":
        images = []
        for P in points:
            images.append("infinity" if P[0] == "infinity"
                          else (P[0], Ycoord(P)))
        return C.F, images, None

    xR = QQ(R[0])
    SQ = PolynomialRing(QQ, 'u')
    u = SQ.gen()
    shifted = FQ(RQ.gen() + xR)
    G = sum(QQ(c) * u ** (2 * g + 2 - i) for i, c in enumerate(list(shifted)))
    d = lcm([QQ(c).denominator() for c in G] + [1])
    G = (d ** 2 * G).change_ring(ZZ)

    images = []
    for P in points:
        if P[0] == "infinity":
            # x = oo -> u = 0, v = d*(2 Y_oo + h_{g+1})
            images.append((QQ(0), d * (2 * QQ(P[1]) + QQ(C.h[g + 1]))))
        elif QQ(P[0]) == xR:
            images.append("infinity")
        else:
            t = QQ(P[0]) - xR
            images.append((1 / t, d * Ycoord(P) / t ** (g + 1)))
    return G, images, xR

def _saved_finite_model(C, points):
    """A rational model on which all marked points, including ``R``, are affine."""
    from sage.all import PolynomialRing, QQ
    if all(P[0] != "infinity" for P in points):
        H = PolynomialRing(QQ, "x")(C.h)
        return C.F, [(P[0], 2*P[1]+H(P[0])) for P in points], None
    # Put the new infinity fibre at a rational x-value containing none of
    # the marked points or Weierstrass points. Original infinities go to u=0.
    RQ = PolynomialRing(QQ, "x")
    F, H = RQ(C.F), RQ(C.h)
    marked = {QQ(P[0]) for P in points if P[0] != "infinity"}
    shift = QQ(0)
    while shift in marked or F(shift) == 0:
        shift += 1
    u = PolynomialRing(QQ, "u").gen()
    d = 2*C.genus()+2
    G = sum(c*u**(d-i) for i,c in enumerate(F(RQ.gen()+shift).list()))
    images = []
    for P in points:
        if P[0] == "infinity":
            images.append((QQ(0), 2*QQ(P[1])+QQ(C.h[C.genus()+1])))
        else:
            t = QQ(P[0])-shift
            images.append((1/t, (2*P[1]+H(P[0]))/t**(C.genus()+1)))
    return G, images, shift

def _saved_triple_model(C, t):
    points = [t["P"], t["Q"], t["R"]]
    G, images, shift = _saved_model_with_pole(C, t["R"], points)
    if images[0] != "infinity" and images[1] != "infinity":
        return G, images, shift, "infinity"
    G, images, shift = _saved_finite_model(C, points)
    return G, images, shift, images[2]


def algebraic_data(row, curve_row):
    """Reconstruct the exact integration model and marked coordinates."""
    from sage.all import PolynomialRing, QQ, latex
    from hyperell_regulator.components import curve_from_row

    curve = curve_from_row(curve_row)
    if row.get('divisors') is not None:
        from ..motivic import integration_chart

        F, _, shift = integration_chart(curve.F.change_ring(QQ), (), row.get('shift'))
        ring = PolynomialRing(QQ, 'u')
        model = ring(F.list())
        if model != ring(re.sub(r'\bx\b', 'u', row['model'])):
            raise ValueError('Reconstructed integration model differs: ' + row['label'])
        return dict(model=str(model), model_latex=latex(model), original_f=str(curve.f),
                    original_h=str(curve.h), shift=None if shift is None else str(shift),
                    scale='1')
    triple = dict(row["triple"])
    for key in ("P", "Q", "R"):
        triple[key] = _parse_point(triple[key])
    model, images, shift, pole = _saved_triple_model(curve, triple)
    ring = PolynomialRing(QQ, "u")
    if ring(model) != ring(re.sub(r"\bx\b", "u", row["model"])):
        raise ValueError("Reconstructed integration model differs: " + row["label"])
    scale = QQ(1)
    if shift is not None and pole == "infinity":
        from sage.all import lcm
        shifted = curve.F.change_ring(QQ)(curve.F.parent().gen()+shift)
        scale = lcm([c.denominator() for c in shifted] + [1])
    def point_json(point):
        return point if point == "infinity" else [str(v) for v in point]
    return dict(model=str(ring(model)), model_latex=latex(ring(model)),
                original_f=str(curve.f), original_h=str(curve.h),
                image_points=[point_json(p) for p in images],
                shift=None if shift is None else str(shift), scale=str(scale),
                pole=point_json(pole))


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


def point_tex(point):
    if point == "infinity":
        return r"\infty"
    if point[0] == "infinity":
        return r"\infty_{" + rational_tex(point[1]) + "}"
    return "(" + ",".join(rational_tex(v) for v in point) + ")"


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


def integrality_tex(verdict):
    return r"\text{" + ("Yes" if verdict is True else "No" if verdict is False else "Undecided") + "}"


def model_pole_tex(row):
    """The common pole in the integration coordinates, in LaTeX."""
    from sage.all import PolynomialRing, QQ

    data = row["algebraic_data"]
    if data["pole"] != "infinity":
        return point_tex(data["pole"])
    ring = PolynomialRing(QQ, "x")
    model = ring(re.sub(r"\bu\b", "x", data["model"]))
    if model.degree() % 2:
        return r"\infty"
    h = ring(data["original_h"])
    pole = row["triple"]["R"]
    if pole[0] == "infinity":
        eta = 2*QQ(pole[1]) + h[model.degree()//2]
    else:
        eta = QQ(data["scale"])*(2*QQ(pole[1]) + h(QQ(pole[0])))
    return point_tex(("infinity", str(eta)))


def relation_tex(row):
    terms = row.get('linear_comparison', {}).get('terms', [])
    if len(terms) > 1:
        return '+'.join(rational_tex(term['coefficient'])+r'\det(D_{%d},D_{%d})' %
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
             r'\usepackage[hidelinks]{hyperref}', r'\begin{document}', r'\begin{landscape}',
             rf'\noindent\textbf{{Beilinson survey: conductor $<{bound}$, {digits}-digit targets.}}',
             r'We use $q=\pi c^+c^-L^{\prime\prime}(\mathrm{std},1)/D$ and full cycle periods,',
             r'with determinant weights $(1,2\operatorname{Re}(z),|z|^2)$.',
             r'The divisors below are the actual coefficient rows of the regulator columns,',
             r'in the original LMFDB coordinates of Table~\ref{tab:curve-models}.',
             r'Packet order means $N=\mathop{\mathrm{lcm}}_P\operatorname{ord}(P-P_0)$ in the Jacobian,',
             r'not the number of packet points. A dash means the exact order was not recorded.',
             r'The endomorphism column is the LMFDB rational endomorphism algebra label;',
             r'\texttt{CM} does not by itself assert that the Jacobian is of CM type.',
             r'\begin{longtable}{@{}lcrp{16cm}c@{}}',
             r'\caption{Selected torsion-packet divisors and integrality.}\label{tab:torsion-integrality}\\',
             r'\toprule Curve & $\mathrm{End}_{\mathbf Q}^0(J)$ & Packet $N$ & Selected divisors & Integral? \\',
             r'\midrule\endfirsthead',
             r'\toprule Curve & $\mathrm{End}_{\mathbf Q}^0(J)$ & Packet $N$ & Selected divisors & Integral? \\',
             r'\midrule\endhead', r'\bottomrule\endfoot']
    for row in ok:
        data = metadata[row['label']]
        verdicts = (None if len(row.get('linear_comparison', {}).get('terms', [])) > 1
                    else row.get('divisor_integrality'))
        if 'triple' in row:
            verdicts = [row['triple'].get(k) for k in ('verdict_P', 'verdict_Q')]
        integral = stacked([integrality_tex(v) for v in verdicts]) if verdicts else '$'+integrality_tex(row.get('integral'))+'$'
        lines.append(' & '.join([curve_link(row['label']),
            r'\texttt{'+escape_tex(data.get('end_alg') or '--')+'}',
            '$'+rational_tex(data.get('packet_torsion_order'))+'$',
            stacked([D['tex'] for D in data['divisors']]), integral])+r' \\')
    lines += [r'\end{longtable}', r'\clearpage',
              r'\begin{longtable}{@{}lp{22cm}@{}}',
              r'\caption{LMFDB models, integration charts and exact algebraic support definitions.}\label{tab:curve-models}\\',
              r'\toprule Curve & $y^2+h(x)y=f(x)$; integration chart and LMFDB support definitions \\ \midrule\endfirsthead',
              r'\toprule Curve & $y^2+h(x)y=f(x)$; integration chart and LMFDB support definitions \\ \midrule\endhead',
              r'\bottomrule\endfoot']
    for row in ok:
        a, d = row['algebraic_data'], metadata[row['label']]
        model = re.sub(r'\bu\b', 'x', a['model_latex'])
        display_model = d.get('lmfdb_model')
        descriptions = ([display_model['equation_tex'], r'\text{Cut model: }v^2='+re.sub(r'\bu\b', 'u', a['model_latex'])]
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
              r'\toprule Curve & $D$ & $c^+$ & $c^-$ & $L^{\prime\prime}(M_{\mathrm{std}},1)$ & $q$ or comparison & Relative error \\ \midrule\endfirsthead',
              r'\toprule Curve & $D$ & $c^+$ & $c^-$ & $L^{\prime\prime}(M_{\mathrm{std}},1)$ & $q$ or comparison & Relative error \\ \midrule\endhead',
              r'\bottomrule\endfoot']
    for row in ok:
        values = [number_tex(row.get(k)) for k in ('det', 'c_plus', 'c_minus', 'L_second_derivative')]
        values += [relation_tex(row), number_tex(comparison_error(row), 2)]
        lines.append(curve_link(row['label'])+' & '+' & '.join('$'+v+'$' for v in values)+r' \\')
    lines += [r'\end{longtable}', r'\noindent Displayed values have 10 significant digits.',
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


def explorer_inputs(report, output, *, bound=BOUND, digits=DIGITS):
    """Export exact divisor support and every center used on the shared cut."""
    raw = Path(report).read_bytes()
    curves = []
    for row in json.loads(raw):
        if (row.get('status') != 'ok' or not 0 < int(row['cond']) < bound
                or int(row['cond']) in EXCLUDED_DISPLAY_CONDUCTORS):
            continue
        if any(row.get(k) is None for k in ('matrix', 'det', 'c_plus', 'c_minus', 'L_second_derivative')):
            continue
        data = row.get('divisor_data') or divisor_metadata(row)
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
            integrality_verified=row.get('integral', all(row.get('triple', {}).get(k) is True for k in ('verdict_P', 'verdict_Q'))),
            linear_comparison=row.get('linear_comparison'), original_curve=row.get('input_curve')))
    snapshot = dict(schema_version=3, source=dict(title=f'Beilinson survey: conductor < {bound}, {digits}-digit targets',
        input_snapshot=str(output), report_snapshot=str(report), report_sha256=hashlib.sha256(raw).hexdigest(),
        conductor_bound=bound, regulator_digits=digits, lvalue_digits=digits,
        excluded_conductors=sorted(EXCLUDED_DISPLAY_CONDUCTORS)), curves=curves)
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
    parser.add_argument('--bound', type=int, default=BOUND)
    parser.add_argument('--digits', type=int, default=DIGITS)
    args = parser.parse_args(argv)
    export(args.input, args.out, args.cache, bound=args.bound, digits=args.digits)
    if args.explorer_inputs:
        explorer_inputs(args.out/"beilinson_100_digits.json", args.explorer_inputs, bound=args.bound, digits=args.digits)


if __name__ == "__main__":
    main()
