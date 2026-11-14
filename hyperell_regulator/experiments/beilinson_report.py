"""Export the completed 100-digit survey with conductor strictly below 1300.

    sage -python -m hyperell_regulator.experiments.beilinson_report

This reads saved results without repeating numerical integrations or L-values.
Exact integration models and marked coordinates are reconstructed. Failures and
component verdicts remain in the JSON and are described separately in the TeX.
"""

import argparse
from collections import Counter
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re


BOUND = 1300
DIGITS = 100
EXCLUDED_DISPLAY_CONDUCTORS = frozenset({1116})
NORMALISATION = "pi c+ c- L''(std,1) / D"


def select_results(results):
    """Keep the current 100-digit policy, including unsuccessful attempts."""
    from .beilinson import REPORT_VERSION

    selected = []
    for original in results:
        settings = original.get("settings", {})
        if not (0 < int(original["cond"]) < BOUND
                and settings.get("report_version") in (9, 10, 11, REPORT_VERSION)
                and settings.get("regulator_digits") == DIGITS
                and settings.get("lvalue_digits") == DIGITS
                and settings.get("with_lvalue") is True):
            continue
        row = deepcopy(original)
        if row.get("status") == "ok":
            if row.get("regulator_digits") != DIGITS:
                raise ValueError("Successful regulator has the wrong precision: " + row["label"])
            comparison = row.get("comparison", {})
            if comparison and comparison.get("normalisation") != NORMALISATION:
                raise ValueError("Unexpected comparison convention: " + row["label"])
            if comparison and (int(row.get("L_digits", 0)) < DIGITS
                               or not isinstance(row.get("L_second_derivative"), str)):
                raise ValueError("Comparison lacks a 100-digit L-value: " + row["label"])
            ratio, rational = comparison.get("ratio"), comparison.get("rational")
            if ratio is not None and rational is not None:
                with localcontext() as ctx:
                    ctx.prec = 180
                    q = Fraction(rational)
                    difference = abs(Decimal(ratio) - Decimal(q.numerator)/Decimal(q.denominator))
                row["absolute_ratio_error"] = str(difference)
        selected.append(row)
    labels = [row["label"] for row in selected]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate curve labels in the saved results")
    return sorted(selected, key=lambda row: (int(row["cond"]), row["label"]))


def algebraic_data(row, curve_row):
    """Reconstruct exact integration coordinates, without quadrature."""
    from sage.all import PolynomialRing, QQ, latex
    from hyperell_regulator.components import curve_from_row
    from . import beilinson as b

    curve = curve_from_row(curve_row)
    triple = dict(row["triple"])
    for key in ("P", "Q", "R"):
        triple[key] = b._parse_point(triple[key])
    model, images, shift, pole = b._model_for_triple(curve, triple)
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
        text = format(Decimal(str(value)), f".{digits}g")
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
    """Label the common pole on the integration model, rather than guessing a sign."""
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


def render(rows, source_hash):
    """Three publication tables using the saved values and integration coordinates."""
    rows = [row for row in rows if int(row.get("cond", 0)) not in EXCLUDED_DISPLAY_CONDUCTORS]
    ok = [row for row in rows if row.get("status") == "ok"]
    lines = [r"\documentclass[10pt]{article}",
             r"\usepackage[a4paper,margin=12mm]{geometry}",
             r"\usepackage{amsmath,amssymb,booktabs,arydshln}",
             r"\usepackage[hidelinks]{hyperref}",
             r"\begin{document}",
             r"\noindent\textbf{Beilinson survey: conductor $<1300$, 100-digit targets.}",
             r"Conductor $1116$ is retained in the JSON and omitted from these tables.",
             r"We use $q=\pi c^+c^-L^{\prime\prime}(\mathrm{std},1)/D$, full cycle periods,",
             r"and determinant weights $(1,2\operatorname{Re}(z),|z|^2)$.",
             r"Displayed decimals have 10 significant digits; errors use the full saved values.",
             r"\begin{table}[htbp]",
             r"\centering\small",
             r"\renewcommand{\arraystretch}{1.5}",
             r"\setlength{\tabcolsep}{6pt}",
             r"\begin{tabular}{@{}cccccc@{}}",
             r"\toprule",
             r"Curve & Point $P$ & Order of $P$ & $\alpha_P$ integral? & Factor $q$ & Error \\",
             r"\midrule"]
    for index, row in enumerate(ok):
        triple, data = row["triple"], row["algebraic_data"]
        comparison = row.get("comparison", {})
        cells = [
            curve_link(row["label"]),
            stacked([point_tex(point) for point in data["image_points"][:2]]),
            stacked([rational_tex(triple[key]) for key in ("order_P", "order_Q")]),
            stacked([integrality_tex(triple.get(key)) for key in ("verdict_P", "verdict_Q")]),
            "$" + rational_tex(comparison.get("rational")) + "$",
            "$" + number_tex(comparison.get("relative_agreement"), 2) + "$",
        ]
        lines.append(" & ".join(cells) + r" \\")
        if index + 1 < len(ok):
            lines.append(r"\hdashline")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\caption{Torsion point classes for " + str(len(ok)) +
        r" genus 2 curves. Each stacked cell lists the two marked points in the order of the regulator columns,",
        r"in the integration coordinates of Table~\ref{tab:curve-models}.",
        r"The order of a displayed point is the order of $[P]-[R]$, with common pole $R$ given there.",
        r"The integrality column records the component criterion: Yes, No, or Undecided;",
        r"Undecided makes no integrality assertion.",
        r"Factor $q$ is the saved rational guess. The error is",
        r"$|\pi c^+c^-L^{\prime\prime}(M_{\mathrm{std}},1)/(qD)-1|$,",
        r"computed from unrounded values. Numerical targets are not certified error bounds.}",
        r"\label{tab:torsion-integrality}",
        r"\end{table}",
        r"\clearpage",
        r"\begin{table}[htbp]\centering\small",
        r"\renewcommand{\arraystretch}{1.35}",
        r"\setlength{\tabcolsep}{6pt}",
        r"\begin{tabular}{@{}ll@{}}",
        r"\toprule",
        r"Curve & $f(x)$ in $y^2=f(x)$; common pole $R$ \\",
        r"\midrule"]
    for row in ok:
        data = row["algebraic_data"]
        model = re.sub(r"\bu\b", "x", data["model_latex"])
        description = model + r";\quad R=" + model_pole_tex(row)
        changed = data["shift"] is not None
        if changed:
            description += r";\quad(s,d)=(" + rational_tex(data["shift"]) + "," + rational_tex(data["scale"]) + ")"
        lines.append(curve_link(row["label"], changed) + " & $" + description + r"$ \\")
    lines += [
        r"\bottomrule\end{tabular}",
        r"\caption{Integration models for the " + str(len(ok)) +
        r" displayed curves; points in Table~\ref{tab:torsion-integrality} use these coordinates.",
        r"Write the LMFDB equation as $y_{\mathrm L}^2+h(x_{\mathrm L})y_{\mathrm L}=f_{\mathrm L}(x_{\mathrm L})$.",
        r"For unmarked rows, $x=x_{\mathrm L}$ and $y=2y_{\mathrm L}+h(x_{\mathrm L})$.",
        r"Rows marked $\dagger$ use $x=(x_{\mathrm L}-s)^{-1}$ and",
        r"$y=d(2y_{\mathrm L}+h(x_{\mathrm L}))/(x_{\mathrm L}-s)^3$,",
        r"with the displayed $(s,d)$. In odd degree the point at infinity is unique;",
        r"in even degree $\infty_\eta$ denotes $y/x^3\to\eta$.}",
        r"\label{tab:curve-models}",
        r"\end{table}",
        r"\clearpage",
        r"\begin{table}[htbp]",
        r"\centering",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\renewcommand{\arraystretch}{1.3}",
        r"\begin{tabular}{@{}lrrrrr@{}}",
        r"\toprule",
        r"Curve & $D$ & $c^+$ & $c^-$ & $L^{\prime\prime}(M_{\mathrm{std}},1)$ & Relative error \\",
        r"\midrule"]
    for row in ok:
        values = [number_tex(row.get(key)) for key in ("det", "c_plus", "c_minus", "L_second_derivative")]
        values.append(number_tex(row.get("comparison", {}).get("relative_agreement"), 2))
        lines.append(curve_link(row["label"]) + " & " +
                     " & ".join("$" + value + "$" for value in values) + r" \\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\caption{Numerical data for Table~\ref{tab:torsion-integrality}.",
        r"Real values are rounded to 10 significant digits; relative errors are calculated from the unrounded data.",
        r"The columns $c^+$ and $c^-$ use full cycle periods.}",
        r"\label{tab:data}",
        r"\end{table}"]
    issues = [row for row in rows if row.get("status") not in ("ok", "no_admissible_triple")
              or (row.get("status") == "ok" and not row.get("comparison"))]
    if issues:
        lines.append(r"\par\noindent\textbf{Incomplete computations (full records in JSON).}")
        for row in issues:
            reason = row.get("error") or row.get("L_error") or row.get("L_status") or row.get("periods_error") or row["status"]
            reason = str(reason).partition("; state=")[0]
            lines.append(r"\par\noindent\texttt{" + escape_tex(row["label"]) + "}: " + escape_tex(reason))
    lines += [r"\par\noindent Source: \texttt{results/beilinson\_results.json}.",
              r"\par\noindent\texttt{SHA-256: " + source_hash + "}", r"\end{document}"]
    return "\n".join(lines) + "\n"


def export(source, output, cache):
    raw = Path(source).read_bytes()
    rows = select_results(json.loads(raw))
    if not rows:
        raise ValueError("No current 100-digit results below conductor 1300")
    needs_models = any(row.get('status') == 'ok' and not row.get('algebraic_data') for row in rows)
    curves = {r["label"]: r for r in json.loads(Path(cache).read_text())} if needs_models else {}
    for row in rows:
        if row.get("status") == "ok" and not row.get('algebraic_data'):
            row["algebraic_data"] = algebraic_data(row, curves[row["label"]])
    checksum = hashlib.sha256(raw).hexdigest()
    tex = render(rows, checksum)
    counts = Counter(row["status"] for row in rows)
    compared = sum(bool(row.get("comparison")) for row in rows)
    displayed = sum(row.get("status") == "ok" and int(row["cond"]) not in EXCLUDED_DISPLAY_CONDUCTORS
                    for row in rows)
    source_name = Path(source)
    if source_name.is_absolute() and source_name.is_relative_to(Path.cwd()):
        source_name = source_name.relative_to(Path.cwd())
    readme = f"""# Beilinson survey: conductor < 1300, 100 digits

This directory contains the current 100-digit Beilinson survey:

- [Full saved data and integration models](beilinson_100_digits.json).
- [Three LaTeX tables: torsion/integrality, models, and numerical data](beilinson_100_digits.tex).

The snapshot has {len(rows)} curve records, {counts.get('ok', 0)} successful
regulators and {compared} numerical comparisons. Unsuccessful attempts and their
reasons are retained in the JSON. The three tables display completed
comparisons. This is a survey of locally cached geometrically simple
curves, without a claim of completeness over all curves.
The conductor-1116 record is retained in the JSON, but omitted from the LaTeX
tables and explorer.
The LaTeX tables and explorer include {displayed} successful regulators.

Regulators and L-values target 100 decimal digits. The determinant uses weights
`(1, 2 Re(z), |z|^2)`, full cycle periods, and
`q = pi*c_plus*c_minus*L_second_derivative/D`. The comparison error is
`abs(q/rational - 1)`. An undecided component verdict is displayed as
`Undecided` in the integrality column. Precision targets are not certified
error bounds.

Full cycle periods use the MST cycle basis recorded in each successful row.
The current code computes D and the signed period minors at a 100-digit target.
L-values use the precision-validated coefficient cache. The quotient uses
these ingredients in the same integration model as its displayed points.

The numerical table displays 10 significant digits. The JSON preserves the
complete saved decimal strings, matrices, precision settings, component
verdicts and integration models. Table points use the displayed integration
models. Reciprocal coordinate changes are marked by a dagger in the model
table, with their shift and scale given explicitly. The JSON retains both
the original points and their integration-model images.

From the repository root in SageMath, recompute the survey (network access is
needed to fetch curve data; cached L-values are reused after validation):

```sh
sage -python -m hyperell_regulator.experiments.beilinson 1300 --fetch --force --digits 100 --reg-digits 100 --edge-seconds 600
sage -python -m hyperell_regulator.experiments.beilinson_report
```

The bound is strict: conductor < 1300. A fresh cache also recomputes L-values
and can take substantially longer. Component/reduction computations require
the package's Sage dependencies. The command retains failed attempts and
does not claim that the locally cached curves form a complete list.

To rebuild the tables from the versioned JSON without repeating numerical
calculations or requiring local curve caches:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report --input reports/beilinson_100_digits.json
```

To typeset the survey, install the LaTeX packages in its preamble, including
`arydshln`. Compile into a separate directory to keep this folder tidy.
Run the final command twice to resolve cross-references:

```sh
mkdir -p /tmp/beilinson-survey
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/beilinson-survey reports/beilinson_100_digits.tex
```

Source: `{source_name}` (read without modification).
Source SHA-256: `{checksum}`.
"""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for filename, contents in (("beilinson_100_digits.json", json.dumps(rows, indent=2, allow_nan=False)+"\n"),
                               ("beilinson_100_digits.tex", tex), ("README.md", readme)):
        path = output/filename
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(contents)
        temporary.replace(path)
    print(f"Exported {len(rows)} records, {counts.get('ok', 0)} regulators, {compared} comparisons to {output}")


def explorer_inputs(report, output):
    """Use the report's exact models, marked points and 100-digit references."""
    raw = Path(report).read_bytes()
    rows = json.loads(raw)
    curves = []
    x_expression = lambda value: re.sub(r"\bu\b", "x", value)
    for row in rows:
        if row.get("status") != "ok" or int(row["cond"]) in EXCLUDED_DISPLAY_CONDUCTORS:
            continue
        a, t = row["algebraic_data"], row["triple"]
        comparison = row.get("comparison", {})
        agreement = comparison.get("relative_agreement")
        verified = bool(comparison.get("rational")) and agreement is not None and Decimal(str(agreement)) < Decimal('1e-95')
        curves.append(dict(label=row["label"], cond=int(row["cond"]),
                           F=x_expression(a["model"]), points=a["image_points"][:2],
                           orders=[int(t["order_P"]), int(t["order_Q"])],
                           pole=a["pole"],
                           original_triple=t,
                           q=comparison.get("rational") if verified else None,
                           quotient=comparison.get("ratio"), rational_guess=comparison.get("rational"),
                           numerical={key: row.get(key) for key in
                                      ("det", "c_plus", "c_minus", "L_second_derivative")},
                           relative_agreement=comparison.get("relative_agreement"),
                           regulator_digits=row["regulator_digits"], L_digits=row.get("L_digits"),
                           comparison_verified=verified,
                           integrality_verified=all(t.get(key) is True for key in ("verdict_P", "verdict_Q"))))
    snapshot = dict(source=dict(title="Beilinson survey: conductor < 1300, 100-digit targets",
                                input_snapshot=str(output), report_snapshot=str(report),
                                report_sha256=hashlib.sha256(raw).hexdigest(),
                                marked_points=["P", "R"],
                                pole_convention="Q is the report's common pole R; viewer zeros P,R are report zeros P,Q.",
                                conductor_bound=BOUND, regulator_digits=DIGITS, lvalue_digits=DIGITS,
                                excluded_conductors=sorted(EXCLUDED_DISPLAY_CONDUCTORS)),
                    curves=curves)
    Path(output).write_text(json.dumps(snapshot, indent=2, allow_nan=False)+"\n")
    return snapshot


def main(argv=None):
    from .beilinson import CACHE, RESULTS_OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=RESULTS_OUT)
    parser.add_argument("--cache", type=Path, default=CACHE)
    parser.add_argument("--out", type=Path, default=Path("reports"))
    parser.add_argument("--explorer-inputs", type=Path,
                        help="also write an explorer input catalog from the new report")
    args = parser.parse_args(argv)
    export(args.input, args.out, args.cache)
    if args.explorer_inputs:
        explorer_inputs(args.out/"beilinson_100_digits.json", args.explorer_inputs)


if __name__ == "__main__":
    main()
