"""Export the completed 100-digit survey with conductor strictly below 1300.

    sage -python -m hyperell_regulator.experiments.beilinson_report

This reads saved results without repeating numerical integrations or L-values.
Exact functions are reconstructed on the saved integration models. Failures and
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
                and settings.get("report_version") == REPORT_VERSION
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
    """Reconstruct the exact functions used for this saved triple, without quadrature."""
    from sage.all import PolynomialRing, QQ, latex
    from hyperell_regulator.components import curve_from_row
    from hyperell_regulator.regulator import function_with_divisor
    from . import beilinson as b

    curve = curve_from_row(curve_row)
    triple = dict(row["triple"])
    for key in ("P", "Q", "R"):
        triple[key] = b._parse_point(triple[key])
    model, images, shift, pole = b._model_for_triple(curve, triple)
    ring = PolynomialRing(QQ, "u")
    if ring(model) != ring(re.sub(r"\bx\b", "u", row["model"])):
        raise ValueError("Reconstructed integration model differs: " + row["label"])
    functions = {}
    for key, point in zip(("P", "Q"), images[:2]):
        a, c = function_with_divisor(model, point, pole, int(triple["order_" + key]))
        a, c = ring.fraction_field()(a), ring.fraction_field()(c)
        expression = latex(a) if a else ""
        if c:
            expression += ("+" if a else "") + r"\left(" + latex(c) + r"\right)v"
        functions[key] = dict(A=str(a), B=str(c), latex=expression or "0")
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
                pole=point_json(pole), functions=functions)


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


def render(rows, source_hash):
    """A standalone document; full decimals and unsuccessful rows stay in JSON."""
    rows = [row for row in rows if int(row.get("cond", 0)) not in EXCLUDED_DISPLAY_CONDUCTORS]
    ok = [row for row in rows if row.get("status") == "ok"]
    lines = [r"\documentclass[10pt]{article}",
             r"\usepackage[a4paper,landscape,margin=12mm]{geometry}",
             r"\usepackage{amsmath,amssymb,booktabs,longtable}",
             r"\begin{document}",
             r"\noindent\textbf{Beilinson survey: conductor $<1300$, 100-digit targets.}",
             r"Conductor $1116$ is retained in the JSON and omitted from these tables.",
             r"The determinant $D$ uses weights $(1,2\operatorname{Re}(z),|z|^2)$ and full cycle periods.",
             r"We compare $q=\pi c^+c^-L^{\prime\prime}(\mathrm{std},1)/D$ with the rational guess $r$.",
             r"The error is $|q/r-1|$. Numerical targets are not certified error bounds.",
             r"Displayed decimals have 10 significant digits; the JSON retains all saved digits.",
             r"A dagger marks an undecided component verdict; no integrality assertion is made for that row.",
             r"\scriptsize", r"\begin{longtable}{lrrrrrrr}", r"\toprule",
             r"Curve & $D$ & $c^+$ & $c^-$ & $L^{\prime\prime}(\mathrm{std},1)$ & $q$ & $r$ & Relative error \\",
             r"\midrule\endhead"]
    for row in ok:
        comparison = row.get("comparison", {})
        t = row["triple"]
        dagger = r"$^{\dagger}$" if any(t.get(k) is not True for k in ("verdict_P", "verdict_Q")) else ""
        values = [number_tex(row.get(k)) for k in ("det", "c_plus", "c_minus", "L_second_derivative")]
        values += [number_tex(comparison.get("ratio")), rational_tex(comparison.get("rational")),
                   number_tex(comparison.get("relative_agreement"), 3)]
        lines.append(escape_tex(row["label"]) + dagger + " & " + " & ".join("$"+v+"$" for v in values) + r" \\")
    lines += [r"\bottomrule\end{longtable}",
              r"\noindent The following table uses original coordinates for $P,Q,R$.",
              r"The saved keys $P,Q$ denote the two zeros and $R$ the common pole.",
              r"On the integration model $v^2=G(u)$, $(\phi_P)=n_P(P-R)$ and $(\phi_Q)=n_Q(Q-R)$.",
              r"For $s=\text{--}$, $(u,v)=(x,2y+h(x))$; otherwise $(u,v)=((x-s)^{-1},d(2y+h(x))/(x-s)^3)$.",
              r"Exact functions use the package's scalar normalization; $D$ includes its norm correction.",
              r"\begin{longtable}{llllrr}", r"\toprule",
              r"Curve & $P$ & $Q$ & $R$ & $n_P$ & $n_Q$ \\", r"\midrule\endhead"]
    for row in ok:
        t, a = row["triple"], row["algebraic_data"]
        lines.append(escape_tex(row["label"]) + " & " + " & ".join("$"+point_tex(t[k])+"$" for k in ("P", "Q", "R")) +
                     " & " + str(t["order_P"]) + " & " + str(t["order_Q"]) + r" \\*")
        shift = r"\text{--}" if a["shift"] is None else rational_tex(a["shift"])
        lines.append(r"\multicolumn{6}{l}{$G(u)=" + a["model_latex"] + r";\quad s=" + shift +
                     r",\ d=" + rational_tex(a["scale"]) + r",\quad \text{pole}=" + point_tex(a["pole"]) + r"$} \\*")
        for key in ("P", "Q"):
            lines.append(r"\multicolumn{6}{p{0.95\linewidth}}{$\phi_" + key + "=" + a["functions"][key]["latex"] + r"$} \\")
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule\end{longtable}", r"\normalsize",
              r"\noindent\textbf{Survey status.}", r"\begin{longtable}{lr}", r"\toprule",
              r"Status & Curves \\", r"\midrule"]
    for status, count in sorted(Counter(row["status"] for row in rows).items()):
        lines.append(escape_tex(status) + " & " + str(count) + r" \\")
    lines += [r"\bottomrule\end{longtable}"]
    issues = [row for row in rows if row.get("status") not in ("ok", "no_admissible_triple")
              or (row.get("status") == "ok" and not row.get("comparison"))]
    if issues:
        lines += [r"\begin{longtable}{lp{0.7\linewidth}}", r"\toprule", r"Curve & Incomplete computation \\", r"\midrule"]
        for row in issues:
            reason = row.get("error") or row.get("L_error") or row.get("L_status") or row.get("periods_error") or row["status"]
            # The JSON retains the complete integrator state, including long
            # decimal coordinates; the table needs the failure reason.
            reason = str(reason).partition("; state=")[0]
            lines.append(escape_tex(row["label"]) + " & " + escape_tex(reason) + r" \\")
        lines += [r"\bottomrule\end{longtable}"]
    lines += [r"\noindent Source: \texttt{results/beilinson\_results.json}.",
              r"\par\noindent\texttt{SHA-256: " + source_hash + "}", r"\end{document}"]
    return "\n".join(lines) + "\n"


def export(source, output, cache):
    raw = Path(source).read_bytes()
    rows = select_results(json.loads(raw))
    if not rows:
        raise ValueError("No current 100-digit results below conductor 1300")
    curves = {r["label"]: r for r in json.loads(Path(cache).read_text())}
    for row in rows:
        if row.get("status") == "ok":
            row["algebraic_data"] = algebraic_data(row, curves[row["label"]])
    checksum = hashlib.sha256(raw).hexdigest()
    tex = render(rows, checksum)
    counts = Counter(row["status"] for row in rows)
    compared = sum(bool(row.get("comparison")) for row in rows)
    displayed = sum(row.get("status") == "ok" and int(row["cond"]) not in EXCLUDED_DISPLAY_CONDUCTORS
                    for row in rows)
    readme = f"""# Beilinson survey: conductor < 1300, 100 digits

This directory contains only the current 100-digit Beilinson survey:

- [Full saved data and exact functions](beilinson_100_digits.json).
- [LaTeX numerical and algebraic tables](beilinson_100_digits.tex).

The snapshot has {len(rows)} curve records, {counts.get('ok', 0)} successful
regulators and {compared} numerical comparisons. Unsuccessful attempts and their
reasons are retained in the JSON; the tables list their status counts and any
numerical failures. This is a survey of locally cached geometrically simple
curves, without a claim of completeness over all curves.
The conductor-1116 record is retained in the JSON, but omitted from the LaTeX
tables and explorer. Table status counts use the same exclusion.
The LaTeX tables and explorer include {displayed} successful regulators.

Regulators and L-values target 100 decimal digits. The determinant uses weights
`(1, 2 Re(z), |z|^2)`, full cycle periods, and
`q = pi*c_plus*c_minus*L_second_derivative/D`. The comparison error is
`abs(q/rational - 1)`. An undecided component verdict remains undecided and is
marked by a dagger. Precision targets are not certified error bounds.

The numerical table displays 10 significant digits. The JSON preserves the
complete saved decimal strings, matrices, precision settings, component
verdicts and exact functions. Points are in the original cached model;
functions are on the displayed integration model, with its coordinate change.

After the calculation has finished, regenerate without repeating integrations:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report
```

The standalone LaTeX document uses `geometry`, `amsmath`, `amssymb`, `booktabs`
and `longtable`. No PDF is generated.

Source: `results/beilinson_results.json` (read without modification).
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
    """Use the report's exact models, functions and 100-digit references."""
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
                           exact_functions=[{key: x_expression(a["functions"][name][key])
                                             for key in ("A", "B")} for name in ("P", "Q")],
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
