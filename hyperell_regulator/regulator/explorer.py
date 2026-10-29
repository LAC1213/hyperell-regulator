"""Export the integrator's actual cuts and line integrals for the web viewer.

Run with ``sage -python -m hyperell_regulator.regulator.explorer --help``.
Numerical values are decimal strings; only display geometry uses float64.
"""

import argparse
import ast
import json
from pathlib import Path


def parse_expression(source, field):
    """Parse arithmetic in QQ(x), without evaluating Python supplied by a visitor."""
    if not isinstance(source, str) or not source.strip() or len(source) > 2000:
        raise ValueError("Enter an expression of at most 2000 characters")
    try:
        tree = ast.parse(source.replace("^", "**"), mode="eval")
    except (SyntaxError, RecursionError) as exc:
        raise ValueError("Use rational arithmetic in x, for example x^5 + (x+1)^2") from exc
    if sum(1 for _ in ast.walk(tree)) > 300:
        raise ValueError("Expression is too large")

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            if abs(node.value).bit_length() > 256:
                raise ValueError("Integer coefficient is too large")
            return field(node.value)
        if isinstance(node, ast.Name) and node.id == "x":
            return field.gen()
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = visit(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp):
            left = visit(node.left)
            if isinstance(node.op, ast.Pow):
                exponent = node.right
                sign = 1
                if isinstance(exponent, ast.UnaryOp) and isinstance(exponent.op, ast.USub):
                    sign, exponent = -1, exponent.operand
                if (not isinstance(exponent, ast.Constant) or type(exponent.value) is not int
                        or not 0 <= exponent.value <= 64):
                    raise ValueError("Powers must be integers from -64 to 64")
                result = left ** (sign * exponent.value)
            else:
                right = visit(node.right)
                if isinstance(node.op, ast.Add):
                    result = left + right
                elif isinstance(node.op, ast.Sub):
                    result = left - right
                elif isinstance(node.op, ast.Mult):
                    result = left * right
                elif isinstance(node.op, ast.Div):
                    result = left / right
                else:
                    raise ValueError("Only +, -, *, / and powers are supported")
            if max(result.numerator().degree(), result.denominator().degree()) > 64:
                raise ValueError("Expression degree exceeds 64")
            coefficients = list(result.numerator()) + list(result.denominator())
            if any(max(c.numerator().nbits(), c.denominator().nbits()) > 4096
                   for c in coefficients):
                raise ValueError("Expression coefficients are too large")
            return result
        raise ValueError("Only integers, rational fractions, x and arithmetic are supported")

    try:
        return visit(tree.body)
    except ZeroDivisionError as exc:
        raise ValueError("Expression divides by zero") from exc


def parse_request(request):
    """Validate a JSON request and return exact Sage objects and precision."""
    from sage.all import PolynomialRing, QQ

    if not isinstance(request, dict):
        raise ValueError("Expected a JSON object with F, A, B and prec")
    ring = PolynomialRing(QQ, "x")
    field = ring.fraction_field()
    F, A, B = (parse_expression(request.get(key), field) for key in ("F", "A", "B"))
    if F.denominator().degree() > 0:
        raise ValueError("F must be a polynomial")
    F = ring(F)
    if not 5 <= F.degree() <= 12:
        raise ValueError("The explorer supports degrees 5 through 12 (genera 2 through 5)")
    if F.gcd(F.derivative()).degree() > 0:
        raise ValueError("F must be squarefree: the curve is singular")
    if not B and A.numerator().degree() <= 0 and A.denominator().degree() <= 0:
        raise ValueError("phi must be nonconstant")
    prec = request.get("prec", 15)
    if type(prec) is not int or not 10 <= prec <= 100:
        raise ValueError("prec must be an integer from 10 through 100")
    return F, (A, B), prec


def _sample_piece(cov, a, b, reference, labels, count):
    """Display samples, continued from the same labelled midpoint as quadrature.

    Finite endpoints are approached rather than evaluated; the last ray is
    truncated. These samples never enter the numerical integration.
    """
    import numpy as np
    from .stokes import follow
    from ._legs import leg_in_w

    t = np.linspace(0.002, 0.998, count)
    w = a + (b - a) * t
    xx = np.empty((count, cov.N), complex)
    yy = np.empty_like(xx)
    rx, ry = cov.fibre(reference)
    for indices in (np.flatnonzero(t <= 0.5)[::-1], np.flatnonzero(t > 0.5)):
        last, x, y = reference, rx.copy(), ry.copy()
        for k in indices:
            try:
                x, y = follow(cov, last, w[k], x, y)
            except RuntimeError:
                # The polynomial-root matcher cannot cross a finite-value
                # point at infinity. Use the integrator's chart-aware
                # continuation there, without doing another quadrature.
                x, y, _, _ = leg_in_w(cov, last, w[k], x, y,
                                      np.zeros(cov.N, complex), integrate=False)
            xx[k], yy[k] = x[labels], y[labels]
            last = w[k]
    return w, xx, yy


def refresh_geometry(data, samples=40):
    """Refresh display paths without recomputing any saved integral."""
    from .cover import Cover
    from .stokes import cut_data

    F, phi, _ = parse_request(data["input"])
    cov = Cover(F, *phi)
    nodes, kinds, direction = cov.cut_path()
    upper, lower, _, _ = cut_data(cov, nodes, direction)
    return _fill_geometry(data, cov, nodes, kinds, direction, upper, lower, samples)


def _fill_geometry(data, cov, nodes, kinds, direction, upper, lower, samples):
    import numpy as np
    from .stokes import cut_scale

    starts = [j for j, kind in enumerate(kinds) if kind != "none"]
    groups = [range(a, b) for a, b in zip(starts, starts[1:] + [len(nodes)])]
    if len(groups) != len(data["edges"]):
        raise RuntimeError("The saved cut topology changed")
    span = cut_scale(nodes)[0]
    warnings = []
    curve_error = function_error = 0.0
    for group, edge in zip(groups, data["edges"]):
        for lift in edge["lifts"]:
            lift["xyz"] = []
        for j in group:
            if any(lower[j].index(upper[j][i]) + 1 != edge["lifts"][i]["minus_sheet"] for i in range(cov.N)):
                raise RuntimeError("The saved sheet labels changed")
            a = complex(nodes[j])
            b = complex(nodes[j+1]) if j+1 < len(nodes) else complex(*edge["base"][-1])
            reference = (a+b)/2 if j+1 < len(nodes) else a+span*direction
            try:
                w, x, y = _sample_piece(cov, a, b, reference, upper[j], samples)
                fx = np.polyval(cov.F, x)
                ce = float(np.max(abs(y*y-fx)/(1+abs(fx))))
                fe = float(np.max(abs(cov.phi_stable(x, y)-w[:, None])/(1+abs(w[:, None]))))
                if not np.isfinite(ce+fe) or max(ce, fe) > 1e-6:
                    raise RuntimeError("Display continuation residuals too large (%g, %g)" % (ce, fe))
                curve_error, function_error = max(curve_error, ce), max(function_error, fe)
                height = np.arcsinh((np.exp(-0.55j)*(y if cov.pure_x else x)).real)
                for i, lift in enumerate(edge["lifts"]):
                    lift["xyz"].extend([[float(z.real), float(z.imag), float(h)] for z, h in zip(w, height[:, i])])
            except (RuntimeError, ValueError, FloatingPointError) as exc:
                warnings.append("Display samples on segment %d unavailable: %s" % (j+1, exc))
                for lift in edge["lifts"]:
                    lift["xyz"].append(None)
    data.update(warnings=warnings, geometry_version=2,
                validation=dict(max_relative_curve_residual=curve_error,
                                max_relative_function_residual=function_error))
    return data


def cut_system_data(F, phi, prec=15, *, samples=40):
    """Compute checked pairings and all lifted-edge integrals in one run.

    The two logarithmic banks satisfy L_plus = L_minus + jump * I.
    Edge labels use the positive-bank sheet order, and tau maps a positive
    bank's sheet to the corresponding negative-bank sheet. Ordinary and
    logarithmic integrals are oriented from zero toward infinity on both
    banks. The viewer reverses signs for boundary traversal when requested.
    """
    import numpy as np
    from ._mp_stokes import compute
    from .cover import Cover
    from .stokes import cut_scale

    cov = Cover(F, *phi)
    if cov.N > 24:
        raise ValueError("The explorer supports covering degree N <= 24")
    places = [cov.places_over(z) for z in (0, None)]
    if any(p is None or len(p) != 1 or p[0]["degree"] != 1
           or p[0]["e"] != cov.N for p in places):
        raise ValueError("phi must have divisor N P - N Q with single rational points P and Q")
    details = {}
    periods, regulator = compute(F, phi, prec, details=details)
    ar, cov = details["ar"], details["cover"]
    nodes, kinds = details["nodes"], details["kinds"]
    direction, jump = details["direction"], details["jump"]
    u, ell = details["u"], details["ell"]
    upper, lower = details["upper"], details["lower"]
    n, g, r = cov.N, cov.g, len(nodes)

    def number(v):
        return str(ar.output(v))

    def cnumber(v):
        v = ar.C(v)
        return [number(v.real()), number(v.imag())]

    def point(v):
        return [float(v.real), float(v.imag)]

    start_indices = [j for j, kind in enumerate(kinds) if kind != "none"]
    groups = [range(a, b) for a, b in zip(start_indices, start_indices[1:] + [r])]
    span = cut_scale(np.array([complex(z) for z in nodes]))[0]
    edges = []
    for edge_index, group in enumerate(groups):
        # A geometric lift has the same pair of boundary sheets along all
        # unramified corners. Refuse to silently combine differing labels.
        tau = [lower[group.start].index(upper[group.start][i]) for i in range(n)]
        for j in group:
            if [lower[j].index(upper[j][i]) for i in range(n)] != tau:
                raise RuntimeError("Sheet labels change across an unramified cut corner")
        base = [point(complex(nodes[j])) for j in group]
        final = group.stop == r
        endpoint = (complex(nodes[-1]) + 1.5 * span * complex(direction)
                    if final else complex(nodes[group.stop]))
        base.append(point(endpoint))
        lifts = []
        for sheet in range(n):
            I = [sum((u[a][j][upper[j][sheet]] for j in group), ar.C(0)) for a in range(g)]
            L = [sum((ell[a][j][upper[j][sheet]] for j in group), ar.C(0)) for a in range(g)]
            lifts.append(dict(sheet=sheet + 1, minus_sheet=tau[sheet] + 1,
                              ordinary=[cnumber(z) for z in I],
                              log_minus=[cnumber(z) for z in L],
                              log_plus=[cnumber(z + jump * v) for z, v in zip(L, I)],
                              xyz=[]))
        edges.append(dict(index=edge_index + 1, base=base, infinite=final,
                          integration_segments=len(group), lifts=lifts))
    pairings = []
    index = 0
    for a in range(g):
        for b in range(a, g):
            pairings.append(dict(forms=[a + 1, b + 1], period=str(periods[index]),
                                 regulator=str(regulator[index]),
                                 pairing=number(-4 * ar.R(regulator[index]) / ar.R.pi())))
            index += 1
    data = dict(schema_version=1, input=dict(F=str(F), A=str(phi[0]), B=str(phi[1]), prec=prec),
                genus=g, degree=n, precision=prec, jump=cnumber(jump),
                divisor=dict(P=str(places[0][0]["place"]), Q=str(places[1][0]["place"]), N=n),
                pairings=pairings, edges=edges,
                conventions=dict(forms="omega_a = x^(a-1) dx/y",
                                 pairing="pairing = -4/pi * package_regulator",
                                 geometry="(Re w, Im w, asinh(Re(exp(-0.55i) q))); q=y for pure x maps, otherwise x. "
                                          "Finite endpoints are approached; the ray is truncated for display. "
                                          "Integrals include the exact endpoints and the entire infinite ray."))
    return _fill_geometry(data, cov, np.array([complex(z) for z in nodes]), kinds,
                          complex(direction), upper, lower, samples)


def compute_request(request):
    """Compute a validated web/CLI request."""
    return cut_system_data(*parse_request(request))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--F", required=True, help="Polynomial in x over QQ")
    parser.add_argument("--A", required=True, help="phi = A(x) + B(x)y")
    parser.add_argument("--B", required=True)
    parser.add_argument("--prec", type=int, default=15)
    parser.add_argument("--out", type=Path, default=Path("cut-system.json"))
    args = parser.parse_args(argv)
    result = compute_request(dict(F=args.F, A=args.A, B=args.B, prec=args.prec))
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print("Saved %s (%d sheets, %d cut edges)" % (args.out, result["degree"], len(result["edges"])))


if __name__ == "__main__":
    main()
