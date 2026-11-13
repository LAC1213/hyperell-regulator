"""Export one x-plane branch diagram per curve, with every logarithmic weight."""
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
    """Validate a polynomial, rational logarithmic centers and precision."""
    from sage.all import PolynomialRing, QQ
    if not isinstance(request, dict):
        raise ValueError('Expected an object with F, centers and prec')
    S = PolynomialRing(QQ, 'x')
    value = parse_expression(request.get('F'), S.fraction_field())
    if value.denominator().degree() > 0:
        raise ValueError('F must be a polynomial')
    F = S(value)
    if not 5 <= F.degree() <= 12 or F.gcd(F.derivative()).degree() > 0:
        raise ValueError('F must be squarefree with degree 5 through 12')
    centers = request.get('centers', ['0'])
    if not isinstance(centers, list) or not 1 <= len(centers) <= 8:
        raise ValueError('Provide one through eight rational centers')
    centers = [QQ(parse_expression(str(a), S.fraction_field())) for a in centers]
    prec = request.get('prec', 15)
    if type(prec) is not int or not 10 <= prec <= 100:
        raise ValueError('prec must be an integer from 10 through 100')
    return F, centers, prec


def cut_system_data(F, centers=(0,), prec=15, *, samples=40):
    """Export shared edges, bank gluing and all ordinary/logarithmic integrals."""
    import numpy as np
    from .stokes import compute, prod
    details = {}
    periods, regulators = compute(F, centers, prec, details=details)
    cut, u, ell = details['cut'], details['u'], details['ell']
    ar, g = cut.ar, cut.g
    number = lambda z: str(ar.output(z))
    cnumber = lambda z: [number(z.real()), number(z.imag())]
    point = lambda z: [float(z.real), float(z.imag)]
    edges = []
    for i, a in enumerate(cut.nodes):
        b = cut.nodes[i+1] if i+1 < len(cut.nodes) else a+ar.R('1.5')*cut.span*cut.direction
        ref, yr = complex(cut.references[i]), complex(cut.reference_y[i])
        roots = [complex(r) for r in cut.roots]
        base = np.linspace(complex(a), complex(b), samples)
        # Sampling follows the same square-root factors as integration. The
        # endpoint inset is display-only; all integrals include endpoints.
        near = np.linspace(complex(a), complex(b), samples+2)[1:-1]
        ys = [yr*prod(np.sqrt((z-r)/(ref-r)) for r in roots) for z in near]
        lifts = []
        for sheet in (0, 1):
            s = cut.upper[i][sheet]
            sign = 1 if s == 0 else -1
            I = [u[a][i][s] for a in range(g)]
            logs = []
            for k in range(len(cut.centers)):
                L = [ell[k][a][i][s] for a in range(g)]
                logs.append(dict(minus=[cnumber(z) for z in L],
                                 plus=[cnumber(z+cut.jumps[k][i]*v) for z, v in zip(L, I)]))
            xyz = [[float(z.real), float(z.imag), float(np.arcsinh((np.exp(-.55j)*sign*y).real))]
                   for z, y in zip(near, ys)]
            lifts.append(dict(sheet=sheet+1, minus_sheet=cut.lower[i].index(s)+1,
                              ordinary=[cnumber(z) for z in I], logs=logs, xyz=xyz))
        edges.append(dict(index=i+1, base=[point(z) for z in base], lifts=lifts,
                          infinite=i+1 == len(cut.nodes),
                          jumps=[cnumber(j[i]) for j in cut.jumps]))
    weights = []
    for a, reg in zip(cut.centers, regulators):
        pairings, k = [], 0
        for i in range(g):
            for j in range(i, g):
                pairings.append(dict(forms=[i+1, j+1], period=str(periods[k]),
                                     regulator=str(reg[k]), pairing=number(-4*ar.R(reg[k])/ar.R.pi())))
                k += 1
        weights.append(dict(center=str(a), label='log|x - (%s)|' % a, pairings=pairings))
    return dict(schema_version=2, input=dict(F=str(cut.F), centers=[str(a) for a in cut.centers], prec=prec),
                degree=2, genus=g, precision=prec, weights=weights, edges=edges,
                nodes=[point(complex(z)) for z in cut.nodes],
                roots=[point(complex(z)) for z in cut.roots],
                direction=point(complex(cut.direction)), warnings=[],
                conventions=dict(forms='omega_a = x^(a-1) dx/y',
                    geometry='(Re x, Im x, asinh(Re(exp(-0.55i) y))). Display rays are truncated; integrals include infinity.',
                    normalization='Vectors use (1, Re(x), |x|^2); double the middle entry in genus-two determinants.'))


def compute_request(request):
    return cut_system_data(*parse_request(request))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--F', required=True)
    parser.add_argument('--centers', nargs='+', default=['0'])
    parser.add_argument('--prec', type=int, default=15)
    parser.add_argument('--out', type=Path, default=Path('cut-system.json'))
    args = parser.parse_args(argv)
    result = compute_request(dict(F=args.F, centers=args.centers, prec=args.prec))
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Saved %s: two sheets, %d shared edges' % (args.out, len(result['edges'])))


if __name__ == '__main__':
    main()
