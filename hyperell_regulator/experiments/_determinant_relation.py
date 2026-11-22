"""Integer-relation searches with the normalized L-value as a required term."""

from itertools import combinations
from math import comb

from sage.all import QQ, ZZ, RealField, gcd, identity_matrix, matrix

from ..certified import decode_ball


def determinant_relation(pairs, c_plus, c_minus, lvalue, *, digits, lvalue_digits,
                         bound, certified):
    r"""Search for A = sum(q_ij det(i,j)) with A = pi c+ c- L''.

    The relation lattice keeps all determinant coordinates, including
    those with dependencies among themselves. Only rows involving A
    can become candidates. Acceptance uses the original inputs.
    """
    if not pairs or lvalue is None or c_plus is None or c_minus is None:
        return {}
    precision = min(digits or 14, int(lvalue_digits))
    P = RealField(max(53, int(3.33*precision)+24))
    if certified:
        C = pairs[0][0].parent()
        L = decode_ball(lvalue) if isinstance(lvalue, dict) else lvalue
        target = C.pi()*C(c_plus)*C(c_minus)*C(L)
        if target.contains_zero():
            return dict(error='Normalized L-value enclosure contains zero')
        values = [P(v.mid()) for v in [target, *(p[0] for p in pairs)]]
    else:
        target = P.pi()*P(c_plus)*P(c_minus)*P(lvalue)
        if not target:
            return dict(error='Normalized L-value is zero')
        values = [target, *(P(p[0]) for p in pairs)]

    maximum = max(map(abs, values))
    cut = P(10)**precision / maximum
    approximations = [ZZ((cut*v).round()) for v in values]
    # The numerical survey reserves three digits for accumulated error in
    # the periods, L-value and determinants. Certified inputs carry that
    # error in their enclosures and use a tighter threshold.
    tolerance = QQ(10)**(-max(3, precision-(2 if certified else 3)))
    max_height = 10**max(1, int(bound))
    complexity_limit = 10**max(3, precision-5)
    best = None
    seen = set()

    def screen(indices):
        nonlocal best
        lattice = identity_matrix(ZZ, len(indices)).augment(
            matrix(ZZ, len(indices), 1, [approximations[i] for i in indices]))
        for row in lattice.LLL().rows():
            if not row[0] or not any(row[1:-1]):
                continue
            divisor = gcd(row[:-1])
            c = tuple(row[:-1]/(divisor if row[0] > 0 else -divisor))
            sparse = tuple((i-1, a) for i, a in zip(indices[1:], c[1:]) if a)
            signature = c[0], sparse
            if signature in seen:
                continue
            seen.add(signature)
            height = max(map(abs, c))
            terms = len(sparse)
            # Require enough precision for the number and height of the
            # coefficients being guessed, including the choice of support.
            if (height > max_height or
                    (2*height+1)**(terms+1)*comb(len(pairs), terms) > complexity_limit):
                continue
            coefficients = [(i, -QQ(a)/c[0]) for i, a in sparse]
            approx = (sum(q*pairs[i][0] for i, q in coefficients) if certified else
                      sum(P(q)*values[i+1] for i, q in coefficients))
            residual = target-approx
            relative = abs(residual/target)
            if not relative < tolerance or (certified and not residual.contains_zero()):
                continue
            score = terms, height, sum(map(abs, c)), float(relative)
            if best is None or score < best[0]:
                best = score, dict(
                    terms=[dict(basis_indices=pairs[i][2:], coefficient=str(q))
                           for i, q in coefficients],
                    normalized_lvalue=target, residual=abs(residual), relative_residual=relative,
                    certified=certified, exact_relation_proved=False)

    for i in range(1, len(values)):
        screen((0, i))
    if best is None and len(values) > 2:
        screen(tuple(range(len(values))))
    if best is None and len(values) > 3:
        for i, j in combinations(range(1, len(values)), 2):
            screen((0, i, j))
    return best[1] if best is not None else {}


def single_determinant_comparison(linear, det):
    """Present an accepted one-term relation in the usual ratio format."""
    q = QQ(linear['terms'][0]['coefficient'])
    norm = linear['normalized_lvalue']
    ratio = norm/det
    agreement = abs(ratio/q-1)
    result = dict(
        normalisation="pi c+ c- L''(std,1) / D", norm_value=norm, ratio=ratio,
        rational=q, rational_is_integer=q.denominator() == 1, relative_agreement=agreement,
        certified=linear['certified'], exact_relation_proved=False)
    if linear['certified']:
        result['relative_agreement_upper'] = str(agreement.upper().exact_rational())
    return result
