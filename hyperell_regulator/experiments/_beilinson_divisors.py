"""Beilinson comparisons on the integral divisor kernel of torsion packets."""

from itertools import combinations
import json
import time

from sage.all import QQ, RealField, ZZ, lcm, matrix

from ..certified import display_float, encode_ball, encode_result
from ..components import curve_points
from ..integrality import _integer_matrix, torsion_packets
from ..motivic import MotivicPacket, _divisor_evaluator
from ..regulator import c_invariants


def _packets(C, row, record, selection):
    if selection == 'weierstrass':
        return [MotivicPacket(C)], []
    if selection != 'rational':
        raise ValueError('packet must be rational, weierstrass or geometric')

    from .beilinson import _parse_point, lmfdb_rank

    points = dict.fromkeys(p for p in curve_points(C, row) if p is not None)
    for pair in record.get('pairs', []):
        points.update(dict.fromkeys(_parse_point(pair[key]) for key in ('P', 'Q')))
    points.pop(None, None)
    if len({x for x, _ in points}) < 3:
        return [], []
    if lmfdb_rank(row) == 0:
        invs = row['mw_invs']
        N = lcm([ZZ(n) for n in (json.loads(invs) if isinstance(invs, str) else invs)] or [1])
        # The arithmetic constructor verifies this exponent on every point.
        return [MotivicPacket(C, points, torsion_multiples=N)], []
    packets, unknown = torsion_packets(C, points)
    return [MotivicPacket(C, P) for P in packets if len({x for x, _ in P}) >= 3], unknown


def _number(value, digits, certified):
    if certified:
        return encode_ball(value)
    if digits is None:
        return float(value)
    return RealField(int(3.33*digits)+24)(value).str(truncate=False)


def _nonzero(det, rows, digits, certified):
    if certified:
        return not det.contains_zero()
    scale = max(abs(x) for x in rows[0])*max(abs(x) for x in rows[1])*max(abs(x) for x in rows[2])
    return abs(det) > scale*(10**(-digits+2) if digits else 1e-12)


def _integral_kernel(packet):
    from .geometric_torsion import encode_point

    basis, complete = packet.integral_basis()
    rows = _integer_matrix(matrix(QQ, len(basis), len(packet.points),
                                  [d.coefficients for d in basis]))
    return dict(points=[encode_point(P) for P in packet.points],
                basis=[[str(a) for a in row] for row in rows], complete=complete,
                places={str(q): dict(method=data.method, unresolved=data.unresolved, exact=data.exact)
                        for q, data in packet.integrality.projected_local_data.items()})


def analyse_divisors(record, C, row, selection, *, verbose, with_periods,
                     with_lvalue, digits, reg_digits, lindep_bound, tries, certified,
                     geometric_seconds):
    from .beilinson import beilinson_comparison, standard_lvalue_cached
    from .geometric_torsion import decode_point

    if certified and reg_digits is None:
        reg_digits = digits
    initial_packet = 'rational' if selection == 'geometric' else selection
    result = dict(label=record['label'], cond=record['cond'], packet_selection=selection,
                  packet_used=initial_packet, unresolved_torsion_pairs=0, integral_kernels=[])
    kernels = []
    try:
        packets, unknown = _packets(C, row, record, initial_packet)
        result['unresolved_torsion_pairs'] = len(unknown)
        for packet in packets:
            kernels.append(_integral_kernel(packet))
    except (ArithmeticError, NotImplementedError, ValueError) as error:
        if selection != 'geometric':
            raise
        result['rational_kernel_error'] = str(error)
    result['integral_kernels'].extend(kernels)

    if selection == 'geometric' and not any(len(k['basis']) >= 2 for k in kernels):
        from ._geometric_kernel import geometric_kernel

        def progress(message):
            print('  %s %s' % (record['label'], message), flush=True)

        progress('Fewer than two integral directions found on any rational packet')
        kernel, info = geometric_kernel(C, row, seconds=geometric_seconds,
                                        certified=certified, progress=progress)
        result.update(info, packet_used='geometric')
        if kernel is None:
            return result
        kernels = (kernel,)
        result['integral_kernels'].append(kernel)

    full_support = result.get('packet_support_complete', True)
    options = {'certified': True} if certified else {}
    attempts, errors = [], []
    for kernel in kernels:
        basis = kernel['basis']
        if len(basis) < 2:
            continue
        if verbose:
            print('  %s integral directions, integrating on a common cut' % len(basis), flush=True)
        clock = time.monotonic()
        try:
            # The kernel calculation verified torsion and integrality.
            # Only the individual centres are needed for the regulator.
            support = [i for i in range(len(kernel['points'])) if any(QQ(d[i]) for d in basis)]
            abscissae = [decode_point(kernel['points'][i])[0] for i in support]
            rows = [[QQ(d[i]) for i in support] for d in basis]
            evaluator, rows, shift = _divisor_evaluator(
                C.F.change_ring(QQ), abscissae, rows, digits=reg_digits, **options)
            period, columns = evaluator.columns(rows)
            pairs = []
            for i, j in combinations(range(len(basis)), 2):
                det, mat = evaluator.determinant_of_columns(columns[i], columns[j], period)
                if _nonzero(det, mat, reg_digits, certified):
                    pairs.append((det, mat, i, j))
        except (ArithmeticError, ValueError, NotImplementedError, TimeoutError) as exc:
            errors.append(str(exc))
            continue
        if pairs:
            pairs.sort(key=lambda p: abs(display_float(p[0])), reverse=True)
            attempts.append((kernel, evaluator.polynomial, shift, pairs,
                             time.monotonic()-clock))
    if not attempts:
        if errors:
            return dict(result, status='numerical_failure', error='; '.join(errors))
        if any(len(k['basis']) >= 2 for k in kernels):
            return dict(result, status='no_nonzero_determinant')
        complete = (full_support and all(k['complete'] for k in kernels)
                    and (result['packet_used'] == 'geometric' or not result['unresolved_torsion_pairs']))
        return dict(result, status='insufficient_integral_rank' if complete else 'integral_rank_undecided')

    if with_lvalue:
        result.update(standard_lvalue_cached(record['label'], C, digits=digits,
                                            h1_conductor=record['cond'], **options))
    choices = []
    numeric = lambda x: _number(x, reg_digits, certified)
    for kernel, F, shift, pairs, seconds in attempts:
        cp = cm = None
        periods_error = None
        if with_periods:
            try:
                cp, cm = c_invariants(F, prec=reg_digits, **options)
            except Exception as exc:
                periods_error = str(exc)
        for det, mat, i, j in pairs[:max(1, int(tries))]:
            comparison = beilinson_comparison(
                det, cp, cm, result.get('L_second_derivative'), bound=lindep_bound,
                digits=reg_digits, **options)
            choices.append(dict(
                status='ok', integral=True, basis_complete=kernel['complete'] and full_support,
                basis_indices=(i, j), points=kernel['points'],
                divisors=(kernel['basis'][i], kernel['basis'][j]),
                model=str(F), shift=shift, det=numeric(det),
                matrix=[[numeric(x) for x in r] for r in mat], comparison=comparison,
                c_plus=None if cp is None else numeric(cp), c_minus=None if cm is None else numeric(cm),
                periods_error=periods_error, regulator_seconds=seconds,
                regulator_digits=reg_digits, regulator_method='stokes',
                regulator_cut='tree' if certified else 'mst'))
    best = min(choices, key=lambda c: display_float(c['comparison']['relative_agreement'])
               if 'relative_agreement' in c['comparison'] else float('inf'))
    result.update(best, n_integral_pairs_attempted=len(choices))
    if certified:
        result.update(regulator_certified=True, certified=(
            (not with_periods or best['c_plus'] is not None)
            and (not with_lvalue or result.get('L_certified', False))))
        return encode_result(result)
    return result
