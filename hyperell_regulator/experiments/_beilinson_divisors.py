"""Beilinson comparisons on the integral divisor kernel of torsion packets."""

from itertools import combinations
import json
import time

from sage.all import QQ, RealField, ZZ, lcm, matrix

from ..certified import display_float, encode_ball, encode_result
from ..components import curve_points
from ..integrality import _integer_matrix, torsion_packets
from .._galois import galois_trace_basis, galois_traces
from .._algebraic import encode_algebraic
from ..motivic import MotivicPacket, _divisor_evaluator
from ..regulator import c_invariants
from ._determinant_relation import determinant_relation, single_determinant_comparison
from ._kernel_cache import cached_kernel


def _packets(C, row, record, selection):
    if selection == 'weierstrass':
        return [MotivicPacket(C)], []
    if selection != 'rational':
        raise ValueError('packet must be rational, weierstrass or geometric')

    from .beilinson import lmfdb_rank

    points = dict.fromkeys(p for p in curve_points(C, row) if p is not None)
    points.pop(None, None)
    if len({x for x, _ in points}) < 2:
        return [], []
    if lmfdb_rank(row) == 0:
        invs = row['mw_invs']
        N = lcm([ZZ(n) for n in (json.loads(invs) if isinstance(invs, str) else invs)] or [1])
        # The arithmetic constructor verifies this exponent on every point.
        return [MotivicPacket(C, points, torsion_multiples=N)], []
    packets, unknown = torsion_packets(C, points)
    return [MotivicPacket(C, P) for P in packets if len({x for x, _ in P}) >= 2], unknown


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


def _integral_kernel(packet, *, cache=True):
    from .geometric_torsion import encode_point

    points = [encode_point(P) for P in packet.points]
    support = ['points', points, packet._torsion_multiples, packet._max_model_degree]

    def compute():
        basis, complete = packet.integral_basis()
        solver = packet.integrality
        places = {str(q): dict(method=data.method, unresolved=list(data.unresolved), exact=data.exact)
                  for q, data in solver.projected_local_data.items()}
        record = _kernel_record(points, [d.coefficients for d in basis], complete,
                              solver._orbit_sum, solver.abscissa_divisor,
                              lambda D: solver.check(D)['q_integral'], places,
                              galois_trace_basis(solver.abscissae))
        record.update(packet_torsion_order=str(lcm(solver.orders)),
                      packet_torsion_order_exact=True, orders_to_base=list(map(str, solver.orders)))
        return record

    return cached_kernel(packet.f, packet.h, support, compute, incomplete_revision=3) if cache else compute()


def _kernel_record(points, integral, complete, orbit_sum, lift, check, places, traces):
    rows = list(integral)
    verdicts = [True]*len(rows)
    space = matrix(QQ, len(rows), len(points), rows)*orbit_sum.transpose()
    span = space.row_space()
    if any(row not in traces.row_space() for row in space):
        raise ValueError('integral directions must lie in the rational Galois-trace space')
    rational_rank = space.rank()
    for abscissa_row in traces:
        if span.dimension() == traces.nrows():
            break
        if abscissa_row not in span:
            divisor = lift(abscissa_row)
            rows.append(tuple(divisor))
            verdicts.append(check(divisor))
            span = span + matrix(QQ, [abscissa_row]).row_space()
    normalized = _integer_matrix(matrix(QQ, len(rows), len(points), rows))
    full_basis = [[str(a) for a in row] for row in normalized]
    return dict(points=points, basis=full_basis[:len(integral)], complete=complete,
                divisor_basis=full_basis, divisor_integrality=verdicts, places=places,
                divisor_descent='galois_trace', base_field='Q',
                rank_over_Q=rational_rank, trace_space_rank_over_Q=traces.rank())


def _kernel_trace_rank(kernel, basis_key='basis'):
    """Rank over QQ after the Galois trace and symmetric-ordinate projection."""
    from .geometric_torsion import decode_point

    points = [decode_point(P) for P in kernel['points']]
    xs = tuple(dict.fromkeys(x for x, _ in points))
    groups = [[i for i, (x, _) in enumerate(points) if x == a] for a in xs]
    rows = matrix(QQ, len(kernel[basis_key]), len(xs),
                  [[sum(QQ(D[i]) for i in group) for group in groups] for D in kernel[basis_key]])
    traced = galois_traces(xs, rows)
    if rows != traced:
        raise ValueError('cached divisor directions are not rational Galois traces')
    return traced.rank()


def _real_logarithm_rank(abscissae, rows):
    """Upper bound using equality of the real moments at conjugate centres."""
    groups = {}
    for i, x in enumerate(abscissae):
        centre = x.conjugate() if x != 'infinity' and x.imag() < 0 else x
        groups.setdefault(centre, []).append(i)
    return matrix(QQ, [[sum(row[i] for i in g) for g in groups.values()] for row in rows]).rank()


def analyse_divisors(record, C, row, selection, **options):
    result = _analyse_packet(record, C, row, selection, **options)
    if (selection != 'geometric' or result['packet_used'] != 'rational'
            or not options.get('ignore_integrality')
            or not (options['with_lvalue'] and options['with_periods'])
            or result['integral_rank'] >= 2 or result.get('linear_comparison', {}).get('terms')
            or (result.get('L_status') and result.get('L_second_derivative') is None)):
        return result
    geometric = _analyse_packet(record, C, row, selection, geometric_only=True, **options)
    if geometric.get('linear_comparison', {}).get('terms'):
        return dict(geometric, rational_attempt=result)
    return dict(result, geometric_attempt=geometric)


def _analyse_packet(record, C, row, selection, *, verbose, with_periods,
                    with_lvalue, digits, reg_digits, lindep_bound, tries, certified,
                    geometric_seconds, ignore_integrality=False, geometric_only=False):
    from .beilinson import beilinson_comparison, standard_lvalue_cached
    from .geometric_torsion import decode_point

    if certified and reg_digits is None:
        reg_digits = digits
    initial_packet = 'rational' if selection == 'geometric' else selection
    basis_key = 'divisor_basis' if ignore_integrality else 'basis'
    result = dict(label=record['label'], cond=record['cond'], packet_selection=selection,
                  packet_used=initial_packet, unresolved_torsion_pairs=0, integral_kernels=[])
    kernels = []
    if verbose and not geometric_only:
        print('  %s checking the %s packet' % (record['label'], initial_packet), flush=True)
    try:
        packets, unknown = ([], []) if geometric_only else _packets(C, row, record, initial_packet)
        result['unresolved_torsion_pairs'] = len(unknown)
        for packet in packets:
            kernels.append(_integral_kernel(packet))
    except (ArithmeticError, NotImplementedError, ValueError) as error:
        if selection != 'geometric':
            raise
        result['rational_kernel_error'] = str(error)
    result['integral_kernels'].extend(kernels)

    if selection == 'geometric' and not any(len(k[basis_key]) >= 2 for k in kernels):
        from ._geometric_kernel import geometric_kernel

        def progress(message):
            print('  %s %s' % (record['label'], message), flush=True)

        if geometric_only:
            progress('No rational-packet relation and fewer than two rational integral directions')
        else:
            qualifier = '' if ignore_integrality else 'integral '
            progress('Fewer than two %sdirections found on any rational packet' % qualifier)
        kernel, info = geometric_kernel(C, row, seconds=geometric_seconds, progress=progress)
        result.update(info, packet_used='geometric')
        if kernel is None:
            return result
        kernels = (kernel,)
        result['integral_kernels'].append(kernel)

    full_support = result.get('packet_support_complete', True)
    complete_support = full_support and (result['packet_used'] == 'geometric'
                                         or not result['unresolved_torsion_pairs'])
    for kernel in kernels:
        kernel['rank_over_Q'] = _kernel_trace_rank(kernel)
    result.update(integral_rank=max((k['rank_over_Q'] for k in kernels), default=0),
                  integral_rank_field='Q', integral_rank_space='galois_traces',
                  integral_rank_complete=complete_support and all(k['complete'] for k in kernels))
    if complete_support and all(k['complete'] or 'rank_upper_bound' in k for k in kernels):
        result['integral_rank_upper_bound'] = max(
            (k.get('rank_upper_bound', k['rank_over_Q']) for k in kernels), default=0)
    options = {'certified': True} if certified else {}
    attempts, errors, symmetry_bounds = [], [], []
    for kernel in kernels:
        basis = kernel[basis_key]
        rank = _kernel_trace_rank(kernel, basis_key)
        if rank < 2:
            continue
        if verbose:
            print('  %s rational Galois-trace directions, integrating on a common cut' % rank, flush=True)
        clock = time.monotonic()
        try:
            # The kernel calculation verified torsion and integrality.
            # Only the individual centres are needed for the regulator.
            support = [i for i in range(len(kernel['points'])) if any(QQ(d[i]) for d in basis)]
            abscissae = [decode_point(kernel['points'][i])[0] for i in support]
            rows = [[QQ(d[i]) for i in support] for d in basis]
            bound = _real_logarithm_rank(abscissae, rows)
            if bound < 2:
                symmetry_bounds.append(bound)
                continue
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
            integration_input = dict(
                centers=[encode_algebraic(x) for x in evaluator.logarithmic_centers],
                abscissae=['infinity' if x == 'infinity' else encode_algebraic(x)
                           for x in evaluator.abscissae],
                divisors=[list(map(str, D)) for D in rows])
            attempts.append((kernel, evaluator.polynomial, shift, pairs,
                             time.monotonic()-clock, integration_input))
    if not attempts:
        if errors:
            return dict(result, status='numerical_failure', error='; '.join(errors))
        eligible = sum(len(k[basis_key]) >= 2 for k in kernels)
        if eligible:
            if len(symmetry_bounds) == eligible:
                result.update(regulator_rank_upper_bound=max(symmetry_bounds),
                              regulator_issue='Conjugate logarithmic centres give dependent real regulator columns')
            return dict(result, status='no_nonzero_determinant')
        if ignore_integrality:
            complete = full_support and (result['packet_used'] == 'geometric'
                                         or not result['unresolved_torsion_pairs'])
            return dict(result, status='insufficient_torsion_rank' if complete else 'torsion_rank_undecided')
        return dict(result, status=('insufficient_integral_rank' if result['integral_rank_complete']
                                   else 'integral_rank_undecided'))

    if with_lvalue:
        if verbose:
            print('  %s evaluating the L-value' % record['label'], flush=True)
        result.update(standard_lvalue_cached(record['label'], C, digits=digits,
                                            h1_conductor=record['cond'], **options))
    choices, pair_searches = [], []
    n_pairs = n_integral = 0
    numeric = lambda x: _number(x, reg_digits, certified)
    for kernel, F, shift, pairs, seconds, integration_input in attempts:
        cp = cm = None
        periods_error = None
        if with_periods:
            try:
                if verbose:
                    print('  %s evaluating the cycle periods' % record['label'], flush=True)
                cp, cm = c_invariants(F, prec=reg_digits, **options)
            except Exception as exc:
                periods_error = str(exc)
        selected_pairs = pairs[:max(1, int(tries))]
        linear, terms = {}, []
        if ignore_integrality:
            linear = determinant_relation(
                pairs, cp, cm, result.get('L_second_derivative'), digits=reg_digits,
                lvalue_digits=result.get('L_digits') or digits, bound=lindep_bound, certified=certified)
            terms = linear.get('terms', [])
            representative = tuple(terms[0]['basis_indices']) if terms else pairs[0][2:]
            selected_pairs = [pair for pair in pairs if pair[2:] == representative]
            pair_searches.append(dict(
                kernel_index=result['integral_kernels'].index(kernel),
                model=str(F), shift=shift, linear_comparison=linear,
                determinants=[dict(basis_indices=(i, j), det=numeric(det))
                              for det, _, i, j in pairs]))
            n_pairs += len(pairs)
            n_integral += sum(all(kernel['divisor_integrality'][k] is True for k in (i, j))
                              for _, _, i, j in pairs)
        for det, mat, i, j in selected_pairs:
            comparison = {}
            if not ignore_integrality:
                comparison = beilinson_comparison(
                    det, cp, cm, result.get('L_second_derivative'), bound=lindep_bound,
                    digits=reg_digits, **options)
            elif len(terms) == 1:
                comparison = single_determinant_comparison(linear, det)
            verdicts = ([kernel['divisor_integrality'][k] for k in (i, j)]
                        if ignore_integrality else [True, True])
            used_verdicts = [kernel['divisor_integrality'][k]
                             for term in terms for k in term['basis_indices']] or verdicts
            integral = False if False in used_verdicts else None if None in used_verdicts else True
            if len(terms) > 1 and integral is False:
                # Individual boundaries may cancel in a linear combination.
                integral = None
            choice = dict(
                status='ok', integral=integral, divisor_integrality=verdicts,
                basis_complete=kernel['complete'] and full_support,
                basis_indices=(i, j), points=kernel['points'],
                divisors=(kernel[basis_key][i], kernel[basis_key][j]),
                packet_torsion_order=kernel.get('packet_torsion_order'),
                packet_torsion_order_exact=kernel.get('packet_torsion_order_exact', False),
                end_alg=row.get('end_alg'), regulator_input=integration_input,
                model=str(F), shift=shift, det=numeric(det),
                matrix=[[numeric(x) for x in r] for r in mat], comparison=comparison,
                c_plus=None if cp is None else numeric(cp), c_minus=None if cm is None else numeric(cm),
                periods_error=periods_error, regulator_seconds=seconds,
                regulator_digits=reg_digits, regulator_method='stokes',
                regulator_cut='tree' if certified else 'mst')
            if ignore_integrality:
                choice.update(linear_comparison=linear, pair_search_index=len(pair_searches)-1)
            choices.append(choice)
    if ignore_integrality:
        best = min(choices, key=lambda c: display_float(c['linear_comparison']['relative_residual'])
                   if c['linear_comparison'].get('terms') else float('inf'))
        result.update(best, pair_searches=pair_searches, n_divisor_pairs_attempted=n_pairs,
                      n_integral_pairs_attempted=n_integral)
    else:
        best = min(choices, key=lambda c: display_float(c['comparison']['relative_agreement'])
                   if 'relative_agreement' in c['comparison'] else float('inf'))
        result.update(best, n_divisor_pairs_attempted=len(choices),
                      n_integral_pairs_attempted=sum(c['integral'] is True for c in choices))
    if certified:
        result.update(regulator_certified=True, certified=(
            (not with_periods or best['c_plus'] is not None)
            and (not with_lvalue or result.get('L_certified', False))))
        return encode_result(result)
    return result
