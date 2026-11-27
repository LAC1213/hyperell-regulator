r"""Beilinson experiments on integral Galois-trace torsion-packet divisors.

Load a curve by LMFDB label with ``hyperell_regulator.load_curve`` and call
``run_curve(C)`` for one curve, or ``run(N)`` for a conductor-bounded survey.
"""

import hashlib
import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
import time
from contextlib import contextmanager, closing, nullcontext
from math import isfinite

from sage.all import RealField, Integer, QQ, RR, ZZ, gp, pi

from hyperell_regulator.components import CACHE, curve_from_row, load_curve, load_curves
from hyperell_regulator.lfunction import StandardLFunction
from hyperell_regulator.experiments.beilinson_local_data import lfunction_options, local_data_metadata
from hyperell_regulator.paths import data_path, ensure_parent
from hyperell_regulator.certified import (
    ANALYTIC_ASSUMPTIONS, BACKEND_REVISION, encode_ball, decode_ball, display_float,
)

DEFAULT_CONDUCTOR_BOUND = 1000
DEFAULT_GEOMETRIC_SECONDS = 30
DEFAULT_CURVE_SECONDS = 180
REPORT_VERSION = 13
BEILINSON_END_ALGS = ('Q', 'CM')
REPORT_COLUMNS = "%-18s %5s %-26s %9s %14s %12s %10s %10s"
LVALUES = data_path('results', 'lvalue_cache.json')
RESULTS_OUT = data_path('results', 'beilinson_results.json')


def _fmt(point):
    if point[0] == 'infinity':
        return 'oo(Y=%s)' % point[1]
    return '(%s, %s)' % tuple(point)


def lmfdb_rank(row):
    r"""The Mordell--Weil rank from ``mw_invs``, or ``None`` if absent."""
    if row is None:
        return None
    inv = row.get("mw_invs")
    if inv is None:
        return None
    if isinstance(inv, str):
        try:
            inv = json.loads(inv)
        except ValueError:
            return None
    return sum(1 for t in inv if ZZ(t) == 0)


def _sage_real(v, digits):
    r"""Convert ``v`` to a Sage real with guard digits."""
    from sage.all import RealField

    R = RealField(int(3.33 * int(digits)) + 24)
    return R(v)


def load_lvalues(path=LVALUES):
    r"""Saved standard `L`-values, keyed by curve label."""
    try:
        with open(path) as fh:
            return json.load(fh)
    except (IOError, OSError, ValueError):
        return {}


def save_lvalue(label, data, path=LVALUES):
    r"""Merge one curve's `L`-data into ``path``, keyed by label."""
    have = load_lvalues(path)
    have[label] = data
    _save_json(path, have, sort_keys=True)
    return have


def _save_json(path, data, *, sort_keys=False):
    """Publish complete cache files even when a calculation is interrupted."""
    destination = Path(ensure_parent(path))
    with NamedTemporaryFile(mode='w', dir=destination.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(data, stream, indent=2, sort_keys=sort_keys, default=str)
            stream.write('\n')
            stream.close()
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)


def standard_lvalue_cached(label, C, digits=15, path=LVALUES, force=False,
                           h1_conductor=None, certified=False, progress=None):
    r"""Compute or reuse `L^{\prime\prime}(M_{\mathrm{std}}(C),1)`.

    Reuse requires sufficient precision and coefficient count. Return the value,
    conductor, functional equation residual, precision and coefficient count;
    ``L_status`` records failure.
    """
    if certified:
        return _certified_lvalue_cached(label, C, digits, path, force, h1_conductor, progress)
    have = load_lvalues(path)
    local_data = local_data_metadata(label)
    try:
        LF = StandardLFunction(C, h1_conductor=h1_conductor,
                              **lfunction_options(label, C))
        if not force and label in have:
            rec = have[label]
            # Records before the PARI precision fix contain at most sixteen
            # correct digits. Also reject a stale conductor or coefficient
            # count before reusing an otherwise sufficiently precise value.
            if rec.get("L_second_derivative") is not None and \
                    rec.get("L_local_data") == local_data and \
                    rec.get("L_pari_precision") and \
                    (digits <= 15 or isinstance(rec.get("L_second_derivative"), str)) and \
                    ZZ(rec.get("L_digits", 0)) >= digits and \
                    rec.get("standard_conductor") == LF.conductor() and \
                    ZZ(rec.get("L_num_coeffs", 0)) >= LF._coeffs_needed(digits):
                rec = dict(rec)
                rec["L_cached"] = True
                return rec
        num_coeffs = LF._coeffs_needed(digits)
        # Build and cache the coefficients once for both evaluations.
        LF.lfun(num_coeffs, digits)
        rec = {"standard_conductor": int(LF.conductor()),
               "feq_check": float(LF.check_functional_equation(num_coeffs,
                                                               digits)),
               "L_second_derivative": LF.value(1, 2, num_coeffs,
                                                      digits).str(truncate=False),
               "L_digits": int(digits), "L_num_coeffs": int(num_coeffs),
               "L_pari_precision": True}
        # the functional equation is the only independent check on the value
        # here, and it is a sharp one: with too few coefficients it fails
        # first.  A ``digits`` that the coefficients do not support is worth
        # saying out loud rather than recording as achieved.
        if rec["feq_check"] > 10.0 ** (-int(digits)):
            rec["L_precision_warning"] = (
                "the functional equation holds only to %.1e, short of the %d "
                "digits asked for" % (rec["feq_check"], int(digits)))
    except NotImplementedError as exc:
        rec = {"L_second_derivative": None, "L_status": str(exc),
               "L_digits": int(digits)}
    except Exception as exc:
        rec = {"L_second_derivative": None,
               "L_status": "%s: %s" % (type(exc).__name__, exc),
               "L_digits": int(digits)}
    if local_data is not None:
        rec['L_local_data'] = local_data
    save_lvalue(label, rec, path)
    rec = dict(rec)
    rec["L_cached"] = False
    return rec


def _certified_lvalue_cached(label, C, digits, path, force, h1_conductor, progress=None):
    """Numerical cache records can never satisfy a certified request."""
    from hyperell_regulator.regulator.balls import narrow
    key = label + ':certified'
    local_data = local_data_metadata(label)
    try:
        LF = StandardLFunction(C, h1_conductor=h1_conductor,
                              **lfunction_options(label, C))
        saved = load_lvalues(path).get(key, {})
        if not force and saved.get('L_certified') and saved.get('L_backend') == BACKEND_REVISION \
                and saved.get('L_local_data') == local_data \
                and saved.get('L_analytic_assumptions') == ANALYTIC_ASSUMPTIONS \
                and saved.get('L_curve') == str(C) \
                and saved.get('standard_conductor') == LF.conductor() \
                and saved.get('L_digits', 0) >= digits:
            try:
                if narrow(decode_ball(saved['L_second_derivative']), digits):
                    return dict(saved, L_cached=True)
            except (ValueError, KeyError, TypeError):
                pass
        result = LF.value(1, 2, digits=digits, certified=True,
                          **({'progress':progress} if progress is not None else {}))
        metadata = LF._certified_metadata
        rec = dict(standard_conductor=int(LF.conductor()), feq_check=None,
                   L_second_derivative=encode_ball(result), L_digits=int(digits),
                   L_num_coeffs=metadata['coefficient_cutoff'], L_certified=True,
                   L_backend=BACKEND_REVISION, L_analytic_assumptions=ANALYTIC_ASSUMPTIONS)
        rec['L_metadata'] = dict(metadata)
        rec['L_curve'] = str(C)
    except Exception as exc:
        rec = dict(L_second_derivative=None, L_status='%s: %s' % (type(exc).__name__, exc),
                   L_digits=int(digits), L_certified=False)
        if 'LF' in locals() and getattr(LF, '_certified_metadata', {}).get('attempts'):
            rec['L_metadata'] = LF._certified_metadata
    if local_data is not None:
        rec['L_local_data'] = local_data
    save_lvalue(key, rec, path)
    return dict(rec, L_cached=False)


def beilinson_comparison(det, c_plus, c_minus, lval, bound=8, digits=None, certified=False):
    r"""A rational approximation to `q=\pi c^+c^-L^{\prime\prime}/\det R`.

    Return ``{}`` for missing data or `\det R=0`. A relation ``[m,n]`` means
    ``m*norm+n*det=0``, giving ``q=-n/m``.
    """
    if lval is None or c_plus is None or c_minus is None:
        return {}
    if certified:
        from hyperell_regulator.regulator.balls import CertificationError
        values = [decode_ball(v) if isinstance(v, dict) else v
                  for v in (det, c_plus, c_minus, lval)]
        if any(not hasattr(v, 'rad') for v in values):
            raise CertificationError('certified comparison requires balls for every input')
        det, c_plus, c_minus, lval = values
        if det.contains_zero():
            return dict(comparison_error='determinant enclosure contains zero')
        F = det.parent()
        norm = F.pi()*F(c_plus)*F(c_minus)*F(lval)
        ratio = norm/det
        out = dict(normalisation="pi c+ c- L''(std,1) / D", norm_value=norm,
                   ratio=ratio, certified=True, exact_relation_proved=False,
                   analytic_assumptions=ANALYTIC_ASSUMPTIONS)
        # Relation discovery uses the midpoint. Only the ball residual and
        # agreement bounds are certified; containment of zero is not equality.
        try:
            P = RealField(F.precision())
            rel = [Integer(t) for t in gp.lindep([P(ratio.mid()), P(1)], bound)]
            if rel and rel[0]:
                q = -QQ(rel[1])/QQ(rel[0])
                out.update(relation=rel, rational=q, rational_is_integer=q.denominator() == 1,
                           residual=abs(rel[0]*norm+rel[1]*det))
                if q:
                    agreement = abs(ratio/q-1)
                    out.update(relative_agreement=agreement,
                               relative_agreement_upper=str(agreement.upper().exact_rational()))
        except Exception as exc:
            out['relation_error'] = str(exc)
        return out
    # in RR the residual could not be smaller than about 1e-14 whatever the
    # ingredients are worth, so the field follows ``digits``
    F = RR if digits is None else RealField(int(3.33 * int(digits)) + 24)
    cv = (lambda v: RR(v)) if digits is None else (lambda v: _sage_real(v, digits))
    det, c_plus, c_minus, lval = (cv(det), cv(c_plus), cv(c_minus), cv(lval))
    norm = F(pi) * c_plus * c_minus * lval
    if det == 0 or norm == 0:
        return {}
    out = {"normalisation": "pi c+ c- L''(std,1) / D",
           "norm_value": norm, "ratio": norm / det}
    try:
        rel = [Integer(t) for t in gp.lindep([out["ratio"], F(1)], bound)]
        out["relation"] = rel
        if rel and rel[0] != 0:
            out["rational"] = -QQ(rel[1]) / QQ(rel[0])
            out["rational_is_integer"] = out["rational"].denominator() == 1
            out["residual"] = abs(rel[0] * norm + rel[1] * det)
            if out["rational"]:
                out["relative_agreement"] = abs(out["ratio"] / out["rational"] - 1)
    except Exception as exc:
        out["relation_error"] = str(exc)
    return out


def load_stored(out=RESULTS_OUT):
    """Saved reports keyed by curve label."""
    if out is None:
        return {}
    try:
        with open(out) as fh:
            return {r["label"]: r for r in json.load(fh)}
    except (OSError, ValueError):
        return {}


def _fetch_curves(max_cond, cache, end_algs):
    """Yield ready curves from one I/O thread; keep all Sage work in the caller."""
    from queue import Queue, Empty, Full
    from threading import Event, Thread
    from hyperell_regulator.components import fetch_lmfdb, load_curves

    events, cancelled = Queue(maxsize=64), Event()

    def send(kind, value=None):
        while not cancelled.is_set():
            try:
                events.put((kind, value), timeout=.1)
                return
            except Full:
                continue
        raise InterruptedError('LMFDB fetch cancelled')

    def fetch():
        try:
            fetch_lmfdb(max_cond=max_cond-1, path=cache, geom_simple=False,
                        end_algs=end_algs, on_curve=lambda row: send('curve', row),
                        _cancel_event=cancelled)
            send('done')
        except BaseException as exc:
            if not cancelled.is_set():
                try:
                    send('error', exc)
                except InterruptedError:
                    pass

    worker = Thread(target=fetch, name='beilinson-fetch', daemon=True)
    worker.start()
    seen = set()
    try:
        while True:
            try:
                kind, value = events.get(timeout=.5)
            except Empty:
                continue
            if kind == 'error':
                raise value
            if kind == 'done':
                break
            if value['label'] not in seen:
                seen.add(value['label'])
                yield value
        # Include previously cached rows outside the refreshed catalogue,
        # preserving the usual cache selection without evaluating duplicates.
        rows = load_curves(path=cache, max_cond=max_cond-1, simple_only=False)
        if end_algs is not None and any(r.get('end_alg') is None for r in rows):
            raise ValueError('End_Q(J) data missing from the curve cache; rerun with --fetch')
        for row in rows:
            if row['label'] not in seen and (end_algs is None or row.get('end_alg') in end_algs):
                yield row
    finally:
        cancelled.set()
        worker.join(timeout=1)


def _format_divisor(coefficients, labels):
    terms = []
    for coefficient, label in zip(map(QQ, coefficients), labels):
        if not coefficient:
            continue
        magnitude = abs(coefficient)
        term = ('%s*' % magnitude if magnitude != 1 else '') + label
        sign = ' - ' if coefficient < 0 else ' + '
        terms.append(sign + term if terms else ('-' if coefficient < 0 else '') + term)
    return ''.join(terms) or '0'


@contextmanager
def _monitored_result(result, verbose, max_edge_seconds, max_edge_steps):
    """Record progress, failures and integration costs."""
    from hyperell_regulator.regulator import IntegrationMonitor

    def progress(event):
        if verbose and event["event"] != "start":
            print("  %s edge %d %s: %.1fs, %s" %
                  (result["label"], event["index"], event["event"],
                   event["seconds"], event["phase"]), flush=True)

    monitor = IntegrationMonitor(max_edge_seconds=max_edge_seconds,
                                 max_edge_steps=max_edge_steps,
                                 callback=progress, progress_interval=5)
    started = time.monotonic()
    try:
        with monitor:
            yield
    except Exception as exc:
        result.update(status="preparation_error",
                      error="%s: %s" % (type(exc).__name__, exc))
    result.update(seconds=time.monotonic()-started,
                  edges_attempted=len(monitor.records),
                  slowest_edge_seconds=max((e["seconds"] for e in monitor.records), default=0))


def _report_path(default, certified):
    path = Path(default)
    return str(path.with_name(path.stem+'_certified'+path.suffix)) if certified else default


def _endomorphism_filter(end_algs):
    if isinstance(end_algs, str):
        end_algs = (end_algs,)
    return None if end_algs is None else tuple(sorted(set(end_algs)))


def _endomorphism_condition(end_algs):
    return None if end_algs is None else 'end_alg in {%s}' % ','.join(end_algs)


def survey_curves(max_cond=DEFAULT_CONDUCTOR_BOUND, *, cache=CACHE, end_algs=BEILINSON_END_ALGS):
    """Cached curve metadata below the strict conductor bound, without pair records."""
    end_algs = _endomorphism_filter(end_algs)
    rows = load_curves(path=cache, max_cond=int(max_cond)-1, simple_only=False)
    if end_algs is not None and any(r.get('end_alg') is None for r in rows):
        raise ValueError('End_Q(J) data missing from the curve cache; rerun with --fetch')
    return [r for r in rows if end_algs is None or r.get('end_alg') in end_algs]


def _curve_record(curve):
    """Metadata for an actual curve, with a model-specific identity if unlabelled."""
    if not all(hasattr(curve, name) for name in ('f', 'h', 'F', 'bad_primes')):
        raise TypeError('expected a genus-two curve with reduction, or an LMFDB label')
    if curve.genus() != 2:
        raise ValueError('the Beilinson experiment requires a genus-two curve')
    equation = [list(map(str, curve.f)), list(map(str, curve.h))]
    row = dict(getattr(curve, 'lmfdb_record', {}))
    if 'eqn' in row:
        reference = curve_from_row(row)
        if reference.f != curve.f or reference.h != curve.h:
            raise ValueError('curve no longer matches its LMFDB metadata')
    row['eqn'] = equation
    if not row.get('label'):
        fingerprint = hashlib.sha256(json.dumps(equation).encode()).hexdigest()
        row['label'] = 'custom.'+fingerprint
    if row.get('cond') is None:
        row['cond'] = int(curve.conductor())
    return row


def analyse(curve, *, cache=CACHE, verbose=True, with_periods=True, with_lvalue=True,
            digits=15, lindep_bound=8, reg_digits=None, tries=1,
            ignore_integrality=False, certified=False, end_algs=None,
            packet='geometric', geometric_seconds=DEFAULT_GEOMETRIC_SECONDS):
    """Compute divisor bases, regulators and the L-value for one curve.

    ``curve`` is a curve with reduction, an LMFDB label, or a complete
    LMFDB metadata row. Explicit curves are not subject to the survey's
    endomorphism filter unless ``end_algs`` is supplied. Divisor bases
    contain only Galois traces over QQ.
    """
    if packet not in ('rational', 'weierstrass', 'geometric'):
        raise ValueError('packet must be rational, weierstrass or geometric')
    if not (isfinite(geometric_seconds) and geometric_seconds > 0):
        raise ValueError('geometric seconds must be positive and finite')
    if isinstance(curve, str):
        curve = load_curve(curve, cache=cache)
    if isinstance(curve, dict):
        curve = curve_from_row(curve)
    row = _curve_record(curve)
    end_algs = _endomorphism_filter(end_algs)
    if end_algs is not None and row.get('end_alg') not in end_algs:
        status = ('End_Q(J) data missing; refresh the curve cache with --fetch'
                  if row.get('end_alg') is None else 'end_alg not in {%s}' % ','.join(end_algs))
        return dict(label=row['label'], cond=row['cond'], status=status)
    from ._beilinson_divisors import analyse_divisors
    return analyse_divisors(row, curve, row, packet, verbose=verbose,
                           with_periods=with_periods, with_lvalue=with_lvalue,
                           digits=digits, reg_digits=reg_digits, lindep_bound=lindep_bound,
                           tries=tries, certified=certified, geometric_seconds=geometric_seconds,
                           ignore_integrality=ignore_integrality)


def _experiment_settings(*, digits, reg_digits, certified, with_lvalue, with_periods,
                         tries, packet, geometric_seconds, max_edge_seconds, max_edge_steps,
                         ignore_integrality, end_algs, curve_seconds=None):
    from ._kernel_cache import KERNEL_CACHE_VERSION
    if int(digits) < 1 or int(tries) < 1:
        raise ValueError('digits and tries must be positive')
    if certified and int(digits) < 5:
        raise ValueError('certified precision must be at least 5 decimal digits')
    if certified and reg_digits is None:
        reg_digits = digits
    if reg_digits is not None and int(reg_digits) < 5:
        raise ValueError('regulator precision must be at least 5 decimal digits')
    if packet not in ('rational', 'weierstrass', 'geometric'):
        raise ValueError('packet must be rational, weierstrass or geometric')
    if not (isfinite(geometric_seconds) and geometric_seconds > 0):
        raise ValueError('geometric seconds must be positive and finite')
    if not (isfinite(max_edge_seconds) and max_edge_seconds > 0) or int(max_edge_steps) < 1:
        raise ValueError('edge budgets must be positive and finite')
    if curve_seconds is None:
        geometry_budget = 3*geometric_seconds if packet == 'geometric' else 0
        curve_seconds = max(DEFAULT_CURVE_SECONDS, geometry_budget+2*max_edge_seconds, 6*int(digits))
    if not isfinite(curve_seconds) or curve_seconds < 0:
        raise ValueError('curve seconds must be nonnegative and finite')
    settings = dict(report_version=REPORT_VERSION, regulator_digits=reg_digits,
                    lvalue_digits=int(digits), with_lvalue=with_lvalue, with_periods=with_periods,
                    endomorphism_condition=_endomorphism_condition(end_algs),
                    ignore_integrality=ignore_integrality, tries=int(tries),
                    max_edge_seconds=max_edge_seconds, max_edge_steps=int(max_edge_steps),
                    curve_seconds=curve_seconds,
                    packet=packet, divisor_normalization='integer',
                    integral_basis_cache_version=KERNEL_CACHE_VERSION)
    if certified:
        settings['certified'] = True
    if packet == 'geometric':
        settings['geometric_seconds'] = geometric_seconds
    if ignore_integrality:
        settings['comparison_method'] = 'linear_determinant_relation_v2'
    return settings


def _stored_result_satisfies(old, settings):
    """Reuse successful records only under matching policies and precision."""
    if old.get('status') != 'ok' or not isinstance(old.get('settings'), dict):
        return False
    if old.get('regulator_method') == 'isogeny_transport':
        return False
    if settings.get('with_lvalue') and old.get('L_local_data') != local_data_metadata(old.get('label')):
        return False
    stored = old['settings']
    if bool(stored.get('certified', False)) != bool(settings.get('certified', False)):
        return False
    if settings.get('certified'):
        from hyperell_regulator.regulator.balls import narrow
        try:
            if (not old.get('certified') or old.get('regulator_cut') != 'tree'
                    or not narrow(decode_ball(old['det']), settings['regulator_digits'])):
                return False
            if settings.get('with_lvalue') and (not old.get('L_certified')
                    or old.get('L_backend') != BACKEND_REVISION
                    or not narrow(decode_ball(old['L_second_derivative']), settings['lvalue_digits'])):
                return False
            if settings.get('with_periods', settings.get('with_lvalue')) and any(
                    not narrow(decode_ball(old[key]), settings['regulator_digits'])
                    for key in ('c_plus', 'c_minus')):
                return False
        except (KeyError, ValueError, TypeError):
            return False
    precision = {'regulator_digits', 'lvalue_digits'}
    execution = {'max_edge_seconds', 'max_edge_steps', 'tries', 'geometric_seconds', 'curve_seconds'}
    if (settings.get('packet') == 'geometric' and old.get('packet_support_complete') is False
            and stored.get('geometric_seconds', 0) < settings.get('geometric_seconds', 0)):
        return False
    if not settings.get('ignore_integrality') and stored.get('tries', 0) < settings.get('tries', 0):
        return False
    if any(stored.get(k, False if k in ('ignore_integrality', 'certified') else None) != v
           for k, v in settings.items() if k not in precision | execution):
        return False
    if not precision <= stored.keys():
        return False
    have, want = stored['regulator_digits'], settings['regulator_digits']
    if have is None:
        if want is not None:
            return False
    elif int(have) < (16 if want is None else int(want)):
        return False
    if settings.get('with_lvalue') and (int(stored['lvalue_digits']) < int(settings['lvalue_digits'])
            or old.get('L_second_derivative') is None):
        return False
    return True


def _report_header():
    print((REPORT_COLUMNS+'  %s') % ('curve', 'cond', 'status', 'integral',
          'determinant', 'q', 'rational', 'rel. error', 'torsion divisors'))


def _report_divisors(result):
    if result.get('status') != 'ok':
        return '-'
    terms = (result.get('linear_comparison') or {}).get('terms')
    if terms:
        search = result['pair_searches'][result['pair_search_index']]
        kernel = result['integral_kernels'][search['kernel_index']]
        labels = [_fmt(P['coordinates']) if P['field_polynomial'] is None else 'P%d' % i
                  for i, P in enumerate(kernel['points'], 1)]
        indices = sorted({i for term in terms for i in term['basis_indices']})
        return ' | '.join('D%d = %s' % (i+1, _format_divisor(kernel['divisor_basis'][i], labels))
                          for i in indices)
    if result.get('divisors') and result.get('points'):
        labels = [_fmt(P['coordinates']) if P['field_polynomial'] is None else 'P%d' % i
                  for i, P in enumerate(result['points'], 1)]
        rows = result['divisors']
    else:
        return '-'
    return ' | '.join('D%d = %s' % (i, _format_divisor(row, labels))
                      for i, row in enumerate(rows, 1))


def _report_row(result, stored=False):
    comparison = result.get("comparison") or {}
    component = "-"
    if 'integral' in result:
        component = {True: 'integral', False: 'nonintegral', None: 'undecided'}[result['integral']]
    def number(value):
        return "-" if value is None else "%.8g" % display_float(value)
    print((REPORT_COLUMNS+"  %s%s") % (
        result["label"], result.get("cond", "?"), result["status"],
        component,
        number(result.get("det")), number(comparison.get("ratio")),
        str(comparison.get("rational", "-")),
        number(comparison.get("relative_agreement")),
        _report_divisors(result),
        " [stored]" if stored else ""), flush=True)
    if result.get('packet_used') == 'geometric' and 'integral_rank' in result:
        rank = result['integral_rank']
        upper = result.get('integral_rank_upper_bound')
        detail = (str(rank) if result['integral_rank_complete'] else
                  '%d to %d' % (rank, upper) if upper is not None else 'at least %d' % rank)
        print('  Geometric packet integral rank ' + detail, flush=True)
    linear = result.get('linear_comparison') or {}
    if linear.get('terms'):
        terms = linear['terms']
        expression = _format_divisor(
            [term['coefficient'] for term in terms],
            ['det(D%d,D%d)' % tuple(i+1 for i in term['basis_indices']) for term in terms])
        print("  pi*c+*c-*L'' ~= %s   relative residual %s" %
              (expression, number(linear['relative_residual'])), flush=True)
    elif result.get('pair_searches') and result.get('L_second_derivative') is not None:
        print('  ' + linear.get('error', 'No linear relation found'), flush=True)
    for key in ("rational_kernel_error", "geometric_issue", "regulator_issue", "error",
                "periods_error", "L_status", "L_precision_warning"):
        if result.get(key):
            print("  %s: %s" % (key, result[key]), flush=True)


def _analyse_monitored(curve, row, settings, *, verbose, end_algs):
    result = dict(label=row['label'], cond=row['cond'], input_curve=row)
    with _monitored_result(result, verbose, settings['max_edge_seconds'], settings['max_edge_steps']):
        result.update(analyse(curve, verbose=verbose, with_periods=settings['with_periods'],
                              with_lvalue=settings['with_lvalue'], digits=settings['lvalue_digits'],
                              reg_digits=settings['regulator_digits'], tries=settings['tries'],
                              lindep_bound=max(1, settings['lvalue_digits']//2),
                              packet=settings['packet'], end_algs=end_algs,
                              geometric_seconds=settings.get('geometric_seconds', DEFAULT_GEOMETRIC_SECONDS),
                              ignore_integrality=settings['ignore_integrality'],
                              certified=settings.get('certified', False)))
    return result


def _run_one(curve, row, settings, *, out, force, verbose, end_algs, have):
    label = row['label']
    old = have.get(label, {})
    if not force and _stored_result_satisfies(old, settings):
        _report_row(old, stored=True)
        return old
    if verbose:
        print('Computing %s ...' % label, flush=True)
    if settings['curve_seconds']:
        from ._curve_worker import run as run_worker

        result = dict(label=label, cond=row['cond'], input_curve=row)
        result.update(run_worker(row, settings, verbose=verbose, end_algs=end_algs))
    else:
        result = _analyse_monitored(curve, row, settings, verbose=verbose, end_algs=end_algs)
    result['settings'] = settings
    have[label] = result
    if out is not None:
        _save_json(out, sorted(have.values(), key=lambda r: (int(r['cond']), r['label'])))
    _report_row(result)
    return result


def run_curve(curve, *, cache=None, fetch=True, refresh=False, out=None, force=False,
              digits=15, reg_digits=None, with_lvalue=True, with_periods=None,
              tries=1, verbose=False, max_edge_seconds=15, max_edge_steps=100000,
              ignore_integrality=False, certified=False, packet='geometric',
              geometric_seconds=DEFAULT_GEOMETRIC_SECONDS, curve_seconds=0):
    """Run the Beilinson experiment on one curve or LMFDB label.

    Return a result dictionary with timings and precision settings; ``out``
    merges it by label into a JSON report. Uncached labels are fetched
    individually. ``curve_seconds`` sets a wall-clock budget.
    """
    settings = _experiment_settings(digits=digits, reg_digits=reg_digits, certified=certified,
        with_lvalue=with_lvalue, with_periods=with_lvalue if with_periods is None else with_periods,
        tries=tries, packet=packet, geometric_seconds=geometric_seconds,
        max_edge_seconds=max_edge_seconds, max_edge_steps=max_edge_steps,
        ignore_integrality=ignore_integrality, end_algs=None, curve_seconds=curve_seconds)
    if isinstance(curve, str):
        curve = load_curve(curve, cache=cache, fetch=fetch, refresh=refresh)
    if isinstance(curve, dict):
        curve = curve_from_row(curve)
    row = _curve_record(curve)
    _report_header()
    return _run_one(curve, row, settings, out=out, force=force, verbose=verbose,
                    end_algs=None, have=load_stored(out))


def run(max_cond=DEFAULT_CONDUCTOR_BOUND, *, cache=CACHE, out=RESULTS_OUT,
        force=False, digits=15, reg_digits=None, fetch=False, with_lvalue=True,
        tries=1, verbose=False, max_edge_seconds=15, max_edge_steps=100000,
        ignore_integrality=False, certified=False, end_algs=BEILINSON_END_ALGS,
        packet='geometric', geometric_seconds=DEFAULT_GEOMETRIC_SECONDS, curve_seconds=None):
    """Run the divisor-basis survey below a strict conductor bound."""
    from collections import Counter
    end_algs = _endomorphism_filter(end_algs)
    if int(max_cond) < 1:
        raise ValueError('conductor bound must be positive')
    settings = _experiment_settings(digits=digits, reg_digits=reg_digits, certified=certified,
        with_lvalue=with_lvalue, with_periods=with_lvalue, tries=tries, packet=packet,
        geometric_seconds=geometric_seconds, max_edge_seconds=max_edge_seconds,
        max_edge_steps=max_edge_steps, ignore_integrality=ignore_integrality, end_algs=end_algs,
        curve_seconds=curve_seconds)
    rows = (_fetch_curves(int(max_cond), cache, end_algs) if fetch else
            survey_curves(int(max_cond), cache=cache, end_algs=end_algs))
    print('%s curves with conductor < %d; %s packet, Galois traces over Q' %
          ('Fetching and analysing' if fetch else '%d cached' % len(rows), int(max_cond), packet))
    _report_header()
    have, results = load_stored(out), []
    with closing(rows) if fetch else nullcontext():
        for row in rows:
            results.append(_run_one(row, row, settings, out=out, force=force, verbose=verbose,
                                    end_algs=end_algs, have=have))
    print('\n'+'; '.join('%s: %d' % item for item in sorted(Counter(r['status'] for r in results).items())))
    print("q = pi*c+*c-*L''(std,1) / D; rel. error = abs(q/rational - 1).")
    if certified:
        print('Certified interval agreement; exact rational equality is not proved.')
        print(ANALYTIC_ASSUMPTIONS)
    if out is not None:
        print('Results: %s' % out)
    return results


def main(argv=None):
    """Run a conductor-bounded survey or a single LMFDB label."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('N', nargs='?', type=int, help='exclusive survey conductor bound (default: 1000)')
    parser.add_argument('--label', help='run only this LMFDB genus-two curve label')
    parser.add_argument('--cache', default=None, help='LMFDB curve cache path')
    parser.add_argument('--fetch', action='store_true', help='refresh LMFDB metadata and points')
    parser.add_argument('--cm-only', action='store_const', dest='end_algs', const=('CM',),
                        default=BEILINSON_END_ALGS, help='survey only: select rational end_alg=CM')
    parser.add_argument('--certified', action='store_true', help='certified regulators and conditional L-values')
    parser.add_argument('--force', action='store_true', help='recompute saved successful experiments')
    parser.add_argument('--no-lvalue', action='store_true', help='only compute regulator determinants')
    parser.add_argument('--verbose', action='store_true', help='show divisor and per-edge progress')
    parser.add_argument('--digits', type=int, default=15, help='L-value precision (default: 15)')
    parser.add_argument('--reg-digits', '--prec', type=int, default=None, help='regulator precision (default: float64)')
    parser.add_argument('--edge-seconds', type=float, default=15, help='time limit per integration edge')
    parser.add_argument('--edge-steps', type=int, default=100000, help='continuation steps per edge')
    parser.add_argument('--out', default=None, help='JSON report destination')
    parser.add_argument('--ignore-integrality', action='store_true', help='compare all Galois-trace divisor pairs')
    parser.add_argument('--packet', choices=('rational', 'weierstrass', 'geometric'), default='geometric')
    parser.add_argument('--geometric-seconds', type=float, default=DEFAULT_GEOMETRIC_SECONDS)
    parser.add_argument('--curve-seconds', type=float, default=None,
                        help='wall-clock budget per curve (default: at least 180s; 0 disables it)')
    parser.add_argument('--tries', type=int, default=1, help='number of integral regulator pairs to compare')
    args = parser.parse_args(argv)
    if args.label and args.N is not None:
        parser.error('choose either --label or a survey conductor bound')
    if args.label and args.end_algs != BEILINSON_END_ALGS:
        parser.error('--cm-only applies to surveys; --label selects one curve directly')
    options = dict(force=args.force, with_lvalue=not args.no_lvalue, verbose=args.verbose,
        digits=args.digits, reg_digits=args.reg_digits, tries=args.tries,
        ignore_integrality=args.ignore_integrality, packet=args.packet, certified=args.certified,
        geometric_seconds=args.geometric_seconds, max_edge_seconds=args.edge_seconds,
        max_edge_steps=args.edge_steps, curve_seconds=args.curve_seconds)
    default = data_path('results', 'beilinson_'+args.label+'.json') if args.label else RESULTS_OUT
    options['out'] = args.out if args.out is not None else _report_path(default, args.certified)
    try:
        _experiment_settings(digits=args.digits, reg_digits=args.reg_digits, certified=args.certified,
            with_lvalue=not args.no_lvalue, with_periods=not args.no_lvalue, tries=args.tries,
            packet=args.packet, geometric_seconds=args.geometric_seconds,
            max_edge_seconds=args.edge_seconds, max_edge_steps=args.edge_steps,
            ignore_integrality=args.ignore_integrality, end_algs=_endomorphism_filter(args.end_algs),
            curve_seconds=args.curve_seconds)
        if args.N is not None and args.N < 1:
            raise ValueError('conductor bound must be positive')
    except ValueError as exc:
        parser.error(str(exc))
    if args.label:
        return run_curve(args.label, cache=args.cache, refresh=args.fetch, **options)
    return run(DEFAULT_CONDUCTOR_BOUND if args.N is None else args.N,
               cache=CACHE if args.cache is None else args.cache,
               fetch=args.fetch, end_algs=args.end_algs, **options)


if __name__ == '__main__':
    main()
