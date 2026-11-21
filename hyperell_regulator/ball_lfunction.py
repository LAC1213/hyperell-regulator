"""Conditional, rigorously bounded L-values using the pinned native backend."""
import math
import os
import hashlib
from pathlib import Path
import subprocess
import tempfile

from sage.all import QQ, ZZ, RealBallField, prime_range

from .certified import ANALYTIC_ASSUMPTIONS, BACKEND_REVISION
from ._subprocess import run as _native_run
from .regulator.balls import CertificationError, digits_value, narrow


def executable():
    path = Path(os.environ.get('HYPERELL_CERTIFIED_LVALUE',
                 Path(__file__).parent/'bin/certified-lvalue'))
    if not path.is_file() or not os.access(path, os.X_OK):
        raise CertificationError('certified L-value backend unavailable; build with '
                                 "Sage's Python: python -m hyperell_regulator.build_certified")
    return str(path.resolve())


def _run(args, progress=None, stage='ball evaluation'):
    timeout = float(os.environ.get('HYPERELL_CERTIFIED_TIMEOUT', '900'))
    try:
        result = _native_run(args, timeout=timeout, progress=progress, label=stage)
    except subprocess.TimeoutExpired as exc:
        raise CertificationError('native L-value computation exceeded %.0f seconds; '
            'adjust HYPERELL_CERTIFIED_TIMEOUT if needed' % timeout) from exc
    if result.returncode:
        raise CertificationError('native L-value backend failed (exit %d): %s' %
                                 (result.returncode, result.stderr.strip()))
    return {line.split()[0]: line.split()[1:] for line in result.stdout.splitlines() if line.strip()}


def _interval(row, bits):
    lo, hi, exponent = map(ZZ, row)
    scale = QQ(2)**exponent
    R = RealBallField(bits)
    return R((lo+hi)*scale/2).add_error(R((hi-lo)*scale/2))


def value(L, s=1, derivative=0, num_coeffs=None, digits=15, progress=None):
    """Return a real ball of absolute radius < 10^-digits, or raise.

    Supported: real, exactly binary64-representable s; derivatives 0 and 1
    away from gamma poles; standard L''(1) via the exact functional-equation
    identity sqrt(N)/(4*pi^2)*L(2). No heuristic zero/rank detection is used.
    """
    from .lfunction import StandardLFunction
    report = progress if progress is not None else lambda message: None
    digits = digits_value(digits)
    if isinstance(s, float):
        s = QQ(s.as_integer_ratio()[0])/s.as_integer_ratio()[1]
    elif hasattr(s, 'exact_rational'):
        s = s.exact_rational()
    else:
        s = QQ(s)
    point_float = float(s)
    exact_float = (QQ(point_float.as_integer_ratio()[0])/point_float.as_integer_ratio()[1]
                   if math.isfinite(point_float) else None)
    if exact_float != s:
        raise NotImplementedError('certified evaluation requires an exactly binary64-representable real argument')
    if int(derivative) != derivative or derivative < 0:
        raise ValueError('derivative must be a nonnegative integer')
    standard = isinstance(L, StandardLFunction)
    bits = math.ceil((digits+40)*math.log2(10))
    R = RealBallField(bits)
    if standard and s == 1 and derivative in (0, 1):
        L._certified_metadata = dict(backend=BACKEND_REVISION, certified=True,
            analytic_assumptions=ANALYTIC_ASSUMPTIONS, root_number=1,
            coefficient_cutoff=0, working_bits=bits)
        return R(0)
    second = standard and s == 1 and derivative == 2
    if derivative > 1 and not second:
        raise NotImplementedError('certified derivatives above order one are implemented only for standard L\'\'(1)')
    point, order = (QQ(2), 0) if second else (s, derivative)
    gamma = [-1, -1, 0, 0, 1] if standard else [0, 0, 1, 1]
    if order and any(point+mu <= 0 and (point+mu)/2 in ZZ for mu in gamma):
        raise NotImplementedError('certified first derivatives at gamma poles are not implemented')
    N = ZZ(L.conductor())
    if not 0 < N < 2**64:
        raise ValueError('native backend requires a positive conductor below 2^64')
    if num_coeffs is not None and (int(num_coeffs) != num_coeffs or num_coeffs < 2):
        raise ValueError('num_coeffs must be an integer at least two')
    key = ('ball', s, int(derivative), num_coeffs, digits)
    if key in L._cache:
        result, metadata = L._cache[key]
        L._certified_metadata = dict(metadata)
        return result
    binary = executable()
    binary_hash = hashlib.sha256(Path(binary).read_bytes()).hexdigest()
    required = math.ceil(digits*math.log2(10))
    attempts = []
    L._certified_metadata = dict(backend=BACKEND_REVISION, certified=False, attempts=[])
    # target_prec is not achieved precision. In particular, planning at a
    # lower precision can leave a coefficient tail wider than requested.
    # Refine using the measured, fully normalized radius, not a fixed count
    # of three small precision increments.
    target = max(160, 2*required+64)
    plan_bits = max(80, required+32)
    max_attempts = 20
    for attempt in range(max_attempts):
        args = [binary, str(target), str(N), 'info', str(required), '0',
                str(len(gamma)), str(float(point)), str(order)]
        planning = list(args)
        if num_coeffs is None:
            # Plan the initial coefficient range at lower precision. Using
            # the much larger internal precision here needlessly stores far
            # more prime matrices. The corrected native tail bound and final
            # radius test determine whether this range actually suffices.
            planning[1] = str(plan_bits)
        report('planning the certified coefficient range')
        native_progress = {'progress': progress} if progress is not None else {}
        if progress is None:
            info = _run(planning)
        else:
            info = _run(planning, progress=progress, stage='coefficient planning')
        full = int(info['coefficient_cutoff'][0])
        cutoff = full if num_coeffs is None else min(int(num_coeffs), full)
        primes = list(prime_range(cutoff+1))
        # Distribute independent reduction forests across the worker pool.
        # Serialising the largest batches made high-conductor, 100-digit
        # evaluations take much longer. Eight is the default; an explicit
        # environment setting overrides it in either direction.
        workers = max(1, int(os.environ.get('HYPELLFROB_THREADS', '8')))
        report('computing Euler factors at %d primes through %d (%d worker%s)' %
               (len(primes),cutoff,workers,'' if workers==1 else 's'))
        L.precompute_h1(primes, nthreads=workers, **native_progress)
        report('assembling exact Euler-factor input')
        with tempfile.TemporaryDirectory(prefix='hyperell-certified-') as directory:
            factors = Path(directory)/'factors.txt'
            with factors.open('w') as out:
                for p in primes:
                    P = L.euler_factor(p)
                    coefficients = [ZZ(P[j]) for j in range(len(gamma)+1)]
                    if P.degree() > len(gamma) or coefficients[0] != 1:
                        raise CertificationError('invalid exact Euler polynomial at %s' % p)
                    out.write(' '.join(map(str, [p]+coefficients))+'\n')
            euler_hash = hashlib.sha256(factors.read_bytes()).hexdigest()
            args[3], args[5] = str(factors), str(cutoff)
            report('native ball evaluation at %d bits (attempt %d)' % (target,attempt+1))
            data = _run(args, **native_progress)
        result = _interval(data['value'], int(data['working_bits'][0]))
        if second:
            R = result.parent()
            result *= R(N).sqrt()/(4*R.pi()**2)
        radius = QQ(result.rad())
        achieved = (int(radius.denominator().nbits()-radius.numerator().nbits())
                    if radius else None)
        attempts.append(dict(coefficient_cutoff=cutoff, target_bits=target,
                             coefficient_plan_bits=int(planning[1]),
                             working_bits=int(data['working_bits'][0]),
                             radius=str(radius), radius_display=str(result.rad()),
                             achieved_absolute_bits=achieved))
        # Preserve completed attempts even if a later native call times out
        # or coefficient precomputation fails.
        L._certified_metadata = dict(backend=BACKEND_REVISION, certified=False,
                                    attempts=list(attempts))
        report('L-value enclosure radius %s (absolute target 10^-%d)' %
               (result.rad(), digits))
        if narrow(result, digits):
            report('certified L-value radius verified')
            metadata = dict(backend=BACKEND_REVISION, certified=True,
                analytic_assumptions=ANALYTIC_ASSUMPTIONS,
                root_number=int(data['root_number'][0]), coefficient_cutoff=cutoff,
                target_bits=target, working_bits=int(data['working_bits'][0]),
                attempts=attempts,
                euler_input_sha256=euler_hash, backend_binary_sha256=binary_hash,
                flint_version=data['flint_version'][0])
            L._certified_metadata = metadata
            L._cache[key] = result, dict(metadata)
            return result
        if num_coeffs is not None:
            break
        if attempt+1 < max_attempts:
            shortfall = max(0, required-achieved)
            increment = max(64, shortfall+32)
            plan_bits += increment
            target = max(target+max(80, increment), plan_bits+64)
            report('enclosure needs refinement (%d bits short); increasing the '
                   'coefficient plan to %d bits and internal target to %d bits' %
                   (shortfall, plan_bits, target))
    last = attempts[-1]
    raise CertificationError('L-value enclosure did not achieve absolute radius < 10^-%d '
        'after %d attempts; last radius %s, coefficient cutoff %d, working precision %d bits' %
        (digits, len(attempts), last['radius_display'], last['coefficient_cutoff'], last['working_bits']))
