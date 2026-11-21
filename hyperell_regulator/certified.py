"""Lossless storage of rigorous balls and common certification metadata."""

ANALYTIC_ASSUMPTIONS = (
    'The supplied Euler factors, conductor and gamma factors are correct; '
    'the completed L-function is holomorphic, has the growth required by '
    'the generic L-function algorithm, and satisfies its functional equation. '
    'Ball arithmetic certifies numerical errors conditional on these hypotheses.'
)
BACKEND_REVISION = 'lfunctions-42014e1-hyperell-2'


def encode_ball(value):
    """Exact rational endpoints, including outward rounding on conversion."""
    return dict(type='real_ball', bits=int(value.parent().precision()),
                lower=str(value.lower().exact_rational()),
                upper=str(value.upper().exact_rational()))


def decode_ball(record):
    from sage.all import QQ, RealBallField
    if not isinstance(record, dict) or record.get('type') != 'real_ball':
        raise ValueError('expected an exact real-ball record')
    lo, hi = QQ(record['lower']), QQ(record['upper'])
    if lo > hi:
        raise ValueError('reversed ball endpoints')
    R = RealBallField(int(record['bits']))
    return R((lo+hi)/2).add_error(R((hi-lo)/2))


def display_float(value):
    """Approximate display only; never use in certified calculations."""
    if isinstance(value, dict):
        value = decode_ball(value)
    return float(value)


def encode_result(value):
    """Recursively preserve ball fields in experiment JSON results."""
    if hasattr(value, 'rad') and hasattr(value, 'lower'):
        return encode_ball(value)
    if isinstance(value, dict):
        return {key: encode_result(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode_result(v) for v in value]
    return value
