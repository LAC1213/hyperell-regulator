"""Run the exact geometric kernel calculation in a time-limited Sage process."""

import json
import subprocess
import sys

from .. import _subprocess
from ..components import curve_from_row
from ..geometric_torsion import GeometricTorsionPacket, geometric_torsion_packet
from ..motivic import MotivicPacket
from ._kernel_cache import cached_kernel


def _kernel_in_worker(row, packet, rational_only, seconds, progress):
    data = packet.metadata['data']

    def compute():
        payload = json.dumps(dict(row=row, packet=data, rational_only=rational_only,
                                  verified=packet.complete), default=str)
        kind = 'rational-abscissa integral kernel' if rational_only else 'full geometric integral kernel'
        progress('Computing the '+kind)
        completed = _subprocess.run(
            [sys.executable, '-m', __spec__.name], input=payload, timeout=seconds,
            progress=progress, label=kind.capitalize(), interval=5)
        if completed.returncode:
            detail = completed.stderr.strip().splitlines()
            reason = detail[-1] if detail else 'worker exited with code %s' % completed.returncode
            raise ArithmeticError(reason)
        return json.loads(completed.stdout)

    # Revisit older incomplete results once under the reduced-support solver.
    return cached_kernel(packet.f, packet.h, ['geometric', data, packet.complete, rational_only],
                         compute, incomplete_revision=None if rational_only else 1)


def geometric_kernel(C, row, *, seconds, progress):
    """Try the full packet before a smaller exact arithmetic fallback."""
    progress('Loading or computing the Poonen torsion packet')
    try:
        packet = geometric_torsion_packet(C, timeout=seconds, progress=progress)
    except subprocess.TimeoutExpired:
        return None, dict(status='packet_timeout', error=f'Poonen computation exceeded {seconds:g}s')
    except (ArithmeticError, OSError, ValueError) as error:
        return None, dict(status='packet_unavailable', error=str(error))

    progress('%s geometric torsion points in the %s packet' %
             (packet.size, 'cached' if packet.metadata.get('cached') else 'verified'))
    info = dict(geometric_points=packet.size)
    has_algebraic_abscissae = any(g is not None and g.degree() > 1 for g, _ in packet.factors)
    try:
        kernel = _kernel_in_worker(row, packet, False, seconds, progress)
    except (ArithmeticError, subprocess.TimeoutExpired) as error:
        reason = (f'The integral kernel calculation exceeded {seconds:g}s'
                  if isinstance(error, subprocess.TimeoutExpired) else str(error))
        if not has_algebraic_abscissae:
            status = 'geometric_timeout' if isinstance(error, subprocess.TimeoutExpired) else 'preparation_error'
            return None, dict(info, status=status, error=reason)
        info['geometric_issue'] = reason
        progress(reason+'. Trying the rational-abscissa support')
        try:
            kernel = _kernel_in_worker(row, packet, True, seconds, progress)
        except (ArithmeticError, subprocess.TimeoutExpired) as fallback_error:
            if isinstance(fallback_error, subprocess.TimeoutExpired):
                return None, dict(info, status='geometric_timeout',
                                  error=f'The rational-abscissa kernel also exceeded {seconds:g}s')
            return None, dict(info, status='preparation_error', error=str(fallback_error))
    info['packet_support_complete'] = len(kernel['points']) == packet.size
    if not info['packet_support_complete']:
        info['geometric_issue'] = 'Used rational-abscissa support. ' + info['geometric_issue']
    return kernel, info


def _main():
    from ._beilinson_divisors import _integral_kernel
    from ._simple_geometric_kernel import simple_nodal_kernel

    data = json.load(sys.stdin)
    packet = GeometricTorsionPacket(curve_from_row(data['row']), data['packet'],
                                    metadata={'complete': data['verified']})
    kernel = None if data['rational_only'] else simple_nodal_kernel(packet)
    if data['rational_only'] and not any(g is None or g.degree() == 1 for g, _ in packet.factors):
        kernel = dict(points=[], basis=[], divisor_basis=[], divisor_integrality=[],
                      complete=True, places={})
    elif kernel is None:
        motivic = MotivicPacket.from_geometric_packet(packet, rational_abscissae=data['rational_only'])
        kernel = _integral_kernel(motivic, cache=False)
    json.dump(kernel, sys.stdout, default=str)


if __name__ == '__main__':
    _main()
