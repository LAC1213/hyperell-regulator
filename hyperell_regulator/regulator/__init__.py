r"""Direct two-sheet periods and logarithmic regulators on y^2=F(x).

Use ``regulator_stokes(F, a, prec=30)`` for ``log|x-a|``, and
``regulators_stokes(F, [a,b], prec=30)`` to integrate several weights on
one shared cut. No torsion functions are needed for numerical integration.
Returned genus-two vectors use weights ``(1, Re(x), |x|^2)``.
"""
from .api import (c_invariants, periods_plane, periods_stokes,
                  periods_and_regulator_stokes, regulators_stokes,
                  regulator_plane, regulator_stokes)
from .diagnostics import IntegrationLimitError, IntegrationMonitor
from .torsion import torsion_order
from .plots import plot_cut_system

__all__ = ['c_invariants', 'periods_plane', 'periods_stokes',
           'periods_and_regulator_stokes', 'regulators_stokes',
           'regulator_plane', 'regulator_stokes', 'plot_cut_system',
           'torsion_order', 'IntegrationLimitError', 'IntegrationMonitor']
