r"""
Beilinson's conjecture for `K_2` of hyperelliptic curves, numerically.

Three things, for a hyperelliptic (in practice genus 2) curve `C/\QQ`:

* :mod:`hyperell_regulator.reduction` -- semistable reduction and the
  *component question*: given `P, Q \in C(\QQ)` and a prime `p`, do `P` and
  `Q` reduce to a smooth point of the same irreducible component of the
  special fibre at `p`?  Cluster pictures ([DDMM]_, through
  ``sage_cluster_pictures``) at odd primes, Liu's algorithm at 2.
  :mod:`hyperell_regulator.components` runs this over the LMFDB list.

* :mod:`hyperell_regulator.lfunction` -- the standard `L`-function
  `L(\mathrm{std}, s)`, the degree 5 piece of `\wedge^2 h^1(C)`, with its bad
  Euler factors read off the cluster pictures, and `L''(\mathrm{std}, 1)`.

* :mod:`hyperell_regulator.regulator` -- the Beilinson regulator of the
  `K_2` element `\{f, g\}` attached to a pair of 2-torsion (Weierstrass)
  points, as a numerical integral, in four independent ways.

:mod:`hyperell_regulator.experiments` puts the three together and compares
the regulator determinant with the `L`-value.

EXAMPLES:

The component question::

    sage: from hyperell_regulator import HyperellipticCurveWithReduction
    sage: R.<x> = ZZ[]
    sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
    sage: C.bad_primes()
    [2, 53, 61]
    sage: C.same_component((0, 1), (0, -1), 53)["answer"]
    True

.. [DDMM] T. Dokchitser, V. Dokchitser, C. Maistret, A. Morgan, *Arithmetic
   of hyperelliptic curves over local fields*, Math. Ann. 385 (2023).
"""

from hyperell_regulator.reduction import (
    HyperellipticCurveWithReduction,
    HyperellipticReduction,
    genus2_semistable_report,
    self_test,
)
from hyperell_regulator.lfunction import (
    CurveLFunction,
    StandardLFunction,
    curve_lvalue,
    standard_lvalue,
    wedge2_euler_factor,
)

__version__ = "0.1.0"

__all__ = [
    "HyperellipticCurveWithReduction",
    "HyperellipticReduction",
    "CurveLFunction",
    "StandardLFunction",
    "curve_lvalue",
    "genus2_semistable_report",
    "self_test",
    "standard_lvalue",
    "wedge2_euler_factor",
]
