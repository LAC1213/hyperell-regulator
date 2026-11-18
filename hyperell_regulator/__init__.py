r"""Periods, regulators and standard `L`-functions of genus two curves.

The regulators arise from classes in `H^3_{\mathcal M}(C^2, \QQ(2))`.

EXAMPLES::

    sage: from hyperell_regulator import HyperellipticCurveWithReduction
    sage: R.<x> = ZZ[]
    sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
    sage: C.bad_primes()
    [2, 53, 61]
    sage: C.same_component((0, 1), (0, -1), 53)["answer"]
    True
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
