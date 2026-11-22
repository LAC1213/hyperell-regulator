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
from hyperell_regulator.integrality import TorsionPacketIntegrality, torsion_packets
from hyperell_regulator.geometric_torsion import (
    GeometricTorsionPacket, geometric_torsion_packet, install_poonen_backend,
)
from hyperell_regulator.number_field_integrality import NumberFieldTorsionPacketIntegrality
from hyperell_regulator.projected_integrality import ProjectedTorsionPacketIntegrality
from hyperell_regulator.motivic import MotivicDivisor, MotivicPacket
from hyperell_regulator.graph_boundary import FibreGraph, GraphBoundary
from hyperell_regulator.stable_model import StableModel, stable_model
from hyperell_regulator.components import load_curve


def beilinson_experiment(curve, **options):
    """Run the Beilinson experiment for one curve with reduction or LMFDB label."""
    from hyperell_regulator.experiments.beilinson import run_curve
    return run_curve(curve, **options)

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
    "TorsionPacketIntegrality",
    "torsion_packets",
    "GeometricTorsionPacket",
    "geometric_torsion_packet",
    "install_poonen_backend",
    "NumberFieldTorsionPacketIntegrality",
    "ProjectedTorsionPacketIntegrality",
    "MotivicDivisor",
    "MotivicPacket",
    "FibreGraph",
    "GraphBoundary",
    "StableModel",
    "stable_model",
    "load_curve",
    "beilinson_experiment",
]
