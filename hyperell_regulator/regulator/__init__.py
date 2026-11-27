r"""Periods and logarithmic regulator integrals on `y^2=F(x)`.

Genus 2 vectors use the weights `(1,\Re(x),|x|^2)`.
"""
from .api import (c_invariants, periods_plane, periods_stokes,
                  periods_and_regulator_stokes, regulators_stokes,
                  regulator_plane, regulator_stokes)
from .diagnostics import IntegrationLimitError, IntegrationMonitor
from .torsion import torsion_order
from .plots import plot_cut_system
from .balls import CertificationError

__all__ = ['c_invariants', 'periods_plane', 'periods_stokes',
           'periods_and_regulator_stokes', 'regulators_stokes',
           'regulator_plane', 'regulator_stokes', 'plot_cut_system',
           'torsion_order', 'IntegrationLimitError', 'IntegrationMonitor', 'CertificationError']
