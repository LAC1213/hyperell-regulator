Periods and regulators
======================

Run the following in Sage or with ``sage -python``. Use Python syntax in
scripts: exponentiation is ``**``.

.. code-block:: python

   from sage.all import QQ, PolynomialRing
   from hyperell_regulator.regulator import (
       periods_and_regulator_stokes, c_invariants,
   )

   x = PolynomialRing(QQ, "x").gen()
   F = x**5 + (x + 1)**2
   phi = (x + 1, 1)

   periods, regulator = periods_and_regulator_stokes(F, phi)
   periods, regulator = periods_and_regulator_stokes(F, phi, prec=30)
   c_plus, c_minus = c_invariants(F, prec=30)

The curve is :math:`y^2 = F(x)`. A pair ``phi = (A, B)`` represents the
function :math:`\varphi = A(x) + B(x)y`; here it is :math:`x + 1 + y`.

Precision
---------

Without ``prec``, the regulator routines use float64 arithmetic. With
``prec=30``, they target 30 decimal digits. Use exact rational coefficients
for arbitrary precision. Convergence failures raise exceptions.

The returned period and regulator vectors have entries indexed by
:math:`(1,1), (1,2), (2,2)` in genus 2. The :math:`c^+` and :math:`c^-`
invariants come from minors of the period matrix.

Period conventions
------------------

For ``c_invariants``, finite branch points are ordered by increasing argument
about a perturbed centroid, starting at the root of largest absolute value
(ties by real, then imaginary part). Period paths use the regulator's
semilunar detours and embeddedness checks. The period matrix contains full
cycle integrals, twice the integrals along the lifted paths; ``c_plus`` and
``c_minus`` are signed pivoted minors of its real and imaginary parts.

For genus two, full periods multiply each minor by four relative to
half-periods on the same paths. Different paths can also change signs and
rational normalizations. Saved reports retain their original conventions;
experiment caches are invalidated automatically when rerun.

Integration methods
-------------------

The Stokes method lifts a cut system to the sheets of the function
:math:`\varphi`. It supports float64 and arbitrary precision. The plane
integrator uses adaptive quadrature over :math:`\mathbb{C}` and provides
an independent float64 computation.

See :mod:`hyperell_regulator.regulator.api` for the numerical entry points,
:mod:`hyperell_regulator.regulator.functions` for constructing functions
with prescribed divisors, and :mod:`hyperell_regulator.regulator.plots`
for plotting cut systems.
