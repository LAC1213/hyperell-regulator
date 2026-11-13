Periods and regulators
======================

The numerical API takes a curve :math:`y^2=F(x)` and rational logarithmic
centers. It integrates directly on the two sheets of the hyperelliptic map.
Run examples in Sage or with ``sage -python``:

.. code-block:: python

   from sage.all import QQ, PolynomialRing
   from hyperell_regulator.regulator import (
       periods_and_regulator_stokes, regulators_stokes, c_invariants,
   )

   x = PolynomialRing(QQ, "x").gen()
   F = x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 1
   periods, log_x = periods_and_regulator_stokes(F, a=0, prec=30)
   periods, columns = regulators_stokes(F, [-1, 0], prec=100)
   c_plus, c_minus = c_invariants(F, prec=30)

``columns[0]`` uses :math:`\log|x+1|`, and ``columns[1]`` uses
:math:`\log|x|`. Both columns and the ordinary periods share one set of cut
edges. No torsion functions or divisors enter numerical integration.

Pass ``pole=b`` to use :math:`\log|(x-a)/(x-b)|` instead. A finite pole adds
its logarithmic center to the shared cut and its column is subtracted.
The default ``pole=None`` means infinity.

Precision and normalization
---------------------------

Without ``prec``, results are float64 arrays; guarded quadrature is still
used internally. With ``prec=30``, results retain 30 decimal digits. Exact
rational coefficients are required, and failed convergence checks raise
exceptions. Precision targets are not certified interval bounds.

The returned vectors have entries indexed by :math:`a\le b`, for the
weights :math:`\operatorname{Re}(x^a\bar{x}^b)/|F|`. In genus two their
order is :math:`(1,\operatorname{Re}(x),|x|^2)`. Double the middle entry for
the Beilinson determinant's convention :math:`(1,2\operatorname{Re}(x),|x|^2)`.

The Stokes algorithm
--------------------

A semilunar cut connects the branch values and logarithmic centers and ends
with an outward ray to infinity. Its two lifted boundaries give the Stokes
sum. Different logarithms can have different jumps on an edge, all read
from the same cut. Second moments reduce to products of ordinary edge
periods after taking the symmetric imaginary part.

Finite edges use :math:`x=a+(b-a)\sin^2(\pi t/2)`. Endpoint root factors
cancel analytically. The infinite ray uses a fixed reciprocal substitution.
Square roots are evaluated as continuous products of linear factors.
There is no general cover inversion, chart selection, or Newton tracking.

``regulator_plane(F, a)`` and ``periods_plane(F)`` retain independent adaptive
quadrature over :math:`\mathbb C` as a slower float64 reference.

Cycle periods
-------------

``c_invariants`` uses full cycle periods, twice the integrals on lifted
branch-to-branch paths. The finite branch chain starts at the root of
largest absolute value, breaking ties by real and imaginary parts.
The invariants are signed pivoted minors of the real and imaginary parts
of the period matrix. Changing a cycle basis can change their signs.

Migration from the function-cover API
-------------------------------------

Replace ``periods_and_regulator_stokes(F, (x-a, 0))`` with
``periods_and_regulator_stokes(F, a)``. For a genus-two determinant use
``regulators_stokes(F, [a, b])`` to compute both columns together.
The general :math:`A(x)+B(x)y` numerical API, fast/slow moment switch,
cover classes and torsion-function constructors have been removed.
``torsion_order`` and ``divisor_class_order`` still certify divisor classes
for the algebraic-cycle and component computations.
