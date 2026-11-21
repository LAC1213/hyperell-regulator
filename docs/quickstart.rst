Periods and regulators
======================

For :math:`C:y^2=F(x)`, the numerical API integrates periods and logarithmic
weights on a common cut. Run in Sage or with ``sage -python``:

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

The columns use :math:`\log|x+1|` and :math:`\log|x|`, respectively.
Pass ``pole=b`` for :math:`\log|(x-a)/(x-b)|`; ``pole=None`` means infinity.
Coefficients and logarithmic centers are rational.

Precision and normalisation
---------------------------

``prec`` is the decimal accuracy target. Omitting it returns float64
arrays. Failed convergence checks raise exceptions; the accuracy target
is not a certified error bound.

The vectors contain plane integrals with weights
:math:`\operatorname{Re}(x^i\bar{x}^j)/|F(x)|`, for :math:`i\le j`.
In genus two the order is :math:`(1,\operatorname{Re}(x),|x|^2)`.
Double the middle entry for the determinant convention
:math:`(1,2\operatorname{Re}(x),|x|^2)`.

Stokes integration
------------------

A Euclidean minimum spanning tree joins the branch values and logarithmic
centers. A ray from its rightmost vertex extends the cut to infinity.
The two lifted boundaries give the Stokes sum. Each logarithmic jump is
supported on the tree path from its center to infinity.

Finite edges use :math:`x=a+(b-a)\sin^2(\pi t/2)`, with endpoint root
factors cancelled analytically. The ray uses a reciprocal substitution.
Second moments reduce to products of edge periods after taking the
symmetric imaginary part.

``regulator_plane(F, a)`` and ``periods_plane(F)`` compute the corresponding
integrals by independent float64 quadrature over :math:`\mathbb C`.

Cycle periods
-------------

``c_invariants`` returns signed pivoted minors of the real and imaginary
parts of the full period matrix. Full periods are twice the integrals on
lifted branch-to-branch tree edges.

In odd degree all :math:`2g` edges give an integral cycle basis. In even
degree omit the leaf edge at the branch value of largest modulus, breaking
ties by real and imaginary parts. These choices fix the period normalisation;
a change of cycle basis can change the chosen minors by rational factors.

``plane.period_matrix`` returns a complex128 array. Its tolerance parameters
select an internal accuracy target of at least 15 decimal digits.

Validated ball arithmetic
-------------------------

Supply ``certified=True`` to ``regulators_stokes``, ``periods_stokes``,
``regulator_stokes``, ``periods_plane``, ``regulator_plane``, ``c_invariants``
or ``plane.period_matrix`` for rigorous Arb enclosures. ``prec`` specifies
decimal absolute accuracy (15 when omitted); the final ball radius is checked.
Certified plane integrals use a tree with one infinity ray; set
``certified_cut="parallel"`` to select the independent comparison backend.
The L-value methods accept the same flag with their existing ``digits`` argument.
See :doc:`ball-arithmetic` for the optional native L backend and its analytic
assumptions.
