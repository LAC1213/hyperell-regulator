Beilinson experiments
=====================

The experiments compare the regulator determinant :math:`D` with

.. math::

   \pi c^+c^-L''(M_{\mathrm{std}}(C),1)=qD,\qquad q\in\mathbb Q^\times.

The determinant uses weights :math:`(1,2\operatorname{Re}(z),|z|^2)`.
``comparison.rational_is_integer`` tests whether the guessed rational
factor has denominator one.

Survey
------

The survey selects LMFDB records with ``end_alg="Q"`` or ``end_alg="CM"``.
Geometric simplicity is not required; other algebra labels, including
``RM``, are excluded. ``--fetch`` refreshes this metadata in older
caches and includes newly eligible curves. Saved reports record the
endomorphism condition so that former geometric-simplicity or Q-only policies
cannot satisfy a new cache request.

Add ``--cm-only`` to restrict the survey to ``end_alg="CM"`` records, for
example with certified 100-digit computations:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 1300 \
       --cm-only --certified --digits 100 --reg-digits 100 \
       --edge-seconds 300 --edge-steps 1000000

Fetch curves with conductor strictly below 1000 and test them:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 1000 --fetch

With ``--fetch``, downloading runs in the background. The experiment starts
the integral-kernel and regulator calculations as soon as each curve's metadata and rational
points are cached. Curves are analysed one at a time with a bounded prefetch
queue, and each completed result is saved immediately. Existing point data
can be used as soon as its endomorphism metadata is refreshed. Results may
appear in availability order; the saved JSON remains sorted by conductor.

Use cached data with 30-digit regulators and :math:`L`-values:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 250 \
       --prec 30 --digits 30

Regulators default to float64, :math:`L`-values to 15 decimal digits, and
edge limits to 15 seconds. ``--prec`` sets regulator precision and
``--digits`` sets :math:`L`-value precision. ``--no-lvalue`` skips
:math:`L`-values. Use ``--help`` for all options.

Add ``--certified`` to compute Arb enclosures for the regulator, periods,
determinant and conditional L-value. In this mode ``--reg-digits`` defaults
to ``--digits``. Allow a larger edge budget, for example
``--certified --digits 100 --edge-seconds 600``. Exact interval endpoints and
analytic assumptions are saved in JSON; numerical cache entries are kept
separate. The rational relation remains a conjecture, with a certified error
bound. See :doc:`ball-arithmetic` for installation and supported evaluations.

Integrality criterion
---------------------

The survey solves for degree-zero divisors with vanishing projected
boundary at all finite places. It computes the regulators of a basis on
one cut per packet and compares pairs of the resulting columns.
Each result records the exact points and basis coefficients in
``integral_kernels``, the two chosen rows in ``divisors``, and whether
the basis spans the full integral kernel. Denominators are cleared in each
basis row before computing the regulator. Scaling the two divisors by
:math:`m` and :math:`n` multiplies the determinant by :math:`mn` and divides
the reported factor :math:`q` by :math:`mn`. The terminal's ``torsion divisors``
column prints both chosen divisors with those integer coefficients in
the original curve coordinates. An algebraic point printed as ``P3``
is the third entry of the result's ``points`` array.
``--tries`` controls how many
nonzero pairs per packet are compared, in decreasing determinant magnitude.
An undecided kernel can still contain certified integral classes.

The default ``--packet rational`` partitions known rational points by torsion
differences. ``--packet weierstrass`` uses the six geometric Weierstrass
points. ``--packet geometric`` searches the rational packets first. It calls
the Poonen backend only if every rational packet has fewer than two
certified integral directions. The JSON field ``packet_used`` records the
packet selection used for the regulator. The survey reports progress during packet
discovery and the exact integral-kernel calculation. Each of these
attempts has a time budget of 30 seconds, adjustable with
``--geometric-seconds``. If arithmetic on the full packet exceeds the
budget or fails, the survey tries the points with rational abscissae
from that verified packet. Certified integration uses this support
directly. The JSON field ``packet_support_complete`` records whether all
points were included. Columns use the individual marked points.
See :doc:`integrality` for the local model and the graph pairing.

The historical triple search is available with ``--packet legacy``.
In that mode ``--exclude-undecided`` requires affirmative component
verdicts. To evaluate all torsion triples up to permutation:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 1000 \
       --ignore-component-criterion --out results/beilinson_all_triples.json

For each unordered triple, the common pole is chosen by cut conditioning.
All outcomes are saved in ``triple_results``; ``--tries`` does not limit
this mode. The top-level fields give the first success, or the first failure.
``same_component_everywhere`` is ``true``, ``false`` or ``null`` (undecided).

Numerical helpers
-----------------

.. code-block:: python

   from hyperell_regulator.experiments.beilinson import (
       beilinson_determinant_stokes, regulator_column_stokes,
   )

   det, matrix, errors = beilinson_determinant_stokes(F, P, Q, digits=30)
   column = regulator_column_stokes(F, P, digits=30)

Here ``P`` and ``Q`` are point coordinates on :math:`y^2=F(x)`.
``pole`` specifies a common finite pole; the default is infinity.
Torsion orders certify the algebraic classes independently of integration.

For computations at 2, ``h1_conductor`` supplies the exact Jacobian
conductor when PARI omits its local exponent. A trivial local
:math:`h^1` factor requires additional inertia data. High precision and
large coefficient counts can require several GiB of PARI stack memory.

The CM-labelled survey through conductor 10000 automatically supplies
manual bad factors for seven otherwise unsupported curves. Their wild
standard conductor exponents are included separately, through the
``StandardLFunction(..., conductor_exponents={p: exponent})`` argument.
The equations, polynomials and local-representation derivation are recorded
in ``reports/beilinson_cm_bad_factors.rst``; the experiment saves this local
data with each L-value and invalidates caches when it changes.

Saved results
-------------

Caches are stored in ``lmfdb_cache/`` and ``results/`` relative to
``HYPERELL_REGULATOR_DATA`` or the working directory. Successful results
are reused when their precision and packet selection suffice. Failed
computations are retried. ``--force`` recomputes successful regulators.

The `saved 100-digit survey
<https://github.com/LAC1213/hyperell-regulator/tree/main/reports>`_
contains tables, full-precision values and integration models.
