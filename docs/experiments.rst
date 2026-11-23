Beilinson experiments
=====================

The experiments compare the regulator determinant :math:`D` with

.. math::

   \pi c^+c^-L''(M_{\mathrm{std}}(C),1)=qD,\qquad q\in\mathbb Q^\times.

The determinant uses weights :math:`(1,2\operatorname{Re}(z),|z|^2)`.
``comparison.rational_is_integer`` tests whether the guessed rational
factor has denominator one.

One curve at a time
-------------------

Load a genus-two curve by its full LMFDB label and run its experiment::

    from hyperell_regulator import load_curve, beilinson_experiment

    C = load_curve('394.a.394.1')
    assert C.bad_primes() == [2, 197]
    result = beilinson_experiment(C, digits=20, reg_digits=20)

``load_curve`` returns a Sage curve with the package's reduction methods.
The attributes ``lmfdb_label`` and ``lmfdb_record`` retain the exact label
and metadata. A cache hit needs no network request. A missing curve is
queried by label in the curve and rational-point tables, then merged into
the cache. Use ``cache=...`` to select a cache, ``fetch=False`` for cached
data only, ``refresh=True`` to refresh one label, and ``prec=...`` to set
the reduction working precision.

``beilinson_experiment`` also accepts an LMFDB label directly or an
unlabelled curve constructed with ``HyperellipticCurveWithReduction``.
Unlabelled curves have model-specific identifiers for safe cache reuse.
Explicitly requested curves are independent of the survey's conductor
bound and rational endomorphism filter. The returned dictionary contains
the divisor bases, numerical results or failure status, precision settings,
timings, and integration work. ``out=...`` optionally saves it in a JSON
report; no report is written by default. The equivalent experiment-module
entry point is ``run_curve(C, ...)``.

The command-line equivalent is:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson \
       --label 394.a.394.1 --digits 20 --reg-digits 20

This writes ``results/beilinson_394.a.394.1.json`` by default. Certified
runs use a separate ``_certified.json`` file. ``--out`` selects another
destination, and ``--fetch`` refreshes only the selected label. A survey
bound cannot be combined with ``--label``.

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

The CLI and ``run()`` run each curve in a separate process with live
elapsed-time updates. The default deadline is at least 180 seconds and
increases with the precision target and the geometric and integration
budgets. Use ``--curve-seconds`` to set a deadline explicitly, or
``--curve-seconds 0`` to disable it. The matching Python keyword is
``curve_seconds``. A timeout is saved as ``curve_timeout`` and the next
curve starts. A later run retries that curve and can reuse exact kernels
already saved. Report and L-value cache entries are written atomically.

Add ``--certified`` to compute Arb enclosures for the regulator, periods,
determinant and conditional L-value. In this mode ``--reg-digits`` defaults
to ``--digits``. Allow a larger edge budget, for example
``--certified --digits 100 --edge-seconds 600``. Exact interval endpoints and
analytic assumptions are saved in JSON; numerical cache entries are kept
separate. The rational relation remains a conjecture, with a certified error
bound. See :doc:`ball-arithmetic` for installation and supported evaluations.

Integrality criterion
---------------------

The survey solves for degree-zero Galois traces over :math:`\mathbf Q`
with vanishing projected boundary at all finite places. Algebraic points
contribute through complete conjugate orbit sums; anti-invariant
differences with zero trace are excluded. It computes the regulators of a basis on
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

``--packet rational`` partitions known rational points by torsion differences.
``--packet weierstrass`` uses the six geometric Weierstrass points.
The default ``--packet geometric`` searches the rational packets first. It calls
the Poonen backend only if every rational packet has fewer than two
eligible directions under the current policy. The JSON field ``packet_used`` records the
packet selection used for the regulator. The survey reports progress during packet
discovery and the exact integral-kernel calculation. Each of these
attempts has a time budget of 30 seconds, adjustable with
``--geometric-seconds``. If arithmetic on the full packet exceeds the
budget or fails, the survey tries the points with rational abscissae
from that verified packet. Both numerical and certified runs attempt the
full geometric kernel. The JSON fields ``integral_rank`` and
``integral_rank_complete`` distinguish an exact rank from a lower bound
in this Galois-trace space. Cached kernels record
``rank_over_Q`` from exact rational matrix rank after the Galois trace
and symmetric-ordinate projections; dependent rows contribute only once.
Saved experiments label this with ``integral_rank_field="Q"`` and
``integral_rank_space="galois_traces"``. Integer scaling of representatives
does not change this rank. Cached kernels also record
``divisor_descent="galois_trace"`` and ``base_field="Q"``; older
number-field bases and the survey results using them are recomputed.
``packet_support_complete`` records whether all points were included.
At an odd prime with a single node, the survey can compute the pair-sum
constraints using a branch-point field for an odd-thickness collision
within one factor. When a rational branch point collides with a root of
another factor, its coefficient pairs with each root of that factor.
This handles the thickness-two node at 7 for ``427.a.2989.1``.
If the computed local constraints already force the kernel to be zero,
the survey records exact rank zero and skips the remaining places.
For ``604.a.9664.1`` the constraints at 151 suffice, so no dyadic model
is needed.
If some directions survive, the remaining local calculations can use
just the points that appear in those directions. This can avoid a large
compositum while still determining the kernel on the full packet. When
available, ``integral_rank_upper_bound`` records the bound from necessary
local conditions.
Conjugate logarithmic centres give equal real moments. The survey uses
this equality to detect a vanishing determinant even when the integral
divisor rank is two or higher. Columns use the individual marked points.
See :doc:`integrality` for the local model and the graph pairing.

Use ``--ignore-integrality`` to search for a rational linear relation
between the normalized :math:`L`-value and all nonzero determinants of
pairs from a basis :math:`D_1,\ldots,D_r` of the degree-zero Galois-trace
torsion divisor space modulo the hyperelliptic identifications.

.. math::

   \pi c^+c^-L^{\prime\prime}(\mathrm{std},1)
   \approx \sum_{i<j} q_{ij}\det(D_i,D_j),
   \qquad q_{ij}\in\mathbb Q.

The survey searches separately on each packet using one integration chart
for all its basis pairs. At 20 digits it requires relative residual below
:math:`10^{-17}` for numerical candidates. The terminal prints the divisors
and the whole relation.
The JSON field ``pair_searches`` lists the determinants and basis indices
for each packet tested. The chosen relation is saved in
``linear_comparison``. The integrality calculation still caches the
integral basis and records the individual divisor verdicts in
``integral_kernels``. Any degree-zero trace basis direction can participate.
For a relation with several determinants, the integrality verdict remains
undecided if any contributing divisor lacks an affirmative verdict.
The packet options apply and ``--tries`` only limits the ordinary
integral-basis comparison.
With ``--packet geometric``, a rational-packet search without a match
also tries the full geometric packet if the rational integral rank is
below two. JSON retains both attempts.

Numerical helpers
-----------------

.. code-block:: python

   from hyperell_regulator.motivic import DivisorRegulators

   evaluator = DivisorRegulators(F, [-1, 0, 'infinity'], digits=30)
   D1, D2 = (1, 0, -1), (0, 1, -1)
   det, matrix = evaluator.determinant(D1, D2)
   period, columns = evaluator.columns((D1, D2))

Each coefficient row has degree zero and refers to the supplied abscissae
on :math:`y^2=F(x)`. All columns use a shared integration cut.
Use ``MotivicPacket`` for exact curve points with torsion verification
and integral Galois-trace divisor generators.

For computations at 2, ``h1_conductor`` supplies the exact Jacobian
conductor when PARI omits its local exponent. A trivial local
:math:`h^1` factor requires additional inertia data. High precision and
large coefficient counts can require several GiB of PARI stack memory.

The experiments automatically supply manual bad factors for known
otherwise unsupported curves, including curves with geometric CM. Their wild
standard conductor exponents are included separately, through the
``StandardLFunction(..., conductor_exponents={p: exponent})`` argument.
The equations, polynomials and derivation summaries are recorded
in ``hyperell_regulator.experiments.beilinson_local_data``; the experiment saves this local
data with each L-value and invalidates caches when it changes.

For ``3125.a.3125.1`` the standard factor at 5 is 1 and the standard
conductor is :math:`5^6`. For ``4096.b.65536.1`` the standard factor at 2
is :math:`1+2T` and the conductor is :math:`2^{12}`. The constructor accepts
exterior-square factors, so these correspond respectively to
``bad_factors={5: 1-5*T}`` and ``bad_factors={2: 1-4*T**2}``, together with
the separately supplied conductor exponents 6 and 12.

The latter curve has two elliptic components after extension at 2.
Its integral-basis computation uses both verified elliptic discs.
Full rational Galois orbits are evaluated using exact Newton-polygon
averages; individual nonrational points can remain unresolved without
preventing a complete basis over :math:`\mathbb Q`.

Saved results
-------------

Caches are stored in ``lmfdb_cache/`` and ``results/`` relative to
``HYPERELL_REGULATOR_DATA`` or the working directory. Successful results
are reused when their precision and packet selection suffice. Failed
computations are retried. Exact bases are saved per curve model and
packet in ``results/integral_bases/`` before regulator integration. This
includes rank-one and incomplete kernels and works independently of the
survey report. Numerical retries and ``--force`` reuse these bases.

The `saved 100-digit survey
<https://github.com/LAC1213/hyperell-regulator/tree/main/reports>`_
contains tables, full-precision values and integration models.

Survey report and explorer
--------------------------

Generate the survey tables from a completed run without recomputing
regulators or L-values:

.. code-block:: sh

   sage -python -m hyperell_regulator.experiments.beilinson_report \
       --input results/beilinson_results_certified.json \
       --bound 1000 --digits 100 --explorer-inputs docs/explorer_inputs.json \
       --extra-results docs/explorer_extra_results.json

The three multipage tables list the selected divisors with their exact
coefficients in the original LMFDB coordinates, the LMFDB models,
integration charts and algebraic support definitions,
and the numerical comparisons. They also record the packet torsion
exponent and the rational endomorphism algebra label. Unknown packet
orders are left blank. Saved certified intervals remain in the JSON.
The explorer regeneration and support display are described in
:doc:`explorer-guide`.
