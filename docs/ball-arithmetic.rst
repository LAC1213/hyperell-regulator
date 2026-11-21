Certified ball arithmetic
==========================

The public regulator, period and L-value methods accept ``certified=True``.
They return Sage Arb balls and check the actual absolute radius against
``10^-prec`` (regulators/periods) or ``10^-digits`` (L-values). Omitting
``prec`` requests 15 decimal digits. Existing calls retain their ordinary
numerical behavior. Failed integration, uncertain nonzero pivots, missing
native dependencies and insufficient final accuracy raise an exception;
there is no fallback to an uncertified answer.

Usage and installation
-----------------------

Regulators and periods need only Sage's ball fields::

   from sage.all import QQ, PolynomialRing
   from hyperell_regulator.regulator import regulators_stokes, c_invariants
   from hyperell_regulator.regulator.plane import period_matrix

   x = PolynomialRing(QQ, 'x').gen()
   F = x**5 + (x+1)**2
   periods, columns = regulators_stokes(F, [0, -1], prec=100, certified=True)
   c_plus, c_minus = c_invariants(F, prec=100, certified=True)
   cycles = period_matrix(F, prec=100, certified=True)

The same option is available on ``periods_stokes``, ``regulator_stokes``,
``periods_and_regulator_stokes``, ``periods_plane``, ``regulator_plane``,
``period_matrix_mp``, ``c_invariants_mp``, ``plane.cplus`` and ``plane.cminus``.
Plane methods select the validated Stokes algorithm when this flag is set.
The input polynomial must be exact, squarefree, rational and of degree at
least five. Logarithmic centres and finite poles must be rational.
Floating-point inputs denote their exact stored binary values; use ``QQ``
inputs to specify exact decimal fractions such as ``QQ(1)/10``.

For L-values the package bundles a pinned, MIT-licensed snapshot of the
Costa--Platt generic L-function sources, with the corrections described below.
The optional executable builds during installation when FLINT 3 and primesieve
headers and libraries are available in the Sage/Conda prefix. It can also be
built explicitly, without network access::

   sage -python -m hyperell_regulator.build_certified

Set ``HYPERELL_REGULATOR_REQUIRE_CERTIFIED=1`` to make its absence an
installation error. ``HYPERELL_CERTIFIED_LVALUE`` can select an explicitly
built executable. Certified coefficient precomputation defaults to eight
workers. ``HYPELLFROB_THREADS`` overrides that default, including larger
settings such as ``HYPELLFROB_THREADS=12``. Each L-value
runs in a fresh native process, with no upstream disk kernel cache.

Large Frobenius batches can require several GiB across the worker pool;
set ``HYPELLFROB_THREADS`` to reduce it. The native FFT multiplication streams
matrix rows to reduce temporary storage. Missing higher primes on a precision
retry use balanced modular products across the cached range. If the native
Frobenius accelerator fails or declines a bulk batch, the computation raises
an error instead of silently starting thousands of per-prime Sage fallback
computations. Native Frobenius and L-value subprocesses have 900-second
timeouts, configurable via ``HYPELLFROB_TIMEOUT`` and
``HYPERELL_CERTIFIED_TIMEOUT``. The certified experiment reports regulator,
period and L-value stages even without ``--verbose``. Native computations
report elapsed time every 30 seconds. ``--verbose`` additionally reports
triples and integration edges; ``--edge-seconds`` applies only to quadrature.

After building::

   from hyperell_regulator import HyperellipticCurveWithReduction, StandardLFunction
   C = HyperellipticCurveWithReduction(x**2+x, x**3+1)
   L = StandardLFunction(C)
   enclosure = L.value(digits=100, certified=True)

``CurveLFunction.value``, ``standard_lvalue`` and ``curve_lvalue`` also accept
the flag. Supported evaluations are real arguments exactly representable in
binary64, derivatives of order zero or one away from gamma poles, and the
standard second derivative at 1. The standard value and first derivative at 1
are exactly zero under the holomorphy assumptions. Other higher derivatives
and non-dyadic arguments raise ``NotImplementedError``. An explicit
``num_coeffs`` limits the coefficient cutoff; its sufficiency is verified
from the resulting error bounds. A returned L-ball certifies **numerical
errors conditional on correct Euler factors, conductor, gamma factors,
holomorphic continuation, the functional equation and the algorithm's growth
hypotheses**. It does not prove those analytic inputs. The wrapper dictionaries
include ``analytic_assumptions`` and backend metadata; ``feq_check=None``
avoids presenting a numerical residual as an analytic proof.

Without an explicit ``num_coeffs``, an insufficient enclosure triggers
adaptive refinement of the coefficient range and working precision, using
the measured radius after normalization. Up to six attempts are allowed.
Progress output and ``L_metadata.attempts`` record each cutoff, precision
and achieved radius; a failure retains these diagnostics in the experiment
cache. The final radius requirement is unchanged.


The Beilinson experiment
------------------------

Run validated regulators, periods, determinants and conditional L-values::

   sage -python -m hyperell_regulator.experiments.beilinson 250 \
       --certified --digits 100 --reg-digits 100 --edge-seconds 600 \
       --out results/beilinson_certified.json

``--reg-digits`` defaults to ``--digits`` in certified mode. The family
experiment also accepts ``--certified``. Edge time and work limits apply to
ball integration as well; validated rays can require more than the default
15 seconds. Successful numerical cache records cannot satisfy a certified
request. Certified L-values use a separate cache key, and saved enclosures
are checked against the requested radius before reuse.
The CLI also uses a separate default report filename ending in
``_certified.json``; an explicit ``--out`` takes precedence.

JSON records store each real ball as exact rational lower and upper endpoints
plus its field precision, using ``certified.encode_ball`` and
``certified.decode_ball``. Console decimals are for display. All determinant,
period normalization and quotient calculations use balls. A determinant ball
must exclude zero before division. Integer relation detection uses a midpoint
to propose a rational; ``relative_agreement_upper`` is an exact rational upper
bound for its relative error. ``exact_relation_proved`` remains false: proving
an error below ``10^-100`` proves an inequality, not rational equality.
The existing saved numerical survey is not retroactively certified.

Certified regulator construction
---------------------------------

Plane moments and logarithmic integrals use a certified embedded tree joining
all branch values and logarithmic centres, with a single horizontal ray from
the rightmost vertex to infinity. Floating-point distances propose a short
tree; certified geometry checks reject crossings and interior cut vertices.
A higher-precision proposal is used when float64 cannot distinguish vertices.
Minimality is unnecessary for the Stokes identity.

Cyclic orders at vertices are determined by exact algebraic half-plane and
cross-product comparisons. The resulting rotation system visits each directed
bank once. Logarithms are continued analytically along the straight segments
and the clockwise exterior sector at each finite vertex. The infinity
connector makes one counterclockwise turn. Ball arguments isolate the exact
integer turns; they are not estimated by sampled phases or rounded without
an enclosure. Every jump is also checked against the component separated by
that edge. The parity of the branch-root turns fixes the sheet gluing.
Local integration sheets are matched by isolating the exact ratio ``+1`` or
``-1`` between two square roots of the same value.

Finite edges use the validated endpoint-regularized kernel below. The single
ray is split at ``B=z+D``, with positive rational ``D``. A regularized finite
segment reaches ``B``; the reciprocal tail uses ``x=z+D/u^8``. This makes
the tail differential vanish at least to order three at ``u=0`` and isolates
its logarithmic singularity as ``-8*log(u)``. The normalized root factors
stay in the right half-plane on the integration interval because the initial
vertex is rightmost. Arb checks holomorphy on complex quadrature
neighbourhoods. Explicit endpoint bounds enclose the omitted pieces, including
logarithmic tails. The Stokes assembly and its final radius check use balls.
A closure check is an additional diagnostic; it does not replace the topology
certificate. An uncertifiable tree raises an exception.

The independent parallel-ray construction remains available explicitly via
``certified_cut="parallel"`` on the plane/Stokes APIs. It orders disjoint rays
using exact transverse coordinates and is useful as a cross-check; the tree
backend does not silently fall back to it. ``certified=True`` and the
experiment's ``--certified`` select the tree by default. Saved parallel-ray
regulators are recomputed when a tree report is requested; their L-values
remain reusable.

Cycle periods retain the existing tree basis. Numerical root approximations
propose edges, then interval geometry certifies their disjoint embedding and
exact algebraic checks exclude interior branch points. Minimality of the tree
is unnecessary. The substitution ``x=a+(b-a)*sin(pi*t/2)^2`` cancels endpoint
branch factors; the same analytic quadrature and endpoint estimates apply.
Signed period minors use certified nonzero ball pivots. Midpoints only choose
among possible pivots; every arithmetic operation and division remains bounded.

Historical integration-kernel checks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The original finite-edge prototype is now ``regulator.balls.edge_integrals``;
its experimental import remains available. Analytic checks of
``integral_0^1 dx/sqrt(x*(1-x))=pi`` and the corresponding logarithmic beta
integrals pass at 100-digit absolute accuracy. Their tests are in
``tests/test_ball_regulator_probe.py``. Public-API tests additionally cover
parallel cuts, several centres, finite poles, Weierstrass endpoints, real
branch points, signed period minors, resource limits and exact JSON storage.

L-values: avoid the second derivative
-------------------------------------

For this repository's standard L-function let

.. math::

   \Gamma_{\mathbb R}(s)&=\pi^{-s/2}\Gamma(s/2),\\
   G(s)&=\Gamma_{\mathbb R}(s-1)^2
          \Gamma_{\mathbb R}(s)^2\Gamma_{\mathbb R}(s+1),\\
   \Lambda(s)&=N^{s/2}G(s)L(s).

Assume that ``Lambda`` is holomorphic at 1 and 2 and satisfies
``Lambda(s) = epsilon*Lambda(3-s)``, with ``epsilon`` equal to ``+1`` or ``-1``.
Then ``L`` has a zero of order at least two at 1. Since

.. math::

   G(s)\sim\frac{4}{\pi(s-1)^2},\qquad G(2)=\frac{1}{2\pi^3},

we obtain

.. math::

   \Lambda(1)=\frac{2\sqrt N}{\pi}L''(1),\qquad
   \Lambda(2)=\frac{N}{2\pi^3}L(2),\qquad
   \boxed{L''(1)=\epsilon\frac{\sqrt N}{4\pi^2}L(2)}.

The repo uses ``epsilon=+1``. For ``249.a.249.1``, ``N=62001`` and the
35-digit PARI evaluations of both sides agreed to about ``5.9e-39`` relative
error. This validates the normalization numerically, not the analytic
assumptions. The production backend applies this identity entirely in Arb; the historical
helper ``standard_second_from_value_at_two`` propagates the same bound.

This removes the need for a general second-derivative routine, or uncertified
finite differences. It does not make the Dirichlet series at 2 absolutely
convergent: the weight-two degree-five coefficients still require analytic
continuation or a smoothed functional-equation method here.

Historical backend investigation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Costa--Platt's generic lfunctions library
<https://github.com/edgarcosta/lfunctions>`_ was downloaded into ``/tmp`` and
built against the FLINT and primesieve supplied by the local Sage environment.
The inspected revision was
``42014e1204f91e93f72b46945c4f5701b3ed7b1b``. It accepts our parameters:

* degree ``5``;
* algebraic-to-analytic normalization ``1``;
* algebraic gamma shifts ``[-1,-1,0,0,1]``; internally these become
  ``[0,0,1,1,2]``, satisfying its nonnegative half-integer restriction;
* exact prime Euler factors, including the bad primes, obtained from
  ``StandardLFunction``.

The C probe prints the coefficient cutoff before requiring input. At
conductor ``62001`` the cutoffs were ``53443`` for an 80-bit target and
``517454`` for a 360-bit target. These are coefficient-index bounds, so Euler
factors are needed at every prime up to the cutoff, including the bad primes.

The original all-in-one computation returned ``ERR_RH_ERROR`` with the
diagnostic ``Found too many zeros``. Investigation established that this is
an inconclusive zero-count bound, not a contradiction in the computed values.
The initial values-only 80-bit run reproduced the same narrow interval as
the full pipeline. That interval has since been **withdrawn**: a subsequent
audit found separate upstream errors in the Fourier error-bound propagation.
The historical record ``reports/ball_lvalue_249_80.json`` is retained with
status ``withdrawn_error_bounds`` and must not be used as a certificate.
With the corrections below, the bridge rejects that run with exit code 9
because its achieved accuracy is below 80 bits.

Fourier error-bound corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``lfunctions_error_bounds.patch`` fixes two issues in upstream ``src/error.c``:

* In a loop over the retained Fourier entries indexed by ``j``, the code
  added the kernel-approximation and Fourier-aliasing errors repeatedly to
  ``res[i]``, where ``i`` was the end of the retained range. It now adds them
  to each ``res[j]`` before the inverse FFT.
* The code computed a frequency-dependent coefficient cutoff but did not use
  it in ``M_error``. The corrected code obtains the included prefix using the
  same rounded bucket index as the actual convolution and passes that cutoff
  into the tail bound. It verifies that these indices are monotone. If the
  partial-sum bound's sufficient conditions fail, subsequent Fourier entries
  are bounded in full with the existing whole-tail bound instead.

The regression program ``ball_error_bounds_test.c`` checks the computed
special value of the product of the quadratic Dirichlet L-functions of
conductors 3 and 4 against an independent FLINT Dirichlet-L enclosure. It also
deliberately enlarges the kernel error and checks that every tested Fourier
entry receives it. The latter test fails on the original upstream error code
and passes on the correction. At an internal target of 400 bits, the corrected
Dirichlet example achieved 210 bits, demonstrating why achieved precision
must be checked separately from requested precision.

The C probe requires a revision symbol provided by this patch, so an
accidentally unpatched shared library cannot silently supply its results.

Completed 100-digit run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The corrected values-only computation for ``249.a.249.1`` succeeded with
42861 prime Euler factors through coefficient index 517454, an internal
600-bit target and 704-bit working arithmetic. Its final enclosure for
``L''(std,1)`` has 398 bits of relative accuracy and radius less than
``5.228e-120``. Both endpoints round to the same 100-decimal-place value.
These bounds are conditional on the stated analytic inputs, including the
functional equation and holomorphy assumptions; the computation does not
prove those inputs.

The result and provenance are in ``reports/ball_lvalue_249_100.json``; the
human-readable report is ``reports/ball_lvalue_249_100.md``. The complete
exact input is saved as ``reports/ball_lvalue_249_100_euler.txt.gz`` and the
backend output as ``reports/ball_lvalue_249_100_run.txt``.

The probe's optional fourth and fifth arguments select the required final
accuracy and the actual coefficient cutoff independently of its internal
target. A reduced cutoff is passed to the corrected tail-bound calculation;
it is not assumed sufficient. The final radius must pass the requested
acceptance check. This run used::

   gzip -dc reports/ball_lvalue_249_100_euler.txt.gz > /tmp/ball-249-euler.txt
   /tmp/ball-lvalue-probe 600 62001 /tmp/ball-249-euler.txt 333 517454
   sage -python -m pytest -q tests/test_saved_ball_lvalue_100.py

The saved-artifact tests use exact rational arithmetic to check the interval
radius and rounding and use Machin's identity to check the normalization
against an independent rigorous pi enclosure. They also verify every input
prime position and the recorded input and source hashes.

Why zero counting failed
^^^^^^^^^^^^^^^^^^^^^^^^

``Lfunc_compute`` computes the sampled completed L-function and root number,
then also determines the apparent rank, searches for zeros, and runs Turing's
method to certify that its zero list is complete. The last step is useful
for zero lists and finite-height RH statements, but is not an input to the
special-value interpolation routine.

Instrumenting ``src/turing.c`` gave, for this degree-five example:

* The output height is ``t0=12.8`` and the counting window is ``h=6.4``.
* The code found 40 zeros, counting both sides of the real axis.
* Its analytic count bound is approximately ``39.9643290043 +/- 3.1365``.
* To accept the count, the code requires this entire interval to be strictly
  inside ``(39,41)``. The lower-endpoint test fails, producing the misleading
  message ``Found too many zeros``. The actual interval includes 40 and
  several other integers.

The main uncertainty is the analytic bound for the integral of ``S(t)`` in
``St_int``, not floating-point rounding. Merely increasing the bit precision
does not make this fixed-window bound narrow enough. The upstream high-degree
regression suite already explicitly tolerates ``ERR_RH_ERROR`` for examples
affected by this limitation.

A diagnostic experiment changing ``TURING_RATIO`` from 16 to 4 increased
``h`` to 25.6 and reduced the count radius to about 0.827. However, at the
same working precision the zero samples became unusable before the end of
the longer window, yielding ``ERR_SOME_DATA`` and an incomplete zero list.
That experiment was not accepted as a successful zero certificate. The
original Turing settings were restored.

Values-only entry point
^^^^^^^^^^^^^^^^^^^^^^^

The adjacent ``lfunctions_values_only.patch`` adds ``Lfunc_compute_values``
to the optional upstream library. It returns immediately after the common
sample computation and root-number normalization, before rank and zero
processing. The original ``Lfunc_compute`` retains all its original behavior.
Neither entry point may be called twice, or after the other, on one object:
the coefficient storage is normalized destructively.

With the error-bound patch applied, the values-only route propagates the
coefficient-tail, convolution, FFT and special-value interpolation errors.
It does not catch and discard
``ERR_RH_ERROR``; the zero-counting step is never requested. All error and
warning flags from the stages that are requested remain rejection conditions
in our bridge. No zero-list completeness or RH claim is made.

The control build, selected by ``-DBALL_PROBE_COMPUTE_ALL``, invokes the
unmodified all-in-one entry point and rejects the zero-count warning with
exit code 6. The special-value computation does not depend on the rank search
or zero-counting verdict. The earlier equality of full and values-only
outputs established this separation, but did not validate the shared
Fourier-error propagation, which required the separate correction above.

The bridge uses ``Lfunc_special_value`` at algebraic ``s=2`` and propagates
the resulting ball through the formula above. It requires the complete
prime-factor range, rejects all backend error/warning flags, checks finiteness
and the achieved relative accuracy, and exports intervals as exact dyadic
endpoints: a row ``label lo hi exponent`` means
``[lo*2^exponent, hi*2^exponent]``. Decimal output is for display only.

Apply the patch to a checkout of the pinned upstream revision, and rebuild
the optional library before compiling the bridge::

   cd /path/to/lfunctions
   git apply /path/to/hyperell-regulator/hyperell_regulator/experiments/lfunctions_values_only.patch
   git apply /path/to/hyperell-regulator/hyperell_regulator/experiments/lfunctions_error_bounds.patch
   make lib

Then, from this repository, adapt the following paths to the local installation::

   cc -O2 hyperell_regulator/experiments/ball_lvalue_probe.c \
     -I/path/to/lfunctions/include -I/path/to/flint/include \
     -L/path/to/lfunctions/build -L/path/to/flint/lib \
     -Wl,-rpath,/path/to/lfunctions/build -Wl,-rpath,/path/to/flint/lib \
     -llfun -lflint -lprimesieve -o /tmp/ball-lvalue-probe

The probe sets no ``cache_dir``, so it recomputes the G kernel without
consulting the upstream on-disk cache. An isolated working directory keeps
the experiment's inputs and outputs together::

   mkdir /tmp/ball-249-80
   cd /tmp/ball-249-80
   /tmp/ball-lvalue-probe 80 62001
   /tmp/ball-lvalue-probe 80 62001 factors.txt

The input file has one row per prime in increasing order:
``p c0 c1 c2 c3 c4 c5``. These are exact integer coefficients of the local
Euler polynomial in ascending order, padded by zeros. For example, run this
in Sage from the repository root to generate the 80-bit input::

   from sage.all import PolynomialRing, QQ, prime_range
   from hyperell_regulator import HyperellipticCurveWithReduction, StandardLFunction

   x = PolynomialRing(QQ, 'x').gen()
   curve = HyperellipticCurveWithReduction(x**2+x, x**3+1)
   L = StandardLFunction(curve)
   primes = list(prime_range(53443 + 1))
   L.precompute_h1(primes)
   with open('/tmp/ball-249-80/factors.txt', 'w') as stream:
       for p in primes:
           factor = L.euler_factor(p)
           stream.write(' '.join(map(str, [p]+[factor[i] for i in range(6)]))+'\n')

The historical probe is separate from the bundled production backend described
above; it remains useful for reproducing the saved certificate.

The optional regression tests can be enabled with::

   HYPERELL_BALL_LVALUE_PROBE=/tmp/ball-lvalue-probe \
   HYPERELL_BALL_LVALUE_FACTORS=/tmp/ball-249-80/factors.txt \
   sage -python -m pytest -q tests/test_ball_lvalue_probe.py

To run the full-path control as well, compile a second bridge with
``-DBALL_PROBE_COMPUTE_ALL`` and set ``HYPERELL_BALL_LVALUE_FULL_CONTROL``
to that binary. The tests check rejection of insufficient achieved precision,
rejection of incomplete prime input, and the full-path zero-count warning.
Set ``HYPERELL_BALL_ERROR_TEST`` to the compiled ``ball_error_bounds_test.c``
program to run the independent Dirichlet-value, first-derivative and
error-propagation checks. First derivatives, including at the central
interpolation node, are checked against independent FLINT Dirichlet-L jets.

Relevant interfaces and references
-----------------------------------

* `Sage complex balls and rigorous integration
  <https://doc.sagemath.org/html/en/reference/rings_numerical/sage/rings/complex_arb.html>`_.
* `Sage certified complex root isolation
  <https://doc.sagemath.org/html/en/reference/polynomial_rings/sage/rings/polynomial/complex_roots.html>`_.
* `FLINT integration algorithm and analyticity contract
  <https://flintlib.org/doc/acb_calc/integrate.html>`_.
* `Generic L-function library API
  <https://github.com/edgarcosta/lfunctions/blob/42014e1204f91e93f72b46945c4f5701b3ed7b1b/doc/api.md>`_.
* `Computing Classical Modular Forms, section 9.2
  <https://research-information.bris.ac.uk/ws/portalfiles/portal/253332734/cmfs_mf_101120_scrubbed.pdf>`_,
  for previous rigorous use of the generic L-function library.
