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

The exact beta-integral identities in ``tests/test_certified_edges.py``
check the certified finite-edge kernel, including both logarithmic
endpoint tails at 100-digit absolute accuracy.

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
assumptions. The production backend applies this identity entirely in Arb.

This removes the need for a general second-derivative routine, or uncertified
finite differences. It does not make the Dirichlet series at 2 absolutely
convergent: the weight-two degree-five coefficients still require analytic
continuation or a smoothed functional-equation method here.

Archived 100-digit L-value certificate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The certificate for 249.a.249.1 used 42861 exact prime Euler factors
through coefficient index 517454. The corrected values-only backend
achieved 398 bits of relative accuracy with radius below
:math:`5.228\times10^{-120}`. Both endpoints round to the same
100-decimal-place value under the stated analytic assumptions.

The enclosure, source checksums, and provenance are in
``reports/ball_lvalue_249_100.json``. The exact input and backend
output are in ``reports/ball_lvalue_249_100_euler.txt.gz`` and
``reports/ball_lvalue_249_100_run.txt``. The archived source files
include the two upstream patches used for this certificate. The
maintained backend vendors the corrected code.

The tests in ``tests/test_saved_ball_lvalue_100.py`` verify the
recorded hashes, the rounding, and the normalization using
independent exact rational bounds for :math:`\pi`.

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
