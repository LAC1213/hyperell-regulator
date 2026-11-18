Beilinson experiments
=====================

The experiments compare the regulator determinant :math:`D` with

.. math::

   \pi c^+c^-L''(M_{\mathrm{std}}(C),1)=qD,\qquad q\in\mathbb Q^\times.

The determinant uses weights :math:`(1,2\operatorname{Re}(z),|z|^2)`.
``comparison.rational_is_integer`` tests whether the guessed rational
factor has denominator one.

Survey and family
-----------------

Fetch curves with conductor strictly below 1000 and test them:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 1000 --fetch

Use cached data with 30-digit regulators and :math:`L`-values:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 250 \
       --prec 30 --digits 30

For the family
:math:`y^2+x(x-1)y=x(x-1)(x^3+nx^2-(n+1)x+1)`, test :math:`n=2`:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.family_n 2 --no-lvalue

Regulators default to float64, :math:`L`-values to 15 decimal digits, and
edge limits to 15 seconds. ``--prec`` sets regulator precision and
``--digits`` sets :math:`L`-value precision. ``--no-lvalue`` skips
:math:`L`-values. For the family, ``--max-conductor`` bounds the standard
conductor inclusively. Use ``--help`` for all options.

Integrality criterion
---------------------

Marked points are distinct and pairwise nonconjugate, with torsion
differences from a common pole. The component criterion requires each
marked pair to reduce smoothly on the same stable component at every prime.
Undecided verdicts are accepted and recorded by default; use
``--exclude-undecided`` to require affirmative verdicts.

To evaluate all torsion triples up to permutation without component tests:

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

Saved results
-------------

Caches are stored in ``lmfdb_cache/`` and ``results/`` relative to
``HYPERELL_REGULATOR_DATA`` or the working directory. Successful results
are reused when their precision and component policy suffice; failed
computations are retried. ``--force`` recomputes successful regulators.

The `saved 100-digit survey
<https://github.com/LAC1213/hyperell-regulator/tree/main/reports>`_
contains tables, full-precision values and integration models.
