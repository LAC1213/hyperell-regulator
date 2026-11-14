Beilinson experiments
=====================

The experiment commands combine regulator determinants, periods and standard
:math:`L`-values, and compare

.. math::

   q = \frac{\pi\,c^+c^-\,L''(\mathrm{std},1)}{D}

with a rational guess. The determinant :math:`D` uses weights
:math:`d_{11}=1`, :math:`d_{12}=2\operatorname{Re}(z)`, and
:math:`d_{22}=|z|^2`.
``comparison.rational_is_integer`` reports whether the guessed rational
has denominator one.

Survey curves
-------------

Fetch curve data and test conductors below 1000 using float64:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 1000 --fetch

Reuse cached data and request 30-digit regulators and L-values:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.beilinson 250 \
       --prec 30 --digits 30 --edge-seconds 120

One-parameter family
--------------------

Test parameter ``n=2`` without computing an L-value:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.family_n 2 --no-lvalue

Both commands support ``--prec``, ``--digits``, ``--no-lvalue``, ``--verbose``,
``--edge-seconds``, ``--force`` and ``--out``. Use ``--help`` for all options.
The default edge time limit is 15 seconds. For ``family_n``,
``--max-conductor`` bounds the standard motive's conductor.

PARI's stack allowance grows with the number of L-function coefficients and
the requested precision. Memory is committed as needed; if the workspace
exceeds the initial estimate, the failed PARI operation is retried with a
larger allowance, reusing the computed coefficients. Large coefficient
counts can require several GiB of memory at high precision.

The Beilinson survey includes triples with undecided component verdicts by
default, so ``388.a.776.1`` is attempted despite the undecided condition at 2.
The saved triple retains those verdicts. Use ``--exclude-undecided`` (or
``allow_undecided=False`` in Python) to require affirmative verdicts; known
negative verdicts are always excluded.

Numerical helpers
-----------------

The determinant helper takes the curve and two point coordinates:

.. code-block:: python

   from hyperell_regulator.experiments.beilinson import (
       beilinson_determinant_stokes, regulator_column_stokes,
   )

   det, matrix, errors = beilinson_determinant_stokes(F, P, Q, digits=30)
   column = regulator_column_stokes(F, P, digits=30)

Torsion orders are no longer arguments to these numerical helpers. Orders
remain part of exact cycle certification and report metadata. Pass numerical
options such as ``digits`` and ``pole`` by keyword. ``triple_conditioning``
and ``rank_triples`` no longer accept the unused ``mp`` flag; ranking uses
the shared cut geometry.

Saved results
-------------

Data and results are cached under ``lmfdb_cache/`` and ``results/`` in the
working directory. Set ``HYPERELL_REGULATOR_DATA`` to choose another location.
A saved regulator whose L-value failed is retried on the next run requesting
L-values.

The `saved 100-digit survey
<https://github.com/LAC1213/hyperell-regulator/tree/main/reports>`_ includes
LaTeX tables, full-precision values, integration models and regeneration
instructions.
