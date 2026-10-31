Beilinson experiments
=====================

The experiment commands combine regulator determinants, periods and standard
:math:`L`-values, and compare

.. math::

   q = \frac{\pi\,c^+c^-\,L''(\mathrm{std},1)}{D}

with a rational guess. The determinant :math:`D` uses weights
:math:`d_{11}=1`, :math:`d_{12}=2\operatorname{Re}(z)`, and
:math:`d_{22}=|z|^2`. The periods :math:`c^\pm` are full cycle minors;
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

The Beilinson survey includes triples with undecided component verdicts by
default, so ``388.a.776.1`` is attempted despite the undecided condition at 2.
The saved triple retains those verdicts. Use ``--exclude-undecided`` (or
``allow_undecided=False`` in Python) to require affirmative verdicts; known
negative verdicts are always excluded.

Saved results
-------------

Data and results are cached under ``lmfdb_cache/`` and ``results/`` in the
working directory. Set ``HYPERELL_REGULATOR_DATA`` to choose another location.

The `saved 100-digit survey
<https://github.com/LAC1213/hyperell-regulator/tree/main/reports>`_ includes
LaTeX tables, full-precision values, exact functions and regeneration
instructions.
