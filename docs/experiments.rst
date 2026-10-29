Beilinson experiments
=====================

The experiment commands combine regulator determinants, periods and standard
:math:`L`-values, and compare

.. math::

   q = \frac{\det}{2\pi\,c^+c^-\,L''(\mathrm{std},1)}

with a rational guess.

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

Saved results
-------------

Data and results are cached under ``lmfdb_cache/`` and ``results/`` in the
working directory. Set ``HYPERELL_REGULATOR_DATA`` to choose another location.

The `saved 100-digit survey
<https://github.com/LAC1213/hyperell-regulator/tree/main/reports>`_ includes
LaTeX tables, full-precision values, exact functions and regeneration
instructions.
