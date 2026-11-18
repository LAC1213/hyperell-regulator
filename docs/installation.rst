Installation
============

Install in an existing SageMath environment:

.. code-block:: bash

   sage -pip install git+https://github.com/LAC1213/hyperell-regulator.git

NumPy, SciPy and the patched ``sage_cluster_pictures`` are installed as
dependencies. From a checkout, use ``sage -pip install .``.
For development and tests:

.. code-block:: bash

   sage -pip install -e '.[test]'
   sage -python -m pytest
   sage -t hyperell_regulator

Frobenius accelerator
---------------------

Installation builds the bundled rforest executable with a private,
thread-safe NTL 11.6.0 and Sage's GMP. The build requires network access,
C11 and C++17 compilers, Perl, Make, and a 64-bit platform with 64-bit GMP
limbs. If the accelerator is unavailable, computations use Sage.

Set ``HYPERELL_REGULATOR_BUILD_JOBS`` to control build parallelism and
``GMP_INC`` and ``GMP_LIB`` for GMP outside Sage's prefix.
To require the accelerator and check it:

.. code-block:: bash

   HYPERELL_REGULATOR_REQUIRE_EULER=1 sage -pip install .
   sage -python -m hyperell_regulator.frobenius --check

The check compares odd and even degree matrices against reference values
with one, four and CPU-count workers.

Frobenius computations use ``HYPELLFROB_THREADS``, defaulting to the logical
CPU count. An explicit ``nthreads`` argument takes precedence:

.. code-block:: python

   frobenius_polynomials(f, primes, nthreads=4)

See the `native build instructions
<https://github.com/LAC1213/hyperell-regulator/tree/main/hyperell-frobenius>`_
for manual builds and the driver format.
