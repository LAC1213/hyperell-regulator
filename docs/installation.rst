Installation
============

Requirements
------------

Install into an existing SageMath environment. SageMath itself is not
installed by pip. The package also uses NumPy, SciPy and a patched version
of ``sage_cluster_pictures``; these dependencies are installed automatically.

.. code-block:: bash

   sage -pip install git+https://github.com/LAC1213/hyperell-regulator.git

From a checkout:

.. code-block:: bash

   sage -pip install .

For development and tests:

.. code-block:: bash

   sage -pip install -e '.[test]'
   sage -python -m pytest
   sage -t hyperell_regulator

Optional Frobenius accelerator
------------------------------

Installation downloads and verifies NTL 11.6.0, then builds a private static
copy for the bundled Frobenius accelerator. This requires network access,
C11 and C++17 compilers, Perl and Make. Sage's own NTL installation is untouched;
GMP comes from Sage. If the build is unavailable, the package uses Sage's
Frobenius implementation.

Build parallelism defaults to the logical CPU count. Set
``HYPERELL_REGULATOR_BUILD_JOBS`` to override it, and
``GMP_INC`` and ``GMP_LIB`` if GMP is outside Sage's prefix. To require the
accelerator and check the installed binary:

.. code-block:: bash

   HYPERELL_REGULATOR_REQUIRE_EULER=1 sage -pip install .
   sage -python -m hyperell_regulator.frobenius --check

See the `native build instructions
<https://github.com/LAC1213/hyperell-regulator/tree/main/hyperell-frobenius>`_
for more detail.

Euler-factor computation also defaults to one worker per logical CPU. Use
``frobenius_polynomials(f, primes, nthreads=4)`` to choose a worker count.
The native build checks odd- and even-degree matrices against fixed verified
references with serial, four-worker and CPU-count runs before installation.

The native accelerator uses the bundled rforest C library for the Python
API, both L-function classes and the native driver. Rforest is the only native
backend; no separate rforest installation or backend setting is needed.
Sources live in ``hyperell-frobenius`` with
the bundled library under ``hyperell-frobenius/rforest``. The simultaneous algorithm benefits dense batches of
primes, including even degree models without a rational branch point.
Python monicizes nonmonic even input using a common quadratic twist and
removes its local character from the returned characteristic polynomial.
Set ``HYPELLFROB_THREADS`` to choose workers for L-function computations;
an explicit ``nthreads`` argument takes precedence. Unsupported jobs retain
Sage fallback. The vendored
FFT arithmetic requires a 64-bit platform with 64-bit GMP limbs.
