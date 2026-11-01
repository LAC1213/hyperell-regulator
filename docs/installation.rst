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
<https://github.com/LAC1213/hyperell-regulator/tree/main/hypellfrob-threaded>`_
for more detail.

Euler-factor computation also defaults to one worker per logical CPU. Use
``frobenius_polynomials(f, primes, nthreads=4)`` to choose a worker count.
The native build checks serial, four-worker and CPU-count results against
each other before installing the binary.

The experimental rforest branch also compiles the vendored rforest C
library into this binary. No separate rforest installation is needed.
Select it with ``HYPELLFROB_BACKEND=rforest`` for L-function computations,
or ``frobenius_polynomials(f, primes, backend="rforest")``. The default
remains ``legacy``. The simultaneous algorithm benefits dense batches of
primes, including even degree models without a rational branch point.
Python monicizes nonmonic even input using a common quadratic twist and
removes its local character from the returned characteristic polynomial.
Set ``HYPELLFROB_THREADS`` to choose workers for L-function computations;
an explicit ``nthreads`` argument takes precedence. Unsupported jobs retain
native and Sage fallback. The vendored
FFT arithmetic requires a 64-bit platform with 64-bit GMP limbs.
