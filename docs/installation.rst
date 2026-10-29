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
a C++17 compiler, Perl and Make. Sage's own NTL installation is untouched;
GMP comes from Sage. If the build is unavailable, the package uses Sage's
Frobenius implementation.

Set ``HYPERELL_REGULATOR_BUILD_JOBS`` to control build parallelism, and
``GMP_INC`` and ``GMP_LIB`` if GMP is outside Sage's prefix. To require the
accelerator and check the installed binary:

.. code-block:: bash

   HYPERELL_REGULATOR_REQUIRE_EULER=1 sage -pip install .
   sage -python -m hyperell_regulator.frobenius --check

See the `native build instructions
<https://github.com/LAC1213/hyperell-regulator/tree/main/hypellfrob-threaded>`_
for more detail.
