r"""
The numerical experiments driving the package.

* :mod:`~hyperell_regulator.experiments.beilinson` -- for each curve in the
  LMFDB list, find a triple of rational points whose differences are torsion
  and which satisfy the component condition, compute the regulator
  determinant, and compare it with `L''(\mathrm{std}, 1)`.

* :mod:`~hyperell_regulator.experiments.family_n` -- the same test for the
  one-parameter family, whose members satisfy the conditions by construction.

Each is runnable, e.g. ``sage -python -m
hyperell_regulator.experiments.beilinson``, and each writes its results to
JSON under :func:`hyperell_regulator.paths.data_dir`.
"""
