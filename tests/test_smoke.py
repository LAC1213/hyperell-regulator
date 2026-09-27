"""Minimal smoke tests that don't require network access or long computation."""

from sage.all import ZZ, PolynomialRing

from hyperell_regulator import HyperellipticCurveWithReduction


def test_import():
    import hyperell_regulator  # noqa: F401
    import hyperell_regulator.components  # noqa: F401
    import hyperell_regulator.lfunction  # noqa: F401
    import hyperell_regulator.paths  # noqa: F401
    import hyperell_regulator.regulator  # noqa: F401
    import hyperell_regulator.regulator.plane  # noqa: F401
    import hyperell_regulator.regulator.stokes  # noqa: F401
    import hyperell_regulator.experiments.beilinson  # noqa: F401
    import hyperell_regulator.experiments.family_n  # noqa: F401


def test_regulator_exports():
    from hyperell_regulator import regulator
    for name in regulator.__all__:
        assert callable(getattr(regulator, name)), name


def test_bad_primes_and_component_question():
    R = PolynomialRing(ZZ, 'x')
    x = R.gen()
    C = HyperellipticCurveWithReduction(x**5 + x**3 + 1)
    assert C.bad_primes() == [2, 53, 61]
    assert C.is_semistable_at(53) is True
    assert C.same_component((0, 1), (0, -1), 53)["answer"] is True
