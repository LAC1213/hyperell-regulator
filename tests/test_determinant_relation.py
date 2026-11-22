"""Relations with an L-value term and potentially dependent determinants."""

import pytest
from sage.all import QQ, RealBallField, RealField

from hyperell_regulator.experiments._determinant_relation import determinant_relation


@pytest.mark.parametrize("certified", [False, True])
def test_relation_requires_nonzero_lvalue_coefficient(certified):
    field = (RealBallField if certified else RealField)(120)
    a, b = field(2).sqrt(), field(3).sqrt()
    pairs = [(a, None, 0, 1), (2*a, None, 0, 2), (b, None, 1, 2)]
    options = dict(digits=25, lvalue_digits=25, bound=8, certified=certified)
    unrelated = field(5).sqrt()/field.pi()
    assert determinant_relation(pairs, field(1), field(1), unrelated, **options) == {}

    lvalue = (3*a-2*b)/field.pi()
    relation = determinant_relation(pairs, field(1), field(1), lvalue, **options)
    coefficients = {term['basis_indices']: QQ(term['coefficient']) for term in relation['terms']}
    assert coefficients.get((0, 1), 0) + 2*coefficients.get((0, 2), 0) == 3
    assert coefficients[(1, 2)] == -2
    assert relation['relative_residual'] < QQ(10)**-23
    assert not relation['exact_relation_proved']
    if certified:
        assert relation['residual'].contains_zero()
