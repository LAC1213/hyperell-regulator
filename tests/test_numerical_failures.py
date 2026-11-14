"""The direct API refuses singular curves and invalid precision."""
import pytest
from hyperell_regulator.regulator import regulator_stokes
from hyperell_regulator.regulator.polynomials import ring


@pytest.mark.parametrize('prec',[0,4,12.5])
def test_invalid_accuracy_is_rejected(prec):
    x=ring().gen()
    with pytest.raises(ValueError,match='prec'):
        regulator_stokes(x**5+x+1,prec=prec)


def test_singular_curve_is_rejected():
    with pytest.raises(ValueError,match='squarefree'):
        regulator_stokes(ring().gen()**5)
