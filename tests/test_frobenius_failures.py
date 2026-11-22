"""A failed bulk accelerator must not silently start hours of Sage fallback."""
import subprocess
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from sage.all import PolynomialRing,ZZ,prime_range

from hyperell_regulator import frobenius
from hyperell_regulator.ball_lfunction import _run
from hyperell_regulator.regulator.balls import CertificationError


def test_bulk_native_crash_is_reported():
    with patch.object(frobenius,'euler_binary',return_value='/euler'), \
         patch('hyperell_regulator._subprocess.run',side_effect=subprocess.CalledProcessError(-11,['/euler'],stderr='failed')):
        with pytest.raises(frobenius.FrobeniusBackendError,match='refusing bulk Sage fallback'):
            frobenius._run([1,1],[(101,2)]*1001,1)
        assert frobenius._run([1,1],[(101,2)],1)=={}


def test_missing_bulk_accelerator_is_reported():
    with patch.object(frobenius,'euler_binary',return_value=None):
        with pytest.raises(frobenius.FrobeniusBackendError,match='unavailable'):
            frobenius._run([1,1],[(101,2)]*1001,1)


def test_mass_declines_do_not_enter_sage_frobenius():
    x=PolynomialRing(ZZ,'x').gen()
    with patch.object(frobenius,'_run',return_value={}), \
         patch('sage.all.HyperellipticCurve',side_effect=AssertionError('unexpected bulk fallback')):
        with pytest.raises(frobenius.FrobeniusBackendError,match='declined'):
            frobenius.frobenius_polynomials(x**5+x+1,list(prime_range(10000))[:1001])


def test_native_lvalue_timeout_is_a_certification_failure(monkeypatch):
    monkeypatch.setenv('HYPERELL_CERTIFIED_TIMEOUT','12')
    with patch('hyperell_regulator.ball_lfunction._native_run',
               side_effect=subprocess.TimeoutExpired(['/backend'],12)):
        with pytest.raises(CertificationError,match='exceeded 12 seconds'):
            _run(['/backend'])
