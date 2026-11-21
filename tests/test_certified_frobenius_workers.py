"""Large certified coefficient ranges must use the configured worker pool."""
from pathlib import Path

import pytest

from hyperell_regulator import ball_lfunction


@pytest.mark.parametrize('configured,expected',
                         [(None, 8), ('2', 2), ('8', 8), ('12', 12), ('32', 32), ('0', 1)])
def test_large_batch_worker_pool(monkeypatch, configured, expected):
    if configured is None:
        monkeypatch.delenv('HYPELLFROB_THREADS', raising=False)
    else:
        monkeypatch.setenv('HYPELLFROB_THREADS', configured)
    monkeypatch.setattr(ball_lfunction, 'executable', lambda: str(Path(__file__)))
    monkeypatch.setattr(ball_lfunction, '_run', lambda args: {'coefficient_cutoff': ['130000']})

    class PrecomputationReached(Exception):
        pass

    class LFunction:
        _cache = {}

        def conductor(self):
            return 1

        def precompute_h1(self, primes, nthreads):
            assert len(primes) > 10000
            assert nthreads == expected
            raise PrecomputationReached

    with pytest.raises(PrecomputationReached):
        ball_lfunction.value(LFunction(), s=2, digits=100)
