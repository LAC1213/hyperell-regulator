"""Supervise a full curve calculation, including its native descendants."""

import json
import os
from pathlib import Path
import pickle
import signal
import subprocess
import sys
import tempfile
import time


def run(row, settings, *, verbose, end_algs):
    """Run one curve with live progress and a wall-clock budget."""
    limit = settings['curve_seconds']
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='beilinson-') as work:
        request = Path(work)/'request.json'
        result = Path(work)/'result.pickle'
        request.write_text(json.dumps(dict(row=row, settings=settings, verbose=verbose,
                                           end_algs=end_algs), default=str))
        # The worker and its unmodified Sage/GP/native subprocesses share a
        # process group so a timeout cannot leave a descendant running.
        grouped = os.name == 'posix'
        try:
            process = subprocess.Popen(
                [sys.executable, '-u', '-m', __spec__.name, str(request), str(result)],
                start_new_session=grouped)
        except OSError as exc:
            return dict(status='preparation_error', error=str(exc), seconds=time.monotonic()-started)
        try:
            while True:
                elapsed = time.monotonic()-started
                if elapsed >= limit:
                    return dict(status='curve_timeout', seconds=elapsed,
                                error='Curve calculation exceeded %g seconds' % limit)
                try:
                    code = process.wait(timeout=min(10, limit-elapsed))
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic()-started
                    if elapsed < limit:
                        print('  %s curve calculation: %.0fs elapsed' % (row['label'], elapsed),
                              flush=True)
                    continue
                if code:
                    return dict(status='worker_failure', seconds=time.monotonic()-started,
                                error='Curve worker exited with code %d' % code)
                try:
                    # The private temporary directory contains only this worker's result.
                    with result.open('rb') as stream:
                        return pickle.load(stream)
                except (OSError, EOFError, pickle.UnpicklingError) as exc:
                    return dict(status='worker_failure', seconds=time.monotonic()-started,
                                error='Could not read curve result: %s' % exc)
        finally:
            try:
                if grouped:
                    os.killpg(process.pid, signal.SIGKILL)
                elif process.poll() is None:
                    process.kill()
            except ProcessLookupError:
                pass
            process.wait()


def _main():
    from .beilinson import _analyse_monitored

    request, result = map(Path, sys.argv[1:])
    args = json.loads(request.read_text())
    args['curve'] = args['row']
    with result.open('wb') as stream:
        pickle.dump(_analyse_monitored(**args), stream)


if __name__ == '__main__':
    _main()
