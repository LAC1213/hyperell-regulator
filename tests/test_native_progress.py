"""Progress reporting must preserve native pipe data, errors and deadlines."""
import subprocess
import sys

import pytest

from hyperell_regulator._subprocess import run


def test_progress_drains_pipes_and_sends_input_once():
    messages = []
    payload = 'abc'*100000
    command = [sys.executable, '-c',
               'import sys,time; data=sys.stdin.read(); '
               'sys.stderr.write("e"*200000); sys.stderr.flush(); '
               'sys.stdout.write(data); sys.stdout.flush(); time.sleep(.15)']
    result = run(command, input=payload, timeout=5, interval=.02,
                 progress=messages.append, label='Euler factors', check=True)
    assert result.stdout == payload
    assert result.stderr == 'e'*200000
    assert messages and all(m.startswith('Euler factors: ') for m in messages)


def test_progress_keeps_deadline_and_captured_errors():
    command = [sys.executable, '-c',
               'import sys,time; print("started",flush=True); time.sleep(10)']
    with pytest.raises(subprocess.TimeoutExpired) as error:
        run(command, timeout=.2, interval=.03, progress=lambda message: None)
    assert 'started' in error.value.output
    assert error.value.timeout == .2


def test_progress_preserves_failed_exit_status():
    with pytest.raises(subprocess.CalledProcessError) as error:
        run([sys.executable, '-c', 'import sys; sys.stderr.write("failed"); sys.exit(7)'],
            timeout=5, progress=lambda message: None, check=True)
    assert error.value.returncode == 7
    assert error.value.stderr == 'failed'
