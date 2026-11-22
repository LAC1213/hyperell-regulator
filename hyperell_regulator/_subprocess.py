"""Native subprocesses with elapsed-time reports and bounded waits."""
import subprocess
import tempfile
import time
from contextlib import ExitStack


def run(args, *, input=None, timeout, progress=None, label='native computation',
        check=False, interval=30):
    started = time.monotonic()
    with ExitStack() as stack:
        data = None
        if input is not None:
            # A file gives the child a stable input stream and EOF. Retrying
            # communicate() after a timeout cannot reliably resume a large
            # partially written stdin pipe on all supported Python versions.
            data = stack.enter_context(tempfile.TemporaryFile(mode='w+t'))
            data.write(input)
            data.seek(0)
        # A descendant can inherit a pipe and keep communicate blocked after
        # the timed-out child exits. Files allow us to wait for the child alone.
        stdout = stack.enter_context(tempfile.TemporaryFile(mode='w+t'))
        stderr = stack.enter_context(tempfile.TemporaryFile(mode='w+t'))
        process = stack.enter_context(subprocess.Popen(args, stdin=data, stdout=stdout, stderr=stderr))
        expired = False
        try:
            while True:
                remaining = timeout-(time.monotonic()-started)
                if remaining <= 0:
                    expired = True
                    break
                try:
                    process.wait(timeout=remaining if progress is None else min(interval, remaining))
                    break
                except subprocess.TimeoutExpired:
                    if progress is not None and time.monotonic()-started < timeout:
                        progress('%s: %.0fs elapsed' % (label, time.monotonic()-started))
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
        stdout.seek(0)
        stderr.seek(0)
        out, err = stdout.read(), stderr.read()
        if expired:
            raise subprocess.TimeoutExpired(args, timeout, output=out, stderr=err)
        result = subprocess.CompletedProcess(args, process.returncode, out, err)
        if check:
            result.check_returncode()
        return result
