"""Native subprocesses with elapsed-time reports while their pipes are drained."""
import subprocess
import tempfile
import time
from contextlib import ExitStack


def run(args, *, input=None, timeout, progress=None, label='native computation',
        check=False, interval=30):
    if progress is None:
        return subprocess.run(args, input=input, text=True, capture_output=True,
                              timeout=timeout, check=check)
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
        process = stack.enter_context(subprocess.Popen(
            args, stdin=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
        try:
            while True:
                remaining = timeout-(time.monotonic()-started)
                if remaining <= 0:
                    process.kill()
                    out, err = process.communicate()
                    raise subprocess.TimeoutExpired(args, timeout, output=out, stderr=err)
                try:
                    out, err = process.communicate(timeout=min(interval, remaining))
                    break
                except subprocess.TimeoutExpired:
                    progress('%s: %.0fs elapsed' % (label, time.monotonic()-started))
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()
        result = subprocess.CompletedProcess(args, process.returncode, out, err)
        if check:
            result.check_returncode()
        return result
