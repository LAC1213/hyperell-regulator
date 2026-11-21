"""Build the optional pinned L-value backend, without network access.

Run with Sage's Python: ``python -m hyperell_regulator.build_certified``.
Requires FLINT >= 3 and primesieve development headers/libraries.
This module also works without importing Sage during package installation.
"""
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


def build(output=None, *, driver=None):
    root = Path(__file__).resolve().parent
    output = Path(output or root / 'bin' / 'certified-lvalue').resolve()
    prefixes = [os.environ.get(k) for k in ('SAGE_LOCAL', 'CONDA_PREFIX')]
    prefixes += [sys.prefix, sys.base_prefix, '/usr/local', '/usr']
    if os.environ.get('SAGE_ROOT'):
        prefixes.insert(0, str(Path(os.environ['SAGE_ROOT'])/'local'))
    for p in (sys.prefix, sys.base_prefix):
        prefixes += [str(q) for q in Path(p).parents if q.name == 'local']
    prefix = next((Path(p) for p in prefixes if p and
                   (Path(p)/'include/flint/acb.h').exists() and
                   (Path(p)/'include/primesieve.h').exists()), None)
    if prefix is None:
        raise RuntimeError('FLINT 3 and primesieve headers not found in Sage/Conda prefix')
    vendor = root / '_vendor' / 'lfunctions'
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='certified-build-', dir=output.parent) as tmp:
        binary = Path(tmp)/'certified-lvalue'
        cmd = shlex.split(os.environ.get('CC', 'cc')) + ['-O2', '-std=gnu11']
        cmd += ['-I'+str(prefix/'include'), '-I'+str(vendor/'include')]
        cmd += [str(p) for p in sorted((vendor/'src').glob('*.c'))]
        cmd += [str(driver or root/'_native/lvalue.c'), '-L'+str(prefix/'lib'),
                '-Wl,-rpath,'+str(prefix/'lib'), '-lflint', '-lprimesieve',
                '-lmpfr', '-lgmp', '-lm', '-o', str(binary)]
        subprocess.run(cmd, check=True)
        os.replace(binary, output)
    return str(output)


if __name__ == '__main__':
    print(build(sys.argv[1] if len(sys.argv) > 1 else None))
