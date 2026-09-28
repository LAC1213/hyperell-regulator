"""Packaging must use Sage's library prefix and the same Make output path."""
import os
import runpy
import tempfile
from pathlib import Path
from unittest.mock import patch
from setuptools import Distribution


def build_module():
    with patch('setuptools.setup'):
        return runpy.run_path(str(Path(__file__).resolve().parents[1] / 'setup.py'))


def test_sage_macos_venv_prefix():
    module = build_module()
    with patch.dict(os.environ, {}, clear=True), patch.object(
            module['sys'], 'prefix',
            '/private/var/tmp/sage/local/var/lib/sage/venv-python3.13'):
        assert '/private/var/tmp/sage/local' in module['_prefixes']()
        assert '/opt/homebrew' in module['_prefixes']()



def test_make_paths_and_smoke_test():
    module = build_module()
    cmd = module['build_with_euler'](Distribution())
    with tempfile.TemporaryDirectory() as directory, patch.dict(
            os.environ, {'NTL_INC': '/sage/include', 'NTL_LIB': '/sage/lib'},
            clear=True), patch('subprocess.run') as run:
        cmd.build_lib = str(Path(directory) / 'lib')
        run.return_value.stdout = '101 1 4 matrix entries'
        binary = cmd._make()
    args = run.call_args_list[0].args[0]
    assert 'CONFIG_MK=' in args
    assert 'BUILD=' + str(Path(binary).parent) in args
    assert 'NTL_LINK=-lntl' in args
    assert 'NTL_LIB=/sage/lib' in args
    assert run.call_args_list[1].args[0] == [binary, '1']
    assert run.call_args_list[1].kwargs['timeout'] == 30
