r"""
Build the bundled ``hypellfrob`` driver as part of installing the package.

``sage -pip install git+<url>`` should leave nothing to do by hand, so the
C++ in ``hypellfrob-threaded`` is compiled here and the resulting ``euler``
binary is put inside the package, where
:func:`hyperell_regulator.frobenius.euler_binary` looks for it first.

NTL and GMP come with Sage, so their headers and libraries are taken from the
prefix of the interpreter doing the install -- which is Sage's own when this
is run as ``sage -pip``.  ``NTL_INC``, ``NTL_LIB`` and ``NTL_LINK`` in the
environment override that.

The build is *not* required: without a compiler, or without NTL, the install
still succeeds and :mod:`hyperell_regulator.frobenius` falls back to Sage's
own Frobenius matrices, which are correct and slower.  Only the speed is
lost, so a failure here warns rather than raises.
"""

import os
import shutil
import subprocess
import sys

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hypellfrob-threaded")
INSIDE = os.path.join("hyperell_regulator", "bin")


def _prefixes():
    """Where to look for NTL and GMP, best guess first."""
    out = []
    for var in ("SAGE_LOCAL", "CONDA_PREFIX"):
        if os.environ.get(var):
            out.append(os.environ[var])
    if os.environ.get("SAGE_ROOT"):
        out.append(os.path.join(os.environ["SAGE_ROOT"], "local"))
    for prefix in (sys.prefix, sys.base_prefix):
        out.append(prefix)
        # Sage's Python can live in local/var/lib/sage/venv-python*.
        parent = prefix
        while parent != os.path.dirname(parent):
            parent = os.path.dirname(parent)
            if os.path.basename(parent) == "local":
                out.append(parent)
    out.extend(("/opt/homebrew", "/usr/local", "/usr"))
    return list(dict.fromkeys(out))



class build_with_euler(build_py):
    def run(self):
        build_py.run(self)
        target = os.path.join(self.build_lib, INSIDE)
        try:
            binary = self._make()
        except Exception as exc:                       # noqa: BLE001
            if os.environ.get("HYPERELL_REGULATOR_REQUIRE_EULER") == "1":
                raise
            print("hyperell-regulator: could not build the hypellfrob driver "
                  "(%s); falling back to Sage's Frobenius matrices" % exc)
            return
        self.mkpath(target)
        shutil.copy(binary, os.path.join(target, "euler"))
        os.chmod(os.path.join(target, "euler"), 0o755)
        print("hyperell-regulator: built %s" % os.path.join(target, "euler"))

    def _make(self):
        inc = os.environ.get("NTL_INC")
        lib = os.environ.get("NTL_LIB")
        if not inc or not lib:
            for p in _prefixes():
                if os.path.exists(os.path.join(p, "include", "NTL", "ZZ.h")):
                    inc = inc or os.path.join(p, "include")
                    lib = lib or os.path.join(p, "lib")
                    break
            else:
                raise RuntimeError("NTL headers not found in %s"
                                   % ", ".join(_prefixes()))
        # Ignore checkout-specific config.mk and never reuse its object files.
        build = os.path.abspath(os.path.join(self.build_lib, "..", "euler-build"))
        args = ["make", "-C", SRC, "CONFIG_MK=", "BUILD=" + build,
                "NTL_INC=" + inc, "NTL_LIB=" + lib,
                "NTL_LINK=" + os.environ.get("NTL_LINK", "-lntl"),
                "GMP_INC=" + os.environ.get("GMP_INC", inc),
                "GMP_LIB=" + os.environ.get("GMP_LIB", lib)]
        subprocess.run(args, check=True)
        binary = os.path.join(build, "euler")
        result = subprocess.run([binary, "1"], input="1 0 0 1 0 1\n101 3\n",
                                text=True, capture_output=True, check=True, timeout=30)
        if result.stdout.split()[:2] != ["101", "1"]:
            raise RuntimeError("Frobenius driver failed its smoke test: " + result.stdout)
        return binary



class binary_distribution(Distribution):
    """The wheel holds a compiled program, so it is not ``py3-none-any``."""

    def has_ext_modules(self):
        return True


setup(cmdclass={"build_py": build_with_euler},
      distclass=binary_distribution,
      package_data={"hyperell_regulator": ["bin/euler"]})
