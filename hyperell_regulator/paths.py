r"""
Where the cached data files live.

The modules of this package keep their results in JSON files: the LMFDB rows
they work from, the component answers, the `L`-values, the Beilinson
determinants.  As scripts these paths were relative to the current directory;
as a package they are relative to a *data directory*, which is

* ``$HYPERELL_REGULATOR_DATA`` if that is set, else
* the current working directory,

so the default behaviour is the one the scripts had.  Writers call
:func:`ensure_parent` first, so a missing ``lmfdb_cache/`` is created rather
than raising, which the scripts did not do.

EXAMPLES::

    sage: import os
    sage: from hyperell_regulator.paths import data_dir, data_path
    sage: os.environ["HYPERELL_REGULATOR_DATA"] = tmp_dir()
    sage: data_path("lmfdb_cache", "g2c.json").startswith(str(data_dir()))
    True
"""

import os
from pathlib import Path

ENV_VAR = "HYPERELL_REGULATOR_DATA"


def data_dir():
    r"""
    The directory the cached JSON files are read from and written to.
    """
    return Path(os.environ.get(ENV_VAR) or Path.cwd())


def data_path(*parts):
    r"""
    The path of a data file inside :func:`data_dir`, as a string.

    Nothing is created here; call :func:`ensure_parent` before writing.
    """
    return str(data_dir().joinpath(*parts))


def ensure_parent(path):
    r"""
    Create the directory ``path`` lives in, if it does not exist, and return
    ``path``.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    return path
