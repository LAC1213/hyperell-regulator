r"""
Where the cached data files live.

LMFDB rows, component answers, `L`-values and Beilinson determinants are
cached in JSON files relative to

* ``$HYPERELL_REGULATOR_DATA`` if that is set, else
* the current working directory,

Writers call :func:`ensure_parent` to create missing directories.

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
