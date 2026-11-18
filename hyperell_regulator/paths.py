r"""Cache paths relative to ``HYPERELL_REGULATOR_DATA`` or the working directory.

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
    r"""The cache directory."""
    return Path(os.environ.get(ENV_VAR) or Path.cwd())


def data_path(*parts):
    r"""A cache path as a string; no directories are created."""
    return str(data_dir().joinpath(*parts))


def ensure_parent(path):
    r"""Create the parent directory and return ``path``."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    return path
