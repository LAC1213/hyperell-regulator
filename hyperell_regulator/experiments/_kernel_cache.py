"""Exact packet kernels cached independently of numerical survey reports."""

import hashlib
import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile

from ..paths import data_path

KERNEL_CACHE_VERSION = 4  # Galois traces over QQ with exact packet torsion orders.


def cached_kernel(f, h, support, compute, *, incomplete_revision=None):
    """Reuse the kernel for an exact model and ordered support."""
    identity = [KERNEL_CACHE_VERSION, list(map(str, f)), list(map(str, h)), support]
    key = hashlib.sha256(json.dumps(identity, sort_keys=True, default=str).encode()).hexdigest()
    path = Path(data_path('results', 'integral_bases', key + '.json'))
    try:
        saved = json.loads(path.read_text())
        if (incomplete_revision is None or saved['complete']
                or saved.get('incomplete_revision') == incomplete_revision):
            return saved
    except (FileNotFoundError, json.JSONDecodeError):
        pass

    kernel = compute()
    if incomplete_revision is not None and not kernel['complete']:
        kernel['incomplete_revision'] = incomplete_revision
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(kernel, stream, default=str)
            stream.write('\n')
            stream.close()
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
    return kernel
