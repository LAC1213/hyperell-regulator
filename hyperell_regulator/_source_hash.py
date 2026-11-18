"""Fingerprint Python code without comments, formatting, or docstrings."""

import ast
import hashlib
import json


SOURCE_HASH_METHOD = 'python-ast-no-docstrings-v1'


def source_sha256(source):
    """Hash source bytes or text, retaining code and non-docstring literals."""
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)) and ast.get_docstring(node) is not None:
            node.body.pop(0)

    def canonical(node):
        if isinstance(node, ast.AST):
            # Omit empty optional fields (e.g. type_params added in Python 3.12)
            # so the same code has the same fingerprint across Python versions.
            return [type(node).__name__, [
                [name, canonical(value)] for name, value in ast.iter_fields(node)
                if value is not None and value != []]]
        if isinstance(node, list):
            return [canonical(value) for value in node]
        return repr(node)

    payload = json.dumps(canonical(tree), ensure_ascii=True, separators=(',', ':'))
    return hashlib.sha256(payload.encode('utf-8')).hexdigest()
