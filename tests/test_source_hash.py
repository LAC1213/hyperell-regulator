"""Source provenance checks that can run without SageMath."""

import ast
import json
from pathlib import Path
import runpy
import textwrap
import unittest


ROOT = Path(__file__).resolve().parents[1]
# Loading this standard-library helper directly avoids the Sage package imports.
source_sha256 = runpy.run_path(str(ROOT/'hyperell_regulator/_source_hash.py'))['source_sha256']


class SourceHashTests(unittest.TestCase):
    def test_module_class_function_async_and_nested_docstrings_are_ignored(self):
        source = textwrap.dedent('''\
            """Original module documentation."""
            class Example:
                """Original class documentation."""
                def method(self):
                    """Original method documentation."""
                    def nested(x):
                        """Original nested documentation."""
                        return x * 2
                    return nested(3)
            async def async_function():
                """Original async documentation."""
                return 7
        ''')
        self.assertEqual(source_sha256(source), source_sha256(
            source.replace('Original', 'Rewritten\nwith more details about')))
        self.assertEqual(source_sha256(source), source_sha256('\n'.join(
            line for line in source.splitlines() if '"""' not in line)))

    def test_comments_and_formatting_are_ignored(self):
        self.assertEqual(
            source_sha256('def f(x):\n    # double x\n    return x * 2\n'),
            source_sha256('\n# revised comment\ndef f( x ): return (x*2)\n'))

    def test_code_and_runtime_strings_still_change_the_hash(self):
        source = 'def f(x):\n    result = x * 2\n    return result, "label"\n'
        for changed in (source.replace('* 2', '* 3'), source.replace('*', '+'),
                        source.replace('"label"', '"other label"'),
                        source.replace('return result', 'return x')):
            with self.subTest(source=changed):
                self.assertNotEqual(source_sha256(source), source_sha256(changed))

    def test_non_docstring_expression_is_retained(self):
        source = 'def f():\n    x = 1\n    "a later string"\n    return x\n'
        self.assertNotEqual(source_sha256(source), source_sha256(
            source.replace('a later string', 'a changed string')))

    def test_indentation_changes_are_retained(self):
        source = 'def f(x):\n    if x:\n        x += 1\n    return x\n'
        self.assertNotEqual(source_sha256(source), source_sha256(
            source.replace('    return x', '        return x')))

    def test_bytes_and_text_have_the_same_hash(self):
        source = '# UTF-8\nvalue = "π"\n'
        self.assertEqual(source_sha256(source), source_sha256(source.encode('utf-8')))

    def test_saved_catalog_matches_current_code_despite_docstring_edits(self):
        catalog = json.loads((ROOT/'docs/_extra/explorer/catalog.json').read_text())
        self.assertEqual(catalog['numerical_source_hash_method'], 'python-ast-no-docstrings-v1')
        for name, checksum in catalog['numerical_source_sha256'].items():
            with self.subTest(path=name):
                source = (ROOT/name).read_bytes()
                self.assertEqual(source_sha256(source), checksum)
                tree = ast.parse(source)
                self.assertIsNotNone(ast.get_docstring(tree))
                tree.body[0].value.value = 'New module documentation.'
                self.assertEqual(source_sha256(ast.unparse(tree)), checksum)


if __name__ == '__main__':
    unittest.main()
