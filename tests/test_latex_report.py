"""Table errors use full precision, and exported algebraic data are reproducible."""
import json
from pathlib import Path
from decimal import Decimal
from tempfile import TemporaryDirectory
from hyperell_regulator.experiments import latex_report as report


def test_table_error_is_absolute_and_not_rounded_to_float():
    assert Decimal(report.ratio_error('0.12500000000000000000000000000000000000001', '1/8')) == Decimal('1e-41')
    assert report.ratio_error(None, None) is None
    assert report.scientific('1.23456e-101', 3) == r'1.23\times 10^{-101}'


def test_saved_table_snapshot_is_complete_and_self_contained():
    path = Path(__file__).resolve().parents[1]/'reports'/'beilinson_100_digits.json'
    rows = json.loads(path.read_text())
    assert len(rows) == 16 and all(r['status'] == 'ok' for r in rows)
    assert sum(r.get('L_second_derivative') is None for r in rows) == 3
    with TemporaryDirectory() as d:
        tex = Path(d)/'table.tex'
        exported = report.export(path, tex, cache=Path(d)/'nonexistent.json')
        assert len(exported) == 16
        text = tex.read_text()
        assert text.count(r'\phi_P=') == text.count(r'\phi_Q=') == 16
        assert r'\delta=|q-r|' in text
        assert text.count(r'\begin{longtable}') == text.count(r'\end{longtable}') == 2
        assert text.count('{') == text.count('}')
        assert all(r['label'] in text for r in rows)
        assert json.loads(tex.with_suffix('.json').read_text()) == rows
