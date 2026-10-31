"""Saved manuscript inputs are transported to the current determinant basis."""

from decimal import Decimal, localcontext
from pathlib import Path
import runpy


def test_table_regeneration_rescales_determinant_and_quotient(tmp_path):
    root = Path(__file__).resolve().parents[1]
    report = runpy.run_path(str(root/'reports/remake_paper_tables.py'))
    source = report['load'](root/'reports/paper_full_periods_input.json')
    original = source['rows'][0]
    rows = report['compute'](dict(rows=[original]), 40, 25)
    row = rows[0]
    with localcontext() as context:
        context.prec = 180
        assert Decimal(row['determinant']) == 2*Decimal(original['determinant'])
    assert row['new_paper_q'] == '-14'
    assert Decimal(row['relative_error']) < Decimal('1e-37')
    target = tmp_path/'tables.tex'
    report['render'](rows, target, 10)
    text = target.read_text()
    assert r'\pi c^+c^-L^{\prime\prime}=qD' in text
    assert r'2\operatorname{Re}(z)' in text
    assert r'2\pi' not in text
