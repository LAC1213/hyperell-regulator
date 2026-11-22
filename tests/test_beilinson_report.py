"""Divisor survey cache policies, formatting and numerical comparisons."""
import contextlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from sage.all import RealField
from hyperell_regulator.experiments import beilinson as b


def test_completed_regulator_with_failed_lvalue_is_retried():
    settings = dict(report_version=b.REPORT_VERSION, regulator_digits=100,
                    lvalue_digits=100, with_lvalue=True)
    old = dict(status='ok', settings=settings, det='5.071007759569784',
               L_second_derivative=None,
               L_status='PariError: the PARI stack overflows')
    assert not b._stored_result_satisfies(old, settings)
    recovered = dict(old, L_status=None, L_second_derivative='0.123456789')
    assert b._stored_result_satisfies(recovered, settings)


def test_cache_policy_retries_failures_and_invalidates_old_reports():
    records=[dict(label='a',cond=249),dict(label='b',cond=295)]
    calls=[]
    def analyse(record,**kwargs):
        calls.append((record['label'],kwargs))
        return dict(record,status='ok' if record['label']=='a' else 'numerical_failure',
                    det=1.25 if record['label']=='a' else None,comparison={})
    with TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()):
        out=Path(d)/'report.json'
        out.write_text(json.dumps([dict(records[0],status='ok',det=999)]))
        with patch.object(b,'survey_curves',return_value=records),patch.object(b,'analyse',side_effect=analyse):
            first=b.run(curve_seconds=0, out=out,with_lvalue=False)
            assert len(calls)==2 and first[0]['det']==1.25
            assert all(k['reg_digits'] is None and k['digits']==15 for _,k in calls)
            calls.clear();b.run(curve_seconds=0, out=out,with_lvalue=False)
            assert [label for label,_ in calls]==['b']
            stored = json.loads(out.read_text())
            stored[0]['regulator_method'] = 'isogeny_transport'
            out.write_text(json.dumps(stored))
            calls.clear();b.run(curve_seconds=0, out=out,with_lvalue=False)
            assert [label for label,_ in calls]==['a','b']
            calls.clear();b.run(curve_seconds=0, out=out,with_lvalue=False,force=True)
            assert len(calls)==2
        assert len(json.loads(out.read_text()))==2


def test_higher_precision_report_is_reused_without_rounding():
    record = dict(label='a', cond=249)
    def analyse(record, **kwargs):
        return dict(record, status='ok', det='1.2345678901234567890123456789', comparison={}, L_second_derivative='1')
    with TemporaryDirectory() as d, contextlib.redirect_stdout(io.StringIO()):
        out = Path(d)/'report.json'
        with patch.object(b, 'survey_curves', return_value=[record]), patch.object(b, 'analyse', side_effect=analyse) as calc:
            b.run(curve_seconds=0, out=out, reg_digits=50, digits=50, max_edge_seconds=120)
            saved = out.read_text()
            calc.reset_mock()
            b.run(curve_seconds=0, out=out, reg_digits=30, digits=30)
            b.run(curve_seconds=0, out=out)
            calc.assert_not_called()
            assert out.read_text() == saved
            b.run(curve_seconds=0, out=out, reg_digits=60, digits=50)
            assert calc.call_count == 1
            b.run(curve_seconds=0, out=out, reg_digits=30, digits=30, force=True)
            assert calc.call_count == 2


def test_precision_cache_boundaries():
    settings = dict(report_version=b.REPORT_VERSION, regulator_digits=30,
                    lvalue_digits=30, with_lvalue=True)
    old = dict(status='ok', settings=dict(settings, regulator_digits=50, lvalue_digits=50), L_second_derivative='1')
    assert b._stored_result_satisfies(old, settings)
    assert not b._stored_result_satisfies(old, dict(settings, lvalue_digits=60))
    assert not b._stored_result_satisfies(old, dict(settings, report_version=-1))
    old['settings']['regulator_digits'] = None
    assert not b._stored_result_satisfies(old, settings)

    old['settings']['regulator_digits'] = 10
    assert not b._stored_result_satisfies(old, dict(settings, regulator_digits=None))
    old['status'] = 'numerical_failure'
    assert not b._stored_result_satisfies(old, settings)


def test_report_prints_exact_torsion_divisors_from_saved_packet(capsys):
    result = dict(label='test', cond=1, status='ok', integral=True,
                  points=[dict(field_polynomial=None, coordinates=['0', '1']),
                          dict(field_polynomial=None, coordinates=['infinity', '-1']),
                          dict(field_polynomial=['-2', '0', '1'], coordinates=[['0', '1'], ['0']])],
                  divisors=[['1', '1', '-2'], ['0', '-1', '1']])
    b._report_header()
    b._report_row(result, stored=True)
    header, line = capsys.readouterr().out.splitlines()
    assert header.endswith('torsion divisors')
    assert line.endswith('D1 = (0, 1) + oo(Y=-1) - 2*P3 | '
                         'D2 = -oo(Y=-1) + P3 [stored]')


def test_comparison_uses_pi_over_determinant_and_reports_integer_guesses():
    R = RealField(180)
    for q in (R(7), R(-28), R(3)/2):
        result = b.beilinson_comparison(R.pi()*6/q, 2, 3, 1, digits=40)
        assert abs(result['ratio']-q) < R('1e-39')
        assert abs(R(result['rational'])-q) < R('1e-39')
        assert result['rational_is_integer'] == (q in (7, -28))
        assert result['relative_agreement'] < R('1e-39')
    assert b.beilinson_comparison(0, 2, 3, 1) == {}
    assert b.beilinson_comparison(1, 2, 3, None) == {}


def test_old_determinant_convention_is_not_reused():
    settings = dict(report_version=b.REPORT_VERSION, regulator_digits=30,
                    lvalue_digits=30, with_lvalue=True)
    old = dict(status='ok', settings=dict(settings, report_version=8))
    assert not b._stored_result_satisfies(old, settings)


def test_integer_relation_does_not_discard_tiny_determinants():
    R = RealField(180)
    det = R('1e-40')
    result = b.beilinson_comparison(det, 1, 1, 40960*det/R.pi(), digits=40)
    assert result['rational'] == 40960
    assert result['rational_is_integer']
    assert result['relative_agreement'] < R('1e-39')
