"""Exact support and full integration cuts in saved survey presentation."""

from copy import deepcopy
from decimal import Decimal

import pytest
from sage.all import QQ, QQbar

from hyperell_regulator._algebraic import decode_algebraic, encode_algebraic
from hyperell_regulator.experiments._survey_divisors import (
    comparison_error, decimal_value, divisor_metadata, packet_order,
)
from hyperell_regulator.experiments.geometric_torsion import encode_point, decode_point
from hyperell_regulator.regulator.explorer import parse_request


def example():
    # The LMFDB model is y^2=(x^6+1)/4; the cut model doubles y.
    return dict(label='test', points=[encode_point(P) for P in
                [(QQbar(0), QQbar(1)/2), (QQbar(0), -QQbar(1)/2), ('infinity', QQbar(1)/2)]],
                divisors=[['2', '0', '-2'], ['0', '1', '-1']], basis_indices=[2, 4],
                algebraic_data=dict(original_f='(x^6+1)/4', original_h='0', model='x^6+1', shift=None, scale='1'),
                regulator_input=dict(centers=[encode_algebraic(QQbar(0)), encode_algebraic(QQbar(2))]),
                packet_torsion_order='6', packet_torsion_order_exact=True, end_alg='Q')


def test_selected_rows_keep_coefficients_infinity_and_extra_cut_centers():
    row = example()
    data = divisor_metadata(row)
    assert [D['id'] for D in data['divisors']] == ['D3', 'D5']
    assert data['matrix_divisor_ids'] == ['D3', 'D5']
    assert data['divisors'][0]['tex'] == r'D_{3}=2\,\left(0,\frac{1}{2}\right)-2\,\infty_{Y=\frac{1}{2}}'
    assert data['support'][0]['coefficients'] == dict(D3='2', D5='0')
    assert data['support'][0]['center_index'] == data['support'][1]['center_index'] == 0
    assert data['support'][2]['center_index'] is None
    assert [a['exact'] for a in data['logarithmic_centers']] == row['regulator_input']['centers']
    assert data['packet_torsion_order'] == '6' and data['end_alg'] == 'Q'
    assert row == example()


def test_recorded_chart_transforms_infinities_and_finite_points_exactly():
    row = example()
    row['shift'] = row['algebraic_data']['shift'] = '0'
    row['regulator_input']['centers'] = [encode_algebraic(0)]
    data = divisor_metadata(row)
    assert [P['infinity'] for P in data['support']] == [True, True, False]
    assert decode_point(data['support'][2]['point']) == (QQbar(0), QQbar(1))
    assert [P['lmfdb']['infinity'] for P in data['support']] == [False, False, True]
    assert decode_point(data['support'][2]['lmfdb']['point']) == ('infinity', QQbar(1)/2)
    assert data['coordinate_system'] == 'lmfdb'
    assert data['lmfdb_model']['f'] == '1/4*x^6 + 1/4'
    row['algebraic_data']['model'] = 'x^6+2'
    with pytest.raises(ValueError, match='off the integration model'):
        divisor_metadata(row)


def test_old_basis_cut_reconstruction_cancels_ordinates_and_includes_unused_rows():
    row = example()
    del row['regulator_input']
    row['divisors'] = [['1', '-1', '0'], ['-1', '1', '0']]
    row['integral_kernels'] = [dict(points=deepcopy(row['points']), basis=row['divisors'])]
    data = divisor_metadata(row)
    assert data['logarithmic_centers'] == []
    row['integral_kernels'][0]['basis'].append(['1', '0', '-1'])
    assert divisor_metadata(row)['logarithmic_centers'][0]['exact'] == {'rational': '0'}


def test_complex_support_is_exact_and_never_projected_into_real_picture():
    row = example()
    i = QQbar(-1).sqrt()
    row['points'][0] = encode_point((i, QQbar(0)))
    row['regulator_input']['centers'].append(encode_algebraic(i))
    data = divisor_metadata(row)
    point = data['support'][0]
    assert not point['real'] and point['x'] == [0., 1.]
    assert point['tex'] == 'P_{1}'
    assert r'\theta_{1}' in point['definition_tex']
    assert decode_point(point['point']) == (i, QQbar(0))
    assert data['logarithmic_centers'][2]['position'] == [0., 1.]


def test_algebraic_center_roundtrip_and_more_than_eight_centers():
    values = [QQbar(2).sqrt(), QQbar(-1).sqrt(), -QQbar(-1).sqrt(), QQbar(QQ(3)/7)]
    for value in values:
        assert decode_algebraic(encode_algebraic(value)) == value
    centers = [encode_algebraic(n) for n in range(10)] + [encode_algebraic(values[1])]
    F, decoded, prec = parse_request(dict(F='x^6+1', centers=centers, prec=20))
    assert len(decoded) == 11 and decoded[-1] == values[1] and prec == 20


def test_packet_exponent_is_not_inferred_from_ambient_jacobian_torsion():
    row = example()
    row['packet_torsion_order_exact'] = False
    row['input_curve'] = dict(mw_invs=[100])
    assert packet_order(row) is None


def test_interval_midpoints_preserve_endpoints_and_error_uses_upper_bound():
    ball = dict(type='real_ball', lower='1/3', upper='2/3')
    assert Decimal(decimal_value(ball)) == Decimal('0.5')
    row = dict(comparison=dict(relative_agreement=ball))
    assert Decimal(comparison_error(row)) > Decimal('.66')
    assert ball == dict(type='real_ball', lower='1/3', upper='2/3')


def test_multi_determinant_comparison_exports_every_contributing_divisor():
    row = example()
    basis = [['2', '0', '-2'], ['0', '1', '-1'], ['1', '-1', '0']]
    row.update(basis_indices=[0, 2], divisors=[basis[0], basis[2]],
               pair_search_index=0, pair_searches=[dict(kernel_index=0)],
               integral_kernels=[dict(points=deepcopy(row['points']), divisor_basis=basis)],
               linear_comparison=dict(terms=[dict(basis_indices=[0, 2], coefficient='2'),
                                             dict(basis_indices=[1, 2], coefficient='-3')]))
    data = divisor_metadata(row)
    assert [D['id'] for D in data['divisors']] == ['D1', 'D2', 'D3']
    assert data['matrix_divisor_ids'] == ['D1', 'D3']
    assert data['support'][0]['coefficients'] == dict(D1='2', D2='0', D3='1')


def test_original_lmfdb_ordinates_and_both_infinity_labels():
    row = example()
    # y^2+(x^3+1)y=x^2+x has two infinities Y=0,-1. Its
    # square model has ordinates 2*y+h(x), which must not leak into labels.
    row['points'] = [encode_point(P) for P in [(QQbar(1), QQbar(-1)+QQbar(3).sqrt()),
                                             ('infinity', QQbar(0)), ('infinity', QQbar(-1))]]
    row['algebraic_data'].update(original_f='x^2+x', original_h='x^3+1',
                                 model='x^6+2*x^3+4*x^2+4*x+1')
    row['regulator_input']['centers'] = [encode_algebraic(1)]
    data = divisor_metadata(row)
    assert data['support'][1]['tex'] == r'\infty_{Y=0}'
    assert data['support'][2]['tex'] == r'\infty_{Y=-1}'
    assert decode_point(data['support'][0]['lmfdb']['point'])[1] == -1+QQbar(3).sqrt()
    assert decode_point(data['support'][0]['point'])[1] == 2*QQbar(3).sqrt()
    assert data['lmfdb_model']['h'] == 'x^3 + 1'
