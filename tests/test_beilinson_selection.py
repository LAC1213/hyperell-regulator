"""The survey selects rational endomorphism labels Q or CM independently of Qbar."""
import io
import json
from urllib.parse import urlparse,parse_qs
from unittest.mock import patch

import pytest

from hyperell_regulator.components import fetch_lmfdb,load_curves
from hyperell_regulator.experiments import beilinson as b


def test_selection_uses_rational_endomorphisms_not_geometric_simplicity(tmp_path):
    rows=[dict(label='q',cond=1,is_simple_base=True,is_simple_geom=False,end_alg='Q'),
          dict(label='rm',cond=2,is_simple_base=True,is_simple_geom=True,end_alg='RM'),
          dict(label='cm',cond=3,is_simple_base=True,is_simple_geom=False,end_alg='CM')]
    cache=tmp_path/'curves.json';cache.write_text(json.dumps(rows))
    assert [r['label'] for r in load_curves(cache,end_q_z=True)] == ['q']
    assert [r['label'] for r in load_curves(cache,end_algs=('Q','CM'))] == ['q','cm']
    assert [r['label'] for r in b.survey_curves(4,cache=cache)] == ['q','cm']
    assert [r['label'] for r in b.survey_curves(4,cache=cache,end_algs=('CM',))] == ['cm']


def test_cm_only_cli_preserves_default_and_certified_options():
    with patch.object(b, 'run', return_value=[]) as run:
        b.main(['1300'])
        assert run.call_args.kwargs['end_algs'] == ('Q', 'CM')
        b.main(['1300', '--cm-only', '--certified', '--digits', '100', '--reg-digits', '100'])
        assert run.call_args.kwargs['end_algs'] == ('CM',)
        assert run.call_args.kwargs['certified'] is True
        assert run.call_args.kwargs['digits'] == run.call_args.kwargs['reg_digits'] == 100


def test_missing_endomorphism_evidence_is_not_geometric_simplicity(tmp_path):
    row=dict(label='q',cond=1,is_simple_geom=True)
    cache=tmp_path/'curves.json';cache.write_text(json.dumps([row]))
    with pytest.raises(ValueError,match='--fetch'):
        b.survey_curves(3,cache=cache)
    row['eqn'] = [[0, -1, 0, 0, 0, 1], []]
    assert '--fetch' in b.analyse(row, verbose=False, end_algs=b.BEILINSON_END_ALGS)['status']


def test_refresh_backfills_metadata_and_preserves_cached_points(tmp_path):
    cache=tmp_path/'curves.json'
    old=[dict(label='q',cond=1,is_simple_geom=True,eqn=[[1],[0]],rat_pts=[[0,1,1]],mw_invs=[7]),
         dict(label='rm',cond=2,is_simple_geom=True,eqn=[[1],[0]],rat_pts=[[0,1,1]])]
    cache.write_text(json.dumps(old))
    catalog=[dict(label=label,cond=i,is_simple_base=True,is_simple_geom=geom,end_alg=endo,eqn=[[1],[0]])
             for i,(label,geom,endo) in enumerate([('q',True,'Q'),('rm',True,'RM'),('new',False,'CM')],1)]
    urls=[]
    def api(request,**kwargs):
        urls.append(request.full_url)
        query=parse_qs(urlparse(request.full_url).query)
        if 'g2c_curves' in request.full_url:
            rows=catalog if query.get('_offset',['0'])==['0'] else []
        else:
            rows=[dict(label='new',rat_pts=[[1,0,1]],mw_invs=[2])]
        return io.BytesIO(json.dumps(dict(data=rows)).encode())
    with patch('urllib.request.urlopen',side_effect=api),patch('time.sleep'):
        data=fetch_lmfdb(max_cond=3,path=cache,end_algs=('Q','CM'))
    bylabel={r['label']:r for r in data}
    assert bylabel['q']['rat_pts']==old[0]['rat_pts'] and bylabel['q']['mw_invs']==[7]
    assert bylabel['rm']['end_alg']=='RM' and not bylabel['new']['is_simple_geom']
    assert [r['label'] for r in load_curves(cache,end_algs=('Q','CM'))]==['q','new']
    query=parse_qs(urlparse(urls[0]).query)
    assert 'end_alg' in query['_fields'][0]
    assert 'is_simple_geom' not in query and query['is_simple_base']==['true']


def test_q_only_policy_does_not_satisfy_q_or_cm_request():
    settings=dict(report_version=b.REPORT_VERSION,regulator_digits=20,lvalue_digits=20,
                  with_lvalue=False,endomorphism_condition=b._endomorphism_condition(('CM','Q')))
    old=dict(status='ok',settings=dict(settings,endomorphism_condition='End_Q(J)=Z'))
    assert not b._stored_result_satisfies(old,settings)
