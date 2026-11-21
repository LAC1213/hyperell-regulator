"""Fetching overlaps sequential analysis and publishes only complete cache rows."""
import io
import json
from threading import Event, current_thread, main_thread
from urllib.parse import parse_qs, urlparse
from unittest.mock import patch

import pytest

from hyperell_regulator import components
from hyperell_regulator.experiments import beilinson as b


def curve(label, cond, algebra='CM'):
    return dict(label=label, cond=cond, end_alg=algebra, is_simple_base=True,
                is_simple_geom=False, eqn=[[1, 1], [0]])


@pytest.mark.parametrize('fail_later', [False, True])
def test_fetch_publishes_cached_and_new_point_rows_before_finishing(tmp_path, fail_later):
    cache = tmp_path/'curves.json'
    old = dict(curve('1.q', 1, 'Q'), end_alg=None, rat_pts=[[0, 1, 1]], mw_invs=[7])
    cache.write_text(json.dumps([old]))
    catalog = [curve('1.q', 1, 'Q'), curve('2.cm', 2),
               curve('3.rm', 3, 'RM'), curve('4.cm', 4)]
    ready = []

    def api(request, **kwargs):
        query = parse_qs(urlparse(request.full_url).query)
        offset = int(query.get('_offset', ['0'])[0])
        if 'g2c_curves' in request.full_url:
            if offset == 0:
                rows = catalog[:3]
            else:
                assert ready == ['1.q']
                rows = catalog[3:] if offset == 3 else []
        elif offset == 0:
            rows = [dict(label='2.cm', rat_pts=[[2, 0, 1]], mw_invs=[2])]
        else:
            assert ready == ['1.q', '2.cm']
            if fail_later:
                raise OSError('connection lost')
            rows = [dict(label='4.cm', rat_pts=[[4, 0, 1]], mw_invs=[4])]
        return io.BytesIO(json.dumps(dict(data=rows)).encode())

    def on_curve(row):
        saved = {r['label']: r for r in json.loads(cache.read_text())}
        assert saved[row['label']] == row
        assert row['end_alg'] in ('Q', 'CM') and row['rat_pts']
        if row['label'] == '2.cm':
            assert '4.cm' not in saved  # Never publish incomplete new rows.
        ready.append(row['label'])

    with patch('urllib.request.urlopen', side_effect=api), patch('time.sleep'):
        if fail_later:
            with pytest.raises(OSError, match='connection lost'):
                components.fetch_lmfdb(max_cond=4, path=cache, end_algs=('Q', 'CM'), on_curve=on_curve)
        else:
            components.fetch_lmfdb(max_cond=4, path=cache, end_algs=('Q', 'CM'), on_curve=on_curve)
    assert ready == (['1.q', '2.cm'] if fail_later else ['1.q', '2.cm', '4.cm'])
    saved = {r['label']: r for r in json.loads(cache.read_text())}
    assert saved['1.q']['rat_pts'] == old['rat_pts']
    assert saved['2.cm']['mw_invs'] == [2]
    assert '3.rm' not in saved


@pytest.mark.parametrize('fail_later', [False, True])
def test_experiment_analyses_and_saves_before_fetch_finishes(tmp_path, fail_later):
    cache, records, out = (tmp_path/name for name in ('curves.json', 'pairs.json', 'out.json'))
    rows = [curve('1.cm', 1), curve('2.cm', 2)]
    records.write_text(json.dumps(rows))
    first_saved, fetch_done = Event(), Event()
    analysed = []

    def fetch(**kwargs):
        assert current_thread() is not main_thread()
        assert kwargs['end_algs'] == ('CM',)
        cache.write_text(json.dumps(rows[:1]))
        kwargs['on_curve'](rows[0])
        assert first_saved.wait(3), 'analysis waited for the entire fetch'
        if fail_later:
            raise OSError('connection lost')
        cache.write_text(json.dumps(rows))
        kwargs['on_curve'](rows[1])
        fetch_done.set()
        return rows

    def analyse(record, **kwargs):
        assert current_thread() is main_thread()
        if record['label'] == '1.cm':
            assert not fetch_done.is_set()
        analysed.append(record['label'])
        return dict(record, status='ok', det=1, comparison={})

    def report(result, **kwargs):
        if result['label'] == '1.cm':
            assert json.loads(out.read_text())[0]['label'] == '1.cm'
            first_saved.set()

    with patch.object(components, 'fetch_lmfdb', side_effect=fetch), \
            patch.object(b, 'analyse', side_effect=analyse), \
            patch.object(b, '_report_row', side_effect=report):
        if fail_later:
            with pytest.raises(OSError, match='connection lost'):
                b.run(3, cache=cache, path=records, out=out, fetch=True,
                      with_lvalue=False, end_algs=('CM',))
        else:
            result = b.run(3, cache=cache, path=records, out=out, fetch=True,
                           with_lvalue=False, end_algs=('CM',))
            assert [r['label'] for r in result] == ['1.cm', '2.cm']
    assert analysed == (['1.cm'] if fail_later else ['1.cm', '2.cm'])
    assert len(json.loads(out.read_text())) == len(analysed)


def test_closing_stream_cancels_producer_even_when_queue_is_full(tmp_path):
    stopped = Event()

    def fetch(**kwargs):
        try:
            for i in range(1000):
                kwargs['on_curve'](curve(str(i), i))
        finally:
            stopped.set()

    with patch.object(components, 'fetch_lmfdb', side_effect=fetch):
        stream = b._fetch_curves(1001, tmp_path/'unused.json', ('CM',))
        assert next(stream)['cond'] == 0
        stream.close()
        assert stopped.wait(2)


def test_empty_fetch_finishes_without_analysis(tmp_path):
    cache = tmp_path/'curves.json'
    with patch('urllib.request.urlopen', return_value=io.BytesIO(b'{"data": []}')), \
            patch.object(b, 'analyse') as analyse:
        assert b.run(100, cache=cache, path=tmp_path/'pairs.json', out=None,
                     fetch=True, with_lvalue=False, end_algs=('CM',)) == []
    analyse.assert_not_called()
    assert json.loads(cache.read_text()) == []
