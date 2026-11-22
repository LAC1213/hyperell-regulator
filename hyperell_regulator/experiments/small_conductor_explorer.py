"""Export one shared two-sheet cut per small-conductor curve."""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import QQ, RealField, matrix
from hyperell_regulator._algebraic import decode_algebraic
from hyperell_regulator.experiments.geometric_torsion import decode_point
from hyperell_regulator._source_hash import SOURCE_HASH_METHOD, source_sha256
from hyperell_regulator.regulator.explorer import cut_system_data, parse_expression
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator import IntegrationMonitor


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def numerical_source_hashes():
    package = Path(__file__).resolve().parents[1]
    names = ['_source_hash.py', 'experiments/small_conductor_explorer.py', 'regulator/stokes.py',
             'regulator/_tree.py', 'regulator/_periods.py', 'regulator/quadrature.py',
             'regulator/precision.py', 'regulator/explorer.py', '_algebraic.py',
             'experiments/_survey_divisors.py', 'experiments/beilinson_report.py', 'experiments/beilinson.py',
             'experiments/beilinson_local_data.py', 'lfunction.py', 'frobenius.py']
    return {'hyperell_regulator/'+name: source_sha256((package/name).read_bytes())
            for name in names}


def compute_curve(entry, out, prec=15, edge_seconds=120, *, table_digits=100, curve=None):
    """Validate divisor columns on the full cut against the completed experiment.

    The displayed high-precision values come from the supplied successful
    run. Fresh lower-precision line integrals check the divisor assembly;
    they never replace certified values or request another L-value.
    """
    F = ring()(parse_expression(entry['F'], ring().fraction_field()))
    centers = [decode_algebraic(c['exact']) for c in entry['logarithmic_centers']]
    center_indices = {x: i for i, x in enumerate(centers)}
    coefficients = []
    for D in entry['divisors']:
        weights = [QQ(0)]*len(centers)
        for P in entry['support']:
            c = QQ(P['coefficients'].get(D['id'], 0))
            if not c or P['infinity']:
                continue
            x = decode_point(P['point'])[0]
            if x not in center_indices:
                # A hyperelliptic pair can cancel before integration.
                if sum(QQ(Q['coefficients'].get(D['id'], 0)) for Q in entry['support']
                       if not Q['infinity'] and decode_point(Q['point'])[0] == x):
                    raise ValueError('Selected divisor uses a center absent from its recorded cut')
                continue
            weights[center_indices[x]] += c
        coefficients.append(weights)
    out = Path(out)
    reference = entry['reference_matrix']
    # Close roots and nearly dependent columns can need guard digits for
    # a determinant even though only the picture's integrals are exported.
    for work in (prec, prec+10, prec+20):
        with IntegrationMonitor(max_edge_seconds=edge_seconds, max_edge_steps=1000000):
            data = cut_system_data(F, centers, work)
        R = RealField(int(3.33*(work+20)))
        per = [R(p['period']) for p in data['weights'][0]['pairings']]
        regs = [[R(p['regulator']) for p in w['pairings']] for w in data['weights']]
        columns = {D['id']: [sum(R(c)*regs[j][i] for j, c in enumerate(weights)) for i in range(3)]
                   for D, weights in zip(entry['divisors'], coefficients)}
        left, right = [columns[name] for name in entry['matrix_divisor_ids']]
        rows = [[R(w)*left[i], R(w)*right[i], R(w)*per[i]] for i, w in enumerate((1, 2, 1))]
        determinant = matrix(R, rows).determinant()
        target = R(entry['numerical']['det'])
        error = abs(determinant-target)/max(abs(target), R('1e-100'))
        entry_error = max(abs(rows[i][j]-R(reference[i][j]))/max(1, abs(R(reference[i][j])))
                          for i in range(3) for j in range(3))
        if error < R('1e-8') and entry_error < R('1e-8'):
            break
    else:
        raise ArithmeticError('Fresh divisor matrix disagrees with saved experiment: '+entry['label'])
    collection = {key: entry[key] for key in ('label', 'divisors', 'support', 'logarithmic_centers',
                  'matrix_divisor_ids', 'packet_torsion_order', 'end_alg', 'coordinate_system', 'lmfdb_model', 'integral_rank', 'integral_rank_field',
                  'integral_rank_space', 'trace_rank_certificate')}
    collection['divisor_center_coefficients'] = [list(map(str, D)) for D in coefficients]
    data['collection'] = collection
    filename = entry['label']+'.json'
    clean = dict(entry)
    clean.update(file='data/'+filename, determinant=entry['numerical']['det'],
                 normalized_matrix=reference, period_basis='mst',
                 validation=dict(determinant_relative_error=str(error), matrix_relative_error=str(entry_error),
                     source='saved_experiment', certified=entry.get('certified', False),
                     L_num_coeffs=entry.get('L_num_coeffs'), feq_check=entry.get('feq_check')))
    data['collection']['experiment'] = {key: clean.get(key) for key in (
        'numerical', 'normalized_matrix', 'reference_intervals', 'certified', 'regulator_digits',
        'L_digits', 'L_certified', 'L_analytic_assumptions', 'q', 'quotient', 'rational_guess',
        'comparison_verified', 'relative_agreement', 'linear_comparison')}
    save(out/'data'/filename, data)
    return clean


def _compute_entry(args):
    entry, out, prec, edge_seconds, table_digits = args
    return compute_curve(entry, out, prec, edge_seconds, table_digits=table_digits)


def main(argv=None):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from multiprocessing import get_context
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=Path('docs/explorer_inputs.json'))
    parser.add_argument('--out', type=Path, default=Path('docs/_extra/explorer'))
    parser.add_argument('--prec', type=int, default=15)
    parser.add_argument('--table-digits', type=int, default=100)
    parser.add_argument('--edge-seconds', type=float, default=120)
    parser.add_argument('--jobs', type=int, default=1)
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    source = json.loads(args.inputs.read_text())
    jobs = [(entry, args.out, args.prec, args.edge_seconds, args.table_digits) for entry in source['curves']]
    curves = [None]*len(jobs)
    if args.jobs == 1:
        for i, job in enumerate(jobs):
            curves[i] = _compute_entry(job)
            print('Validated '+curves[i]['label'], flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs, mp_context=get_context('spawn')) as pool:
            pending = {pool.submit(_compute_entry, job): i for i, job in enumerate(jobs)}
            for future in as_completed(pending):
                i = pending[future]
                curves[i] = future.result()
                print('Validated '+curves[i]['label'], flush=True)
    save(args.out/'catalog.json', dict(schema_version=3, display_digits=10,
          computation_digits=args.prec, table_digits=args.table_digits,
          period_basis='mst', source=source['source'],
          numerical_source_hash_method=SOURCE_HASH_METHOD,
          numerical_source_sha256=numerical_source_hashes(),
          input_sha256=hashlib.sha256(args.inputs.read_bytes()).hexdigest(), curves=curves))
    print('Saved %d curves, %d diagrams' % (len(curves), len(curves)))

if __name__ == '__main__':
    main()
