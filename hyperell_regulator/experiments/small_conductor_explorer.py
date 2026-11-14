"""Export one shared two-sheet cut per small-conductor curve."""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import QQ, RealField, matrix
from hyperell_regulator.regulator.explorer import cut_system_data, parse_expression
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator import IntegrationMonitor, c_invariants
from hyperell_regulator.experiments import beilinson as b


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def numerical_source_hashes():
    package = Path(__file__).resolve().parents[1]
    names = ['experiments/small_conductor_explorer.py', 'regulator/stokes.py',
             'regulator/_tree.py', 'regulator/_periods.py', 'regulator/quadrature.py',
             'regulator/precision.py', 'experiments/beilinson.py', 'lfunction.py', 'frobenius.py']
    return {'hyperell_regulator/'+name: hashlib.sha256((package/name).read_bytes()).hexdigest()
            for name in names}


def compute_curve(entry, out, prec=15, edge_seconds=120, *, table_digits=100, curve=None):
    F = ring()(parse_expression(entry['F'], ring().fraction_field()))
    centers = [QQ(point[0]) for point in entry['points']]
    pole = entry.get('pole', 'infinity')
    if pole != 'infinity':
        centers.append(QQ(pole[0]))
    with IntegrationMonitor(max_edge_seconds=edge_seconds, max_edge_steps=1000000):
        data = cut_system_data(F, centers, prec)
    R = RealField(int(3.322*(prec+15)))
    per = [R(p['period']) for p in data['weights'][0]['pairings']]
    regs = [[R(p['regulator']) for p in w['pairings']] for w in data['weights']]
    if pole != 'infinity':
        base = regs[data['input']['centers'].index(str(QQ(pole[0]))) ]
        regs = [[a-b for a,b in zip(col,base)] for col in regs]
    # Equal x-coordinates give identical weights, even for distinct points.
    chosen = [regs[data['input']['centers'].index(str(QQ(p[0])))] for p in entry['points']]
    rows = [[w*chosen[0][i], w*chosen[1][i], w*per[i]] for i,w in enumerate((1,2,1))]
    determinant = matrix(R, rows).determinant()
    reference = R(entry['numerical']['det'])
    error = abs(determinant-reference)/max(abs(reference), R('1e-100'))
    if error > R(10)**(-min(prec-4, 10)):
        raise ArithmeticError('Fresh determinant disagrees with saved reference: '+entry['label'])
    # The survey is a regression reference, not the source of displayed numbers.
    if curve is None:
        raise ValueError('The original curve is required to validate its cached L-value')
    lvalue = b.standard_lvalue_cached(entry['label'], curve, digits=table_digits)
    if lvalue.get('L_second_derivative') is None or lvalue.get('L_precision_warning'):
        raise ArithmeticError('L-value does not support the table precision: '+entry['label'])
    with IntegrationMonitor(max_edge_seconds=edge_seconds, max_edge_steps=1000000):
        table_det, table_rows, _ = b.beilinson_determinant_stokes(
            F, *entry['points'], digits=table_digits, pole=pole)
    cp, cm = c_invariants(F, prec=table_digits)
    comparison = b.beilinson_comparison(table_det, cp, cm, lvalue['L_second_derivative'],
                                      bound=table_digits//2, digits=table_digits)
    agreement = comparison.get('relative_agreement')
    verified = bool(comparison.get('rational')) and agreement is not None and agreement < R('1e-95')
    table_R = b._sage_real(1, table_digits).parent()
    table_error = abs(table_R(table_det)/table_R(entry['numerical']['det'])-1)
    if table_error > table_R(10)**(-table_digits+3):
        raise ArithmeticError('Table determinant disagrees with saved reference: '+entry['label'])
    number = lambda value: table_R(value).str(truncate=False)
    data['collection'] = dict(label=entry['label'], points=entry['points'], pole=pole)
    filename = entry['label']+'.json'
    save(out/'data'/filename, data)
    clean = {k:v for k,v in entry.items() if k != 'exact_functions'}
    clean.update(file='data/'+filename, determinant=number(table_det),
                 normalized_matrix=[[number(v) for v in row] for row in table_rows],
                 numerical=dict(det=number(table_det), c_plus=number(cp), c_minus=number(cm),
                                L_second_derivative=lvalue['L_second_derivative']),
                 q=str(comparison['rational']) if verified else None,
                 quotient=number(comparison['ratio']),
                 rational_guess=str(comparison['rational']) if comparison.get('rational') else None,
                 relative_agreement=number(agreement) if agreement is not None else None,
                 comparison_verified=verified, period_basis='mst', regulator_digits=table_digits,
                 L_digits=lvalue['L_digits'],
                 validation=dict(determinant_relative_error=str(error),
                                 table_determinant_relative_error=str(table_error),
                                 L_num_coeffs=lvalue['L_num_coeffs'],
                                 L_cached=lvalue['L_cached'], feq_check=lvalue['feq_check']))
    return clean


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=Path('docs/explorer_inputs.json'))
    parser.add_argument('--out', type=Path, default=Path('docs/_extra/explorer'))
    parser.add_argument('--prec', type=int, default=15)
    parser.add_argument('--table-digits', type=int, default=100)
    parser.add_argument('--edge-seconds', type=float, default=120)
    args = parser.parse_args(argv)
    source = json.loads(args.inputs.read_text())
    records = {r['label']: r for r in b.load_results(b.CACHE)}
    curves = []
    for entry in source['curves']:
        curves.append(compute_curve(entry, args.out, args.prec, args.edge_seconds,
                                    table_digits=args.table_digits,
                                    curve=b.curve_from_row(records[entry['label']])))
        print('Validated '+entry['label'], flush=True)
    save(args.out/'catalog.json', dict(schema_version=2, display_digits=10,
          computation_digits=args.prec, table_digits=args.table_digits,
          period_basis='mst', source=source['source'],
          numerical_source_sha256=numerical_source_hashes(),
          input_sha256=hashlib.sha256(args.inputs.read_bytes()).hexdigest(), curves=curves))
    print('Saved %d curves, %d diagrams' % (len(curves), len(curves)))


if __name__ == '__main__':
    main()
