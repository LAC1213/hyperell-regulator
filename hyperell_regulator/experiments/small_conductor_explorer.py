"""Export one shared two-sheet cut per small-conductor curve."""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import QQ, RealField, matrix
from hyperell_regulator.regulator.explorer import cut_system_data, parse_expression
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator import IntegrationMonitor


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def compute_curve(entry, out, prec=15, edge_seconds=120):
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
    data['collection'] = dict(label=entry['label'], points=entry['points'], pole=pole)
    filename = entry['label']+'.json'
    save(out/'data'/filename, data)
    clean = {k:v for k,v in entry.items() if k != 'exact_functions'}
    return dict(**clean, file='data/'+filename,
                determinant=str(determinant), normalized_matrix=[[str(v) for v in row] for row in rows],
                validation=dict(determinant_relative_error=str(error)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=Path('docs/explorer_inputs.json'))
    parser.add_argument('--out', type=Path, default=Path('docs/_extra/explorer'))
    parser.add_argument('--prec', type=int, default=15)
    parser.add_argument('--edge-seconds', type=float, default=120)
    args = parser.parse_args(argv)
    source = json.loads(args.inputs.read_text())
    curves = []
    for entry in source['curves']:
        curves.append(compute_curve(entry, args.out, args.prec, args.edge_seconds))
        print('Validated '+entry['label'], flush=True)
    save(args.out/'catalog.json', dict(schema_version=2, display_digits=10,
          computation_digits=args.prec, source=source['source'],
          input_sha256=hashlib.sha256(args.inputs.read_bytes()).hexdigest(), curves=curves))
    print('Saved %d curves, %d diagrams' % (len(curves), len(curves)))


if __name__ == '__main__':
    main()
