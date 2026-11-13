"""Compare direct two-sheet integrals with saved conductor-249 references."""
import argparse
import json
from pathlib import Path
from time import perf_counter
from sage.all import Matrix, RealField, QQ
from hyperell_regulator.regulator import IntegrationMonitor, regulators_stokes
from hyperell_regulator.regulator.polynomials import ring


def compare(reference, digits):
    F=ring()(reference['model'].replace('u','x'))
    centers=[QQ(reference[k]) for k in ('u_P','u_Q')]
    field=RealField(int((digits+15)*3.322))
    details={}
    start=perf_counter()
    with IntegrationMonitor(max_edge_seconds=120,max_edge_steps=1000000) as monitor:
        per,columns=regulators_stokes(F,centers,prec=digits,details=details)
    elapsed=perf_counter()-start
    matrix=[[field(w)*field(columns[0][i]),field(w)*field(columns[1][i]),field(w)*field(per[i])]
            for i,w in enumerate((1,2,1))]
    determinant=Matrix(field,matrix).determinant()
    entry_error=max(abs(a-field(b)) for row,saved in zip(matrix,reference['matrix']) for a,b in zip(row,saved))
    det_error=abs(determinant-field(reference['det']))
    tolerance=field(10)**(-digits+3)
    assert entry_error<tolerance,entry_error
    assert det_error<tolerance,det_error
    numerator=field.pi()*field(reference['c_plus'])*field(reference['c_minus'])*field(reference['L_second_derivative'])
    cut=details['cut']
    print('%s: %d digits, %.3fs, %d shared edges, matrix difference %s' %
          (reference['label'],digits,elapsed,len(monitor.records),entry_error),flush=True)
    return dict(label=reference['label'],digits=digits,model=str(F),degree=2,
                centers=[str(a) for a in centers],matrix=[[str(v) for v in row] for row in matrix],
                determinant=str(determinant),ratio=str(numerator/determinant),
                log_abs_x_wedge_imaginary_part=str(-4*matrix[1][1]),
                max_absolute_matrix_error=str(entry_error),absolute_determinant_error=str(det_error),
                seconds=elapsed,saved_high_degree_seconds=reference['regulator_seconds'],
                edges=len(monitor.records),nodes=[[float(z.real()),float(z.imag())] for z in cut.nodes])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits',type=int,default=100)
    parser.add_argument('--output',type=Path,default=Path('results/249_x_coordinate_comparison.json'))
    args=parser.parse_args(argv)
    rows=json.loads(Path('reports/beilinson_100_digits.json').read_text())
    references=[next(r for r in rows if r['label']=='249.a.249.1'),
                json.loads(Path('results/249_rational_approximations/input.json').read_text())]
    results=[compare(ref,args.digits) for ref in references]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    main()
