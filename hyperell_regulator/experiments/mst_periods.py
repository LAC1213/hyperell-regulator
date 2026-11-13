"""Compute full MST cycle periods and compare their normalization with the survey."""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import RealField, QQ
from hyperell_regulator.regulator.precision import period_matrix_mp
from hyperell_regulator.regulator._periods import pivoted_minor
from hyperell_regulator.regulator.polynomials import ring


def compare(row, digits):
    F=ring()(row['model'].replace('u','x'))
    P=period_matrix_mp(F,digits)
    R=RealField(int((digits+30)*3.322))
    tolerance=R(10)**(-digits+3)
    cp,cm=[R(pivoted_minor([[part(z) for z in rr] for rr in P],R(10)**(-digits-3)))
           for part in (lambda z:z.real(),lambda z:z.imag())]
    factors=[]
    for value,key in ((cp,'c_plus'),(cm,'c_minus')):
        ratio=value/R(row[key])
        rational=ratio.nearby_rational(max_denominator=1024)
        if abs(ratio-R(rational))>tolerance:
            raise ArithmeticError('period normalizations do not have a resolved rational ratio')
        factors.append(rational)
    q=R.pi()*cp*cm*R(row['L_second_derivative'])/R(row['det'])
    old=row.get('comparison',{})
    rational=None
    if old.get('rational') and R(old.get('relative_agreement','1'))<R('1e-95'):
        rational=QQ(old['rational'])*factors[0]*factors[1]
        if abs(q/R(rational)-1)>tolerance:
            raise ArithmeticError('new period normalization disagrees with the regulator comparison')
    print('%s: c+/old=%s, c-/old=%s, q=%s' %
          (row['label'],factors[0],factors[1],rational if rational is not None else q.n(50)),flush=True)
    return dict(label=row['label'],model=str(F),digits=digits,
                c_plus=str(cp),c_minus=str(cm),c_plus_factor=str(factors[0]),c_minus_factor=str(factors[1]),
                old_c_plus=row['c_plus'],old_c_minus=row['c_minus'],
                quotient=str(q),rational=str(rational) if rational is not None else None,
                periods=[[[str(z.real()),str(z.imag())] for z in rr] for rr in P])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits',type=int,default=100)
    parser.add_argument('--output',type=Path,default=Path('results/mst_periods.json'))
    args=parser.parse_args(argv)
    source=Path('reports/beilinson_100_digits.json')
    rows=[r for r in json.loads(source.read_text()) if r.get('status')=='ok']
    extra=Path('results/249_rational_approximations/input.json')
    if extra.exists():rows.append(json.loads(extra.read_text()))
    results=[compare(row,args.digits) for row in rows]
    paths=[Path(__file__),Path('hyperell_regulator/regulator/_periods.py'),Path('hyperell_regulator/regulator/_tree.py'),Path('hyperell_regulator/regulator/precision.py')]
    output=dict(digits=args.digits,reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},curves=results)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':
    main()
