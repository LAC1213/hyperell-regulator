"""Logarithmic integrals and the 100-digit conductor-249 example."""
import json
from pathlib import Path
from sage.all import RealField
from hyperell_regulator.regulator.stokes import Cut
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.motivic import DivisorRegulators
from hyperell_regulator.experiments.geometric_torsion import decode_point
from hyperell_regulator.experiments._survey_divisors import decimal_value


def test_regular_zero_has_identity_gluing_and_piecewise_log_jumps():
    x=ring().gen();cut=Cut(x**6+4*x**5+4*x**4+2*x**3+1,[-1,0])
    k=next(i for i,z in enumerate(cut.nodes) if z==0)
    assert cut.kinds[k]=='regular'
    # No sheet flip occurs around the regular vertex, including when it
    # is a tree junction: compose the incident gluing permutations.
    incident=[i for i,(a,b) in enumerate(cut.edges) if k in (a,b)]
    flips=sum(cut.upper[i]!=cut.lower[i] for i in incident)
    assert flips%2==0
    # Different logarithms share geometry but have different jumps before
    # the second marked center; no uniform-jump assumption is permitted.
    assert cut.jumps[0]!=cut.jumps[1]


def test_conductor_249_at_one_hundred_digits():
    row=next(r for r in json.loads((Path(__file__).parents[1]/'reports/beilinson_100_digits.json').read_text()) if r['label']=='249.a.249.1')
    F=ring()(row['model'].replace('u','x'))
    points=row['divisor_data']['support']
    names=row['divisor_data']['matrix_divisor_ids']
    abscissae=[decode_point(P['point'])[0] for P in points]
    divisors=[[P['coefficients'][name] for P in points] for name in names]
    det,M=DivisorRegulators(F,abscissae,digits=100).determinant(*divisors)
    R=RealField(400)
    assert max(abs(R(a)-R(decimal_value(b))) for rr,ss in zip(M,row['matrix']) for a,b in zip(rr,ss))<R('1e-97')
    assert abs(R(det)-R(decimal_value(row['det'])))<R('1e-96')
