Integrality of combinations
==============================

``TorsionPacketIntegrality`` implements the supplied cocycle criterion for
classes on a rational torsion packet :math:`T`. Modulo the decomposable
Tate part, a combination depends on its total divisor:

.. math::

   \sum_j c_j\alpha_{\varphi_j} \equiv B(D),\qquad
   D=\sum_j c_j N_j(P_j-Q_j).

The sufficient criterion for :math:`q_*B(D)` is that :math:`D` has degree
zero on every component of every bad fibre of a regular proper model.
The implementation uses exact rational arithmetic. Rational coefficients
are allowed; kernel bases consist of integer coefficients.

For example, on LMFDB 277.a.277.1::

    from sage.all import QQ, PolynomialRing
    from hyperell_regulator import (
        HyperellipticCurveWithReduction, TorsionPacketIntegrality,
    )

    x = PolynomialRing(QQ, 'x').gen()
    C = HyperellipticCurveWithReduction(-x*x-x, x**3+x*x+x+1)
    T = [('infinity', 0), ('infinity', -1), (0, 0), (0, -1)]
    test = TorsionPacketIntegrality(C, T)
    result = test.check([1, -1, 0, 0])
    assert result['q_integral'] is True
    lattice = test.kernel()

Point coordinates refer to :math:`y^2+h(x)y=f(x)`. The second infinity
coordinate is :math:`\lim y/x^{g+1}`. The constructor verifies torsion of
every difference from the first point, and internally moves all points to
an affine chart for reduction. ``test.orders`` records their exact orders.
Non-torsion differences are rejected; undetermined orders raise an error.
``torsion_packets(C, points)`` can first partition a larger set of points
by verified torsion differences, retaining undetermined comparisons.

For classes with :math:`\operatorname{div}(\varphi_j)=N_j(P_j-O)`, supply
their divisors as generators::

    generators = []
    for j, N in enumerate(test.orders[1:], 1):
        D = [0] * len(T)
        D[j], D[0] = N, -N
        generators.append(D)
    combinations = test.kernel(generators)
    for coefficients in combinations['basis'].rows():
        D = test.combine(generators, coefficients)
        assert test.check(D)['q_integral'] is True

More general generators are arbitrary degree-zero coefficient vectors or
mappings from points to coefficients. In particular, they can themselves
represent sums :math:`B(\sum_i a_iP_i)`. ``combine`` aggregates all point
coefficients before reduction, so cancellation at unresolved points is
handled correctly. The integer kernel is saturated; independently clearing
denominators of a rational vector-space basis would not suffice.

``kernel()`` returns the basis in point coordinates. With supplied
generators, ``basis`` gives their integer combination coefficients and
``divisors`` gives the resulting point coefficients. Relations among the
generators are retained as zero rows in ``divisors``.

Interpretation and limits
-------------------------

``check`` returns ``criterion_satisfied`` equal to ``True``, ``False`` or
``None``. The value ``q_integral=True`` certifies integrality by the
multidegree criterion. A nonzero multidegree leaves the projected boundary
undetermined by this test; use the projected boundary calculation below.

Smooth points on a semistable stable component lie on its strict transform
in a regular model, so their component labels suffice. At 2, the code also
factors the original integral affine fibre and uses its smooth points:
resolving a proper compactification preserves this smooth open. This works
even when reduction is wild. Integral coordinates and the Jacobian
criterion are checked individually; reducible affine fibres are handled by
labelling the unique component through each smooth point. Points outside
this affine patch, points on contracted chains, and other unavailable
nonsemistable fibres remain unresolved.
Their net coefficients must vanish for the reported certificate. The
returned ``complete`` flag states whether every point was resolved; otherwise
the kernel is a sufficient sublattice of the full multidegree kernel.
Failure at a fully resolved prime can still disprove the sufficient
criterion when another prime is unresolved.

The additional mean-log/Tate condition for :math:`p'_*` is not tested, and
rescaling functions is not performed. A nonzero formal divisor, or a
positive formal kernel rank, does not certify a nonzero or independent
motivic class. In particular, hyperelliptic conjugacy may kill classes
after projection.

Projected boundaries and integral classes
-----------------------------------------

``ProjectedTorsionPacketIntegrality`` computes the local obstruction
after applying the standard projector. It uses the vertical potential
on the graph of a nodal fibre, with the diagonal Tate direction removed.
At 2 it also uses verified auxiliary charts for supported wild equations.
The module computes local constraints from the curve and points.

For 249.a.6723.1 the following combination is integral::

    from sage.all import QQ, PolynomialRing
    from hyperell_regulator import (
        HyperellipticCurveWithReduction, ProjectedTorsionPacketIntegrality,
    )

    R = PolynomialRing(QQ, 'x')
    C = HyperellipticCurveWithReduction(
        R([2, 3, 1, 1, 0, -1]), R([1, 0, 0, 1]))
    solver = ProjectedTorsionPacketIntegrality.from_abscissae(C, [3, -1, None, 0])
    D = solver.abscissa_divisor([3, 0, 1, -4])
    assert solver.check(D)['q_integral'] is True

    generators = [
        solver.abscissa_divisor([1, 0, -1, 0]),
        solver.abscissa_divisor([0, 1, -1, 0]),
        solver.abscissa_divisor([0, 0, -1, 1]),
    ]
    solution = solver.kernel(generators)
    assert solution['complete']
    assert solution['abscissa_basis'].nrows() == 2

Use ``ProjectedTorsionPacketIntegrality(C, points, field=K)`` to supply
point coordinates in an absolute number field. Optional
``torsion_multiples`` are verified before the local calculation. For
a rational curve, ``from_abscissae`` accepts exact algebraic abscissae
and constructs a common field for both ordinates. The module retains
distinct algebraic abscissae. It uses hyperelliptic conjugacy when the
two ordinates are in the verified packet, since they have equal
standard-projected classes.

``check(D)['q_integral']`` is ``True`` for an integral class,
``False`` for a witnessed nonzero projected boundary, and ``None``
when the local constraints leave a place undetermined. Local results
give the method, boundary values and reduction evidence.
``kernel()`` solves on the whole packet. With supplied generators,
``basis`` is the saturated integer lattice of their integral
combinations and ``divisors`` lists the corresponding point divisors.
``abscissa_basis`` gives the independent rational abscissa directions.

When a model is unavailable, the multidegree condition still
certifies integral combinations. The necessary constraints give
an upper bound and the sufficient constraints give a lower bound
on the kernel. ``complete`` holds when those bounds agree in the
requested span. Arbitrary new point supports on wild models can give
incomplete kernels.

Divisor classes and stable models
---------------------------------

Use ``MotivicPacket`` to work directly with degree-zero divisors. Its
default points are the six geometric Weierstrass points. Supply an
explicit list of torsion points or call
``MotivicPacket.from_geometric_packet(packet)`` to include nonbranch
points of a computed Weierstrass-based packet.

.. code-block:: python

   from sage.all import PolynomialRing, QQ
   from hyperell_regulator import HyperellipticCurveWithReduction, MotivicPacket

   x = PolynomialRing(QQ, 'x').gen()
   C = HyperellipticCurveWithReduction(x*(x-1)*(x-2)*(x-3)*(x-4))
   packet = MotivicPacket(C)
   D = packet.divisor({(0, 0): 1, (1, 0): -1})
   E = packet.divisor({(2, 0): 1, (3, 0): -1})
   determinant = D.determinant(E, digits=20)
   graph_boundary = D.boundary(3)
   integral_divisors, complete = packet.integral_basis()

A single divisor gives one logarithmic regulator column. Its weight
is :math:`\sum_i a_i\log|x-x(P_i)|` in the integration chart, with zero
weight for a mark at infinity. Two divisor columns and the period column
give the Beilinson determinant with weights
:math:`(1,2\operatorname{Re}(x),|x|^2)`. The evaluator uses a common chart
and cut for all requested columns. Certified integration accepts exact
algebraic abscissae in that chart.

The public ``stable_model(C, prime, field=K, points=points)`` constructor
accepts an absolute field and a prime ideal of that field. Over
:math:`\mathbf Q` use a rational prime. The returned ``graph`` retains
vertices for specializations, with ``stable_graph`` obtained by
contracting unstable rational components. On split tame charts the
component equations and node coordinates describe the normalizations
and their incidences. For a metric edge of length :math:`n_e` the
conductance is :math:`1/n_e`. On dyadic compact-type charts the lengths
use :math:`v(2)=1`, with the finite local chart extension specified by
the centre, the square completion and the equation for the disc scale.

To plot the stable fibre of a cached LMFDB curve, run::

    sage -python -m hyperell_regulator.dual_graph_plot 249.a.6723.1 3 --out reduction.svg

Choose ``.png``, ``.svg`` or ``.pdf`` for the output. The plot shows the genus
of each component and the thickness of each node in the model's valuation.
For a model already computed in Sage, use
``hyperell_regulator.dual_graph_plot.plot_dual_graph(model)`` to obtain
a Sage graphics object.

For a divisor with component multidegree :math:`d`, the graph solves
:math:`Lu=d` with one potential fixed to zero. The positive-genus
component pairings are :math:`u_i`, and the negative of

.. math::

   Q_u(\gamma,\delta)
     = \sum_{e=(i,j)} n_e (u_i+u_j)\gamma_e\delta_e

is the pairing on tensors of graph cycles in the product blowup
convention. ``GraphBoundary.projected`` gives coordinates after quotienting
by the diagonal Tate class. Inserting a vertex along an edge gives
linear interpolation of the potential and the same pairing on cycles.
Stacking the projected maps for the finite places gives the linear system
for ``integral_basis()``. By restriction and corestriction the rational
vanishing condition can be tested after finite extension. The returned
basis spans independent divisor directions modulo the hyperelliptic
identifications used by the standard projection.

``MotivicPacket.integral_basis()`` restricts these directions to Galois
traces over :math:`\mathbf Q` before solving the local equations. After
identifying the two ordinates over each abscissa, it uses normalized sums
over complete Galois orbits of abscissae. Degree-zero differences of
these sums span the trace space; denominators are cleared for the
experiment's integer divisor representatives. Incomplete orbits in an
explicitly supplied support contribute no direction. The experiment's
optional search ignoring integrality uses this same trace space.

For example, on ``394.a.394.1`` the quadratic conjugates :math:`P_2,P_3`
give :math:`\operatorname{Tr}(P_2-P_3)=0`. This locally integral direction
is therefore excluded from the basis over :math:`\mathbf Q`; the
remaining integral direction is
:math:`2(1,-1)-\infty_0-\infty_{-1}`.
``packet.integrality.kernel()`` retains the unrestricted number-field
kernel for local boundary calculations. The experiment invalidates
cached bases computed before the Galois-trace restriction.

At an odd prime with one thickness-one node on an irreducible fibre,
the two Weierstrass points specializing to that node meet the midpoint
of the loop after base change. The projected graph boundary vanishes
exactly when their coefficients sum to zero, provided all other marked
points specialize to smooth points of the fibre. When the colliding
pair comes from an irreducible branch factor of prime degree, inertia
acts by a transposition. The graph of its Galois conjugates is connected
by primitivity, so these transpositions generate the full symmetric
group on the factor. The experiment imposes the zero-sum condition for
every pair in that factor. When this check applies at every bad prime,
it gives the exact kernel without constructing a global splitting field.

The split marked double-cover construction applies in odd residue
characteristic. At 2 the model constructor handles nodal equations,
smooth disc charts it can verify, and the verified compact-type equations.
Some other wild fibres have exact projected constraints for particular
specializations. When the full boundary is unavailable,
``boundary()`` raises ``NotImplementedError``. The integral-basis
solver uses any available sufficient component test and sets
``complete=False`` when the necessary and sufficient kernels differ.

The odd-characteristic construction follows
`Gehrunger and Pink <https://arxiv.org/abs/2112.05550>`_.
The finite-place integrality test uses
`Scholl's comparison of integral motivic cohomology and unramified classes <https://arxiv.org/abs/0710.5453>`_.
