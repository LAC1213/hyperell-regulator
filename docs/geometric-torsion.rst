Geometric torsion packets
============================

``geometric_torsion_packet`` computes the complete packet over
:math:`\overline{\mathbb Q}` containing the Weierstrass points of a smooth
genus-two curve over :math:`\mathbb Q`. It runs the author's PARI/GP backend
for `Computing torsion points on curves
<https://math.mit.edu/~poonen/papers/curvetorsion.pdf>`_, including the
upstream exact verification and completeness checks. It accepts success
only after the backend's successful completion and validation of the
packet factors and point count. A timeout raises an error, and never
returns an empty or purportedly complete partial packet.

The external `torsion.gp <https://math.mit.edu/~poonen/papers/torsion.gp>`_
is downloaded separately, checked against a pinned SHA256, and kept in
the local results cache. It is not redistributed with this package.
Install it once, then run in Sage::

    from sage.all import QQ, PolynomialRing
    from hyperell_regulator import (
        HyperellipticCurveWithReduction, install_poonen_backend,
        geometric_torsion_packet,
    )

    install_poonen_backend()   # downloads the pinned upstream GP source
    x = PolynomialRing(QQ, 'x').gen()
    C = HyperellipticCurveWithReduction(
        2+3*x+x*x+x**3-x**5, 1+x**3)
    packet = geometric_torsion_packet(C, basepoint=(-1, 0), timeout=300)
    assert packet.complete and packet.size == 10
    for P in packet.points():
        print(P, packet.order(P))

``backend`` can specify an existing GP file, and ``gp`` can specify the
GP executable. Alternatively set ``HYPERELL_POONEN_GP``. An already
downloaded backend can be installed without network access using
``install_poonen_backend(source_path=...)``. Successful computations and
their full transcripts are cached under ``results/poonen/`` within
``HYPERELL_REGULATOR_DATA``. Changing the base point reuses the computed
packet, while checking membership of the new base point exactly.

Supported base points
-------------------------

All Weierstrass points belong to one packet. The method also accepts any
non-Weierstrass point in that packet, including algebraic base points.
The complete packet is unchanged when the base point moves within it.
If a requested base point lies outside this packet, the method raises
``NotImplementedError``. General packets attached to non-Weierstrass base
points outside it are not implemented; the hyperelliptic packet must not
be returned as the answer for such a base point.

``points()`` enumerates exact :math:`\mathrm{QQbar}` coordinate pairs on
the original :math:`y^2+h(x)y=f(x)` model. Infinity is represented by
``('infinity', eta)`` with :math:`\eta=\lim y/x^3`. Polynomial factors in
``packet.factors`` compactly represent all x-coordinates and retain the
verified orders of :math:`2P-\kappa`, with :math:`\kappa` the canonical
degree-two divisor class. For a point :math:`P` and selected base :math:`O`,
twice the lcm of their doubled orders annihilates :math:`P-O`.
``order(P)`` uses this proved multiple and exact Jacobian arithmetic over
a number field to determine the precise order. It does not apply bounds
for rational torsion to algebraic points. A field-degree limit guards
large composita and reports unsupported computations explicitly.

New classes over the algebraic closure
-----------------------------------------

For 249.a.6723.1 with base :math:`O=(-1,0)`, the complete packet has
10 points, including the six previously used rational points. The four
new points are

.. math::

   P_a=\left(a,-\frac{1+a^3}{2}\right),\qquad
   a^4-2a^3-a^2-2a-3=0.

Each :math:`P_a-O` has exact order 2, and the explicit function

.. math::

   \varphi_a=\frac{x-a}{x+1}

satisfies :math:`\operatorname{div}(\varphi_a)=2(P_a-O)`.
These define additional classes :math:`B(P_a-O)` over number fields.
The numerical experiment compares their regulator moments with those
of previously certified rational pair classes, modulo the period vector.

Noncertified Stokes integration now accepts exact real or complex
algebraic logarithmic centers and poles. Algebraic centers that coincide
with branch points are identified exactly before numerical integration.
Certified integration also accepts exact algebraic centers and poles.

Run the geometric experiment::

    sage -python -m hyperell_regulator.experiments.geometric_torsion \
        --label 249.a.6723.1 --label 256.a.512.1 \
        --digits 25 --out results/geometric_torsion.json

The experiment output contains exact number-field
coordinates, specified embeddings, exact orders, backend fingerprints,
and numerical regulator vectors. Regulator comparisons at the chosen
complex embeddings give numerical evidence of new directions. They do
not certify nonvanishing by interval arithmetic, integrality over number
fields, or descent to :math:`\mathbb Q`.
