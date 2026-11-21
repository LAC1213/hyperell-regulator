Manual bad factors for the CM-labelled Beilinson survey
=======================================================

The command is::

   sage -python -m hyperell_regulator.experiments.beilinson 10000 \
       --digits 20 --reg-digits 20 --cm-only --ignore-component-criterion

The initial run selected 19 curves. Eleven had admissible torsion triples;
seven of those were missing at least one standard Euler factor. The other
eight curves had no admissible triple, so their L-values were not requested.

Here ``CM`` is the LMFDB label for the rational endomorphism algebra. It
does not assert that every Jacobian in this selection is of CM type. Six
of the seven missing cases come from non-CM newforms with an inner twist.

The module ``hyperell_regulator.experiments.beilinson_local_data`` records
the following standard Euler polynomials, their local conductor exponents,
the original equations and the corresponding newform labels. Coefficients
in the Python table are ascending. The factors at other bad primes are
still computed by the package.

================== ===== ================================= ======== ============
Curve              Prime Standard Euler polynomial         Exponent Full N_std
================== ===== ================================= ======== ============
324.a.648.1        3     1-3T                              8        26244
784.c.614656.1     2     (1+2T)(1+2T+4T^2)                 2        9604
1296.a.20736.1     2     (1+2T)(1+2T+4T^2)                 2        26244
1296.a.20736.1     3     1-3T                              8        26244
2704.a.43264.1     2     (1+2T)(1-2T+4T^2)                 2        114244
3969.c.35721.1     3     1-3T                              8        15752961
4096.b.65536.1     2     1+2T                              12       4096
8649.c.700569.1    3     1-3T                              4        74805201
================== ===== ================================= ======== ============

The existing ``bad_factors`` argument to ``StandardLFunction`` takes the
wedge-square polynomial, so the experiment multiplies each entry above
by (1-pT). The new ``conductor_exponents`` argument takes standard local
Artin conductor exponents, including the Swan contribution. Supplying the
polynomial without this conductor correction would give a wrong functional
equation in the wild cases.

Local representation calculation
--------------------------------

The corresponding weight-two newforms are, respectively,
``18.2.c.a``, ``28.2.e.a``, ``36.2.e.a``, ``52.2.h.a``, ``63.2.i.a``,
``64.2.b.a`` and ``93.2.g.b``. These associations are listed on the
`LMFDB genus-two curve pages <https://www.lmfdb.org/Genus2Curve/Q/>`_.
For example, the `324 curve page
<https://www.lmfdb.org/Genus2Curve/Q/324/a/648/1>`_ links to
`18.2.c.a <https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/18/2/c/a/>`_.

For the six non-CM newforms let rho denote the two-dimensional
representation, and epsilon its nebentypus. Conjugation of the coefficient
field gives rho_bar = rho_dual tensor Q_l(-1), so

.. math::

   V_{\rm std} =
   (\epsilon\oplus\epsilon^{-1}\oplus\operatorname{Ad}^0\rho)
      \otimes\mathbb Q_\ell(-1).

At 3 the relevant local representations are principal series with one
unramified character and one of conductor exponent r. There is one
unramified adjoint line, giving 1-3T; the other four standard characters
have conductor r. Hence a_3(std)=4r. Sage's exact ``LocalComponent``
calculation gives r=2 for levels 18, 36 and 63, and r=1 for level 93.
This determines both the factor and the conductor, without fitting them
to the functional equation.

At 2 the forms of levels 28, 36 and 52 have tame supercuspidal local
representations induced from the unramified quadratic extension. Inertia
has order three. The unramified line in the adjoint representation is the
unramified quadratic character; its Frobenius eigenvalue after the Tate
twist is -2. The other two adjoint dimensions have no invariants and
contribute conductor 2. The nebentypus is unramified at 2, giving

.. math::

   P^{\rm std}_2(T)=(1+2T)
       \left(1-2(\epsilon(2)+\epsilon(2)^{-1})T+4T^2\right).

The nebentypus value has order three for levels 28 and 36, and order six
for level 52. This gives the two displayed cubic polynomials and local
conductor exponent 2.

For the Bolza curve (4096.b), the CM decomposition used in
``experiments/bolza_triple.py`` is

.. math::

   L({\rm std},s)=L(f_8,s)
     L(\chi_8,s-1)^2 L(\chi_{-8},s-1),

where f_8 is the weight-three CM newform of level 8, with a_2=-2.
All three quadratic characters have conductor 8. Thus P_2(std)=1+2T,
and the conductor is 8 * 8^2 * 8 = 4096.

Reproduction and checks
-----------------------

The principal-series and supercuspidal data can be obtained with Sage's
`local component implementation
<https://doc.sagemath.org/html/en/reference/modfrm/sage/modular/local_comp/local_comp.html>`_::

   from sage.all import *
   from sage.modular.local_comp.local_comp import LocalComponent
   K = CyclotomicField(6)
   chi = next(c for c in DirichletGroup(18, K)
              if c.is_even() and c.conductor() == 9 and c.order() == 3)
   f = CuspForms(chi, 2).newforms('a')[0]
   local = LocalComponent(f, 3)
   print(local.species(), local.conductor(), list(local.characters()))

This yields principal series of conductor 2 with character conductors
(0,2). For levels 28, 36, 52 choose the even nebentypus of conductor
7, 9, 13 and order 3, 3, 6 respectively, and examine the component at 2.
For levels 63 and 93 choose the matching degree-two newform among the
even characters of conductor 63 and 93 and order 6; its local component
at 3 has character conductors (0,2) and (0,1) respectively. Conjugate
newforms give the same standard data. Curve/newform matching was checked
against the h^1 traces at good primes 2, 5 and 11 for these last two cases.

Independent numerical functional equation checks at 12 requested digits
gave residuals below 4.3e-22 for all seven completed L-functions. These
checks validate the implementation numerically; they are not substitutes
for the local representation calculation. The requested command is a
numerical computation, since it does not include ``--certified``.

The completed 20-digit rerun is recorded in
``reports/beilinson_cm_bad_factors.json``. All seven missing-factor cases
succeeded, with functional equation residuals below 2e-36. All eleven
curves having admissible triples have L-values; the other eight have no
admissible triple. The largest calculation used 1,009,887 coefficients.
Frobenius computations used eight workers and BLAS/OpenMP used one thread.

Both numerical and certified experiment paths use these overrides. Saved
L-values and reports include the complete manual data and its fingerprint
in ``L_local_data``. Cache entries without matching data are recomputed.
The helpers refuse to apply a label's factors to a different equation.
