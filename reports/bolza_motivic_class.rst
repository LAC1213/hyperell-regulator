Integral motivic classes for the Bolza curve
===========================================

Let C/Q be y^2=x^5-x, let K=Q(sqrt(2)), put s=sqrt(2), and put
epsilon=1+s. Define two involutions of C_K by

.. math::

   t_1(x,y)=\left(\frac{x+1}{x-1},\frac{2s\,y}{(x-1)^3}\right),
   \qquad
   t_2(x,y)=\left(\frac{1-x}{x+1},\frac{2s\,y}{(x+1)^3}\right).

Their graphs give two classes over Q:

.. math::

   \Xi_j=\operatorname{Cor}_{K/\mathbb Q}
         (\Gamma_{t_j},\epsilon)
       \in\mathrm{CH}^2(C\times C,1)
       =H^3_{\mathcal M}(C\times C,\mathbb Z(2)).

Here (Gamma,epsilon) is the higher Chow symbol consisting of the curve
Gamma and its constant rational function epsilon. Corestriction is
proper pushforward under (C x C)_K -> C x C, so this is an actual class
over Q, not just a class over K. It is decomposable over K; no claim of
indecomposability is made.

Cycle condition and finite-prime integrality
--------------------------------------------

The divisor of epsilon on each graph is zero. Furthermore

.. math::

   (1+\sqrt2)(1-\sqrt2)=-1,
   \qquad (1+\sqrt2)^{-1}=\sqrt2-1,

so epsilon is an algebraic integer and a unit in O_K=Z[sqrt(2)].
It has valuation zero at every finite prime, including the prime above 2.
On the normalization of the graph closure in an arithmetic model it
remains a unit. Thus it has no horizontal or vertical divisor, even
after resolving the model. This gives an integral class over O_K;
corestriction preserves this property and gives an integral class over Z.
In particular this argument does not rely on component tests at 2.

Nonzero primitive regulators
----------------------------

Let sigma send s to -s and let iota be the hyperelliptic involution.
The formulas give t_j^sigma=iota composed with t_j, while
abs(epsilon^sigma)=epsilon^{-1} in the real embedding with s>0.
Consequently the real regulator of Xi_j is, up to the fixed overall
regulator normalization,

.. math::

   \log(1+\sqrt2)
   \bigl([\Gamma_{t_j}]-[\Gamma_{\iota t_j}]\bigr).

On the H^1 tensor H^1 part, iota acts by -1. The difference therefore
acts by twice the correspondence for t_j. On the holomorphic
differentials (dx/y, x dx/y), the pullback matrices are

.. math::

   T_1=\frac1{\sqrt2}\begin{pmatrix}1&-1\\-1&-1\end{pmatrix},
   \qquad
   T_2=\frac1{\sqrt2}\begin{pmatrix}-1&-1\\-1&1\end{pmatrix}.

These are independent nonzero trace-zero matrices. The involutions
preserve the canonical polarization, so their Rosati adjoints are their
inverses, hence themselves. The correspondences are therefore in the
symmetric, trace-free part, giving nonzero components in the degree-five
standard motive. Since log(1+sqrt(2)) is nonzero, the real regulators of
Xi_1 and Xi_2 are nonzero and independent. Both motivic classes are
therefore nontorsion.

Reproduction
-------------

The following verifies the automorphism identities, the unit, the
Galois conjugation, the matrices and their independence with exact Sage
arithmetic. It also encloses log(epsilon) with Arb. No Frobenius batches
or native L-value computations are needed.

.. code-block:: bash

   conda activate sage
   sage -python -m hyperell_regulator.experiments.bolza_motivic

Regulator matrix in the package normalization
---------------------------------------------

The example now also computes the regulator matrix with certified package
periods. It uses the same completed-square model as the Beilinson
experiment, Y^2=4(x^5-x), with differentials dx/Y and x dx/Y.

If G is the unweighted Gram matrix returned by ``periods_stokes``, the
graph regulator pairing matrix for Xi_j is

.. math::

   H_j=2\log(\varepsilon)G T_j^{\mathsf t}.

The package convention is (i/4) times the integral over the curve; on
the diagonal this equals its complex-plane integral. The factor 2 in
H_j comes from corestriction. A pairing matrix H contributes the column
(H_00,H_01+H_10,H_11), matching the weights
(1,2 Re(x),|x|^2).

For this curve G=b*I by symmetry. Put

.. math::

   b=2.24329395360277215016538560614442858646\ldots,
   \qquad k=\sqrt2\,b\log(1+\sqrt2)
     =2.79615482587842960151911256909178058999\ldots.

Thus the two regulator columns and the period column are

.. math::

   M=[r(\Xi_1),r(\Xi_2),P]
    =\begin{pmatrix}k&-k&b\\-2k&-2k&0\\-k&k&b\end{pmatrix},
   \qquad \det M=-8bk^2.

At 20 decimal places this is

.. math::

   \begin{pmatrix}
   2.79615482587842960152&-2.79615482587842960152&2.24329395360277215017\\
   -5.59230965175685920304&-5.59230965175685920304&0\\
   -2.79615482587842960152&2.79615482587842960152&2.24329395360277215017
   \end{pmatrix},

with determinant -140.313223770893026630766112534787318219698498253...
The saved determinant ball has radius less than 7e-46. Exact rational
endpoints are in ``reports/bolza_motivic_regulator_20.json``.

The matrix was also checked against independently integrated logarithmic
columns r_0, r_1 for log|x| and log|x-1|. In the same row convention, the
enclosures overlap 4(r_1-r_0)-2*log(2)*P and
4*r_1-2*log(2)*P, respectively.

To reproduce the matrix and save the enclosures:

.. code-block:: bash

   sage -python -m hyperell_regulator.experiments.bolza_motivic \
       --digits 20 --json reports/bolza_motivic_regulator_20.json

For the original differential basis dx/y, x dx/y on y^2=x^5-x, every
matrix entry is multiplied by 4, and the augmented determinant by 64.

Observed standard L-value relation
-----------------------------------

In this normalization and column order, the computation gives

.. math::

   \frac{\pi c^+c^-L''(\mathrm{std},1)}{\det M}
       \simeq\frac1{12},
   \qquad \det M\simeq12\pi c^+c^-L''(\mathrm{std},1).

The relevant values are

.. math::

   L''(\mathrm{std},1)=0.36979847590703155470833631435221574482\ldots,
   \qquad c^+=-c^-=3.17249673357460100265574838315919522805\ldots.

The certified comparison at 40 requested digits is saved in
``reports/bolza_motivic_lvalue_40.json``. It records the ratio enclosure,
relative agreement bound, signed periods, L-value, matrix, exact ball
endpoints and backend metadata. As elsewhere, L-value bounds are conditional
on the stated analytic assumptions. The record does not claim a proof of
exact rational equality.

Reproduce using two Frobenius workers:

.. code-block:: bash

   HYPELLFROB_THREADS=2 sage -python -m hyperell_regulator.experiments.bolza_motivic \
       --digits 40 --with-lvalue --json reports/bolza_motivic_lvalue_40.json

Torsion-point regulator representatives
---------------------------------------

Put W_a=(a,0), for a=0,1, and O=infinity. The functions x and x-1
have divisors 2W_0-2O and 2W_1-2O. They define genuine Bloch cycles

.. math::

   B_a=(\Delta_C,x-a)-(C\times W_a,x-a)-(O\times C,x-a).

Define also the constant, decomposable cycle

.. math::

   E_2=(\Delta_C,2)-(C\times O,2)-(O\times C,2).

Its mixed-form regulator column is log(2)*P. Consequently the classes

.. math::

   \widetilde\Xi_1=4(B_1-B_0)-2E_2,
   \qquad \widetilde\Xi_2=4B_1-2E_2

have the logarithmic regulator columns checked numerically above.
These checks establish overlapping certified regulator enclosures;
they are not a proof that Xi_j equals tilde-Xi_j in a higher Chow group.
Such a motivic identity would require an explicit higher Chow relation
or a suitable injectivity argument. Integrality of the graph-and-unit
classes must likewise not be transferred to these representatives merely
from agreement of their regulators.

Component separation at 2 and the determinant
---------------------------------------------

The rational L-value agreement does not prove integrality of B_0 and B_1.
The same-component criterion is a sufficient way to prove integrality;
failing it does not by itself compute a nonzero boundary class. The
repository tests smooth stable components at odd primes and currently
returns an undecided verdict at 2. These tests must not be presented as
a computation of the components in the minimal regular model at 2.

Arithmetic integrality is determined by the finite-prime localization
boundary. For the standard motive one must also take the standard
projector into account; a boundary of a representative in the full
cohomology of C x C need not survive that projection. See Scholl,
``Integral elements of K-theory and products of modular curves II``,
section 1: https://arxiv.org/pdf/0710.5453 .

There is a concrete reason why this determinant can miss a correction
at 2. The decomposable class E_2 above has regulator log(2)*P and can
have a finite-prime boundary at 2, since 2 is not a 2-adic unit. But

.. math::

   \det(r_0,r_1+a\log(2)P,P)=\det(r_0,r_1,P).

In particular, replacing r_1 by r_1-log(2)*P/2 leaves the determinant
unchanged. The regulator columns of Xi_1 and Xi_2 match the checked
combinations 4*(r_1-log(2)*P/2-r_0) and
4*(r_1-log(2)*P/2). Their determinant is -16 times the original
torsion-column determinant. This linear-algebra identity explains why
the observed rational L-value relation is insensitive to such period
corrections. It is not a computation of the 2-adic boundary of B_0 or B_1,
nor does it prove equality with the integral graph-and-unit classes.

The exact 2-torsion data behind the graph automorphisms are the partitions

.. math::

   t_1:\quad \{W_0,W_{-1}\},\{W_1,O\},\{W_i,W_{-i}\},
   \qquad
   t_2:\quad \{W_0,W_1\},\{W_{-1},O\},\{W_i,W_{-i}\}.

Each involution interchanges the two points in each pair. The corresponding
subgroups of J[2] are generated by [W_0-W_-1], [W_1-O] and by
[W_0-W_1], [W_-1-O], respectively. Their common third nonzero element
is [W_i-W_-i], represented over Q by the Mumford pair (x^2+1,0).
All these order-two and sum relations were checked exactly in Sage.

Relation to the earlier triple
------------------------------

There is also the Bloch cycle attached just to P=(0,0), O=infinity,
and f=x, since div(x)=2P-2O:

.. math::

   Z=(\Delta_C,x)+(C\times\{P\},1/x)+(\{O\}\times C,1/x).

Its three generic-fibre divisors cancel. This construction is described
in section 4.1 of Sarkar--Sreekantan, ``Extensions of Motives and the
Fundamental Group``:
https://www.isibang.ac.in/~statmath/eprints/2014/8.pdf .

A third point used only to normalize f need not be torsion. In particular,
(2,sqrt(30)) can normalize this single cycle by replacing x with x/2.
It cannot supply the second torsion pair needed for the earlier
two-column determinant. The unit classes above avoid that problem and
have the direct finite-prime integrality proof just given.
