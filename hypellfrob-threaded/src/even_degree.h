/* ============================================================================

   even_degree.h

   Frobenius matrices for hyperelliptic curves given by an EVEN degree model
   y^2 = Q(x), deg Q = 2g+2, which hypellfrob::matrix() does not handle.

   This is NOT David Harvey's code and is not his algorithm's O(sqrt p)
   method. It is a direct, straightforward implementation of Kedlaya's
   algorithm whose cost is linear in p, written to (a) support even degree at
   all and (b) serve as an independent reference the fast path can be checked
   against. It shares hypellfrob's conventions and its GPL v2-or-later
   licence; see COPYING and README.md.

   Kedlaya, "Counting points on hyperelliptic curves using Monsky-Washnitzer
   cohomology", J. Ramanujan Math. Soc. 16 (2001).

============================================================================ */

#ifndef HYPELLFROB_EVEN_DEGREE_H
#define HYPELLFROB_EVEN_DEGREE_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>

namespace hypellfrob {

/*
Frobenius acting on H^1 of the AFFINE curve y^2 = Q(x) with the two points
at infinity removed, in the basis

    x^i dx/y,   0 <= i <= deg(Q) - 2.

Works for either parity of deg(Q):

  deg Q = 2g+1 : dimension 2g. One point at infinity, so this already IS
                 H^1 of the smooth projective curve, and the result should
                 agree with hypellfrob::matrix() exactly. That agreement is
                 what validates this implementation (see check_even).
  deg Q = 2g+2 : dimension 2g+1. Two points at infinity; H^1 of the affine
                 curve carries one extra class, the one with nonvanishing
                 residue there.

`output` is reduced mod p^N.

extra_precision is how many p-adic digits beyond N to carry internally;
reductions divide by integers that can be divisible by p, and each such
division costs precision. Raise it if reduce_precision_failure() reports
trouble. The default is generous.

PRECONDITIONS:
   p prime (not checked), p > (2N-1)*deg(Q), N >= 1.
   Q monic of degree >= 3, and Q mod p squarefree.

RETURN VALUE: 1 on success, 0 if a precondition fails.
*/
int matrix_affine_reference(NTL::mat_ZZ& output, const NTL::ZZ& p, int N,
                            const NTL::ZZX& Q, int extra_precision = 12);

/*
Frobenius on H^1 of the SMOOTH PROJECTIVE curve: a 2g x 2g matrix mod p^N,
for either parity of deg(Q).

For odd deg(Q) this is matrix_affine_reference() unchanged.

For even deg(Q) it is the restriction of that matrix to the subspace of
differentials with vanishing residue at infinity. The residue map is
Frobenius-equivariant with a factor of p, so that subspace is invariant and

    charpoly(Frob | H^1_affine) = (T - p) * charpoly(Frob | H^1_curve).

The residue of x^i dx/y at infinity is the coefficient of u^(i-g) in
P(u)^(-1/2), where P(u) = u^deg(Q) Q(1/u); it vanishes for i < g and is 1
at i = g. So the residue-free basis is

    x^i dx/y                        for i < g,
    x^i dx/y - rho_i * x^g dx/y     for i > g,

and a residue-free vector's coordinates in it are just its coordinates in
the monomial basis with the g-th dropped.
*/
int matrix_reference(NTL::mat_ZZ& output, const NTL::ZZ& p, int N,
                     const NTL::ZZX& Q, int extra_precision = 12);

} // namespace hypellfrob

#endif
