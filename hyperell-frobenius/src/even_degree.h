// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef HYPELLFROB_EVEN_DEGREE_H
#define HYPELLFROB_EVEN_DEGREE_H

#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>

namespace hypellfrob {
// For a monic Q of degree 2g+2, restrict affine Frobenius in the basis
// x^i dx/y (0 <= i <= 2g) to the kernel of the residue at infinity.
// If rho_i = [u^(i-g)] (u^deg(Q) Q(1/u))^(-1/2), the projective basis is
// e_i for i < g and e_i - rho_i e_g for i > g.
// Verify rho F = p rho and return a 2g-by-2g matrix modulo p^N.
void restrict_even_frobenius(NTL::mat_ZZ& output, const NTL::mat_ZZ& affine,
                            const NTL::ZZ& p, int N, const NTL::ZZX& Q);
}
#endif
