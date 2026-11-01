// SPDX-License-Identifier: GPL-2.0-or-later
// Restrict even-degree affine Frobenius to projective cohomology.
// This is the residue projection used by the simultaneous rforest algorithm.
#include "even_degree.h"
#include <NTL/ZZ_pX.h>
#include <stdexcept>
#include <vector>

NTL_CLIENT
namespace hypellfrob {
namespace {
// Residue covector: rho_i = [u^(i-g)] P(u)^(-1/2), P(u) = u^d Q(1/u).
std::vector<ZZ_p> residue_covector(const ZZ_pX& Q, int g, int d)
{
    ZZ_pX P;
    for (long k = 0; k <= d; k++) SetCoeff(P, d - k, coeff(Q, k));
    ZZ_pX A = InvTrunc(P, g + 1);
    std::vector<ZZ_p> c(g + 1);
    c[0] = to_ZZ_p(1);
    ZZ_p half = to_ZZ_p(1) / to_ZZ_p(2);
    for (int m = 1; m <= g; m++) {
        ZZ_p acc = coeff(A, m);
        for (int a = 1; a < m; a++) acc -= c[a] * c[m - a];
        c[m] = acc * half;
    }
    std::vector<ZZ_p> rho(2 * g + 1);
    for (int i = 0; i <= 2 * g; i++) rho[i] = (i >= g) ? c[i - g] : to_ZZ_p(0);
    return rho;
}


} // anonymous namespace

void restrict_even_frobenius(mat_ZZ& output, const mat_ZZ& M,
                            const ZZ& p, int N, const ZZX& Q)
{
    const int d = deg(Q);
    if (d%2 || d < 4 || M.NumRows() != d-1 || M.NumCols() != d-1)
        throw std::invalid_argument("invalid even-degree affine Frobenius matrix");

    int g = (d - 2) / 2;

    ZZ_pContext caller;  caller.save();
    ZZ_p::init(power(p, N));
    ZZ_pX Qp;  conv(Qp, Q);
    std::vector<ZZ_p> rho = residue_covector(Qp, g, d);

    // Restrict to ker(rho). With rho_g = 1, a basis is e_i (i < g) and
    // e_i - rho_i e_g (i > g); a residue-free vector's coordinates in it are
    // its monomial coordinates with the g-th dropped.
    output.SetDims(2 * g, 2 * g);
    ZZ mod = power(p, N);
    for (int j = 0; j <= 2*g; ++j) {
        ZZ_p value;
        for (int i = 0; i <= 2*g; ++i) value += rho[i]*to_ZZ_p(M[i][j]);
        if (value != to_ZZ_p(p)*rho[j]) {
            caller.restore();
            throw std::runtime_error("even-degree Frobenius residue invariant failed");
        }
    }
    for (int jj = 0, jout = 0; jj <= 2 * g; jj++) {
        if (jj == g) continue;
        for (int ff = 0, fout = 0; ff <= 2 * g; ff++) {
            if (ff == g) continue;
            ZZ_p val = to_ZZ_p(M[ff][jj]) - rho[jj] * to_ZZ_p(M[ff][g]);
            output[fout][jout] = rep(val) % mod;
            fout++;
        }
        jout++;
    }

    caller.restore();
}

} // namespace hypellfrob
