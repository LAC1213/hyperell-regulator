// SPDX-License-Identifier: GPL-2.0-or-later
// Harvey, arXiv:1210.8239, Sections 3--5, with an even-degree extension.
#ifndef HYPELLFROB_AVERAGE_H
#define HYPELLFROB_AVERAGE_H

#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <vector>

namespace hypellfrob {
struct AverageJob {
    long p;
    int precision;
    bool ok = false;
    NTL::mat_ZZ matrix;
};

// Unhandled jobs retain ok=false. No single-prime fallback is hidden here.
// Threads distribute independent (a,b) reduction trees, not individual primes.
void matrices_average(std::vector<AverageJob>& jobs, const NTL::ZZX& Q,
                      int nthreads = 1);
// Upstream rforest C library with bounded forests, for both parities.
// Even models use paths avoiding the horizontal residue singularity and
// return the same residue-free projective basis as matrix_reference().
// kappa=-1 selects the automatic forest size; nonnegative values override it.
void matrices_rforest(std::vector<AverageJob>& jobs, const NTL::ZZX& Q,
                      int nthreads = 1, int kappa = -1);
}
#endif
