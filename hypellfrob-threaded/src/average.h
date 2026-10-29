// SPDX-License-Identifier: GPL-2.0-or-later
// Harvey, arXiv:1210.8239, Sections 3--5 (odd monic integral models).
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
}
#endif
