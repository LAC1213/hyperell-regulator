// Frobenius matrices for one curve at many primes, in parallel.
//
// Usage: euler [nthreads]
//
// Reads a job from stdin:
//   line 1:  the coefficients of Q, ascending, whitespace separated
//   then:    one line per job, "p N"
//   (blank lines and lines beginning with '#' are ignored)
//
// Writes one line per job, in the input order:
//   "<p> 1 <dim> <entries...>"   row-major, the 2g x 2g matrix mod p^N
//   "<p> 0"                      matrix() declined: bad reduction, or a
//                                precondition of the fast path failed
//
// Odd degree goes through hypellfrob::matrix(), which is O(sqrt p); even
// degree through matrix_reference(), which is correct for either parity but
// linear in p, so very much slower. The caller only sends an even degree
// model when the curve has no rational Weierstrass point -- with one, it
// moves that point to infinity first and sends the odd model instead.
//
// The point of the batching is that the caller pays NTL's setup and this
// process's startup once for a whole L-function's worth of primes, and the
// primes run concurrently.
#include "hypellfrob.h"
#include "even_degree.h"
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>

NTL_CLIENT
using namespace hypellfrob;

struct Job {
    long p;
    int N;
    int ok = 0;
    mat_ZZ M;
};

int main(int argc, char** argv) {
    int nthreads = argc > 1 ? atoi(argv[1]) : (int)std::thread::hardware_concurrency();
    if (nthreads < 1) nthreads = 1;

    std::string line;
    ZZX Q;
    int deg = -1;
    while (std::getline(std::cin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        std::string tok;
        long i = 0;
        while (is >> tok) {
            SetCoeff(Q, i, to_ZZ(tok.c_str()));
            i++;
        }
        deg = (int)i - 1;
        break;
    }
    if (deg < 3) {
        std::cerr << "euler: expected the coefficients of Q on the first line"
                  << std::endl;
        return 2;
    }

    std::vector<Job> jobs;
    while (std::getline(std::cin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        Job j;
        if (!(is >> j.p)) continue;
        if (!(is >> j.N)) j.N = 1;
        jobs.push_back(j);
    }

    std::atomic<size_t> next_index(0);
    auto worker = [&]() {
        for (;;) {
            size_t k = next_index.fetch_add(1);
            if (k >= jobs.size()) break;
            Job& j = jobs[k];
            ZZ pz = to_ZZ(j.p);
            try {
                j.ok = (deg % 2 == 1) ? matrix(j.M, pz, j.N, Q)
                                      : matrix_reference(j.M, pz, j.N, Q);
            } catch (...) {
                j.ok = 0;
            }
        }
    };

    std::vector<std::thread> pool;
    for (int t = 0; t < nthreads; t++) pool.emplace_back(worker);
    for (auto& t : pool) t.join();

    for (const Job& j : jobs) {
        std::cout << j.p << ' ' << j.ok;
        if (j.ok) {
            std::cout << ' ' << j.M.NumRows();
            for (long r = 0; r < j.M.NumRows(); r++)
                for (long c = 0; c < j.M.NumCols(); c++)
                    std::cout << ' ' << j.M[r][c];
        }
        std::cout << '\n';
    }
    return 0;
}
