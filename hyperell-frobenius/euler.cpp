// Frobenius matrices for one curve at many primes, in parallel.
//
// Usage: euler [nthreads] [rforest|rforest-only] [kappa]
// rforest (default): simultaneous algorithm using the vendored C library.
// rforest-only is an alias: unsupported jobs are always declined.
//
// Reads a job from stdin:
//   line 1:  the coefficients of Q, ascending, whitespace separated
//   then:    one line per job, "p N"
//   (blank lines and lines beginning with '#' are ignored)
//
// Writes one line per job, in the input order:
//   "<p> 1 <dim> <entries...>"   row-major, the 2g x 2g matrix mod p^N
//   "<p> 0"                      declined: bad reduction, or a
//                                precondition of the fast path failed
//
#include "average.h"
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <cstdlib>

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
    std::string backend = argc > 2 ? argv[2] : "rforest";
    if (backend != "rforest" && backend != "rforest-only") {
        std::cerr << "euler: backend must be rforest or rforest-only\n";
        return 2;
    }
    int kappa = -1;
    if (argc > 3) {
        char* end;
        long value = std::strtol(argv[3], &end, 10);
        if (end == argv[3] || *end || value < -1 || value > 30) {
            std::cerr << "euler: rforest kappa must be -1..30\n";
            return 2;
        }
        kappa = value;
    }

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

    std::vector<AverageJob> average;
    for (const auto& j : jobs) average.push_back({j.p, j.N});
    try {
        matrices_rforest(average, Q, nthreads, kappa);
    } catch (const std::exception& e) {
        std::cerr << "euler: rforest failed: " << e.what() << '\n';
        return 1;
    }
    for (size_t k = 0; k < jobs.size(); ++k)
        if (average[k].ok) { jobs[k].ok = 1; jobs[k].M = std::move(average[k].matrix); }

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
