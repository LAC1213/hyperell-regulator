// SPDX-License-Identifier: GPL-2.0-or-later
// An implementation of David Harvey, "Counting points on hyperelliptic
// curves in average polynomial time", arXiv:1210.8239v3, Sections 3--5.
// Original implementation for hyperell-regulator, 2026.
// Even-degree extension: shifted horizontal/diagonal paths avoid the
// residue singularity; restrict the affine answer to the residue kernel.
//
// All matrices act on columns in the basis of W_{s,t}. In particular the
// accumulating products are M_0 M_1 ... M_{(p-1)/2}, in that order.
#include "average.h"
#include "even_degree.h"
#include <gmp.h>
#include "../rforest/rforest.h"
#include <NTL/mat_ZZ.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>

extern "C" {
void hw_mpz_setup(void);
void hw_mpz_clear(void);
extern thread_local int hw_disable_fft;
extern thread_local int rforest_forward_differences;
}

NTL_CLIENT
namespace hypellfrob {
namespace {

struct Fraction {
    mat_ZZ m;
    ZZ d = ZZ(1);
};

// Cancel only exact integer contents, before constructing the modular tree.
void cancel(Fraction& f) {
    ZZ h = abs(f.d);
    for (long i = 0; i < f.m.NumRows() && h != 1; ++i)
        for (long j = 0; j < f.m.NumCols() && h != 1; ++j)
            GCD(h, h, f.m[i][j]);
    if (f.d < 0) h = -h;
    if (h != 1) {
        f.d /= h;
        for (long i = 0; i < f.m.NumRows(); ++i)
            for (long j = 0; j < f.m.NumCols(); ++j) f.m[i][j] /= h;
    }
}

struct Reductions {
    long d;
    ZZX q;
    ZZ delta;
    mat_ZZ R, S;

    explicit Reductions(const ZZX& Q) : d(deg(Q)), q(Q) {
        // Lemma 3: adjugate of (R,S) -> R Q + S Q', with deg R<d-1,
        // deg S<d. Its first d columns give delta*x^i=R_i Q+S_i Q'.
        mat_ZZ syl, adj;
        syl.SetDims(2*d-1, 2*d-1);
        for (long j = 0; j < d-1; ++j)
            for (long i = 0; i <= d; ++i) syl[i+j][j] = coeff(q, i);
        for (long j = 0; j < d; ++j)
            for (long i = 0; i < d; ++i)
                syl[i+j][d-1+j] = (i+1)*coeff(q, i+1);
        inv(delta, adj, syl, 1);
        if (IsZero(delta)) throw std::invalid_argument("Q is not squarefree");
        R.SetDims(d, d); S.SetDims(d, d);
        for (long col = 0; col < d; ++col) {
            for (long row = 0; row < d-1; ++row) R[row][col] = adj[row][col];
            for (long row = 0; row < d; ++row) S[row][col] = adj[d-1+row][col];
        }
    }

    Fraction step(long s, long t, char kind) const {
        Fraction f;
        f.m.SetDims(d, d);
        if (kind == 'H') { // Lemma 5
            f.d = to_ZZ(d)*(2*t-1)-2*s;
            for (long i = 1; i < d; ++i) f.m[i][i-1] = f.d;
            for (long i = 0; i < d; ++i)
                f.m[i][d-1] = (to_ZZ(2*s)-to_ZZ(2*t-1)*i)*coeff(q, i);
        } else if (kind == 'D') { // Lemma 6
            f.d = delta*(2*t-1);
            for (long i = 0; i < d; ++i)
                for (long j = 0; j < d; ++j) {
                    f.m[i][j] = 2*(to_ZZ(s)+i)*S[i][j];
                    if (i) f.m[i][j] += (2*t-1)*R[i-1][j];
                }
        } else { // Lemma 7; only permitted when c_0 != 0
            const ZZ& c = coeff(q, 0);
            if (IsZero(c)) throw std::logic_error("vertical reduction with c0=0");
            f.d = delta*c*(2*t-1);
            for (long i = 0; i < d; ++i)
                for (long j = 0; j < d; ++j) {
                    f.m[i][j] = S[0][j]*coeff(q, i+1)*
                        (to_ZZ(2*t-3)*(i+1)-2*s) + (2*t-1)*c*R[i][j];
                    if (i+1 < d) f.m[i][j] += 2*c*(to_ZZ(s)+i+1)*S[i+1][j];
                }
        }
        return f;
    }

    // Reduce a bounded starting block without horizontal steps at t>0.
    // This remains valid for even degree: their horizontal maps can have
    // a zero denominator on the residue line s=d(t-1/2).
    Fraction initial_block(long a, long b, long r = 0) const {
        Fraction f;
        ident(f.m, d);
        long s = a*(2*r+1)-1, t = (b*(2*r+1)-1)/2;
        auto apply = [&](char kind) {
            Fraction z = step(s, t, kind);
            f.m = z.m*f.m; f.d *= z.d;
            if (kind != 'V') --s;
            if (kind != 'H') --t;
        };
        while (t > s+1) apply('V');
        while (t > 0) apply('D');
        while (s >= 0) apply('H');
        cancel(f);
        return f;
    }

};

struct Pair {
    long a, b;
    // Coefficient of U^{a,b} in column i, scaled by the common power of 2.
    std::vector<ZZ> weights;
    struct PolynomialBlock;
    std::shared_ptr<PolynomialBlock> polynomial;
};

// The matrix product of Proposition 8 is a polynomial matrix in r. Remove
// its exact common polynomial content with the scalar denominator before
// passing it to rforest. This cancels whole affine factors, rather than only
// cancelling their values after evaluating each block separately.
struct Pair::PolynomialBlock {
    long d, degree, start = 0;
    std::vector<ZZX> entries;
    ZZX denominator;
    std::vector<ZZX> denominator_factors;
    Fraction initial;
    long height_bits = 0;

    PolynomialBlock(const Reductions& red, long a, long b)
        : d(red.d), entries(d*d) {
        // In even degree the horizontal denominator can vanish at an
        // integer. Below the residue line, do H before D, so H moves away
        // from that line. Above it, D before H is safe. A bounded prefix
        // puts the low-slope path strictly below the line at its start.
        const long gap = d*b-2*a;
        const bool horizontal_first = d%2 == 0 && gap > 0 && 2*a > b;
        if (horizontal_first) start = std::max(0L, (d-1-gap/2)/gap);
        initial = red.initial_block(a, b, start);
        for (long i = 0; i < d; ++i) SetCoeff(entries[i*d+i], 0, 1);
        SetCoeff(denominator, 0, 1);
        std::vector<ZZX> factors;
        long s = a-1, t = (b-1)/2;
        auto apply = [&](char kind) {
            Fraction f0 = red.step(s, t, kind);
            Fraction f1 = red.step(s+2*a, t+b, kind);
            mat_ZZ slope = f1.m-f0.m;
            std::vector<ZZX> next(d*d);
            for (long i = 0; i < d; ++i)
                for (long k = 0; k < d; ++k) {
                    ZZX linear;
                    SetCoeff(linear, 0, f0.m[i][k]);
                    SetCoeff(linear, 1, slope[i][k]);
                    if (IsZero(linear)) continue;
                    for (long j = 0; j < d; ++j)
                        if (!IsZero(entries[k*d+j]))
                            next[i*d+j] += linear*entries[k*d+j];
                }
            entries.swap(next);
            ZZX linear;
            SetCoeff(linear, 0, f0.d);
            SetCoeff(linear, 1, f1.d-f0.d);
            denominator *= linear;
            ZZ content;
            GCD(content, coeff(linear, 0), coeff(linear, 1));
            if (LeadCoeff(linear) < 0) content = -content;
            linear /= content;
            if (deg(linear) > 0) factors.push_back(std::move(linear));
            if (kind != 'V') --s;
            if (kind != 'H') --t;
        };
        if (horizontal_first)
            for (long j = b; j < 2*a; ++j) apply('H');
        for (long j = 0; j < std::min(2*a, b); ++j) apply('D');
        if (!horizontal_first)
            for (long j = b; j < 2*a; ++j) apply('H');
        for (long j = 2*a; j < b; ++j) apply('V');
        // All nonconstant denominator factors are already known affine
        // polynomials. Exact division by these primitive factors is much
        // cheaper than generic modular polynomial GCDs, especially at the
        // higher precisions used by the genus-three regression cases.
        for (const auto& factor : factors) {
            std::vector<ZZX> quotients(d*d);
            bool common = true;
            for (long ij = 0; ij < d*d; ++ij)
                if (!IsZero(entries[ij]) && !divide(quotients[ij], entries[ij], factor)) {
                    common = false;
                    break;
                }
            if (common) {
                entries.swap(quotients);
                ZZX quotient;
                if (!divide(quotient, denominator, factor))
                    throw std::logic_error("nonexact polynomial denominator cancellation");
                denominator = quotient;
            }
        }
        ZZ common;
        content(common, denominator);
        for (const auto& entry : entries)
            for (long k = 0; k <= deg(entry) && common != 1; ++k)
                GCD(common, common, coeff(entry, k));
        if (common != 1) {
            denominator /= common;
            for (auto& entry : entries) entry /= common;
        }
        // Recover the remaining affine factors to compute the exact p-adic
        // loss after polynomial cancellation. The leftover constant is a
        // unit at every p-admissible prime.
        ZZX rest = denominator, quotient;
        for (const auto& factor : factors)
            if (deg(rest) > 0 && divide(quotient, rest, factor)) {
                denominator_factors.push_back(factor);
                rest = quotient;
            }
        if (deg(rest) != 0) throw std::logic_error("unfactored reduction denominator");
        degree = 0;
        for (const auto& entry : entries) {
            degree = std::max(degree, deg(entry));
            for (long i = 0; i <= deg(entry); ++i)
                height_bits = std::max(height_bits, NumBits(coeff(entry, i)));
        }
    }

    long valuation(long p) const {
        long v = 0;
        const ZZ pz = to_ZZ(p);
        for (const auto& factor : denominator_factors) {
            ZZ root = (-coeff(factor, 0)*InvMod(coeff(factor, 1)%pz, pz))%pz;
            if (root > start && root <= (p-1)/2) ++v;
        }
        return v;
    }
};

class MpzVector {
    mpz_t* values;
    size_t length;
public:
    explicit MpzVector(size_t n) : values(new mpz_t[n]), length(n) {
        for (size_t i = 0; i < length; ++i) mpz_init(values[i]);
    }
    ~MpzVector() {
        for (size_t i = 0; i < length; ++i) mpz_clear(values[i]);
        delete[] values;
    }
    MpzVector(const MpzVector&) = delete;
    MpzVector& operator=(const MpzVector&) = delete;
    mpz_ptr operator[](size_t i) { return values[i]; }
    mpz_t* data() { return values; }
};

struct GmpConversion {
    std::vector<unsigned char> bytes;
    void to_gmp(mpz_ptr out, const ZZ& value) {
        if (NumBits(value) < NTL_BITS_PER_LONG-1) {
            mpz_set_si(out, to_long(value));
            return;
        }
        long n = NumBytes(value);
        bytes.resize(n);
        BytesFromZZ(bytes.data(), value, n);
        mpz_import(out, n, -1, 1, 0, 0, bytes.data());
        if (value < 0) mpz_neg(out, out);
    }
    void from_gmp(ZZ& out, mpz_srcptr value) {
        if (mpz_fits_slong_p(value)) { conv(out, mpz_get_si(value)); return; }
        bytes.resize((mpz_sizeinbase(value, 2)+7)/8);
        size_t n;
        mpz_export(bytes.data(), &n, -1, 1, 0, 0, value);
        ZZFromBytes(out, bytes.data(), n);
        if (mpz_sgn(value) < 0) NTL::negate(out, out);
    }
};

// rforest's FFT plans have thread-local ownership in the vendored source.
// Keep one context alive per worker so numerator/denominator calls and the
// following reduction pair reuse their initialized FFT tables.
struct RforestContext {
    bool active;
    explicit RforestContext(bool enabled) : active(enabled) {
        if (active) {
            const char* fft = std::getenv("HYPELLFROB_RFOREST_FFT");
            hw_disable_fft = fft && std::strcmp(fft, "0") == 0;
            const char* evaluator = std::getenv("HYPELLFROB_RFOREST_EVAL");
            rforest_forward_differences = !(evaluator && std::strcmp(evaluator, "powers") == 0);
            hw_mpz_setup();
        }
    }
    ~RforestContext() { if (active) hw_mpz_clear(); }
};

struct ForestOutput {
    MpzVector matrix, denominator;
    long d, rows;
    int kappa_used;
    std::vector<long> losses;
    explicit ForestOutput(const Pair::PolynomialBlock& block,
                          const std::vector<AverageJob>& jobs,
                          const std::vector<size_t>& indices, int mu, int kappa)
        : matrix(indices.size()*(block.d-1)*block.d), denominator(indices.size()),
          d(block.d), rows(block.d-1) {
        GmpConversion conversion;
        const long n = indices.size();
        MpzVector moduli(n), product(2), initial(rows*d), scalar_initial(1);
        MpzVector polynomial(d*d*(block.degree+1));
        const long scalar_degree = deg(block.denominator);
        MpzVector scalar_polynomial(scalar_degree+1);
        std::vector<long> endpoints(n);
        losses.resize(n);
        for (long i = 0; i < n; ++i) {
            long p = jobs[indices[i]].p;
            losses[i] = block.valuation(p);
            mpz_ui_pow_ui(moduli[i], p, mu+losses[i]);
            endpoints[i] = (p+1)/2;
        }
        mproduct(product[0], moduli.data(), n);
        mpz_set(product[1], product[0]);
        for (long i = 0; i < rows; ++i)
            for (long j = 0; j < d; ++j)
                conversion.to_gmp(initial[i*d+j], block.initial.m[i+1][j]);
        conversion.to_gmp(scalar_initial[0], block.initial.d);
        for (long ij = 0; ij < d*d; ++ij)
            for (long k = 0; k <= block.degree; ++k)
                conversion.to_gmp(polynomial[ij*(block.degree+1)+k], coeff(block.entries[ij], k));
        for (long k = 0; k <= scalar_degree; ++k)
            conversion.to_gmp(scalar_polynomial[k], coeff(block.denominator, k));
        if (kappa < 0) {
            // Use O(log N) subtrees, preserving average polynomial time.
            // A fixed leaf count would instead require O(N) updates of a
            // batch-sized prefix and could reintroduce quadratic work.
            long bits = 0;
            for (long v = n; v; v >>= 1) ++bits;
            long bound_bits = NumBits(to_ZZ(endpoints.back()));
            kappa = std::min(std::max(0L, NumBits(to_ZZ(bound_bits))+1),
                             std::max(0L, bits-3));
        }
        kappa_used = kappa;
        rforest(matrix.data(), initial.data(), rows, polynomial.data(), block.degree,
                d, moduli.data(), block.start+1, endpoints.data(), n, product[0], kappa);
        rforest(denominator.data(), scalar_initial.data(), 1, scalar_polynomial.data(),
                scalar_degree, 1, moduli.data(), block.start+1, endpoints.data(), n, product[1], kappa);
    }

    void column(long i, GmpConversion& conversion, Fraction& out) {
        out.m.SetDims(rows, 1);
        for (long row = 0; row < rows; ++row)
            conversion.from_gmp(out.m[row][0], matrix[(i*rows+row)*d]);
        conversion.from_gmp(out.d, denominator[i]);
    }
};

void compute_group(std::vector<AverageJob>& jobs, const std::vector<size_t>& indices,
                   const Reductions& red, int mu, int nthreads,
                   int kappa) {
    long dim = red.d-1;
    // Proposition 2. alpha_j = sum_{k=j}^{mu-1} (-1)^j binom(2k,k)/4^k binom(k,j).
    // Use a common denominator 4^(mu-1), a unit at every odd prime.
    ZZ alpha_den = power2_ZZ(2*(mu-1));
    std::vector<ZZ> alpha(mu);
    for (long j = 0; j < mu; ++j)
        for (long k = j; k < mu; ++k) {
            ZZ central(1), choose(1);
            for (long v = 1; v <= k; ++v) central = central*(k+v)/v;
            for (long v = 1; v <= j; ++v) choose = choose*(k-v+1)/v;
            ZZ term = central*choose*power2_ZZ(2*(mu-1-k));
            alpha[j] += j%2 ? -term : term;
        }
    std::vector<Pair> pairs;
    ZZX qpower;
    SetCoeff(qpower, 0, 1);
    for (long j = 0; j < mu; ++j) {
        for (long a = 1; a <= dim+deg(qpower); ++a) {
            Pair pair{a, 2*j+1, std::vector<ZZ>(dim)};
            bool nonzero = false;
            for (long i = 0; i < dim; ++i) {
                pair.weights[i] = alpha[j]*coeff(qpower, a-i-1);
                nonzero |= !IsZero(pair.weights[i]);
            }
            if (nonzero) pairs.push_back(std::move(pair));
        }
        qpower *= red.q;
    }
    // p-admissibility (Section 5) is checked for all nonzero terms together.
    std::vector<size_t> eligible;
    for (size_t idx : indices) {
        long p = jobs[idx].p;
        bool ok = p > red.d*(2*mu-1) && red.delta % p != 0;
        for (const auto& pair : pairs)
            if (p <= red.d*pair.b+2*pair.a ||
                (pair.b > 2*pair.a && coeff(red.q, 0)%p == 0)) ok = false;
        if (ok) { eligible.push_back(idx); jobs[idx].matrix.SetDims(dim, dim); }
    }
    if (eligible.empty()) return;
    std::sort(eligible.begin(), eligible.end(), [&](size_t i, size_t j) {
        return jobs[i].p < jobs[j].p;
    });
    for (auto& pair : pairs)
        pair.polynomial = std::make_shared<Pair::PolynomialBlock>(red, pair.a, pair.b);
    // The bounded prefix is exact over Q. Its denominator must be a
    // p-adic unit; exceptional small primes are left to the Sage caller.
    eligible.erase(std::remove_if(eligible.begin(), eligible.end(), [&](size_t idx) {
        for (const auto& pair : pairs)
            if (pair.polynomial->initial.d % jobs[idx].p == 0 ||
                (jobs[idx].p-1)/2 < pair.polynomial->start) return true;
        return false;
    }), eligible.end());
    if (eligible.empty()) return;
    // Start the expensive trees first; bounded forests make this possible
    // without the old full-tree implementation's large memory spike.
    std::stable_sort(pairs.begin(), pairs.end(), [](const Pair& x, const Pair& y) {
        return x.polynomial->degree*64+x.polynomial->height_bits >
               y.polynomial->degree*64+y.polynomial->height_bits;
    });
    std::atomic<size_t> next(0);
    std::mutex output_mutex, error_mutex;
    std::exception_ptr error;
    auto worker = [&]() {
        try {
            RforestContext context(true);
            for (;;) {
                size_t pi = next.fetch_add(1);
                if (pi >= pairs.size()) break;
                const Pair& pair = pairs[pi];
                const auto pair_start = std::chrono::steady_clock::now();
                ForestOutput forest(*pair.polynomial, jobs, eligible, mu, kappa);
                GmpConversion conversion;
                Fraction value;
                std::vector<ZZ> column(dim);
                const auto forest_done = std::chrono::steady_clock::now();
                for (size_t ei = 0; ei < eligible.size(); ++ei) {
                    size_t idx = eligible[ei];
                    long p = jobs[idx].p;
                    long rho = forest.losses[ei];
                    forest.column(ei, conversion, value);
                    ZZ modulus = power(to_ZZ(p), mu), den = value.d;
                    ZZ pvaluation(1);
                    long v = 0;
                    if (IsZero(den)) throw std::runtime_error("denominator precision exhausted");
                    while (den%p == 0) { den /= p; pvaluation *= p; ++v; }
                    if (v > rho) throw std::runtime_error("denominator valuation exceeds bound");
                    ZZ inverse = InvMod((den*alpha_den)%modulus, modulus);
                    for (long row = 0; row < dim; ++row) {
                        ZZ num = p*value.m[row][0];
                        if (num%pvaluation != 0)
                            throw std::runtime_error("nonintegral p*U in Frobenius expansion");
                        column[row] = ((num/pvaluation)*inverse)%modulus;
                    }
                    std::lock_guard<std::mutex> lock(output_mutex);
                    for (long col = 0; col < dim; ++col)
                        if (!IsZero(pair.weights[col]))
                            for (long row = 0; row < dim; ++row)
                                jobs[idx].matrix[row][col] = (jobs[idx].matrix[row][col] +
                                    column[row]*pair.weights[col])%modulus;
                }
                if (std::getenv("HYPELLFROB_PROFILE")) {
                    double fs = std::chrono::duration<double>(forest_done-pair_start).count();
                    double assembly = std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-forest_done).count();
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cerr << "rforest a=" << pair.a << " b=" << pair.b
                              << " degree=" << pair.polynomial->degree
                              << " denominator_degree=" << deg(pair.polynomial->denominator)
                              << " start=" << pair.polynomial->start
                              << " kappa=" << forest.kappa_used
                              << " forest_seconds=" << fs << " assembly_seconds=" << assembly << '\n';
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!error) error = std::current_exception();
        }
    };
    std::vector<std::thread> pool;
    long workers = std::min<long>(std::max(1, nthreads), pairs.size());
    for (long i = 0; i < workers; ++i) pool.emplace_back(worker);
    for (auto& thread : pool) thread.join();
    if (error) std::rethrow_exception(error);
    for (size_t idx : eligible) {
        if (red.d%2 == 0) {
            mat_ZZ projective;
            restrict_even_frobenius(projective, jobs[idx].matrix,
                                   to_ZZ(jobs[idx].p), mu, red.q);
            jobs[idx].matrix = std::move(projective);
        }
        jobs[idx].ok = true;
    }
}
} // namespace

void matrices_rforest(std::vector<AverageJob>& jobs, const ZZX& Q, int nthreads, int kappa) {
    for (auto& job : jobs) { job.ok = false; job.matrix.kill(); }
    if (deg(Q) < 3 || LeadCoeff(Q) != 1) return;
    if (kappa < -1 || kappa > 30) throw std::invalid_argument("rforest kappa must be -1..30");
    std::map<int, std::vector<size_t>> groups;
    for (size_t i = 0; i < jobs.size(); ++i)
        if (jobs[i].precision > 0 && jobs[i].p >= 3 && ProbPrime(jobs[i].p))
            groups[jobs[i].precision].push_back(i);
    if (groups.empty()) return;
    try {
        Reductions red(Q);
        for (const auto& group : groups)
            compute_group(jobs, group.second, red, group.first, nthreads, kappa);
    } catch (const std::invalid_argument&) {
        return;
    }
}
} // namespace hypellfrob
