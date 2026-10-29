// SPDX-License-Identifier: GPL-2.0-or-later
// An implementation of David Harvey, "Counting points on hyperelliptic
// curves in average polynomial time", arXiv:1210.8239v3, Sections 3--5.
// Original implementation for hyperell-regulator, 2026.
//
// All matrices act on columns in the basis of W_{s,t}. In particular the
// accumulating products are M_0 M_1 ... M_{(p-1)/2}, in that order.
#include "average.h"
#include <NTL/mat_ZZ.h>
#include <algorithm>
#include <atomic>
#include <exception>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>

NTL_CLIENT
namespace hypellfrob {
namespace {

void reduce_mod(mat_ZZ& a, const ZZ& m) {
    for (long i = 0; i < a.NumRows(); ++i)
        for (long j = 0; j < a.NumCols(); ++j) a[i][j] %= m;
}

void reduce_centered(mat_ZZ& a, const ZZ& m, const ZZ& half) {
    for (long i = 0; i < a.NumRows(); ++i)
        for (long j = 0; j < a.NumCols(); ++j) {
            ZZ& x = a[i][j];
            // Turning a small negative integer into m-1 at every product
            // node would destroy the size bounds of Proposition 4.
            if (abs(x) <= half) continue;
            x %= m;
            if (x > half) x -= m;
        }
}

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

    Fraction block(long a, long b, long r, bool normalize = true) const {
        Fraction f;
        ident(f.m, d);
        long s = a*(2*r+1)-1, t = (b*(2*r+1)-1)/2;
        auto apply = [&](char kind) {
            Fraction z = step(s, t, kind);
            f.m = z.m*f.m; f.d *= z.d;
            if (kind != 'V') --s;
            if (kind != 'H') --t;
        };
        if (r) { // Proposition 8
            for (long j = 0; j < std::min(2*a, b); ++j) apply('D');
            for (long j = b; j < 2*a; ++j) apply('H');
            for (long j = 2*a; j < b; ++j) apply('V');
        } else { // Proposition 9: vertical steps must precede the final s=0 step
            while (t > s+1) apply('V');
            while (t > 0) apply('D');
            while (s >= 0) apply('H');
        }
        if (normalize) cancel(f);
        return f;
    }
};

// Proposition 8 makes every entry a polynomial in r of degree at most
// max(2a,b). Forward differences evaluate these polynomials using additions,
// avoiding O(N) repetitions of all the small reduction-matrix products.
class BlockSequence {
    std::vector<Fraction> differences;
public:
    BlockSequence(const Reductions& red, long a, long b) {
        long degree = std::max(2*a, b);
        std::vector<Fraction> values;
        for (long r = 1; r <= degree+1; ++r)
            values.push_back(red.block(a, b, r, false));
        for (long order = 0; order <= degree; ++order) {
            differences.push_back(values[0]);
            for (long i = 0; i < degree-order; ++i) {
                values[i].m = values[i+1].m-values[i].m;
                values[i].d = values[i+1].d-values[i].d;
            }
        }
    }
    Fraction next() {
        Fraction f = differences[0];
        for (size_t i = 0; i+1 < differences.size(); ++i) {
            differences[i].m += differences[i+1].m;
            differences[i].d += differences[i+1].d;
        }
        cancel(f);
        return f;
    }
};

struct Target {
    size_t job;
    long k;
    ZZ modulus;
};

// Section 3. Leaves represent indices k and return the prefix THROUGH M_k.
// Store products M_{k+1}, so the right child's prefix includes its left sibling.
// Matrices need not be reduced modulo the product of *unrelated* primes.
// Capping every product at the root modulus bounds coefficient growth without
// changing a residue at any requested leaf.
class RemainderTree {
    const Reductions& red;
    long a, b, base;
    std::vector<Fraction> products;
    std::vector<ZZ> moduli;

    void descend(long node, const Fraction& prefix,
                 const std::vector<Target>& targets,
                 std::vector<Fraction>& answers) const {
        if (moduli[node] == 1) return;
        if (node >= base) {
            long k = node-base;
            auto it = std::lower_bound(targets.begin(), targets.end(), k,
                [](const Target& t, long v) { return t.k < v; });
            for (; it != targets.end() && it->k == k; ++it) {
                Fraction f = prefix;
                reduce_mod(f.m, it->modulus); f.d %= it->modulus;
                answers[it->job] = std::move(f);
            }
            return;
        }
        for (long side = 0; side < 2; ++side) {
            long child = 2*node+side;
            if (moduli[child] == 1) continue;
            Fraction f = prefix;
            if (side) {
                f.m *= products[2*node].m;
                f.d *= products[2*node].d;
            }
            reduce_mod(f.m, moduli[child]); f.d %= moduli[child];
            descend(child, f, targets, answers);
        }
    }

public:
    RemainderTree(const Reductions& red_, long a_, long b_)
        : red(red_), a(a_), b(b_), base(1) {}

    void run(const std::vector<Target>& targets, std::vector<Fraction>& answers) {
        long last = targets.back().k;
        while (base <= last) base *= 2;
        products.resize(2*base); moduli.assign(2*base, ZZ(1));
        // Duplicate jobs at the same p need only the largest requested power.
        for (const auto& t : targets)
            if (t.modulus > moduli[base+t.k]) moduli[base+t.k] = t.modulus;
        for (long i = base-1; i > 0; --i) moduli[i] = moduli[2*i]*moduli[2*i+1];
        // This threshold has O(N) bits. Compute it once: computing m/2
        // at every one of the O(N) product nodes adds quadratic bit-work.
        const ZZ half_modulus = moduli[1]/2;
        BlockSequence sequence(red, a, b);
        for (long k = 0; k < base; ++k) {
            if (k < last) products[base+k] = sequence.next();
            else ident(products[base+k].m, red.d);
        }
        for (long i = base-1; i > 0; --i) {
            // No prefix ever uses products beyond the last target.
            if (i == 1) continue;
            products[i].m = products[2*i].m*products[2*i+1].m;
            products[i].d = products[2*i].d*products[2*i+1].d;
            reduce_centered(products[i].m, moduli[1], half_modulus);
            products[i].d %= moduli[1];
            // Descent uses only left-sibling products. The right child's
            // product has had its last use now; its own left descendants
            // are kept for descent inside that subtree.
            products[2*i+1].m.kill();
            clear(products[2*i+1].d);
        }
        // The root product is never formed, so its right child has no use.
        if (base > 1) {
            products[3].m.kill();
            clear(products[3].d);
        }
        Fraction initial = red.block(a, b, 0);
        // The first coordinate in W_{-1,0} is identically zero. Only its
        // remaining d-1 rows are needed in all accumulating prefixes.
        mat_ZZ rows;
        rows.SetDims(red.d-1, red.d);
        for (long i = 1; i < red.d; ++i) rows[i-1] = initial.m[i];
        initial.m = std::move(rows);
        reduce_mod(initial.m, moduli[1]); initial.d %= moduli[1];
        descend(1, initial, targets, answers);
    }
};

struct Pair {
    long a, b;
    // Coefficient of U^{a,b} in column i, scaled by the common power of 2.
    std::vector<ZZ> weights;
};

long denominator_valuation(long d, long a, long b, long p) {
    // Sharpen Proposition 10's rho: each horizontal linear factor has a
    // single root mod p, but only roots in [1,(p-1)/2] occur in the product.
    // p-admissibility bounds every factor strictly below p^2 and makes the
    // final reduction's denominator a unit. Cancelled contents can only
    // decrease this valuation, so it remains a rigorous upper bound.
    long v = (b-1)/2;
    if (b <= 2*a) {
        ZZ pz = to_ZZ(p);
        ZZ inverse = InvMod((2*(to_ZZ(d)*b-2*a))%pz, pz);
        for (long j = 0; j < 2*a-b; ++j) {
            ZZ constant = to_ZZ(d)*(b+2)+2*(a-b-1-j);
            ZZ root = (constant*inverse)%pz;
            if (root >= 1 && root <= (p-1)/2) ++v;
        }
    }
    return v;
}

void compute_group(std::vector<AverageJob>& jobs, const std::vector<size_t>& indices,
                   const Reductions& red, int mu, int nthreads) {
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
    std::atomic<size_t> next(0);
    std::mutex output_mutex, error_mutex;
    std::exception_ptr error;
    auto worker = [&]() {
        try {
            for (;;) {
                size_t pi = next.fetch_add(1);
                if (pi >= pairs.size()) break;
                const Pair& pair = pairs[pi];
                std::vector<Target> targets;
                for (size_t idx : eligible) {
                    // Need p*U mod p^mu. The denominator needs mu+rho
                    // digits even though the numerator has an extra p:
                    // p*U can be a unit after cancellation of p^rho.
                    long rho = denominator_valuation(red.d, pair.a, pair.b, jobs[idx].p);
                    targets.push_back({idx, (jobs[idx].p-1)/2,
                        power(to_ZZ(jobs[idx].p), mu+rho)});
                }
                std::vector<Fraction> answers(jobs.size());
                RemainderTree(red, pair.a, pair.b).run(targets, answers);
                for (size_t idx : eligible) {
                    long p = jobs[idx].p;
                    long rho = denominator_valuation(red.d, pair.a, pair.b, p);
                    ZZ modulus = power(to_ZZ(p), mu), den = answers[idx].d;
                    ZZ pvaluation(1);
                    long v = 0;
                    if (IsZero(den)) throw std::runtime_error("denominator precision exhausted");
                    while (den%p == 0) { den /= p; pvaluation *= p; ++v; }
                    if (v > rho) throw std::runtime_error("denominator valuation exceeds bound");
                    ZZ inverse = InvMod((den*alpha_den)%modulus, modulus);
                    std::vector<ZZ> column(dim);
                    for (long row = 0; row < dim; ++row) {
                        ZZ num = p*answers[idx].m[row][0];
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
    for (size_t idx : eligible) jobs[idx].ok = true;
}
} // namespace

void matrices_average(std::vector<AverageJob>& jobs, const ZZX& Q, int nthreads) {
    for (auto& job : jobs) { job.ok = false; job.matrix.kill(); }
    if (deg(Q) < 3 || deg(Q)%2 == 0 || LeadCoeff(Q) != 1) return;
    // Group by requested precision, never impose a small prime's precision on
    // every larger prime. Each group remains a simultaneous many-prime job.
    std::map<int, std::vector<size_t>> groups;
    for (size_t i = 0; i < jobs.size(); ++i)
        if (jobs[i].precision > 0 && jobs[i].p >= 3 && ProbPrime(jobs[i].p))
            groups[jobs[i].precision].push_back(i);
    if (groups.empty()) return;
    try {
        Reductions red(Q);
        for (const auto& group : groups)
            compute_group(jobs, group.second, red, group.first, nthreads);
    } catch (const std::invalid_argument&) {
        // Singular over Q: there are no good primes.
        return;
    }
}
} // namespace hypellfrob
