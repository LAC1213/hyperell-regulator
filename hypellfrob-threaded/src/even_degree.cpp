/* ============================================================================

   even_degree.cpp   -- see even_degree.h

   Kedlaya's algorithm, written directly, for either parity of deg(Q).
   Linear in p; this is the reference implementation, not the fast path.

   The whole reason this file exists is the even-degree singularity. Write
   the reduction relation, multiplied through by 2, as

       Lambda(x^k) = 2k x^(k-1) Q(x) + (2 - tt) x^k Q'(x),

   which is exact as a differential against dx/y^tt (it is 2 d(x^k y^(2-tt))
   rewritten using 2y dy = Q' dx). Its degree is k + d - 1 and its leading
   coefficient is

       lam_k = 2k + (2 - tt) d.

   Reducing the top term of f uses Lambda at k = deg(f) - d + 1 and divides
   by lam_k. That vanishes at k = d(tt-2)/2, an integer EXACTLY when d is
   even -- so for even degree the greedy reduction gets stuck, once, at
   absolute degree

       n_stuck = d(tt-2)/2 + d - 1.

   The resolution is that n_stuck is simultaneously three things: the stuck
   degree, the lowest degree at which x^n dx/y^tt has nonvanishing residue at
   infinity ((g+1)tt - 1), and the degree of x^g Q(x)^(t-1) where tt = 2t-1.
   The last of these is the point: since y^(2t-2) = Q^(t-1),

       x^g Q^(t-1) dx / y^tt   IS   x^g dx / y,

   the same cohomology class at every level. So the stuck monomial is the
   leading term of a representative of that one class, and

       c x^(n_stuck) dx/y^tt = c x^g dx/y - c (x^g Q^(t-1) - x^(n_stuck)) dx/y^tt

   splits it into a contribution to basis vector g -- which bypasses the rest
   of the reduction entirely, being level-independent -- plus a remainder of
   degree < n_stuck that reduces normally. That is the whole fix.

============================================================================ */

#include "even_degree.h"

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <vector>
#include <stdexcept>

NTL_CLIENT

namespace hypellfrob {

namespace {

// ---------------------------------------------------------------------------
// p-adic division by an ordinary integer that may be divisible by p.
//
// The divisors here -- lam_k in the horizontal step, tt-2 in the vertical one
// -- are honest integers whose valuation is exact and known, and BOTH are
// genuinely divisible by p for some of the values they take as the reduction
// sweeps (tt-2 == 0 mod p whenever tt == 2 mod p, for instance). So these
// divisions really do lower the valuation; they are not exact divisions.
//
// The caller handles that by pre-multiplying everything by p^scale, which
// gives every intermediate enough valuation to absorb the drops. A failure
// here therefore means the scale was too small, and the caller retries with a
// larger one rather than returning a wrong answer.
// ---------------------------------------------------------------------------
struct PrecisionFailure : std::runtime_error {
    PrecisionFailure() : std::runtime_error("even_degree: out of p-adic precision") {}
};

ZZ_p divide_by_integer(const ZZ_p& a, const ZZ& b, const ZZ& p)
{
    if (IsZero(b)) throw std::runtime_error("even_degree: division by exact zero");

    ZZ bb = b;
    long sign = 1;
    if (bb < 0) { bb = -bb; sign = -1; }

    long v = 0;
    while (bb % p == 0) { bb /= p; v++; }

    ZZ ar = rep(a);
    if (v > 0) {
        ZZ pv = power(p, v);
        if (ar % pv != 0) throw PrecisionFailure();
        ar /= pv;
    }

    ZZ m = ZZ_p::modulus();
    ZZ inv;
    InvMod(inv, bb % m, m);
    ZZ_p r = to_ZZ_p(ar) * to_ZZ_p(inv);
    return (sign > 0) ? r : -r;
}

ZZ_pX divide_by_integer(const ZZ_pX& f, const ZZ& b, const ZZ& p)
{
    ZZ_pX r;
    r.SetLength(deg(f) + 1);
    for (long i = 0; i <= deg(f); i++)
        SetCoeff(r, i, divide_by_integer(coeff(f, i), b, p));
    r.normalize();
    return r;
}

// ---------------------------------------------------------------------------
// Lambda(x^k) = 2k x^(k-1) Q + (2-tt) x^k Q',  degree k + d - 1.
// ---------------------------------------------------------------------------
ZZ_pX lambda_poly(long k, const ZZ& tt, const ZZ_pX& Q, const ZZ_pX& dQ)
{
    ZZ_pX a = Q;                       // 2k x^(k-1) Q
    a *= to_ZZ_p(to_ZZ(2 * k));
    if (k >= 1) a <<= (k - 1);
    else clear(a);                     // k = 0 kills this term

    ZZ_pX b = dQ;                      // (2-tt) x^k Q'
    b *= to_ZZ_p(to_ZZ(2) - tt);
    b <<= k;

    return a + b;
}

// ---------------------------------------------------------------------------
// Reduce f(x) dx / y^tt to the basis x^0..x^(d-2) dx/y.
//
// `lambda` accumulates the coefficient of the residue class x^g dx/y split
// off at the singular step; it deliberately does NOT pass through the
// vertical reduction, because that class is the same at every level.
// ---------------------------------------------------------------------------
void reduce(ZZ_pX& f, ZZ tt, const ZZ_pX& Q, const ZZ_pX& dQ,
            const std::vector<ZZ_pX>& R, const std::vector<ZZ_pX>& S,
            const ZZ& p, int g, int d, ZZ_p& lambda)
{
    for (;;) {
        // ---- horizontal: bring deg f down to d-2 --------------------------
        while (deg(f) > d - 2) {
            long n = deg(f);
            long k = n - (d - 1);
            ZZ lam = to_ZZ(2 * k) + (to_ZZ(2) - tt) * d;

            if (IsZero(lam)) {
                // The singular step. n is simultaneously the stuck degree,
                // the residue threshold, and deg(x^g Q^(t-1)) -- see the file
                // header. Split the class off and carry on.
                ZZ t = (tt + 1) / 2;
                ZZ_pX G = power(Q, to_long(t - 1));
                G <<= g;                                  // x^g Q^(t-1)
                if (deg(G) != n)
                    throw std::runtime_error("even_degree: residue representative "
                                              "has the wrong degree");
                ZZ_p c = coeff(f, n);
                lambda += c;
                f -= c * G;
                continue;
            }

            ZZ_p c = divide_by_integer(coeff(f, n), lam, p);
            f -= c * lambda_poly(k, tt, Q, dQ);
        }

        if (tt == 1) break;

        // ---- vertical: tt -> tt-2 -----------------------------------------
        // x^i dx/y^tt = [R_i + 2 S_i' / (tt-2)] dx/y^(tt-2), from
        // x^i = R_i Q + S_i Q' and 2 d(S_i y^(2-tt)).
        ZZ_pX nf;
        for (long i = 0; i <= d - 2; i++) {
            ZZ_p fi = coeff(f, i);
            if (IsZero(fi)) continue;
            // Divide AFTER multiplying by f_i: S_i is a fixed Bezout
            // polynomial with no p-adic room of its own, whereas f_i carries
            // the caller's p^scale, which is what absorbs v_p(tt-2).
            nf += fi * R[i];
            nf += divide_by_integer((to_ZZ_p(2) * fi) * diff(S[i]), tt - 2, p);
        }
        f = nf;
        tt -= 2;
    }
}

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

// ---------------------------------------------------------------------------

static int matrix_affine_at_scale(mat_ZZ& output, const ZZ& p, int N,
                                  const ZZX& Q, int extra_precision);

// Retry with a larger scale if the reduction runs out of valuation, rather
// than returning a wrong answer or a bare failure.
int matrix_affine_reference(mat_ZZ& output, const ZZ& p, int N,
                            const ZZX& Q, int extra_precision)
{
    for (int attempt = 0; attempt < 6; attempt++) {
        int r = matrix_affine_at_scale(output, p, N, Q, extra_precision);
        if (r) return r;
        extra_precision *= 2;
    }
    return 0;
}

static int matrix_affine_at_scale(mat_ZZ& output, const ZZ& p, int N,
                                  const ZZX& Q, int extra_precision)
{
    if (N < 1 || p < 3) return 0;
    int d = deg(Q);
    if (d < 3 || !IsOne(LeadCoeff(Q))) return 0;
    if (p <= (2 * N - 1) * d) return 0;

    int dim = d - 1;                       // 2g for odd d, 2g+1 for even d
    int g = (d % 2 == 0) ? (d - 2) / 2 : (d - 1) / 2;

    ZZ_pContext caller;  caller.save();

    // Everything is computed p^scale times too large, so that the divisions
    // in the reduction (which genuinely lower the valuation) always have room;
    // the factor is divided out at the end. The modulus carries the scale plus
    // the N digits actually wanted plus a little slack.
    const long scale = extra_precision;
    ZZ_p::init(power(p, N + scale + 4));

    ZZ_pX Qp;  conv(Qp, Q);
    ZZ_pX dQp = diff(Qp);

    // Bezout: x^i = R_i Q + S_i Q'. Needs gcd(Q, Q') = 1, i.e. Q squarefree
    // mod p; check that first in zz_p, where GCD is meaningful.
    {
        zz_pContext sv; sv.save();
        zz_p::init(to_long(p));
        zz_pX q, dq, gg;
        conv(q, Q);
        if (deg(q) != d) { sv.restore(); caller.restore(); return 0; }
        dq = diff(q);
        if (IsZero(dq)) { sv.restore(); caller.restore(); return 0; }
        GCD(gg, q, dq);
        bool ok = (deg(gg) == 0);
        sv.restore();
        if (!ok) { caller.restore(); return 0; }
    }

    // Bezout x^i = R_i Q + S_i Q'. NTL's ZZ_pX XGCD is not valid over a
    // composite modulus, so do it over Z: XGCD gives s*Q + t*Q' = res(Q,Q'),
    // and the resultant is a unit mod p exactly when Q is squarefree mod p
    // (already checked above), so dividing through by it is legitimate.
    ZZ_pX u, v;
    {
        ZZ res;
        ZZX sZ, tZ, dQZ = diff(Q);
        XGCD(res, sZ, tZ, Q, dQZ);
        if (res % p == 0) { caller.restore(); return 0; }
        ZZ_pX sp, tp;
        conv(sp, sZ);  conv(tp, tZ);
        ZZ_p rinv = to_ZZ_p(1) / to_ZZ_p(res);
        u = sp * rinv;  v = tp * rinv;
    }

    std::vector<ZZ_pX> R(d - 1), S(d - 1);
    for (int i = 0; i <= d - 2; i++) {
        // x^i = (x^i u) Q + (x^i v) Q'; reduce so deg S_i < d, folding the
        // quotient into R_i, which keeps both degrees small.
        ZZ_pX xi;  SetCoeff(xi, i, 1);
        ZZ_pX Si = (xi * v) % Qp;
        ZZ_pX Ri = (xi - Si * dQp) / Qp;
        R[i] = Ri;  S[i] = Si;
    }

    // ---- Frobenius ---------------------------------------------------------
    // F*(x^i dx/y) = sum_{j<N} p * C(-1/2, j) * x^(p(i+1)-1) * D^j dx/y^(p(2j+1))
    // with D = Q(x^p) - Q(x)^p, which is divisible by p, so the j-th term
    // carries p^(j+1) and terms with j >= N vanish mod p^N.
    long pl = to_long(p);
    ZZ_pX Qxp;                                     // Q(x^p)
    for (long k = 0; k <= d; k++) SetCoeff(Qxp, k * pl, coeff(Qp, k));
    ZZ_pX D = Qxp - power(Qp, pl);

    // The j-th Frobenius term carries p^(j+1), so terms with j >= N vanish
    // mod p^N -- but the reduction divides by powers of p along the way, so
    // carry a few extra terms rather than cutting exactly at N.
    const int jmax = N + 4;

    // C(-1/2, j) = (-1)^j * binom(2j, j) / 4^j
    std::vector<ZZ_p> binom(jmax);
    for (int j = 0; j < jmax; j++) {
        ZZ bin = to_ZZ(1);                         // binom(2j, j)
        for (int a = 1; a <= j; a++) bin = bin * (2 * j - a + 1) / a;
        ZZ_p val = to_ZZ_p(bin) / power(to_ZZ_p(4), j);
        binom[j] = (j % 2 == 0) ? val : -val;
    }

    output.SetDims(dim, dim);

    bool ok = true;
    for (int i = 0; i < dim && ok; i++) {
        ZZ_pX total;
        ZZ_p lambda = to_ZZ_p(0);

        ZZ_pX Dpow;  set(Dpow);                    // D^0 = 1
        for (int j = 0; j < jmax; j++) {
            ZZ_pX f = Dpow;
            f <<= (pl * (i + 1) - 1);              // x^(p(i+1)-1) * D^j
            f *= to_ZZ_p(p) * binom[j] * to_ZZ_p(power(p, scale));

            ZZ tt = p * (2 * j + 1);
            try {
                reduce(f, tt, Qp, dQp, R, S, p, g, d, lambda);
            } catch (const PrecisionFailure&) {
                ok = false; break;      // caller retries with a larger scale
            }
            total += f;
            Dpow *= D;
        }
        if (!ok) break;

        // The residue class split off at the singular steps is x^g dx/y.
        if (d % 2 == 0) {
            ZZ_pX xg;  SetCoeff(xg, g, 1);
            total += lambda * xg;
        }

        // Undo the p^scale pre-multiplication. The true values are integral,
        // so this must divide exactly; if it does not, the scale was too small.
        ZZ pscale = power(p, scale);
        ZZ mod = power(p, N);
        for (int f = 0; f < dim; f++) {
            ZZ c = rep(coeff(total, f));
            if (c % pscale != 0) { ok = false; break; }   // retry with more scale
            output[f][i] = (c / pscale) % mod;
        }
    }

    caller.restore();
    return ok ? 1 : 0;
}

int matrix_reference(mat_ZZ& output, const ZZ& p, int N,
                     const ZZX& Q, int extra_precision)
{
    int d = deg(Q);
    mat_ZZ M;
    if (!matrix_affine_reference(M, p, N, Q, extra_precision)) return 0;

    if (d % 2 != 0) { output = M; return 1; }      // odd: already H^1 of the curve

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
    return 1;
}

} // namespace hypellfrob
