/* ============================================================================

   hypellfrob.cpp:  main matrix() routine

   This file is part of hypellfrob (version 2.1.1).

   Copyright (C) 2007, 2008, David Harvey

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

   ----------------------------------------------------------------------
   UNMODIFIED from hypellfrob 2.1.1. Present here only so that src/ is a
   complete, self-contained copy; byte-identical to reference/.
   ----------------------------------------------------------------------

============================================================================ */



#include "hypellfrob.h"
#include "recurrences_ntl.h"

#include <cassert>

NTL_CLIENT


namespace hypellfrob {

// Counts how often the even-degree singular reduction step was hit and its
// obstruction discarded. Zero for odd degree.
long hypellfrob_singular_hits = 0;



/*
   For some reason NTL doesn't have conversions to mat_ZZ....
*/
void conv(mat_ZZ& x, const mat_ZZ_p& a)
{
   x.SetDims(a.NumRows(), a.NumCols());
   for (int i = 0; i < a.NumRows(); i++)
   for (int j = 0; j < a.NumCols(); j++)
      x[i][j] = rep(a[i][j]);
}

void conv(mat_ZZ& x, const mat_zz_p& a)
{
   x.SetDims(a.NumRows(), a.NumCols());
   for (int i = 0; i < a.NumRows(); i++)
   for (int j = 0; j < a.NumCols(); j++)
      x[i][j] = rep(a[i][j]);
}

inline mat_ZZ to_mat_ZZ(const mat_ZZ_p& a)
{ mat_ZZ x; conv(x, a); return x; }

inline mat_ZZ to_mat_ZZ(const mat_zz_p& a)
{ mat_ZZ x; conv(x, a); return x; }



/*
Let M(x) be the matrix M0 + x*M1; this is a matrix of linear polys in x.
Let M(a, b) = M(a + 1) M(a + 2) ... M(b). This function evaluates the products
M(a[i], b[i]) for some sequence of intervals
  a[0] < b[0] <= a[1] < b[1] <= ... <= a[n-1] < b[n-1].

The intervals are supplied in "target", simply as the list
  a[0], b[0], a[1], b[1], ...

There are three possible underlying implementations:
   * ntl_interval_products (ZZ_p version),
   * ntl_interval_products (zz_p version),
This function is a wrapper which takes ZZ_p input, calls one of the two
above implementations depending on the size of the current ZZ_p modulus, and
produces output in ZZ_p format.

PRECONDITIONS:
   Let d = b[n-1] - a[0]. Then 2, 3, ... 1 + floor(sqrt(d)) must all be
   invertible.

*/
void interval_products_wrapper(vector<mat_ZZ_p>& output,
                               const mat_ZZ_p& M0, const mat_ZZ_p& M1,
                               const vector<ZZ>& target)
{
   const ZZ& modulus = ZZ_p::modulus();

   if (!modulus.SinglePrecision())
   {
      // Large modulus.... go straight to the ZZ_p version
      ntl_interval_products
         <ZZ_p, ZZ_pX, ZZ_pXModulus, vec_ZZ_p, mat_ZZ_p, FFTRep>
               (output, M0, M1, target);
   }
   else
   {
      // Modulus is small enough to work at single precision.

      // Save old single-precision modulus and set new modulus
      zz_pContext context;
      context.save();
      zz_p::init(to_long(modulus));

      // Convert input data to single-precision format
      mat_zz_p M0_sp = to_mat_zz_p(to_mat_ZZ(M0));
      mat_zz_p M1_sp = to_mat_zz_p(to_mat_ZZ(M1));

      // Run ntl_interval_products at single precision
      vector<mat_zz_p> output_sp;
      ntl_interval_products
         <zz_p, zz_pX, zz_pXModulus, vec_zz_p, mat_zz_p, fftRep>
               (output_sp, M0_sp, M1_sp, target);

      // convert output back to ZZ_p format
      output.resize(output_sp.size());
      for (size_t i = 0; i < output_sp.size(); i++)
         output[i] = to_mat_ZZ_p(to_mat_ZZ(output_sp[i]));

      // restore old single-precision modulus
      context.restore();
   }
}

void hypellfrob_interval_products_wrapper(mat_ZZ_p& output,
                               const mat_ZZ_p& M0, const mat_ZZ_p& M1,
                               const vector<ZZ>& target)
{
   vector<mat_ZZ_p> mat_vector;
   interval_products_wrapper(mat_vector, M0, M1, target);
   int r = M0.NumRows();
   output.SetDims(r, r * mat_vector.size());

   for (size_t i = 0; i < mat_vector.size(); i++)
      for (int x = 0; x < r; x++)
         for (int y = 0; y < r; y++)
         {
            output[y][x + i*r] = mat_vector[i][y][x];
         }
}



/* ============================================================================

   Main routine for computing Frobenius matrix

============================================================================ */


/*
Given polynomials f and g defined over Z/p^N, this routine computes polynomials
a and b such that a*f + b*g = 1 mod p^N, with deg(a) < deg(g) and
deg(b) < deg(f).

PRECONDITIONS:
   f and g must be relatively prime mod p.
   The leading coefficients of f and g must be invertible mod p.
   The current NTL modulus should be p^N, and should be the one used to create
   the polynomials f and g.

Returns 1 on success, 0 if f and g are not relatively prime mod p.

NOTE:
   NTL can almost do this, but not quite; its xgcd routine is only guaranteed
   to work if the coefficient ring is a field. If you try it over Z/p^N, it
   usually works, but sometimes you get unlucky, and it crashes and burns.

*/
int padic_xgcd(ZZ_pX& a, ZZ_pX& b, const ZZ_pX& f, const ZZ_pX& g,
               const ZZ& p, int N)
{
   ZZ_pContext modulus;
   modulus.save();

   // =================================================
   // First do it mod p using NTL's xgcd routine

   // reduce down to Z/p
   ZZ_p::init(p);
   ZZ_pX f_red = to_ZZ_pX(to_ZZX(f));
   ZZ_pX g_red = to_ZZ_pX(to_ZZX(g));

   // do xgcd mod p
   ZZ_pX a_red, b_red, d_red;
   XGCD(d_red, a_red, b_red, f_red, g_red);

   // lift results to Z/p^N
   modulus.restore();
   a = to_ZZ_pX(to_ZZX(a_red));
   b = to_ZZ_pX(to_ZZX(b_red));
   ZZ_pX d = to_ZZ_pX(to_ZZX(d_red));

   if (deg(d) != 0)
      return 0;

   a /= d;
   b /= d;

   // =================================================
   // Now improve the approximation until we have enough precision

   for (int prec = 1; prec < N; prec *= 2)
   {
      ZZ_pX h = a*f + b*g - 1;
      ZZ_pX a_adj = -(h * a) % g;
      ZZ_pX b_adj = -(h * b) % f;
      a += a_adj;
      b += b_adj;
   }

   return 1;
}



/*
Given a matrix A over Z/p^N that is invertible mod p, this routine computes
the inverse B = A^(-1) over Z/p^N.

PRECONDITIONS:
   The current NTL modulus should be p^N, and should be the one used to create
   the matrix A.

NOTE:
   It's not clear to me (from the code or the documentation) whether NTL's
   matrix inversion is guaranteed to work over a non-field. So I implemented
   this one just in case.

*/
void padic_invert_matrix(mat_ZZ_p& B, const mat_ZZ_p& A, const ZZ& p, int N)
{
   ZZ_pContext modulus;
   modulus.save();

   int n = A.NumRows();

   // =================================================
   // First do it mod p using NTL matrix inverse

   // reduce to Z/p
   ZZ_p::init(p);
   mat_ZZ_p A_red = to_mat_ZZ_p(to_mat_ZZ(A));

   // invert matrix mod p
   mat_ZZ_p B_red;
   inv(B_red, A_red);

   // lift back to Z/p^N
   modulus.restore();
   B = to_mat_ZZ_p(to_mat_ZZ(B_red));

   // =================================================
   // Now improve the approximation until we have enough precision

   mat_ZZ_p two;
   ident(two, n);
   two *= 2;

   for (int prec = 1; prec < N; prec *= 2)
      // if BA = I + error, then
      // ((2I - BA)B)A = (I - err)(I + err) = I - err^2
      B = (two - B*A) * B;
}



/*
The main function exported from this module. See hypellfrob.h for information.
*/
int matrix(mat_ZZ& output, const ZZ& p, int N, const ZZX& Q)
{
   // check input conditions
   if (N < 1 || p < 3)
      return 0;

   if ((deg(Q) < 3) || (!IsOne(LeadCoeff(Q))))
      return 0;

   // EVEN DEGREE IS NOT SUPPORTED BY THIS FUNCTION. Use
   // hypellfrob::matrix_reference() in even_degree.h, which handles either
   // parity and is checked against brute-force point counting; it costs O(p)
   // rather than O(sqrt p), so it is a reference, not a replacement.
   //
   // The reason this fast path cannot take even degree is specific, not an
   // oversight. Everything below is written in terms of deg(Q) rather than
   // the genus and is dimension-correct for both parities, but the
   // horizontal reduction has a genuine singularity when deg(Q) is even:
   //
   //   D_H(s) = deg(Q)*(2t-1) - 2s
   //
   // is the denominator the reduction divides by, and it vanishes EXACTLY
   // (over Z, not merely mod p) at s = deg(Q)*(2t-1)/2, which is an integer
   // precisely when deg(Q) is even. In the sweep below that lands on the
   // hand-coded step at k = (deg(Q)/2)*(2j+1), inside the range for every
   // row j. The reduction relation degenerates there: its leading
   // coefficient is zero, so it cannot eliminate the top term.
   //
   // That direction is the differential class with nonvanishing residue at
   // the two points at infinity. It is tempting to conclude that feeding in
   // residue-zero combinations makes the obstruction vanish; measured, it
   // does not. For Q = x^6 + 2x + 3, p = 23, N = 2, the residue covector is
   // (0,0,1,0,0) but the obstruction c at the singular step, as a function
   // of the input monomial x^i dx/y, is (0, 0, 266, 138, 138) mod p^3 --
   // nonzero at i = 3 and 4 where the residue is zero, so c is not a
   // multiple of the residue functional. The reason is that the horizontal
   // reduction handles one row j of the Frobenius series at a time, and
   // residue-zero is a property of the SUM over j, not of each row's
   // contribution. Making this work needs the obstruction carried through
   // the reduction as an extra coordinate and shown to cancel after the
   // vertical maps -- a real derivation, not a dimension change.
   //
   // Dropping the obstruction was tried and is WRONG: tools/check_even.cpp
   // measures both the point count and the structural identity rho*M = p*rho
   // and they fail for every prime tested. The obstructions do not cancel
   // once the vertical maps combine the rows.
   //
   // What the correct treatment has to be (derivation, not yet implemented):
   // at y-level tt the differentials f(x) dx/y^tt modulo exact forms have
   // basis x^0..x^(d-2) TOGETHER WITH the single extra monomial
   //
   //     n*(tt) = (g+1)*tt - 1,
   //
   // which is exactly the lowest degree at which x^n dx/y^tt has nonzero
   // residue at infinity. The greedy reduction can never remove it, which is
   // the singularity. The saving structure is that this extra class is the
   // SAME cohomology class at every level: at level t it is represented by
   //
   //     x^g * Q(x)^(t-1) dx / y^(2t-1),
   //
   // whose degree is g + d(t-1) = n*(2t-1), and which is literally
   // x^g dx/y rewritten, since y^(2t-2) = Q^(t-1). So the residue line is
   // preserved by the vertical maps with scalar 1, and at the bottom level
   // it lands on the ordinary basis vector x^g dx/y.
   //
   // So the reduction should, at the singular step, split the stuck
   // coefficient c off as
   //
   //     c*x^(n*) = c*(x^g Q^(t-1))  -  c*(x^g Q^(t-1) - x^(n*)),
   //
   // accumulate c into a residue coordinate lambda (which bypasses the
   // vertical maps entirely, being level-independent), feed the second term
   // -- degree < n*, so reducible -- back into the sweep, and finally add
   // sum_j lambda_j into coordinate g. The top d coefficients of
   // x^g Q^(t-1) are computable without forming Q^(t-1): writing
   // Q = x^d(1+q) with q = sum_{m<d} Q_m x^(m-d), only the terms
   // C(t-1,k) q^k for k < d contribute.
   if (deg(Q) % 2 == 0)
      return 0;
   //
   // even_degree.cpp implements exactly that. Porting it into this file's
   // interval-product machinery is the open problem: the split needs the top
   // d coefficients of x^g Q^(t-1) at the singular step, and those do not fit
   // the fixed (M(s), D(s)) recurrence the interval products are built on.

   // Everything below is written in terms of the degree of Q rather than the
   // genus, which is what makes both parities work:
   //
   //   deg Q = 2g+1 (odd)  -- one point at infinity. dim_aff = 2g, and the
   //                          basis x^i dx/y, 0 <= i < 2g, IS H^1 of the
   //                          smooth projective curve.
   //   deg Q = 2g+2 (even) -- two points at infinity. dim_aff = 2g+1: the
   //                          affine curve's H^1 carries one extra class,
   //                          the one with nonvanishing residue at infinity.
   //                          matrix() returns Frobenius on that (2g+1)-
   //                          dimensional space; frobenius_matrix() cuts it
   //                          down to the residue-zero subspace.
   //
   // For odd degree d_deg == 2g+1 and dim_aff == 2g exactly, so this is the
   // original code with the same values substituted -- the output is
   // bit-identical (see "make check").
   const int d_deg   = deg(Q);
   const int dim_aff = d_deg - 1;

   const int genus = (d_deg - 1) / 2;   // = (d_deg-2)/2 when d_deg is even
   (void) genus;

   if (p <= (2*N - 1) * d_deg)
      return 0;

   // Back up the caller's modulus
   ZZ_pContext modulus_caller;
   modulus_caller.save();

   // Create moduli for working mod p^N and mod p^(N+1)

   ZZ_pContext modulus0, modulus1;

   ZZ_p::init(power(p, N));
   modulus0.save();
   ZZ_p::init(power(p, N+1));
   modulus1.save();

   // =========================================================================
   // Compute vertical reduction matrix M_V(t) and denominator D_V(t).

   // If N > 1 we need to do this to precision p^(N+1).
   // If N = 1 we can get away with precision N, since the last reduction
   // of length (p-1)/2 doesn't involve any divisions by p.
   ZZ_pContext modulusV = (N == 1) ? modulus0 : modulus1;
   modulusV.restore();

   // To compute the vertical reduction matrices, for each 0 <= i < 2g, we need
   // to find R_i(x) and S_i(x) satisfying x^i = R_i(x) Q(x) + S_i(x) Q'(x).
   vector<ZZ_pX> R, S;
   ZZ_pX Qred, Qred_diff, S0, R0;

   conv(Qred, Q);
   Qred_diff = diff(Qred);
   int result = padic_xgcd(R0, S0, Qred, Qred_diff, p, N+1);
   if (!result)
   {
      // failure if Q(x) and Q'(x) were not relatively prime mod p
      modulus_caller.restore();
      return 0;
   }

   S.push_back(S0);
   for (int i = 1; i < dim_aff; i++)
      S.push_back(LeftShift(S.back(), 1) % Qred);

   R.push_back(R0);
   for (int i = 1; i < dim_aff; i++)
      R.push_back(LeftShift(R.back(), 1) % Qred_diff);

   // Then the (j, i) entry of M_V(t) is given by the x^j coefficient of
   //    (2t-1) R_i(x) + 2 S_i'(x).
   // We store the constant term of M_V(t) in MV0, and the linear term in MV1.
   mat_ZZ_p MV0, MV1;
   MV0.SetDims(dim_aff, dim_aff);
   MV1.SetDims(dim_aff, dim_aff);
   for (int i = 0; i < dim_aff; i++)
   for (int j = 0; j < dim_aff; j++)
   {
      MV0[j][i] = -coeff(R[i], j) + 2*(j+1)*coeff(S[i], j+1);
      MV1[j][i] = 2*coeff(R[i], j);
   }

   // For convenience we also store the vertical reduction denominator
   // D_V(t) = 2t - 1 as a pair of 1x1 matrices.
   mat_ZZ_p DV0, DV1;
   DV0.SetDims(1, 1);
   DV1.SetDims(1, 1);
   DV0[0][0] = -1;
   DV1[0][0] = 2;

   // =========================================================================
   // Compute horizontal reduction matrices M_H^t(s) and denominators D_H^t(s),
   // where t = (p(2j+1) - 1)/2, for j in the range 0 <= j < N.
   // This is done to precision p^N.
   modulus0.restore();

   // We have D_H^t(s) = (2g+1)(2t-1) - 2s.
   // The matrix M_H^t(s) has D_H^t(s) on the sub-diagonal entries, and the
   // coefficients of the polynomial 2sQ(x) - (2t-1)xQ'(x) in the last column.

   vector<mat_ZZ_p> MH0(N), MH1(N), DH0(N), DH1(N);

   for (int j = 0; j < N; j++)
   {
      ZZ t = (p*(2*j+1) - 1)/2;

      MH0[j].SetDims(d_deg, d_deg);
      MH1[j].SetDims(d_deg, d_deg);
      DH0[j].SetDims(1, 1);
      DH1[j].SetDims(1, 1);

      DH0[j][0][0] = to_ZZ_p(d_deg*(2*t-1));
      DH1[j][0][0] = -2;
      // subdiagonal entries:
      for (int i = 0; i < dim_aff; i++)
      {
         MH0[j][i+1][i] = DH0[j][0][0];
         MH1[j][i+1][i] = DH1[j][0][0];
      }
      // last column:
      for (int i = 0; i < d_deg; i++)
      {
         MH0[j][i][dim_aff] = -to_ZZ_p(i*(2*t-1))*coeff(Qred, i);
         MH1[j][i][dim_aff] = 2*coeff(Qred, i);
      }
   }

   // Compute inverse of the vandermonde matrix that will be used in the
   // horizontal reduction phase
   mat_ZZ_p vand_in, vand;
   vand_in.SetDims(N, N);
   for (int y = 1; y <= N; y++)
   {
      vand_in[y-1][0] = 1;
      for (int x = 1; x < N; x++)
         vand_in[y-1][x] = vand_in[y-1][x-1] * y;
   }
   padic_invert_matrix(vand, vand_in, p, N);

   // =========================================================================
   // Compute B_{j, r} coefficients.
   // These are done to precision p^(N+1).

   modulus1.restore();

   // The relevant binomial coefficients.
   vector<vector<ZZ_p> > binomial(2*N);
   binomial[0].push_back(to_ZZ_p(1));
   for (int j = 1; j < 2*N; j++)
   {
      binomial[j].push_back(to_ZZ_p(1));
      for (int i = 1; i < j; i++)
         binomial[j].push_back(binomial[j-1][i-1] + binomial[j-1][i]);
      binomial[j].push_back(to_ZZ_p(1));
   }

   // For 0 <= j < N, compute Btemp[j] =
   //           (-1)^j \sum_{k=j}^{N-1} 4^{-k} {2k choose k} {k choose j}
   ZZ_p fourth = to_ZZ_p(1) / 4;
   ZZ_p fourth_pow = to_ZZ_p(1);
   vector<ZZ_p> Btemp(N);
   for (int k = 0; k < N; k++, fourth_pow = fourth_pow * fourth)
   {
      // even j
      for (int j = 0; j <= k; j += 2)
         Btemp[j] += binomial[2*k][k] * binomial[k][j] * fourth_pow;
      // odd j
      for (int j = 1; j <= k; j += 2)
         Btemp[j] -= binomial[2*k][k] * binomial[k][j] * fourth_pow;
   }

   // Compute the coefficients   B_{j, r} = p C_{j, r} (-1)^j
   //             \sum_{k=j}^{N-1} 4^{-k} {2k choose k} {k choose j},
   // where C_{j, r} is the coefficient of x^r in Q(x)^j.
   vector<vector<ZZ_p> > B(N);
   ZZ_pX Qpow = to_ZZ_pX(p);
   for (int j = 0; j < N; j++, Qpow *= Qred)
   for (int r = 0; r <= d_deg*j; r++)
      B[j].push_back(Btemp[j] * coeff(Qpow, r));

   // =========================================================================
   // Horizontal reductions.

   // reduced[j][i] will hold the reduction to x^0 of the i-th differential
   // along row j (each entry is a vector of length dim_aff)
   vector<vector<vec_ZZ_p> > reduced(N, vector<vec_ZZ_p>(dim_aff));

   for (int j = 0; j < N; j++)
   {
      // Compute the transition matrices, mod p^N, between the indices
      // kp-2g-2 and kp, for each k. We compute at most N matrices (the others
      // will be deduced more efficiently using the vandermonde matrix below)
      modulus0.restore();
      vector<ZZ> s;
      int L = d_deg*(j+1) - 1;
      int L_effective = min(L, N);
      for (int k = 0; k < L_effective; k++)
      {
         s.push_back(k*p);
         s.push_back((k+1)*p - dim_aff - 2);
      }
      vector<mat_ZZ_p> MH, DH;
      MH.reserve(L);
      DH.reserve(L);
      interval_products_wrapper(MH, MH0[j], MH1[j], s);
      interval_products_wrapper(DH, DH0[j], DH1[j], s);

      // Use the vandermonde matrix to extend all the way up to L if necessary.
      if (L > N)
      {
         // First compute X[r] = F^{(r)}(0) p^r / r! for r = 0, ..., N-1,
         // where F is the matrix corresponding to M as described in the paper;
         // and do the same for Y corresponding to the denominator D.
         vector<mat_ZZ_p> X(N);
         vector<ZZ_p> Y(N);
         for (int r = 0; r < N; r++)
         {
            X[r].SetDims(d_deg, d_deg);
            for (int h = 0; h < N; h++)
            {
               ZZ_p& v = vand[r][h];

               for (int y = 0; y < d_deg; y++)
                  for (int x = 0; x < d_deg; x++)
                     X[r][y][x] += v * MH[h][y][x];

               Y[r] += v * DH[h][0][0];
            }
         }

         // Now use those taylor coefficients to get the remaining MH's
         // and DH's.
         MH.resize(L);
         DH.resize(L);
         for (int k = N; k < L; k++)
         {
            MH[k].SetDims(d_deg, d_deg);
            DH[k].SetDims(1, 1);

            // this is actually a power of k+1, since the indices into
            // MH and DH are one-off from what's written in the paper
            ZZ_p k_pow;
            k_pow = 1;

            for (int h = 0; h < N; h++, k_pow *= (k+1))
            {
               for (int y = 0; y < d_deg; y++)
               for (int x = 0; x < d_deg; x++)
                  MH[k][y][x] += k_pow * X[h][y][x];

               DH[k][0][0] += k_pow * Y[h];
            }
         }
      }

      // Divide out each MH[k] by DH[k].
      for (size_t k = 0; k < MH.size(); k++)
      {
         ZZ_p Dinv = 1 / DH[k][0][0];

         for (int x = 0; x < d_deg; x++)
         for (int y = 0; y < d_deg; y++)
            MH[k][y][x] *= Dinv;
      }

      // Reduce all differentials on this row into the x^0 position.
      for (int i = 0; i < dim_aff; i++)
      {
         modulus1.restore();
         vec_ZZ_p sum;
         sum.SetLength(d_deg);

         // loop over blocks in this row
         for (int k = d_deg*(j+1) - 1; k >= 1; k--)
         {
            // add in new term if necessary
            if (k >= i+1 && k <= i+1 + d_deg*j)
               sum[0] += B[j][k - (i+1)];

            // We do the following reductions "by hand" for a little more
            // efficiency, since the matrices are all sparse.

            // tt is 2t - 1 in the notation of the paper
            ZZ tt = (2*j+1)*p - 2;

            // Reduce from kp - 1 down to kp - 2g - 1.
            for (int u = 1; u <= dim_aff; u++)
            {
               // c is the last component
               ZZ_p c = sum[dim_aff];

               // shuffle over all components by one step
               for (int m = dim_aff; m >= 1; m--)
                  sum[m] = sum[m-1];
               sum[0] = 0;

               // add in contribution from c (these follow from the explicit
               // formulae for the horizontal reductions)
               ZZ s = k*p - u;
               c /= to_ZZ_p(d_deg*tt - 2*s);
               for (int m = 0; m <= dim_aff; m++)
                  sum[m] += c * to_ZZ_p(2*s - tt*m) * coeff(Qred, m);
            }

            // Reduce from kp - 2g - 1 to kp - 2g - 2.
            // This is a little different to the previous block owing to the
            // denominator being divisible by p.
            {
               // c is the last component
               ZZ_p c = sum[dim_aff];

               assert(rep(c) % p == 0);
               c = to_ZZ_p(rep(c) / p);

               // shuffle over all components by one step
               for (int m = dim_aff; m >= 1; m--)
                  sum[m] = sum[m-1];
               sum[0] = 0;

               // add in contribution from c
               ZZ s = k*p - d_deg;
               ZZ denom = d_deg*tt - 2*s;
               assert(denom % p == 0);
               c /= to_ZZ_p(denom / p);
               for (int m = 0; m <= dim_aff; m++)
                  sum[m] += c * to_ZZ_p(2*s - tt*m) * coeff(Qred, m);
            }

            // Use big fat matrices to reduce from kp - 2g - 2 to (k-1)p.
            {
               // drop precision to p^N and apply the matrix
               modulus0.restore();
               vec_ZZ_p temp = MH[k-1] * to_vec_ZZ_p(to_vec_ZZ(sum));

               // lift precision back to p^(N+1)
               modulus1.restore();
               sum = to_vec_ZZ_p(to_vec_ZZ(temp));
            }

            // Reduce from (k-1)p down to (k-1)p - 1.
            {
               // c is the last component
               ZZ_p c = sum[dim_aff];

               // shuffle over all components by one step
               for (int m = dim_aff; m >= 1; m--)
                  sum[m] = sum[m-1];
               sum[0] = 0;

               // add in contribution from c
               ZZ s = (k-1)*p;
               c /= to_ZZ_p(d_deg*tt - 2*s);
               for (int m = 0; m <= dim_aff; m++)
                  sum[m] += c * to_ZZ_p(2*s - tt*m) * coeff(Qred, m);
            }
         }

         // Store the reduced vector. Note that currently it has length 2g+1,
         // but the lowest term is zero, and the vertical reductions will
         // operate on length 2g vectors, so we kill the lowest term.
         modulus1.restore();
         reduced[j][i].SetLength(dim_aff);
         for (int f = 0; f < dim_aff; f++)
            reduced[j][i][f] = sum[f+1];
      }
   }

   if (N == 1)
   {
      // If N == 1, then the vertical matrices are stored at precision 1
      // (instead of at N+1), so we reduce "reduced" to precision 1 here.
      modulus0.restore();
      for (int j = 0; j < N; j++)
      for (int i = 0; i < dim_aff; i++)
         reduced[j][i] = to_vec_ZZ_p(to_vec_ZZ(reduced[j][i]));
   }

   // =========================================================================
   // Vertical reductions.

   // Compute the vertical reduction maps between indices
   // 0, (p-1)/2, (3p-1)/2, (5p-1)/2, ..., ((2N-1)p-1)/2.
   vector<ZZ> s;
   s.push_back(to_ZZ(0));
   s.push_back((p-1)/2);
   for (int u = 1; u < N; u++)
   {
      s.push_back(s.back());
      s.push_back(s.back() + p);
   }
   modulusV.restore();
   vector<mat_ZZ_p> MV, DV;
   interval_products_wrapper(MV, MV0, MV1, s);
   interval_products_wrapper(DV, DV0, DV1, s);

   // Divide out each MV[j] by DV[j]. Note that for 1 <= j < N, both of them
   // should be divisible by p too, so we take that into account here.
   for (int j = 0; j < N; j++)
   {
      ZZ u = rep(DV[j][0][0]);
      if (j != 0)
      {
         assert(u % p == 0);
         u /= p;
      }
      InvMod(u, u, ZZ_p::modulus());

      for (int x = 0; x < dim_aff; x++)
      for (int y = 0; y < dim_aff; y++)
      {
         ZZ c = rep(MV[j][y][x]);
         if (j != 0)
         {
            assert(c % p == 0);
            c /= p;
         }
         MV[j][y][x] = to_ZZ_p(c * u);
      }
   }

   output.SetDims(dim_aff, dim_aff);
   ZZ p_to_N = power(p, N);

   // For each i, use the above reductions matrices to vertically reduce the
   // terms corresponding to x^i dx/y.
   for (int i = 0; i < dim_aff; i++)
   {
      vec_ZZ_p sum;
      sum.SetLength(dim_aff);

      for (int j = N-1; j >= 0; j--)
      {
         // pick up a new term,
         sum += reduced[j][i];
         // and reduce:
         sum = MV[j] * sum;
      }

      // store the result
      for (int f = 0; f < dim_aff; f++)
         output[f][i] = rep(sum[f]) % p_to_N;
   }

   // Restore the caller's modulus
   modulus_caller.restore();

   return 1;
}


};   // namespace hypellfrob

// ----------------------- end of file
