// Copyright Edgar Costa 2026
// See LICENSE file for license details.
#ifndef _SYM_POWER_INCLUDE
#define _SYM_POWER_INCLUDE

#include <flint/fmpz_poly.h>

#ifdef __cplusplus
extern "C" {
#endif

  /* Form the Sym^k good-prime Euler factor of an elliptic curve from its a_p.
   *
   * At a prime p of GOOD reduction the Satake parameters alpha, beta satisfy
   * alpha + beta = a_p and alpha*beta = p.  The local factor of the k-th
   * symmetric power is the integer polynomial of degree k+1
   *     out(T) = prod_{i=0}^{k} (1 - alpha^i beta^{k-i} T),
   * computed without complex roots from the Lucas sequence
   *     V_0 = 2, V_1 = a_p, V_m = a_p*V_{m-1} - p*V_{m-2}   (V_m = alpha^m+beta^m):
   *     out = prod_{j=0}^{floor((k-1)/2)} (1 - p^j*V_{k-2j}*T + p^k*T^2)
   *           * (1 - p^{k/2}*T)        [the extra linear factor only when k even].
   * With k=1 this is the curve's own factor 1 - a_p*T + p*T^2.
   *
   * Coefficients are exact and routinely exceed 64 bits (the leading one has
   * absolute value p^{k(k+1)/2}), hence the fmpz_poly output.  Requires k >= 1
   * (k < 1 yields the trivial factor 1 rather than reading out of bounds).
   *
   * GOOD PRIMES ONLY: this is a pure function of (a_p, p, k).  Bad/ramified
   * factors are not a function of a_p (e.g. Sym^4 of 27.a at p=3 is 1 - p^2*T,
   * not derivable from a_3 = 0) and are the caller's responsibility.
   *
   * This is the single Sym^k factor implementation: it is consumed by
   * examples/ec_sym.cpp and by the high-degree regression driver
   * (test/highdeg_check.cpp), and is verified against Pari lfunsympow/lfuneuler
   * for k = 1..8 (see test/sympow_lpoly_test.c). */
  void sym_power_lpoly(fmpz_poly_t out, slong a_p, ulong p, int k);

#ifdef __cplusplus
}
#endif
#endif
