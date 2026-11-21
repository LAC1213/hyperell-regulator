// Copyright Edgar Costa 2026
// See LICENSE file for license details.
#include <assert.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include "sym_power.h"

// See include/sym_power.h for the math (exact fmpz arithmetic).
void sym_power_lpoly(fmpz_poly_t out, slong a_p, ulong p, int k)
{
  assert(k >= 1);
  if (k < 1) {  // out of contract; stay defined (no out-of-bounds) when NDEBUG drops the assert
    fmpz_poly_one(out);
    return;
  }

  fmpz_t pp, pk, pj, c1, t;
  fmpz_init_set_ui(pp, p);          // p
  fmpz_init(pk);
  fmpz_pow_ui(pk, pp, (ulong) k);   // p^k
  fmpz_init_set_ui(pj, 1);          // p^j, running power starting at p^0 = 1
  fmpz_init(c1);
  fmpz_init(t);

  // Lucas sequence V_m = alpha^m + beta^m: V_0 = 2, V_1 = a_p,
  // V_m = a_p*V_{m-1} - p*V_{m-2}.  (_fmpz_vec elements are passed as V + i.)
  fmpz *V = _fmpz_vec_init(k + 1);
  fmpz_set_ui(V + 0, 2);
  fmpz_set_si(V + 1, a_p);
  for (int m = 2; m <= k; m++) {
    fmpz_mul_si(V + m, V + (m - 1), a_p);   // a_p*V_{m-1}
    fmpz_mul(t, pp, V + (m - 2));           // p*V_{m-2}
    fmpz_sub(V + m, V + m, t);
  }

  fmpz_poly_t fac;
  fmpz_poly_init(fac);
  fmpz_poly_one(out);                   // out = 1

  // floor((k-1)/2)+1 = (k+1)/2 quadratic factors 1 - p^j*V_{k-2j}*T + p^k*T^2.
  for (int j = 0; j < (k + 1) / 2; j++) {
    fmpz_poly_zero(fac);
    fmpz_poly_set_coeff_ui(fac, 0, 1);
    fmpz_mul(c1, pj, V + (k - 2 * j));  // p^j*V_{k-2j}
    fmpz_neg(c1, c1);
    fmpz_poly_set_coeff_fmpz(fac, 1, c1);
    fmpz_poly_set_coeff_fmpz(fac, 2, pk);
    fmpz_poly_mul(out, out, fac);
    fmpz_mul(pj, pj, pp);               // advance p^j -> p^{j+1}
  }

  // Even k carries one extra linear factor 1 - p^{k/2}*T.
  if (k % 2 == 0) {
    fmpz_poly_zero(fac);
    fmpz_poly_set_coeff_ui(fac, 0, 1);
    fmpz_pow_ui(t, pp, (ulong) (k / 2));   // p^{k/2}
    fmpz_neg(t, t);
    fmpz_poly_set_coeff_fmpz(fac, 1, t);
    fmpz_poly_mul(out, out, fac);
  }

  fmpz_poly_clear(fac);
  _fmpz_vec_clear(V, k + 1);
  fmpz_clear(pp);
  fmpz_clear(pk);
  fmpz_clear(pj);
  fmpz_clear(c1);
  fmpz_clear(t);
}
