// Copyright Dave Platt 2024
// See LICENSE file for license details.
//
// Use Booker's generalisation of Turing's method to confirm RH
// to height B/(OUTPUT_RATIO*degree).
// see Artin’s Conjecture, Turing’s Method, and the Riemann Hypothesis
// B is hard wired to 512 in g.c so height is 64/degree

// Note that RH is only verifiable for rank <=1
// This code will trust larger values of rank and will
// try to verify RH on the assumption rank is correct.

#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C" {
#endif
#ifdef TURING

// We use N(t)=Phi(t)+S(t) (4-2)
//
// with
//
// Pi*Phi(t)=arg(eps)+log(N)*t/2-r*t*log(Pi)/2+Im sum_j=1^r
// loggamma((1/2+iy+mu_j)/2)
//
// then integrate from t0 to t0+h using both the l-function and its conjugate
// so we can ignore arg of the root number

// int im loggamma(sigma+It) t=a..b
// use Im loggamma(z)=Im [(z-1/2)log(z/e)] + theta(1/(8 |Im z|)) (p 393)
bool imint(arb_t res, arb_t sigma, arb_t a, arb_t b, uint64_t prec) {
  if (verbose) {
    printf("Integrating with sig=");
    arb_printd(sigma, 10);
    printf(" a=");
    arb_printd(a, 10);
    printf(" and b=");
    arb_printd(b, 10);
    printf("\n");
  }
  static arb_t a2, b2, s2, b2s2, a2s2, tmp1, tmp2, tmp3, shalf;
  static bool init = false;
  if (!init) {
    init = true;
    arb_init(a2);
    arb_init(b2);
    arb_init(s2);
    arb_init(b2s2);
    arb_init(a2s2);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);
    arb_init(shalf);
  }
  arb_set_d(tmp1, -0.5);
  arb_add(shalf, sigma, tmp1, prec); // sigma-1/2

  arb_mul(s2, sigma, sigma, prec);
  arb_mul(a2, a, a, prec);
  arb_mul(b2, b, b, prec);
  arb_add(tmp1, a2, s2, prec);
  arb_log(a2s2, tmp1, prec);
  arb_add(tmp1, b2, s2, prec);
  arb_log(b2s2, tmp1, prec);
  arb_sub(tmp1, s2, sigma, prec);
  arb_sub(tmp2, tmp1, b2, prec);
  arb_mul(tmp3, tmp2, b2s2, prec); // log(b^2+sig^2)(sig^2-b^2-sig)
  arb_sub(tmp2, tmp1, a2, prec);
  arb_mul(tmp1, tmp2, a2s2, prec); // log(a^2+sig^2)(sig^2-a^2-sig)
  arb_sub(tmp2, tmp3, tmp1, prec);
  arb_sub(tmp1, b2, a2, prec);
  arb_mul_ui(tmp3, tmp1, 3, prec);
  arb_add(res, tmp2, tmp3, prec);
  arb_mul_2exp_si(res, res, -2);
  if (verbose) {
    printf("value 1 = ");
    arb_printd(res, 10);
    printf("\n");
  }

  arb_div(tmp1, b, sigma, prec);
  arb_atan(tmp2, tmp1, prec);
  arb_mul(tmp1, tmp2, b, prec); // b atan(b/sig)
  arb_div(tmp2, a, sigma, prec);
  arb_atan(tmp3, tmp2, prec);
  arb_mul(tmp2, tmp3, a, prec);    // a atan(a/sig)
  arb_sub(tmp3, tmp2, tmp1, prec); // a atan(a/sigma)-b atan(b/sigma)
  arb_mul(tmp2, tmp3, shalf, prec);
  if (verbose) {
    printf("value 2 = ");
    arb_printd(tmp2, 10);
    printf("\n");
  }

  arb_add(res, res, tmp2, prec);
  arb_neg(res, res);
  arb_div(tmp1, b, a, prec);
  arb_log(tmp2, tmp1, prec);
  arb_mul_2exp_si(tmp2, tmp2, -3);
  arb_add_error(res, tmp2);
  if (verbose) {
    printf("im_int returning ");
    arb_printd(res, 10);
    printf("\n");
  }
  return (true);
}

// Q(s) is analytic conductor defined pg 387 col 2
void logQ(arb_t res, acb_t s, Lfunc *L, int64_t prec) {
  // printf("In logQ with s=");acb_printd(s,10);printf("\n");
  static bool init = false;
  static arb_t two_pi, tmp1, tmp2;
  static acb_t stmp1, stmp2, stmp3;
  if (!init) {
    init = true;
    arb_init(two_pi);
    arb_init(tmp1);
    arb_init(tmp2);
    acb_init(stmp1);
    acb_init(stmp2);
    acb_init(stmp3);
    arb_const_pi(two_pi, prec);
    arb_mul_2exp_si(two_pi, two_pi, 1);
  }
  acb_set_ui(stmp1, L->conductor);
  arb_set(acb_imagref(stmp2), acb_imagref(s));
  for (uint64_t j = 0; j < L->degree; j++) {
    arb_set_d(tmp2, L->mus[j]);
    arb_add(acb_realref(stmp2), acb_realref(s), tmp2, prec);
    acb_mul(stmp3, stmp1, stmp2, prec);
    acb_div_arb(stmp1, stmp3, two_pi, prec);
  }
  if (verbose) {
    printf("Analytic conductor = ");
    acb_printd(stmp1, 10);
    printf("\n");
  }
  acb_abs(tmp1, stmp1, prec);
  arb_log(res, tmp1, prec);
  if (verbose) {
    printf("LogQ returning ");
    arb_printd(res, 10);
    printf("\n");
  }
}

// see 4-10.
// Note:
//  1. X must be > 5.
//  2. Im mu_j = 0 here
bool set_X(arb_t res, uint64_t r, double *mus, double one_over_B,
           int64_t prec) {
  double max_mu = mus[0];
  for (uint64_t j = 1; j < r; j++)
    if (mus[j] > max_mu)
      max_mu = mus[j];
  if (verbose)
    printf("max mu_j = %10.8e\n", max_mu);
  max_mu += 2.5;
  arb_t tmp, max_t;
  arb_init(tmp);
  arb_init(max_t);
  arb_set_d(tmp, OUTPUT_RATIO * one_over_B);
  arb_inv(max_t, tmp, prec);
  arb_sqr(max_t, max_t, prec);
  arb_set_d(tmp, max_mu * max_mu);
  arb_sub(max_t, max_t, tmp, prec);
  arb_sub_ui(tmp, max_t, 25, prec);
  if (!(arb_is_positive(tmp))) {
    fprintf(stderr, "Error setting X (eq 4-10). A mu_j was too large doing "
                    "Turing's method.\n");
    arb_clear(tmp);
    arb_clear(max_t);
    return false;
  }
  arb_sqrt(res, max_t, prec);
  arb_clear(tmp);
  arb_clear(max_t);
  return true;
}

// return bound for 2\int\limits_{t_0}^{t0+h} S(t) \dif t
// see Th 4.6
// see Remark 4.2 of https://arxiv.org/pdf/2508.03023 PALOJARVI + ZHAO
// times 2 for l-function and its conjugate.
bool St_int(arb_t res, arb_t h, arb_t t0, Lfunc *L, int64_t prec) {
  static bool init = false;
  static acb_t s;
  static arb_t Q1, Q2, l2_half, pi, rc_theta_etc, tmp, tmp1;
  if (!init) {
    init = true;
    arb_init(tmp);
    arb_init(tmp1);
    acb_init(s);
    arb_init(pi);
    arb_set_d(acb_realref(s), 1.5);
    arb_init(Q1);
    arb_init(Q2);
    arb_init(l2_half);
    arb_log_ui(Q1, 2, prec);
    arb_set_d(Q2, -0.5);
    arb_add(l2_half, Q1, Q2, prec);
    arb_init(rc_theta_etc);
    arb_const_pi(pi, prec);
  }
  // rc_theta_etc depends on L (via L->X and L->degree), so it must be
  // recomputed for every L-function, not cached in the init block.
  // c_theta for theta=0 (Booker Thm 4.6). The paper's printed bound 5.65055 is
  // loose by exactly log zeta(3/2)=0.96026; the displayed-expression value is
  // 4.69028764131 (verified, and equal to Palojarvi-Zhao c_0(eps=1/2)). 4.6902877
  // rounds UP, so it stays a valid, strictly-tighter upper bound.
  arb_set_d(Q2, 4.6902877);
  arb_set_ui(tmp1, 8);
  arb_div_ui(tmp1, tmp1, 10, prec);              // 0.8 (exact)
  arb_sub_ui(tmp, L->X, 5, prec);                  // X-5
  arb_div(Q1, tmp1, tmp, prec);                    // 0.8/(X-5)
  arb_add(Q1, Q1, Q2, prec);                       // c_theta + 0.8/(X-5)
  arb_mul_ui(rc_theta_etc, Q1, L->degree, prec);   // c_\theta r+0.8 r/(X-5)
  if (verbose) {
    printf("c theta bit = ");
    arb_printd(rc_theta_etc, 10);
    printf("\n");
  }
  arb_set(acb_imagref(s), t0); // 3/2+it0
  logQ(Q2, s, L, prec);
  arb_mul(tmp, Q2, l2_half, prec);
  arb_add(acb_imagref(s), acb_imagref(s), h, prec); // 3/2+it1
  logQ(Q1, s, L, prec);
  arb_mul_2exp_si(Q1, Q1, -2);
  arb_add(tmp1, tmp, Q1, prec);
  arb_add(tmp, tmp1, rc_theta_etc, prec);
  arb_div(tmp1, tmp, pi, prec);
  arb_mul_2exp_si(tmp1, tmp1, 1); // this is an upper bound
  arb_neg(tmp, tmp1);             // lower bound
  arb_union(res, tmp, tmp1, prec);
  return true;
}

// because our zeros are held as intervals, this computes N_left and N_right
// simultaneously (see 4 of Booker)
bool turing_int(arb_t res, arb_t t0, arb_t h, Lfunc *L, uint64_t *count,
                uint64_t side, int64_t prec) {
  static bool init = false;
  static arb_t tmp, t0_plus_h;
  if (!init) {
    init = true;
    arb_init(tmp);
    arb_init(t0_plus_h);
  }

  arb_add(t0_plus_h, t0, h, prec);
  arb_zero(res);
  uint64_t z_ptr = 0, N_t0 = 0;
  bool below_t0 = true;
  while (true) {
    if (z_ptr >= MAX_ZEROS) {
      fprintf(stderr,
              "Ran out of room for zeros before end of Turing Region.\n");
      return false;
    }
    if (arb_is_zero(L->zeros[side][z_ptr])) // end of zeros
      break;
    if (below_t0) {
      arb_sub(tmp, L->zeros[side][z_ptr], t0, prec);
      if (arb_is_negative(tmp)) {
        N_t0++;
        z_ptr++;
        continue;
      } else
        below_t0 = false;
    }
    arb_sub(tmp, t0_plus_h, L->zeros[side][z_ptr], prec);
    arb_add(res, res, tmp, prec);
    z_ptr++;
  }
  count[0] = N_t0;
  return true;
}

// computes Turing's integral into res and returns # zeros actually found.
uint64_t turing_count(arb_t res, Lfunc *L, int64_t prec) {
  static bool init = false;
  static arb_t tint, tmp, tmp1, h, t0, sint, pi, rlogpi, t0hbit;
  if (!init) {
    init = true;
    arb_init(tint);
    arb_init(sint);
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(h);
    arb_init(t0);
    arb_init(pi);
    arb_const_pi(pi, prec);
    arb_init(rlogpi);
    arb_init(t0hbit);
  }
  // h, t0, rlogpi and t0hbit all depend on L (via L->B = 512/degree and
  // L->degree), so they must be recomputed for every L-function rather than
  // cached once in the init block.
  arb_div_ui(h, L->B, TURING_RATIO, prec);
  if (verbose) {
    printf("Turing h set to ");
    arb_printd(h, 10);
    printf("\n");
  }
  arb_div_ui(t0, L->B, OUTPUT_RATIO, prec);
  if (verbose) {
    printf("Turing t0 set to ");
    arb_printd(t0, 10);
    printf("\n");
  }
  arb_log(tmp, pi, prec);
  arb_mul_ui(rlogpi, tmp, L->degree, prec);
  arb_mul(tmp, t0, h, prec);
  arb_mul_2exp_si(tmp, tmp, 1);
  arb_mul(tmp1, h, h, prec);
  arb_add(sint, tmp1, tmp, prec);
  arb_div(t0hbit, sint, pi, prec);
  arb_mul_2exp_si(t0hbit, t0hbit, -1); // (2t0h+h^2)/2Pi
  uint64_t zeros_found = 0;
  if (!turing_int(tint, t0, h, L, &zeros_found, 0, prec))
    return 0;
  if (verbose) {
    printf("Turing int [0] = ");
    arb_printd(tint, 20);
    printf("\n");
  }
  if (L->self_dual == YES) {
    zeros_found <<= 1;
    arb_mul_2exp_si(tint, tint, 1);
  } else {
    uint64_t conj_zeros_found;
    if (!turing_int(tmp, t0, h, L, &conj_zeros_found, 1, prec))
      return 0;
    zeros_found += conj_zeros_found;
    arb_add(tint, tint, tmp, prec);
  }
  zeros_found +=
      L->rank; // not rigorous if rank > 1 unless it can be proven otherwise
  if (verbose) {
    printf("Turing int [*] = ");
    arb_printd(tint, 20);
    printf("\n");
  }
  if (!St_int(sint, h, t0, L, prec))
    return 0;
  if (verbose) {
    printf("St_int returned ");
    arb_printd(sint, 10);
    printf("\n");
  }

  arb_div(tmp1, L->imint, pi, prec);
  arb_add(tmp, tmp1, sint, prec);
  arb_sub(sint, tmp, tint, prec);

  arb_log_ui(tmp1, L->conductor, prec);
  arb_sub(tmp, tmp1, rlogpi, prec); // log(N)-tlog Pi
  arb_mul(tmp1, tmp, t0hbit, prec); //

  arb_add(tmp, tmp1, sint, prec);
  arb_div(res, tmp, h, prec);

  if (verbose) {
    printf("turing_count will try to return ");
    arb_printd(res, 10);
    printf("\n");
  }
  return zeros_found;
}

Lerror_t turing_check_RH(Lfunc *L, int64_t prec) {
  static bool init = false;
  static arb_t tmp, sigma, a, b, t0, h, tcount;
  if (!init) {
    init = true;
    arb_init(tmp);
    arb_init(sigma);
    arb_init(a);
    arb_init(b);
    arb_init(t0);
    arb_init(h);
    arb_init(tcount);
  }

  arb_zero(L->imint);
  arb_div_ui(t0, L->B, OUTPUT_RATIO, prec);
  arb_div_ui(h, L->B, TURING_RATIO, prec);
  arb_set(a, t0);
  arb_add(b, t0, h, prec);
  if (verbose) {
    printf("t0=");
    arb_printd(t0, 20);
    printf("\nh=");
    arb_printd(h, 20);
    printf("\n");
  }
  arb_mul_2exp_si(a, a, -1); // /2 for change of variable in integral
  arb_mul_2exp_si(b, b, -1);

  // compute sum over j=1 to degree of int of Im(lngamma((1/2+it+mu_j)/2)) from
  // t0 to t0+h
  for (uint64_t r = 0; r < L->degree; r++) {
    arb_set_d(sigma, L->mus[r] / 2.0 + 0.25);
    if (!imint(tmp, sigma, a, b, prec))
      return ERR_RH_ERROR;
    arb_add(L->imint, L->imint, tmp, prec);
  }
  arb_mul_2exp_si(L->imint, L->imint, 2);
  if (verbose) {
    printf("imint=");
    arb_printd(L->imint, 20);
    printf("\n");
  }

  if (!set_X(L->X, L->degree, L->mus, L->one_over_B, prec))
    return ERR_RH_ERROR;
  if (verbose) {
    printf("X set to ");
    arb_printd(L->X, 20);
    printf("\n");
  }

  uint64_t zeros_found = turing_count(tcount, L, prec);
  if (verbose) {
    printf("Turing Integral returned ");
    arb_printd(tcount, 20);
    printf("\n");
  }
  if (zeros_found == 0)
    return ERR_RH_ERROR;

  if (verbose)
    printf("We found %lu zeros.\n", zeros_found);

  // turing count should now bracket zeros_found

  arb_sub_ui(tmp, tcount, zeros_found - 1, prec);
  if (!arb_is_positive(tmp)) {
    fprintf(stderr, "Found too many zeros.\n");
    return ERR_RH_ERROR;
  }
  arb_sub_ui(tmp, tmp, 2, prec);
  if (!arb_is_negative(tmp)) {
    fprintf(stderr, "Found too few zeros.\n");
    return ERR_RH_ERROR;
  }

  return 0; // success
}

#ifdef __cplusplus
}
#endif

#endif
