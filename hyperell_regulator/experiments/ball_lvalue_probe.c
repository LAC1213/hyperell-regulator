/* Experimental bridge to edgarcosta/lfunctions, tested at commit
 * 42014e1204f91e93f72b46945c4f5701b3ed7b1b with the adjacent
 * lfunctions_values_only.patch and lfunctions_error_bounds.patch applied.
 * Not built by this package.
 *
 * The probe leaves cache_dir NULL, so upstream recomputes its G kernel
 * without consulting an on-disk cache. Use an isolated working directory
 * to keep external-library experiments reproducible.
 * Usage: probe TARGET_BITS STANDARD_CONDUCTOR [EULER_FACTOR_FILE [REQUIRED_BITS [COEFFICIENT_CUTOFF]]]
 * TARGET_BITS controls the algorithm's internal target; REQUIRED_BITS is
 * checked against the final enclosure and defaults to TARGET_BITS.
 * Without the file, print the required coefficient cutoff and exit.
 * File format: one row per prime <= cutoff, "p c0 c1 c2 c3 c4 c5".
 * Coefficients are exact integers, in ascending order, padded with zeros.
 *
 * Enclosures are conditional on correct arithmetic input and the analytic
 * hypotheses of the upstream algorithm (including its functional equation).
 */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/acb.h>
#include <primesieve.h>
#include "glfunc.h"

extern int Lfunc_error_bounds_revision(void);

static void interval(const char *label, const arb_t x)
{
    fmpz_t lo, hi, exponent;
    fmpz_init(lo); fmpz_init(hi); fmpz_init(exponent);
    arb_get_interval_fmpz_2exp(lo, hi, exponent, x);
    printf("%s ", label);
    fmpz_print(lo); printf(" "); fmpz_print(hi); printf(" ");
    fmpz_print(exponent); printf("\n");
    fmpz_clear(lo); fmpz_clear(hi); fmpz_clear(exponent);
}

int main(int argc, char **argv)
{
    if (argc < 3 || argc > 6) return 2;
    if (Lfunc_error_bounds_revision() != 1) return 2;
    char *end;
    long bits = strtol(argv[1], &end, 10);
    if (*end || bits < 32 || bits > 4096) return 2;
    uint64_t N = strtoull(argv[2], &end, 10);
    if (*end || !N || argv[2][0] == '-') return 2;
    long required = bits;
    if (argc >= 5) {
        required = strtol(argv[4], &end, 10);
        if (*end || required < 1 || required > bits) return 2;
    }
    double mus[5] = {-1, -1, 0, 0, 1};
    Lparams_t params = {0};
    params.degree = 5; params.conductor = N; params.normalisation = 1;
    params.mus = mus; params.target_prec = bits;
    params.self_dual = YES; params.rank = DK;
    Lerror_t error = 0;
    Lfunc_t L = Lfunc_init_advanced(&params, &error);
    if (error) { fprint_errors(stderr, error); return 3; }
    uint64_t cutoff = Lfunc_nmax(L);
    if (argc == 6) {
        uint64_t requested = strtoull(argv[5], &end, 10);
        if (*end || !requested || requested > cutoff || argv[5][0] == '-') return 2;
        if (requested < cutoff && !Lfunc_reduce_nmax(L, requested)) return 2;
        cutoff = requested;
    }
    printf("coefficient_cutoff %" PRIu64 "\n", cutoff);
    printf("target_bits %ld\nworking_bits %" PRId64 "\n", bits, Lfunc_wprec(L));
    printf("required_bits %ld\n", required);
    printf("flint_version %s\n", flint_version); fflush(stdout);
    if (argc == 3) { Lfunc_clear(L); return 0; }
    FILE *input = fopen(argv[3], "r");
    if (!input) return 4;
    size_t count;
    uint64_t *primes = primesieve_generate_primes(0, cutoff, &count, UINT64_PRIMES);
    acb_poly_t poly; acb_poly_init(poly);
    fmpz_t coefficient; fmpz_init(coefficient);
    arb_t exact; arb_init(exact);
    acb_t z; acb_init(z);
    char text[8192];
    for (size_t i = 0; i < count; ++i) {
        uint64_t p;
        if (fscanf(input, "%" SCNu64, &p) != 1 || p != primes[i]) return 5;
        acb_poly_zero(poly);
        for (int j = 0; j < 6; ++j) {
            if (fscanf(input, "%8191s", text) != 1 || fmpz_set_str(coefficient, text, 10)) return 5;
            if (j == 0 && !fmpz_is_one(coefficient)) return 5;
            arb_set_fmpz(exact, coefficient); acb_set_arb(z, exact);
            acb_poly_set_coeff_acb(poly, j, z);
        }
        Lfunc_use_lpoly(L, p, poly);
    }
    if (fscanf(input, "%8191s", text) != EOF) return 5;
    fclose(input); primesieve_free(primes);
#ifdef BALL_PROBE_COMPUTE_ALL
    /* Regression control: preserve upstream's complete zero-counting path. */
    printf("computation full\n");
    error = Lfunc_compute(L);
#else
    printf("computation values_only\n");
    error = Lfunc_compute_values(L);
#endif
    if (error) { fprint_errors(stderr, error); return 6; }
    if (!arb_contains_si(acb_realref(Lfunc_sign(L)), 1)
            || !arb_contains_zero(acb_imagref(Lfunc_sign(L)))
            || arb_contains_si(acb_realref(Lfunc_sign(L)), -1)) return 8;
    error = Lfunc_special_value(z, L, 2.0, 0.0);
    if (error) { fprint_errors(stderr, error); return 7; }
    if (!acb_is_finite(z) || !arb_contains_zero(acb_imagref(z))) return 8;
    interval("L_at_2_real", acb_realref(z));
    interval("L_at_2_imag", acb_imagref(z));
    /* A real L(2) lies in its real projection; do all conversion in Arb. */
    slong precision = Lfunc_wprec(L);
    arb_t pi, scale, result;
    arb_init(pi); arb_init(scale); arb_init(result);
    arb_const_pi(pi, precision); arb_mul(scale, pi, pi, precision);
    arb_mul_2exp_si(scale, scale, 2);
    arb_set_ui(result, N); arb_sqrt(result, result, precision);
    arb_div(scale, result, scale, precision);
    arb_mul(result, scale, acb_realref(z), precision);
    if (!arb_is_finite(result) || arb_rel_accuracy_bits(result) < required) return 9;
    printf("achieved_relative_bits %ld\n", arb_rel_accuracy_bits(result));
    interval("normalization_scale", scale);
    interval("L_second_derivative", result);
    printf("L_second_decimal "); arb_printn(result, (slong)(bits * 0.302) + 10, 0); printf("\n");
    arb_clear(pi); arb_clear(scale); arb_clear(result);
    acb_clear(z); arb_clear(exact); fmpz_clear(coefficient);
    acb_poly_clear(poly); Lfunc_clear(L);
    return 0;
}
