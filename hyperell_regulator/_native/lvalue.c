/* Pinned, patched Arb backend. See _vendor/lfunctions/UPSTREAM.md.
 * argv: target_bits conductor factors|info required_bits cutoff degree s derivative
 * Exact integer Euler factors; exact dyadic output intervals. */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    if (argc != 9) return 2;
    int standard = atoi(argv[6]) == 5;
    int degree = standard ? 5 : 4;
    double point = strtod(argv[7], NULL);
    int derivative = atoi(argv[8]);
    if (derivative < 0 || derivative > 1) return 2;
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
    if (!standard) { mus[0]=0; mus[1]=0; mus[2]=1; mus[3]=1; }
    Lparams_t params = {0};
    params.degree = degree; params.conductor = N; params.normalisation = standard ? 1 : 0.5;
    params.mus = mus; params.target_prec = bits;
    params.self_dual = YES; params.rank = DK;
    Lerror_t error = 0;
    Lfunc_t L = Lfunc_init_advanced(&params, &error);
    if (error) { fprint_errors(stderr, error); return 3; }
    uint64_t cutoff = Lfunc_nmax(L);
    if (strtoull(argv[5], NULL, 10)) {
        uint64_t requested = strtoull(argv[5], &end, 10);
        if (*end || !requested || requested > cutoff || argv[5][0] == '-') return 2;
        if (requested < cutoff && !Lfunc_reduce_nmax(L, requested)) return 2;
        cutoff = requested;
    }
    printf("coefficient_cutoff %" PRIu64 "\n", cutoff);
    printf("target_bits %ld\nworking_bits %" PRId64 "\n", bits, Lfunc_wprec(L));
    printf("required_bits %ld\n", required);
    printf("flint_version %s\n", flint_version); fflush(stdout);
    if (!strcmp(argv[3], "info")) { Lfunc_clear(L); return 0; }
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
        for (int j = 0; j <= degree; ++j) {
            if (fscanf(input, "%8191s", text) != 1 || fmpz_set_str(coefficient, text, 10)) return 5;
            if (j == 0 && !fmpz_is_one(coefficient)) return 5;
            arb_set_fmpz(exact, coefficient); acb_set_arb(z, exact);
            acb_poly_set_coeff_acb(poly, j, z);
        }
        Lfunc_use_lpoly(L, p, poly);
    }
    if (fscanf(input, "%8191s", text) != EOF) return 5;
    fclose(input); primesieve_free(primes);
    error = Lfunc_compute_values(L);
    if (error) { fprint_errors(stderr, error); return 6; }
    const acb_struct *sign = Lfunc_sign(L);
    int eps = arb_contains_si(acb_realref(sign), 1) ? 1 : -1;
    if (!arb_contains_zero(acb_imagref(sign)) ||
        !arb_contains_si(acb_realref(sign), eps) ||
        arb_contains_si(acb_realref(sign), -eps) || (standard && eps != 1)) return 8;
    printf("root_number %d\n", eps);
    acb_t dz; acb_init(dz);
    error = Lfunc_special_value_choice(z, dz, L, point, 0.0, false, derivative == 1);
    if (derivative) acb_set(z, dz);
    acb_clear(dz);
    if (error) { fprint_errors(stderr, error); return 7; }
    if (!acb_is_finite(z) || !arb_contains_zero(acb_imagref(z))) return 8;
    interval("value", acb_realref(z));
    printf("achieved_relative_bits %ld\n", arb_rel_accuracy_bits(acb_realref(z)));
    acb_clear(z); arb_clear(exact); fmpz_clear(coefficient);
    acb_poly_clear(poly); Lfunc_clear(L);
    return 0;
}
