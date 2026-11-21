// Regression checks for the FFT matrix products used by the Euler driver.
#include <assert.h>
#include <stdio.h>
#include "rforest/hwmpz.h"
#include "rforest/rforest.h"

static void check_product(int rows, int dim, unsigned left_bits, unsigned right_bits)
{
    mpz_t *a = mpz_vec_alloc_and_init(rows*dim);
    mpz_t *b = mpz_vec_alloc_and_init(dim*dim);
    mpz_t *actual = mpz_vec_alloc_and_init(rows*dim);
    mpz_t *expected = mpz_vec_alloc_and_init(rows*dim);
    mpz_t work;
    mpz_init(work);
    for (int pass = 0; pass < 4; ++pass) {
        for (int i = 0; i < rows*dim; ++i) {
            mpz_set_ui(a[i], 17+i+pass);
            mpz_mul_2exp(a[i], a[i], left_bits+(i%7));
            mpz_add_ui(a[i], a[i], 31+i);
            if ((i+pass)%2) mpz_neg(a[i], a[i]);
            // A zero row followed by a dense row exercises transform reuse.
            if (pass == 1 && i < dim) mpz_set_ui(a[i], 0);
            if (pass == 2 && i%3) mpz_set_ui(a[i], 0);
            if (pass == 3) mpz_set_ui(a[i], 0);
        }
        for (int i = 0; i < dim*dim; ++i) {
            mpz_set_ui(b[i], 43+i);
            mpz_mul_2exp(b[i], b[i], right_bits+(i%5));
            mpz_add_ui(b[i], b[i], 7+i);
            if (i%3) mpz_neg(b[i], b[i]);
            if (pass == 2 && i%2) mpz_set_ui(b[i], 0);
            if (pass == 3) mpz_set_ui(b[i], 0);
        }
        mpz_rmatrix_mult_naive(expected, a, rows, b, dim, work);
        mpz_rmatrix_mult_fft(actual, a, rows, b, dim, work);
        for (int i = 0; i < rows*dim; ++i)
            assert(mpz_cmp(actual[i], expected[i]) == 0);
    }
    mpz_clear(work);
    mpz_vec_clear_and_free(a, rows*dim);
    mpz_vec_clear_and_free(b, dim*dim);
    mpz_vec_clear_and_free(actual, rows*dim);
    mpz_vec_clear_and_free(expected, rows*dim);
}

static void check_sparse_forest(int kappa, int differences)
{
    // Noncommuting affine transition matrices; long gaps and a duplicate
    // endpoint reproduce a precision retry after a cached prime range.
    long endpoints[] = {20000, 20003, 80000, 80000};
    unsigned long primes[] = {101, 103, 107, 109};
    long coefficients[] = {2, 1, -1, 1, 1, 0, 3, 2};
    mpz_t *polynomial = mpz_vec_alloc_and_init(8);
    mpz_t *moduli = mpz_vec_alloc_and_init(4);
    mpz_t *actual = mpz_vec_alloc_and_init(8);
    mpz_t *expected = mpz_vec_alloc_and_init(8);
    mpz_t *v = mpz_vec_alloc_and_init(2);
    mpz_t *reference = mpz_vec_alloc_and_init(2);
    mpz_t *step = mpz_vec_alloc_and_init(4);
    mpz_t *work = mpz_vec_alloc_and_init(4);
    mpz_t z;
    mpz_init(z);
    for (int i = 0; i < 8; ++i) mpz_set_si(polynomial[i], coefficients[i]);
    for (int i = 0; i < 4; ++i) mpz_ui_pow_ui(moduli[i], primes[i], 2);
    mproduct(z, moduli, 4);
    mpz_set_ui(v[0], 1); mpz_set_ui(v[1], 2);
    mpz_vec_set(reference, v, 2);
    for (long i = 0, k = 0; i < 4; ++i) {
        while (k < endpoints[i]) {
            mpz_poly_matrix_eval_si(step, polynomial, 2, 1, k++, work);
            mpz_rmatrix_mult_mod_naive(work, reference, 1, step, 2, z, work[3]);
            mpz_vec_set(reference, work, 2);
        }
        mpz_vec_mod_naive(expected+2*i, reference, 2, moduli[i]);
    }
    rforest_configure(0, differences);
    rforest(actual, v, 1, polynomial, 1, 2, moduli, 0, endpoints, 4, z, kappa);
    for (int i = 0; i < 8; ++i) assert(mpz_cmp(actual[i], expected[i]) == 0);
    assert(mpz_cmp_ui(z, 1) == 0);
    mpz_clear(z);
    mpz_vec_clear_and_free(polynomial, 8);
    mpz_vec_clear_and_free(moduli, 4);
    mpz_vec_clear_and_free(actual, 8);
    mpz_vec_clear_and_free(expected, 8);
    mpz_vec_clear_and_free(v, 2);
    mpz_vec_clear_and_free(reference, 2);
    mpz_vec_clear_and_free(step, 4);
    mpz_vec_clear_and_free(work, 4);
}

int main(void)
{
    hw_mpz_setup();
    check_product(4, 5, 65536, 512);
    check_product(4, 5, 512, 65536);
    check_product(5, 5, 32768, 32768);
    check_product(2, 7, 4096, 8192);
    check_product(1, 3, 1024, 8192);
    check_product(1, 1, 32, 64);
    check_product(3, 1, 512, 256);
    for (int differences = 0; differences <= 1; ++differences) {
        check_sparse_forest(0, differences);
        check_sparse_forest(2, differences);
        check_sparse_forest(20, differences);
    }
    hw_mpz_clear();
    puts("FFT matrix products and sparse forests agree with exact GMP arithmetic");
    return 0;
}
