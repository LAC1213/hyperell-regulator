#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpzfft.h"
#include "zzmem.h"
#include "hwmem.h"
#include "hwmpz.h"

#define HW_ZZ_PRIMES            4

RF_TLS int hw_disable_fft;
RF_TLS int mpzfft_threads = 1; // setting this to a value other than 1 is typically not all that helpful (better to parallelize at a higher level)

static RF_TLS zz_moduli_t zz_moduli;
static RF_TLS int mpzfft_initialized;
static RF_TLS unsigned context_depth;

static inline void hw_mpzfft_setup ()
{
    if ( hw_disable_fft || mpzfft_initialized ) return;
    zz_moduli_init (&zz_moduli, ZZ_MAX_PRIMES);
#if HW_MEM_TRACKING
    extern unsigned long zz_overhead;
    zz_overhead = zz_mem_peak;
#endif
    mpzfft_initialized = 1;
}

static inline void hw_mpzfft_clear ()
{
    if ( hw_disable_fft || ! mpzfft_initialized ) return;
    zz_moduli_clear (&zz_moduli);
    mpzfft_initialized = 0;
}

void hw_mpz_setup (void) { if (context_depth++ == 0) hw_mpzfft_setup (); }
void hw_mpz_clear (void) { assert(context_depth); if (--context_depth == 0) hw_mpzfft_clear (); }

mpz_t *mpz_vec_mod_fft (mpz_t *A, mpz_t *B, long n, mpz_t m)
{
    mpzfft_mod_t mod;
    long b, c;

    assert ( mpzfft_initialized );
    b = mpz_sizeinbase (m, 2);
    for ( long i = 0 ; i < n ; i++ ) { c = mpz_sizeinbase (B[i], 2);  if ( c > b ) b = c; }
    mpzfft_mod_init (&mod, b, m, HW_ZZ_PRIMES, &zz_moduli, mpzfft_threads);
    for ( long i = 0 ; i < n ; i++ ) { mpzfft_mod_mod (&mod, A[i], B[i], 1); }
    mpzfft_mod_clear (&mod);
    return A;
}

mpz_t *mpz_vec_mod_init_fft (mpz_t *A, mpz_t *B, long n, mpz_t m, mpz_t w)
{
    mpzfft_mod_t mod;
    long b, c;

    assert ( mpzfft_initialized );
    b = mpz_sizeinbase (m, 2);
    for ( long i = 0 ; i < n ; i++ ) { c = mpz_sizeinbase (B[i], 2);  if ( c > b ) b = c; }
    mpzfft_mod_init (&mod, b, m, HW_ZZ_PRIMES, &zz_moduli, mpzfft_threads);
    for ( long i = 0 ; i < n ; i++ ) { mpzfft_mod_mod (&mod, w, B[i], 1);  mpz_init_set (A[i], w); }
    mpzfft_mod_clear (&mod);
    return A;
}

static inline size_t mpz_rmatrix_product_max_bits (mpz_t *A, int r, mpz_t *B, int d)
{
    size_t left = 0, right = 0;
    for (int i = 0; i < r*d; ++i) {
        size_t n = mpz_size(A[i]); if (n > left) left = n;
    }
    for (int i = 0; i < d*d; ++i) {
        size_t n = mpz_size(B[i]); if (n > right) right = n;
    }
    // Each dot product has d terms. Bound its absolute value by
    // d * 2^(bits(A)+bits(B)), including signs and the addition carry.
    // Padding both inputs to the larger size wastes large transforms in
    // the unbalanced prefix updates of a remainder forest.
    size_t carry = 0;
    for (unsigned n = d-1; n; n >>= 1) ++carry;
    size_t bits = GMP_NUMB_BITS*(left+right)+carry;
    return bits ? bits : 1; // the FFT parameter builder requires bits >= 1
}

mpz_t *mpz_rmatrix_mult_fft (mpz_t *C, mpz_t *A, int r, mpz_t *B, int d, mpz_t w)
{
    mpzfft_params_t params;
    mpzfft_t *AT, *BT;
    assert ( mpzfft_initialized );
    
    mpzfft_params_init (&params, mpz_rmatrix_product_max_bits (A, r, B, d), d, HW_ZZ_PRIMES, &zz_moduli);

    // Reuse B's transforms and stream A a row at a time. In particular,
    // a 5x5 product retains 30 transforms instead of 50. The public
    // multiplication contract forbids C aliasing either input.
    BT = hw_malloc(d*d*sizeof(mpzfft_t));
    for (int i = 0; i < d*d; ++i) {
        mpzfft_init(BT[i], &params); mpzfft_fft(BT[i], B[i], mpzfft_threads);
    }
    AT = hw_malloc(d*sizeof(mpzfft_t));
    for (int i = 0; i < d; ++i) mpzfft_init(AT[i], &params);
    for (int row = 0; row < r; ++row) {
        for (int i = 0; i < d; ++i) mpzfft_fft(AT[i], A[row*d+i], mpzfft_threads);
        mpzfft_matrix_mul(AT, AT, BT, 1, d, d, mpzfft_threads);
        for (int i = 0; i < d; ++i) mpzfft_ifft(C[row*d+i], AT[i], mpzfft_threads);
    }
    for (int i = 0; i < d; ++i) mpzfft_clear(AT[i]);
    for (int i = 0; i < d*d; ++i) mpzfft_clear(BT[i]);
    hw_free(AT, d*sizeof(mpzfft_t));
    hw_free(BT, d*d*sizeof(mpzfft_t));

    mpzfft_params_clear (&params);
    return C;
}
