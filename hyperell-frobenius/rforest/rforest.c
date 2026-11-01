#include <assert.h>
#include <limits.h>
#include "hwmem.h"
#include "hwmpz.h"
#include "rtree.h"
#include "rforest.h"

#define _max(a,b) ((a)>(b)?(a):(b))

// Local optimization: polynomial matrices are evaluated at consecutive
// integers throughout the forest. Forward differences replace repeated
// multiplications by exact GMP additions. The original evaluator remains
// available for arithmetic and performance comparisons.
_Thread_local int rforest_forward_differences = 1;

void rforest_configure(int disable_fft, int forward_differences)
{
    hw_disable_fft = disable_fft;
    rforest_forward_differences = forward_differences;
}

static mpz_t *transition_differences(mpz_t *M, int deg, int dim,
                                    long kbase, mpz_t *w)
{
    if (!rforest_forward_differences || kbase > LONG_MAX-deg) return NULL;
    long size = (long)dim*dim, count = size*(deg+1);
    mpz_t *values = mpz_vec_alloc_and_init(count);
    mpz_t *diff = mpz_vec_alloc_and_init(count);
    for (int j = 0; j <= deg; ++j)
        mpz_poly_matrix_eval_si(values+j*size, M, dim, deg, kbase+j, w);
    for (int order = 0; order <= deg; ++order) {
        mpz_vec_set(diff+order*size, values, size);
        for (int j = 0; j < deg-order; ++j)
            for (long i = 0; i < size; ++i)
                mpz_sub(values[j*size+i], values[(j+1)*size+i], values[j*size+i]);
    }
    mpz_vec_clear_and_free(values, count);
    return diff;
}

static inline void transition_eval(mpz_t *out, mpz_t *M, int deg, int dim,
                                   long k, mpz_t *w, mpz_t *diff)
{
    if (!diff) { mpz_poly_matrix_eval_si(out, M, dim, deg, k, w); return; }
    long size = (long)dim*dim;
    mpz_vec_set(out, diff, size);
    for (int j = 0; j < deg; ++j)
        for (long i = 0; i < size; ++i)
            if (mpz_sgn(diff[(j+1)*size+i]))
                mpz_add(diff[j*size+i], diff[j*size+i], diff[(j+1)*size+i]);
}

void rforest (mpz_t *A, mpz_t *V, int rows, mpz_t *M, int deg, int dim, mpz_t *m, long kbase, long *ks, long n, mpz_t z, int kappa)
{
    assert ( A && V && rows > 0 && M && deg >= 0 && dim > 0 && m && ks && n >= 0 && kappa >= 0 );
    if ( !n ) return;

    int ell = _ui_len(n) - kappa; // we will use <= 2^kappa trees of height ell.
    if ( ell < 0 ) ell = 0;

    hw_mem_init(0);
    hw_mpz_setup();

    // working space
    mpz_t *w = mpz_vec_alloc_and_init2 (_max(deg,rows*dim)+1, mpz_bits(z));
    mpz_t *diff = transition_differences(M, deg, dim, kbase, w);

    if ( !ell ) {

        mpz_t *Mk = mpz_vec_alloc_and_init (dim*dim);
        for ( long k = kbase, i = 0 ; i < n ; i++ ) {
            long nextk = ks[i];  assert(nextk >= k);
            while ( k < nextk ) {
                transition_eval (Mk, M, deg, dim, k++, w, diff);
                mpz_rmatrix_mult_mod_inplace (V, rows, Mk, dim, z, w);
            }
            mpz_vec_mod_hard (A+i*rows*dim, V, rows*dim, m[i]);
            mpz_divexact (z, z, m[i]);
            mpz_rmatrix_mod (V, V, rows, dim, z);
        }
        mpz_vec_clear_and_free (Mk, dim*dim);

    } else {

        mpz_t *W0 = mpz_matrix_alloc_and_init (dim), *W1 = mpz_matrix_alloc_and_init (dim);
        mpz_t **mtree = rtree_alloc (ell, 1), **Mtree = rtree_alloc (ell, dim*dim), **Rtree = rtree_alloc (ell, rows*dim);
        long i, s, k = kbase, t = 1L << ell;
        for ( s = 0 ; s <= n ; s += t ) {   // we need to insert 1 at the start of the modulus array so we have n+1 moduli
            long i0 = s;
            long i1 = i0+t;

            // set moduli for current subtree, pad with 1's if necessary, note that i shifts by 1
            if ( !i0 ) mpz_set_one(mtree[ell][0]);  // insert 1 at the begining
            for ( long i = i0?i0:1 ; i < i1 ; i++ ) if ( i <= n ) mpz_set (mtree[ell][i-i0],m[i-1]); else mpz_set_one(mtree[ell][i-i0]);

            // set leaves of matrix tree using transition matrix
            // note that i is not shifted by 1
            mpz_t *Mk = Mtree[ell];
            for ( i = i0 ; i < i1 && i < n ; i++, Mk += dim*dim ) {
                long nextk = ks[i];  assert(nextk >= k);
                if ( k == nextk ) mpz_matrix_set_one (Mk, dim); else transition_eval (Mk, M, deg, dim, k++, w, diff);
                while ( k < nextk ) {
                    mpz_matrix_set (W0, Mk, dim);
                    transition_eval (W1, M, deg, dim, k++, w, diff);
                    mpz_matrix_mult (Mk, W0, W1, dim, w[0]);
                }
            }
            for ( ; i < i1 ; i++, Mk += dim*dim ) mpz_matrix_set_one (Mk, dim); // pad with the identity matrix

            // Run the remainder tree algorithm: build/build/reduce
            rtree_build (mtree, ell, 1);  rtree_build (Mtree, ell, dim);
            rtree_reduce_rows (Rtree, V, Mtree, mtree, ell, dim, rows);

            // copy leaves of current tree to output -- use mpz_vec_mod_naive to force a hard mod (leaf values should be small)
            // note that i shifts by 1
            for ( i = i0?i0:1 ; i < i1 && i <= n ; i++ ) {
                mpz_vec_mod_naive (A+(i-1)*rows*dim, Rtree[ell] + (i-i0)*(rows*dim), rows*dim, m[i-1]);
            }
            // update the modulus and transfer vector
            mpz_divexact (z, z, mtree[0][0]);
            mpz_rmatrix_mult_mod_inplace (V, rows, Mtree[0], dim, z, w);
        }
        mpz_matrix_clear_and_free (W0, dim); mpz_matrix_clear_and_free (W1, dim);
        rtree_free (mtree, ell, 1);  rtree_free (Mtree, ell, dim*dim);  rtree_free (Rtree, ell, rows*dim);
    }

    if (diff) mpz_vec_clear_and_free(diff, (long)dim*dim*(deg+1));
    mpz_vec_clear_and_free (w, _max(deg,rows*dim)+1);
    hw_mpz_clear();
    hw_mem_clear();

    // for consistency reduce V = V mod z before returning (note mpz_rmatrix_mod does a soft reduction)
    mpz_vec_mod_naive (V, V, rows*dim, z);
}


void mproduct (mpz_t z, mpz_t *m, long n)
{
    mpz_t *w = mpz_vec_alloc_and_init (n);
    mpz_vec_set (w, m, n);
    mpz_vec_product (z, w, n);
    mpz_vec_clear_and_free (w, n);
}
