/* Regression checks for lfunctions_error_bounds.patch. Run in a fresh process
 * at each precision: error-test [target_bits], default 80. */
#include <assert.h>
#include <stdlib.h>
#include <flint/acb_dirichlet.h>
#include "glfunc.h"
#include "glfunc_internals.h"

static void factors(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *data)
{
    (void)d; (void)prec; (void)data;
    long a = p == 2 ? 0 : (p % 4 == 1 ? 1 : -1);
    long b = p == 3 ? 0 : (p % 3 == 1 ? 1 : -1);
    acb_poly_zero(poly);
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, 1, -a-b);
    acb_poly_set_coeff_si(poly, 2, a*b);
}

int main(int argc, char **argv)
{
    long bits = argc > 1 ? strtol(argv[1], NULL, 10) : 80;
    double mus[2] = {1, 1};
    Lparams_t p = {0};
    p.degree = 2; p.conductor = 12; p.mus = mus;
    p.target_prec = bits; p.self_dual = YES; p.rank = DK;
    Lerror_t error;
    Lfunc_t handle = Lfunc_init_advanced(&p, &error);
    assert(error == 0);
    error = Lfunc_use_all_lpolys(handle, factors, NULL);
    assert(error == 0);
    error = Lfunc_compute_values(handle);
    assert(error == 0);
    acb_t got, expected, s, a, b;
    acb_init(got); acb_init(expected); acb_init(s); acb_init(a); acb_init(b);
    acb_set_ui(s, 2);
    error = Lfunc_special_value(got, handle, 2, 0);
    assert(error == 0);
    dirichlet_group_t g3, g4;
    dirichlet_group_init(g3, 3); dirichlet_group_init(g4, 4);
    dirichlet_char_t c3, c4;
    dirichlet_char_init(c3, g3); dirichlet_char_init(c4, g4);
    dirichlet_char_log(c3, g3, 2); dirichlet_char_log(c4, g4, 3);
    acb_dirichlet_l(a, s, g3, c3, bits + 200);
    acb_dirichlet_l(b, s, g4, c4, bits + 200);
    acb_mul(expected, a, b, bits + 200);
    assert(acb_contains(got, expected));
    printf("PASS: special value contains independent FLINT enclosure (target %ld, achieved %ld bits)\n", bits, acb_rel_accuracy_bits(got));
    fflush(stdout);
    assert(acb_rel_accuracy_bits(got) >= bits/2);
    /* Independently test the completed-to-uncompleted derivative conversion,
       both at an interpolation node (the centre) and away from it. */
    acb_ptr jet3 = _acb_vec_init(2), jet4 = _acb_vec_init(2);
    acb_t derivative, product;
    acb_init(derivative); acb_init(product);
    for (int centre=0; centre<2; centre++) {
        double point = centre ? 0.5 : 2.0;
        acb_set_d(s, point);
        error = Lfunc_special_value_choice(got, derivative, handle, point, 0, false, true);
        assert(error == 0);
        acb_dirichlet_l_jet(jet3, s, g3, c3, 0, 2, bits+200);
        acb_dirichlet_l_jet(jet4, s, g4, c4, 0, 2, bits+200);
        acb_mul(expected, jet3+1, jet4, bits+200);
        acb_mul(product, jet3, jet4+1, bits+200);
        acb_add(expected, expected, product, bits+200);
        assert(acb_contains(derivative, expected));
        assert(acb_rel_accuracy_bits(derivative) >= bits/2);
    }
    printf("PASS: first derivatives contain independent FLINT jets, including central node\n");
    _acb_vec_clear(jet3,2); _acb_vec_clear(jet4,2);
    acb_clear(derivative); acb_clear(product);
    Lfunc_clear(handle);

    // Inflate eq59 deliberately. Every retained Fourier coefficient must
    // receive this uncertainty before the inverse FFT, not just res[i].
    handle = Lfunc_init_advanced(&p, &error);
    assert(error == 0);
    (void)Lfunc_nmax(handle);
    Lfunc *L = (Lfunc *)handle;
    for(uint64_t j=0; j<L->fft_NN; ++j) acb_set_ui(L->res[j], 10);
    arb_one(L->sum_ans); arb_one(L->eq59); arb_mul_2exp_si(L->eq59, L->eq59, -20);
    error = do_pre_iFFT_errors(L);
    assert(error == 0);
    for(uint64_t j=0; j<16; ++j)
        assert(mag_cmp_2exp_si(arb_radref(acb_realref(L->res[j])), -20) >= 0);
    printf("PASS: all checked Fourier entries received eq59 uncertainty\n");
    Lfunc_clear(handle);
    return 0;
}
