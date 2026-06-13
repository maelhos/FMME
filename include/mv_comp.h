#pragma once

#include <NTL/lzz_pX.h>
#include "lzz_pX_CRT.h"
#include <array>
#include "lzz_pXi.h"
#include <assert.h>
#include <lzz_p_extra.h>
#include "lzz_p_Kedlaya_Umans.h"

NTL_CLIENT
PML_CLIENT

constexpr long segmentation_dim(long deg)
{
    double d = (double)deg;
    double ld = log(d + 1);
    constexpr double sqrt7 = 2.6457513110645907f;

    return roundl(sqrt7 * sqrt(ld / log(ld)));
}

template <size_t d>
void CompModKU(zz_pX& x, const zz_pX& g, const zz_pX& h, const zz_pXModulus& f)
{
    assert(deg(g) < deg(f));
    assert(deg(h) < deg(f));
    assert(deg(f) == d);

    constexpr long n = segmentation_dim(d);
    constexpr long li = std::ceil(std::pow(d + 1, (double)1 / (double)n));
    constexpr long ln = std::ceil((d + 1) / std::pow(li, n - 1));

    //cout << li << " " << ln << " " << n << " " << d << endl;
    if constexpr ((li < 2) || (ln < 2))
    {
        cout << "li too small ! d = " << d << ", n = " << n << "\n\n";
        return; 
    }

    array<long, n - 1> degs{};
    for (long i = 0; i < n - 1; i++)
        degs[i] = li - 1;

    zz_pXi<n> G;
    from_kronecker(G, g, degs);

    zz_pX hi(h);
    Vec<zz_pX> pows{};
    pows.SetLength(n);

    long D = 0;

    for (long i = 0; i < n - 1; i++)
    {
        pows[i] = hi;
        PowerMod(hi, hi, degs[i] + 1, f);
        D += degs[i] + 1;
    }
    pows[n - 1] = hi;
    D += ln - n;

    long N = D*(d - 1) + 1;
    zz_p q{};
    element_of_order(q, 2*N);
    zz_pX_Multipoint_Geometric ev(q, N);

    Vec<Vec<zz_p>> alpha_i{};
    ev.evaluate_vector(alpha_i, pows);

    Vec<array<zz_p, n>> arr_alpha_i{};
    arr_alpha_i.SetLength(alpha_i.length());
    for (long i = 0; i < alpha_i.length(); i++)
    {
        assert(alpha_i[i].length() == n);
        arr_alpha_i[i] = *reinterpret_cast<std::array<zz_p, n>*>(alpha_i[i].data());
    }

    zz_p_Kedlaya_Umans ev_mv{};

    Vec<zz_p> rho_i{};
    Vec<ZZ> dummy{};
    ev_mv.multimodEvaluate(dummy, rho_i, G, arr_alpha_i, 1);

    zz_pX rho{};
    ev.interpolate(rho, rho_i);
    rem(x, rho, f);
}