#include "lzz_pX_Multipoint_Field.h"
#include <assert.h>
#include "lzz_pX_extra_extra.h"
 
zz_pX_Multipoint_Field::zz_pX_Multipoint_Field() : r((zz_p::modulus() - 1) / 2)
{
    element_of_order(q, zz_p::modulus() - 1);
    assert(q.LoopHole());
    qs = q*q;

    ev_even = zz_pX_Multipoint_Geometric(q, to_zz_p(1), r);
    ev_odd = zz_pX_Multipoint_Geometric(q, q, r);
}

/*
Complexity : O~(p + deg(f))
*/
void zz_pX_Multipoint_Field::evaluate(Vec<zz_p>& val, const zz_pX& f) const
{
    long p = zz_p::modulus();
    val.SetLength(p);
    long d = deg(f);

    if (d < 1)
    {
        zz_p c = (d == 0) ? f.rep[0] : to_zz_p(0);

        for (long i = 0; i < p; i++)
            val[i] = c;
        return;
    }

    if (d <= 64)
    {
        evaluateWholeFieldNaive(val, f);
        return;
    }

    zz_pX f_even(f);
    zz_pX f_odd(f);

    for (long i = d; i >= r; i--)
    {
        f_even.rep[i - r] += f_even.rep[i];
        f_odd.rep[i - r] -= f_odd.rep[i];
    }

    long nr = min(r, d + 1);
    f_even.SetLength(nr); f_even.normalize();
    f_odd.SetLength(nr); f_odd.normalize();

    vec_zz_p v_even, v_odd;

    val[0] = f.rep[0];
    ev_even.evaluate(v_even, f_even);
    ev_odd.evaluate(v_odd, f_odd);

    zz_p v(1);

    for (long i = 0; i < r; i++)
    {
        val[v.LoopHole()] = v_even[i];
        val[(v * q).LoopHole()] = v_odd[i];
        v *= qs;
    }
}

void evaluateWholeFieldNaive(Vec<zz_p>& val, const zz_pX& pol)
{
    long p = zz_p::modulus();
    val.SetLength(p);

    for (long i = 0; i < p; i++)
        val[i] = eval(pol, to_zz_p(i));
}

void zz_pX_Multipoint_Field::reduceFermat(zz_pX& a, const zz_pX& b)
{
    long p = zz_p::modulus();
    // the NTL makes the asumptions that zz_p are reduced which can be false
    // here since we just changed moduli, we reinstantiate in a safe way
    a.rep.SetLength(b.rep.length());
    for (long i = 0; i < b.rep.length(); i++)
        a.rep[i] = to_zz_p(b.rep[i]._zz_p__rep);
    a.normalize(); // new moduli = maybe lead term is 0...

    long d = deg(a);

    for (long i = d; i >= p; i--)
        a.rep[i - p + 1] += a.rep[i];

    a.rep.SetLength(min(p, d + 1));
    a.normalize();
}

void zz_pX_Multipoint_Field::reduceFermat(zz_pXY& a, const zz_pXY& b)
{
    long p = zz_p::modulus();
    a.rep.SetLength(b.rep.length());
    for (long i = 0; i < a.rep.length(); i++)
        zz_pX_Multipoint_Field::reduceFermat(a.rep[i], b.rep[i]);

    long d = a.degY();

    for (long i = d; i >= p; i--)
        a.rep[i - p + 1] += a.rep[i];

    a.rep.SetLength(min(p, d + 1));
    a.normalize();
}

/*
Complexity : O~(degY(pol) * p + degY(pol)*degX(pol) + p^2)
*/
void zz_pX_Multipoint_Field::evaluateReduced(Vec<Vec<zz_p>>& val, const zz_pXY& pol) const
{
    long p = zz_p::modulus();
    assert(pol.degX() < p && pol.degY() < p);

    Vec<Vec<zz_p>> t_out;

    long l = pol.degY() + 1;
    t_out.SetLength(l);
    val.SetLength(p);

    for (long i = 0; i < l; i++)
        evaluate(t_out[i], pol.rep[i]);

    for (long x = 0; x < p; x++)
    {
        zz_pX tp;
        tp.rep.SetLength(l);
        for (long i = 0; i < l; i++)
            tp.rep[i] = t_out[i][x];
        tp.normalize();

        evaluate(val[x], tp);
    }
}