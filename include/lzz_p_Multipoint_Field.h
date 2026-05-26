#pragma once

#include <NTL/lzz_pX.h>
#include "lzz_pX_CRT.h"
#include "lzz_pXY.h"
#include "lzz_pXi.h"
#include "NVec.h"
#include <array>
#include <assert.h>


NTL_CLIENT
PML_CLIENT

void reduceFermat(zz_pX& a, const zz_pX& b);
void reduceFermat(zz_pXY& a, const zz_pXY& b);
void evaluateWholeFieldNaive(Vec<zz_p>& val, const zz_pX& pol);


class zz_p_Multipoint_Field
{
private:
    long r;
    zz_p q;
    zz_p qs;
    zz_pX_Multipoint_Geometric ev_even;
    zz_pX_Multipoint_Geometric ev_odd;
public:
    zz_p_Multipoint_Field();
    ~zz_p_Multipoint_Field() {}

    void evaluate(Vec<zz_p>& val, const zz_pX& pol) const;
    void evaluateReduced(Vec<Vec<zz_p>>& val, const zz_pXY& pol) const; // assume degX < p, degY < p
    inline void evaluate(Vec<Vec<zz_p>>& val, const zz_pXY& pol) const
    {
        zz_pXY r_pol;
        reduceFermat(r_pol, pol);
        evaluateReduced(val, r_pol);
    }

    template <size_t K>
    void _aggregate_evals(NVec<K+1>& val, Vec<NVec<K>>& t_out) const
    {
        long p = zz_p::modulus();
        long l = t_out.length();
        
        for (long i = 0; i < p; i++)
            val[i].SetLength(p);

        if constexpr (K == 1)
        {
            for (long j = 0; j < p; j++)
            {
                zz_pXi<1> poly;
                poly.SetLength(l);
                for (long i = 0; i < l; i++)
                    poly[i] = t_out[i][j];
                poly.normalize();
                
                NVec<1> eval_res;
                evaluateReduced(eval_res, poly);
                
                for (long xN = 0; xN < p; xN++)
                    val[xN][j] = eval_res[xN];
            }
        }
        else
        {
            for (long xK = 0; xK < p; xK++)
            {
                Vec<NVec<K-1>> t_out_sub;
                t_out_sub.SetLength(l);

                /* Note : the whole problem of RAM model in the complexity of KU
                comes from here, where generally transposition of n-dimention
                array is expensive, but due to RAM model which is basically the fact
                that we can "move" (or swap for NTL) pointers (because they are 
                non-contiguous adresses of the heap) the swap is constant time but 
                actually fully rewriting would be expensive. */

                for (long i = 0; i < l; i++)
                    swap(t_out_sub[i], t_out[i][xK]);
                
                NVec<K> val_sub;
                val_sub.SetLength(p);
                _aggregate_evals<K-1>(val_sub, t_out_sub);
                
                for (long xN = 0; xN < p; xN++)
                    swap(val[xN][xK], val_sub[xN]);
            }
        }
    }

    template <size_t N>
    void evaluateReduced(NVec<N>& val, const zz_pXi<N>& pol) const
    {
        if constexpr (N == 1)
        {
            evaluate(val.rep, static_cast<zz_pX>(pol));
        }
        else
        {
            long p = zz_p::modulus();
            assert(pol.degXn() < p);

            long l = pol.degXn() + 1;
            val.SetLength(p);

            Vec<NVec<N-1>> t_out;
            t_out.SetLength(l);

            for (long i = 0; i < l; i++)
                evaluateReduced(t_out[i], pol.rep[i]);

            _aggregate_evals<N-1>(val, t_out);
        }
    }

    template <size_t N>
    void evaluate(NVec<N>& val, const zz_pXi<N>& pol) const
    {
        zz_pXi<N> r_pol;
        reduceFermat(r_pol, pol);
        evaluateReduced(val, r_pol);
    }
};

template <auto N>
void reduceFermat(zz_pXi<N>& a, const zz_pXi<N>& b)
{
    long p = zz_p::modulus();
    a.rep.SetLength(b.rep.length());
    for (long i = 0; i < a.rep.length(); i++)
        reduceFermat(a.rep[i], b.rep[i]);

    long d = a.degXn();

    for (long i = d; i >= p; i--)
        a.rep[i - p + 1] += a.rep[i];

    a.rep.SetLength(min(p, d + 1));
    a.normalize();
}

template <>
inline void reduceFermat(zz_pXi<1>& a, const zz_pXi<1>& b)
{
    zz_pX r;
    reduceFermat(r, zz_pX(b));
    a = zz_pXi<1>(r);
}
