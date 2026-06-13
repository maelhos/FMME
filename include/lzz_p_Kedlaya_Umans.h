#pragma once

#include <NTL/lzz_pX.h>
#include <NTL/ZZ.h>
#include <NTL/tools.h>
#include <unordered_map>
#include <algorithm>
#include "lzz_pX_CRT.h"
#include "lzz_p_Multipoint_Field.h"
#include "prime_sieve.h"

NTL_CLIENT
PML_CLIENT

#define get_time GetWallTime
#define PRINT_DEBUG (0)

#if (PRINT_DEBUG)
    #define cout_debug std::cout
#else
    #define cout_debug if (true) {} else std::cout
#endif


class zz_p_Kedlaya_Umans
{
private:
    unordered_map<long, zz_p_Multipoint_Field> evaluators;
    prime_sieve ps;
public:
    zz_p_Kedlaya_Umans() : evaluators{}, ps{}
        { }
    ~zz_p_Kedlaya_Umans() {}


    void evaluateLookupNaive(Vec<zz_p>& val, const zz_pXY& pol, const zz_p_Multipoint_Field& ev, const Vec<pair<zz_p, zz_p>>& pts);
    
    long computePhi(const zz_pXY& pol);
    template <size_t N> long computePhi(const zz_pXi<N>& pol)
    {
        array<long, N> degs = degrees(pol);
    
        long ret = 0;
        for (long e : degs)
            ret += e;    
        ret /= 8*N;

        return 5; //max(ret, 5);
    }

    ZZ computeUpper(const zz_pXY& pol);
    template <size_t N> ZZ computeUpper(const zz_pXi<N>& pol)
    {
        ZZ p(zz_p::modulus());
        array<long, N> degs = degrees(pol);

        ZZ ret(p);
        for (long e : degs)
            div(ret, ret*(power(p, e + 1) - 1), p - 1);

        return ret;
    }

    void getEvaluator(zz_p_Multipoint_Field& ev);
    void removeEvaluator(zz_p_Multipoint_Field& ev);

    // asserts that degX,Y pol < p
    void evaluateLookupNaiveReduced(Vec<zz_p>& val, const zz_pXY& pol, const zz_p_Multipoint_Field& ev, const Vec<pair<zz_p, zz_p>>& pts);
    template <size_t N>
    void evaluateLookupNaiveReduced(Vec<zz_p>& val, const zz_pXi<N>& pol, const zz_p_Multipoint_Field& ev, const Vec<array<zz_p, N>>& pts)
    {
        NVec<N> field_vals;
        ev.evaluateReduced(field_vals, pol);
        
        val.SetLength(pts.length());
        for (long i = 0; i < pts.length(); i++)
        {
            val[i] = getEval(field_vals, pts[i]);
        }
    }
    // asserts that degX,Y pol < p
    void multimodEvaluate(Vec<ZZ>& valz, Vec<zz_p>& val, const zz_pXY& pol, const Vec<pair<zz_p, zz_p>>& pts, long t);
    
    template <size_t N>
    void multimodEvaluate(Vec<ZZ>& valz, Vec<zz_p>& val, const zz_pXi<N>& pol, const Vec<array<zz_p, N>>& pts, long t)
    {
        long r = zz_p::modulus();
        array<long, N> degs = degrees(pol);
        
        bool is_const = true;
        for (long d : degs)
        {
            assert(d < r);
            if (d > 0) is_const = false;
        }
            
        if (t == 0 || is_const)
        {
            Vec<long> idx;
            idx.SetLength(pts.length());
            for (long i = 0; i < pts.length(); i++)
                idx[i] = i;

            sort(&idx[0], &idx[0] + pts.length(), [&pts](long a, long b) {
                for (size_t i = 0; i < N; i++)
                {
                    long xa = rep(pts[a][i]);
                    long xb = rep(pts[b][i]);
                    if (xa != xb) return xa < xb;
                }
                return false;
            });

            long diff_pts = 1;
            array<zz_p, N> prev = pts[idx[0]];
            for (long i = 1; i < pts.length(); i++)
            {
                bool eq = true;
                for (size_t j = 0; j < N; j++)
                {
                    if (pts[idx[i]][j] != prev[j])
                    {
                        eq = false;
                        break;
                    }
                }
                if (!eq)
                {
                    diff_pts++;
                    prev = pts[idx[i]];
                }
            }
            
            long p = zz_p::modulus();
            long d_prod = 1;
            for (long d : degs) d_prod *= (d + 1);
            
            long cost_naive = d_prod * diff_pts + pts.length();
            long cost_whole_field = power_long(p, N) + pts.length();

            if (cost_whole_field < cost_naive)
            {
                cout_debug << "Using whole field\n";
                zz_p_Multipoint_Field ev;
                getEvaluator(ev);
                evaluateLookupNaiveReduced(val, pol, ev, pts);
                removeEvaluator(ev);
            }
            else
            {
                cout_debug << "Using naive\n";
                val.SetLength(pts.length());
                zz_p prev_val;
                array<zz_p, N> prev_pair = pts[idx[0]];
                
                prev_val = eval(pol, prev_pair);
                val[idx[0]] = prev_val;

                for (long i = 1; i < pts.length(); i++)
                {
                    bool eq = true;
                    for (size_t j = 0; j < N; j++)
                    {
                        if (pts[idx[i]][j] != prev_pair[j])
                        {
                            eq = false;
                            break;
                        }
                    }
                    
                    if (!eq)
                    {
                        prev_pair = pts[idx[i]];
                        prev_val = eval(pol, prev_pair);
                    }
                    val[idx[i]] = prev_val;
                }
            }
            return;
        }

        long phi = computePhi(pol);
        ZZ Up = computeUpper(pol);

        Vec<long> basis;
        cout_debug << "KU computing basis" << endl;
        double t1 = get_time();
        ps.boundedMultiModularBasis(basis, phi, Up);
        cout_debug << "phi = " << phi << ", Up = " << Up << endl;
        cout_debug << "basis : " << basis << endl;
        ZZ prod_b(1);
        for (long i = 0; i < basis.length(); i++) prod_b *= basis[i];
        
        cout_debug << "prod = " << prod_b << " vs p = " << zz_p::modulus() << endl;
        cout_debug << "KU -> " << get_time()-t1 << "s" << endl;

        Vec<ZZ> val_z;
        Vec<ZZ> _;
        long nb_pts = pts.length();
        val_z.SetLength(nb_pts);

        cout_debug << "KU fist moduli" << endl;
        t1 = get_time();
        {
            zz_pPush push(basis[0]);

            Vec<zz_p> vals_p;
            Vec<array<zz_p, N>> pts_p;
            zz_pXi<N> pol_r;
            reduceFermat(pol_r, pol);

            pts_p.SetLength(nb_pts);
            for (long i = 0; i < nb_pts; i++)
            {
                for (size_t j = 0; j < N; j++)
                    pts_p[i][j] = to_zz_p(rep(pts[i][j]));
            }
            
            multimodEvaluate(_, vals_p, pol_r, pts_p, t - 1);
            for (long i = 0; i < nb_pts; i++)
                val_z[i] = vals_p[i].LoopHole();
        }
        cout_debug << "KU -> " << get_time()-t1 << "s" << endl;
        
        ZZ mod(basis[0]);
        for (long i = 1; i < basis.length(); i++)
        {
            zz_pPush push(basis[i]);

            Vec<zz_p> vals_p;
            Vec<array<zz_p, N>> pts_p;
            zz_pXi<N> pol_r;
            reduceFermat(pol_r, pol);

            pts_p.SetLength(nb_pts);
            for (long k = 0; k < nb_pts; k++)
            {
                for (size_t j = 0; j < N; j++)
                    pts_p[k][j] = to_zz_p(rep(pts[k][j]));
            }
            
            cout_debug << "KU " << i << " / " << basis.length() << " moduli" << endl;
            t1 = get_time();
            multimodEvaluate(_, vals_p, pol_r, pts_p, t - 1);
            cout_debug << "KU -> " << get_time()-t1 << "s" << endl;
            
            cout_debug << "KU CRT" << endl;
            t1 = get_time();
            CRT(val_z, mod, vals_p);
            
            cout_debug << "KU -> " << get_time()-t1 << "s" << endl << endl;
        }

        val.SetLength(nb_pts);
        for (long i = 0; i < nb_pts; i++)
            conv(val[i], val_z[i]);
        valz = val_z;
    }
};

ZZ eval_ZZ_zz_pX(const zz_pX& f, const ZZ a);
ZZ eval_ZZ_zz_pXY(const zz_pXY& f, const ZZ a, const ZZ b);

template <size_t N>
ZZ _eval_ZZ_zz_pXi_impl(const zz_pXi<N>& c, const ZZ* point)
{
    ZZ ret(0);
    for (long d = c.degXn(); d >= 0; d--)
    {
        ret *= point[N - 1];
        ret += _eval_ZZ_zz_pXi_impl(c.rep[d], point);
    }
    return ret;
}

inline ZZ _eval_ZZ_zz_pXi_impl(const zz_pXi<1>& c, const ZZ* point)
{
    return eval_ZZ_zz_pX(static_cast<zz_pX>(c), point[0]);
}

template <unsigned long N>
inline ZZ eval_ZZ_zz_pXi(const zz_pXi<N>& c, const array<ZZ, N>& point)
{
    return _eval_ZZ_zz_pXi_impl(c, point.data());
}