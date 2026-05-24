#include "lzz_pXY_Kedlaya_Umans.h"
#include <assert.h>
#include <bit>
#include "util.h"

#if (1)
    #define cout_debug std::cout
#else
    #define cout_debug if (true) {} else std::cout
#endif

ZZ eval_ZZ_zz_pX(const zz_pX& f, const ZZ a)
{
    ZZ acc;
    long i;
 
    clear(acc);
    for (i = deg(f); i >= 0; i--) {
       mul(acc, acc, a);
       add(acc, acc, ZZ(f.rep[i]._zz_p__rep));
    }
 
    return acc;
}

ZZ eval_ZZ_zz_pXY(const zz_pXY& f, const ZZ a, const ZZ b)
{
    ZZ acc;
    long i;
 
    clear(acc);
    for (i = f.degY(); i >= 0; i--) {
       mul(acc, acc, b);
       add(acc, acc, eval_ZZ_zz_pX(f.rep[i], a));
    }
 
    return acc;
}

ZZ zz_pXY_Kedlaya_Umans::computeUpper(const zz_pXY& pol)
{
    ZZ p(zz_p::modulus());
    ZZ ret{};
    div(ret, (power(p, pol.degX() + 2) - p) * (power(p, pol.degY() + 1) - 1), (p - 1) * (p - 1));

    return ret;
}

void zz_pXY_Kedlaya_Umans::evaluateLookupNaive(Vec<zz_p>& val, const zz_pXY& pol, const zz_pX_Multipoint_Field& ev, const Vec<pair<zz_p, zz_p>>& pts)
{
    //cout_debug << "[+] Call to whole field eval with p = " << zz_p::modulus() << endl;
    Vec<Vec<zz_p>> vals_field;
    ev.evaluate(vals_field, pol);
    val.SetLength(pts.length());
    //cout_debug << "[1/2] Done whole field eval" << endl;

    for (long i = 0; i < pts.length(); i++)
    {
        zz_p x, y;
        x = pts[i].first; y = pts[i].second;
        //cout_debug << "x = " << x << ", y = " << y << endl << flush;
        val[i] = vals_field[x.LoopHole()][y.LoopHole()];
    }

    //cout_debug << "[2/2] Done whole field lookup" << endl;
}

/*
# Complexity : O~(degY(pol) * p + degY(pol)*degX(pol) + p^2 + len(pts))
*/
void zz_pXY_Kedlaya_Umans::evaluateLookupNaiveReduced(Vec<zz_p>& val, const zz_pXY& pol, const zz_pX_Multipoint_Field& ev, const Vec<pair<zz_p, zz_p>>& pts)
{
    //cout_debug << "[+] Call to whole field eval with p = " << zz_p::modulus() << endl;
    assert((pol.degX() < zz_p::modulus()) && (pol.degY() < zz_p::modulus()));
    Vec<Vec<zz_p>> vals_field;
    ev.evaluateReduced(vals_field, pol);
    val.SetLength(pts.length());
    //cout_debug << "[1/2] Done whole field eval" << endl;

    for (long i = 0; i < pts.length(); i++)
    {
        zz_p x, y;
        x = pts[i].first; y = pts[i].second;
        //cout_debug << "x = " << x << ", y = " << y << endl << flush;
        val[i] = vals_field[x.LoopHole()][y.LoopHole()];
    }

    //cout_debug << "[2/2] Done whole field lookup" << endl;
}

long zz_pXY_Kedlaya_Umans::computePhi(const zz_pXY& pol)
{
    return (pol.degX() + pol.degY()) / 8;
}

void zz_pXY_Kedlaya_Umans::getEvaluator(zz_pX_Multipoint_Field& ev)
{
    long r = zz_p::modulus();
    if (!evaluators.contains(r))
        evaluators[r] = zz_pX_Multipoint_Field();
    
    ev = evaluators[r];
}


void zz_pXY_Kedlaya_Umans::multimodEvaluate(Vec<ZZ>& valz, Vec<zz_p>& val, const zz_pXY& pol, const Vec<pair<zz_p, zz_p>>& pts, long t)
{
    long r = zz_p::modulus();
    //cout_debug << "[+] Call to multimod eval with p = " << r << " and t = " << t << endl << flush;

    assert(pol.degX() < r && pol.degY() < r);
    if (t == 0 || pol.tdeg() == 0)
    {
        // we make a quick calculation on which is cheaper 
        Vec<long> idx;
        idx.SetLength(pts.length());
        for (long i = 0; i < pts.length(); i++)
            idx[i] = i;

        sort(&idx[0], &idx[0] + pts.length(), [&pts](long a, long b) {
            long x_a = rep(pts[a].first);
            long x_b = rep(pts[b].first);
            if (x_a != x_b) return x_a < x_b;
            return rep(pts[a].second) < rep(pts[b].second);
        });

        long diff_pts = 1;
        pair<zz_p, zz_p> prev = pts[idx[0]];
        for (long i = 1; i < pts.length(); i++)
        {
            if (pts[idx[i]] != prev)
            {
                diff_pts++;
                prev = pts[idx[i]];
            }
        }
        
        long cost_naive = pol.degX() * pol.degY() * diff_pts + pts.length();
        long p = zz_p::modulus();
        long cost_whole_field = pol.degY() * p + pol.degX() * pol.degY() + p*p + pts.length();

        if (cost_whole_field < cost_naive)
        {
            cout_debug << "Using whole field\n";
            zz_pX_Multipoint_Field ev;
            getEvaluator(ev);
            evaluateLookupNaiveReduced(val, pol, ev, pts);
        }
        else
        {
            cout_debug << "Using naive\n";
            val.SetLength(pts.length());
            zz_p prev_val{};
            pair<zz_p, zz_p> prev_pair{};

            prev_pair = pts[idx[0]];
            eval(prev_val, evaluate(pol, prev_pair.first), prev_pair.second);
            val[idx[0]] = prev_val;

            for (long i = 1; i < pts.length(); i++)
            {
                if (pts[idx[i]] != prev_pair)
                {
                    prev_pair = pts[idx[i]];
                    eval(prev_val, evaluate(pol, prev_pair.first), prev_pair.second);
                }

                val[idx[i]] = prev_val;
            }
        }

        return;
    }

    long phi = computePhi(pol);

    ZZ Up = computeUpper(pol); // we replace power_ZZ(r, phi) by some tighter bound

    Vec<long> basis{};
    cout_debug << "KU computing basis" << endl;
    double t1 = get_time();
    ps.boundedMultiModularBasis(basis, phi, Up);
    cout_debug << "phi = " << phi << ", Up = " << Up << endl;
    cout_debug << "basis : " << basis << endl;
    ZZ prod_b(1);
    for (long i = 0; i < basis.length(); i++) prod_b *= basis[i];
    
    cout_debug << "prod = " << prod_b << " vs p = " << zz_p::modulus() << endl;
    cout_debug << "KU -> " << get_time()-t1 << "s" << endl;
    //cout_debug << "[1/?] basis : [";
    //for (auto &&pr : basis) cout_debug << pr << " ";
    //cout_debug << "]" << endl;

    Vec<ZZ> val_z{};
    Vec<ZZ> _{};
    long nb_pts = pts.length();
    val_z.SetLength(nb_pts);

    // first moduli
    cout_debug << "KU fist moduli" << endl;
    t1 = get_time();
    {
        zz_pPush push(basis[0]);
        //cout_debug << "[1/?] Multimod 1 / " << basis.length() << "prime = " << zz_p::modulus() << endl;

        Vec<zz_p> vals_p{};
        Vec<pair<zz_p, zz_p>> pts_p{};
        zz_pXY pol_r{};
        zz_pX_Multipoint_Field::reduceFermat(pol_r, pol);

        // reinstantialte and reduce pts with new modulus
        pts_p.SetLength(pts.length());
        for (long i = 0; i < pts.length(); i++)
            pts_p[i] = make_pair(to_zz_p(pts[i].first._zz_p__rep), to_zz_p(pts[i].second._zz_p__rep));
        

        multimodEvaluate(_, vals_p, pol_r, pts_p, t - 1);
        for (long i = 0; i < nb_pts; i++)
            val_z[i] = vals_p[i].LoopHole();
    }
    cout_debug << "KU -> " << get_time()-t1 << "s" << endl;
    
    ZZ mod(basis[0]);
    // the rest
    for (long i = 1; i < basis.length(); i++)
    {
        zz_pPush push(basis[i]); // switch the modulus for this frame
        //cout_debug << "[1/?] Multimod " << i + 1 << " / " << basis.length() << "prime = " << zz_p::modulus() << endl;

        Vec<zz_p> vals_p{};
        Vec<pair<zz_p, zz_p>> pts_p{};
        zz_pXY pol_r{};
        zz_pX_Multipoint_Field::reduceFermat(pol_r, pol);

        // reinstantialte and reduce pts with new modulus
        pts_p.SetLength(pts.length());
        for (long i = 0; i < pts.length(); i++)
            pts_p[i] = make_pair(to_zz_p(pts[i].first._zz_p__rep), to_zz_p(pts[i].second._zz_p__rep));
        
        cout_debug << "KU " << i << " / " << basis.length() << " moduli" << endl;
        t1 = get_time();
        multimodEvaluate(_, vals_p, pol_r, pts_p, t - 1);
        cout_debug << "KU -> " << get_time()-t1 << "s" << endl;
        
        cout_debug << "KU CRT" << endl;
        t1 = get_time();
        CRT(val_z, mod, vals_p);
        
        cout_debug << "KU -> " << get_time()-t1 << "s" << endl << endl;;
    }

    val.SetLength(nb_pts);
    for (long i = 0; i < nb_pts; i++)
        conv(val[i], val_z[i]);
    valz = val_z;
    //cout_debug << "exit with p = " << zz_p::modulus() << endl;
}