
#include "lzz_pX_Multipoint_Field.h"
#include "lzz_pX_extra.h"
#include "lzz_pXY.h"
#include "lzz_pXY_Kedlaya_Umans.h"
#include "kedlaya_prime_sieve.h"
#include <NTL/ALL_FEATURES.h>
#include <NTL/lzz_pX.h>
#include "util.h"

NTL_CLIENT

int test_MultipointField2D(){
    //zz_p::init(14559583);
    zz_p::init(173);
    double t1, t2;

    t1 = get_time();
    zz_pX_Multipoint_Field ev1{};
    t1 = get_time()-t1;


    cout << "[+] time geometric trick : " << t1 << endl;
    cout << "[+] time naive enumerate : 0" << endl;

    long p = zz_p::modulus();
    for (long i = 0; i < 100; i += 1)
    {
        zz_pXY a = random_zz_pXY(1 << i, 1 << i);

        Vec<Vec<zz_p>> v1;
        cout << "[+] Test of degree " << i << " in X" << endl;
        cout << "[-----------------------------------------------]" << endl;

        t1 = get_time();
        ev1.evaluate(v1, a);
        t1 = get_time()-t1;

        t2 = get_time();
        bool cond = 1;
        for (long x = 0; x < p; x++)
        {
            for (long y = 0; y < p; y++)
            {
                cond &= (v1[x][y] == eval(evaluate(a, to_zz_p(x)), to_zz_p(y)));
            }
        }
        
        t2 = get_time()-t2;

        cout << "[-----------------------------------------------]" << endl;
        cout << "[+] time geometric trick 2D : " << t1 << endl;
        cout << "[+] time naive enumerate 2D : " << t2 << endl;
        cout << "[+] equal : " << cond << endl;
        cout << "[+] dX : " << (1 << i) << endl;
        cout << "[-----------------------------------------------]" << endl;
    }

    return 0;
}

int test_MultipointField(){
    //zz_p::init(14559583);
    zz_p::init(49529);
    double t1, t2, t3;

    t1 = get_time();
    zz_pX_Multipoint_Field ev1{};
    t1 = get_time()-t1;

    vec_zz_p pts(INIT_SIZE, zz_p::modulus());
    for (long i = 0; i < zz_p::modulus(); i++) pts[i] = i;

    t2 = get_time();
    zz_pX_Multipoint_General ev2(pts);
    t2 = get_time()-t2;

    cout << "[+] time geometric trick : " << t1 << endl;
    cout << "[+] time subproduct tree : " << t2 << endl;
    cout << "[+] time naive enumerate : 0" << endl;

    for (long i = 0; i < 100; i += 1)
    {
        zz_pX a = random_zz_pX(1 << i);

        vec_zz_p v1, v2, v3;
        cout << "[+] Test of degree " << i << " in X" << endl;
        cout << "[-----------------------------------------------]" << endl;

        t1 = get_time();
        ev1.evaluate(v1, a);
        t1 = get_time()-t1;

        t2 = get_time();
        ev2.evaluate(v2, a);
        t2 = get_time()-t2;

        t3 = get_time();
        evaluateWholeFieldNaive(v3, a);
        t3 = get_time()-t3;

        cout << "[-----------------------------------------------]" << endl;
        cout << "[+] time geometric trick : " << t1 << endl;
        cout << "[+] time subproduct tree : " << t2 << endl;
        cout << "[+] time naive enumerate : " << t3 << endl;
        cout << "[+] equal : " << ((v1 == v2) && (v1 == v3)) << endl;
        cout << "[+] dX : " << (1 << i) << endl;
        cout << "[-----------------------------------------------]" << endl;
    }

    return 0;
}


void test_kedlaya_prime_sieve()
{
    kedlaya_prime_sieve ps(100);

    cout << "first 100 primes : " << endl;
    for (long i = 0; i < ps.primes().length(); i++)
    {
        cout << ps.primes()[i] << " ";
    }
    cout << endl;

    cout << "extension to 200 for test" << endl;
    ps.extend(200);
    cout << "first 200 primes : " << endl;
    for (long i = 0; i < ps.primes().length(); i++)
    {
        cout << ps.primes()[i] << " ";
    }
    cout << endl;

}

void test_kedlaya_multimod()
{
    kedlaya_prime_sieve ps{};

    long low = 200;
    ZZ top(INIT_VAL, "111222333444555");
    cout << "out : " << endl;
    Vec<long> basis{};

    ps.boundedMultiModularBasis(basis, low, top);

    ZZ P(1);

    for (long i = 0; i < basis.length(); i++)
    {
        cout << basis[i] << " ";
        P *= basis[i];
    }
    cout << endl;
    cout << "top  = " << top << endl;
    cout << "prod = " << P << endl;

    cout << "befr = " << P / ps.getPrime(ps.primes().length() - 1) << endl;
    cout << "over = " << P * ps.getPrime(ps.primes().length()) << endl;
}


int test_Kedlaya_Umans_Multimod(){
    zz_p::init(0x64ec6dd0392073ULL);
    //zz_p::init(14559583);
    //zz_p::init(49529);
    //zz_p::init(173);
    double t2, t3;
    zz_pXY_Kedlaya_Umans ev1{};

    int T = 1;

    for (long i = 7; i < 100; i += 1) // i = 7
    {
        zz_pXY a = random_zz_pXY(1 << i, 1 << i);
        Vec<zz_p> v1, v2, v3;
        Vec<ZZ> vz1, vz2;
        for (long j = 20; j < 25; j += 1) // j = 20
        {
            long nbr_pts = (1ULL << j);
            Vec<pair<zz_p, zz_p>> pts{}; pts.SetLength(nbr_pts);
            for (long k = 0; k < nbr_pts; k++)
                pts[k] = make_pair(random_zz_p(), random_zz_p());
            
            cout << "[+] Test of degree " << (1 << i) << " in X and Y and " << nbr_pts << " points" << endl;
            cout << "[-----------------------------------------------]" << endl;
            
            
            //t1 = get_time();
            //cout << "KU for t=0" << endl;
            //ev1.multimodEvaluate(vz1, v1, a, pts, 0);
            //t1 = get_time()-t1;

            t2 = get_time();
            cout << "KU for t=" << T << endl;
            ev1.multimodEvaluate(vz2, v2, a, pts, T);
            t2 = get_time()-t2;

            t3 = get_time();
            v3.SetLength(nbr_pts);
            for (long k = 0; k < nbr_pts; k++)
                eval(v3[k], evaluate(a, pts[k].first), pts[k].second);
            t3 = get_time() - t3;

            cout << "[-----------------------------------------------]" << endl;
            //cout << "[+] time multimod t = 0 : " << t1 << endl;
            cout << "[+] time multimod t = " << T << " : " << t2 << endl;
            cout << "[+] time naive eval : " << t3 << endl;
            //cout << "[+] equal v1 / v3: " << (v1 == v3) << endl;
            cout << "[+] equal v2 / v3: " << (v2 == v3) << endl;
            cout << "[+] v3 is the reference" << endl;
            cout << "[+] dX : " << (1 << i) << endl;

            cout << "[-----------------------------------------------]" << endl;
            cout << endl;

        }
    }

    return 0;
}

int main(int argc, char* argv[]){
    //test_MultipointField2D();
    test_Kedlaya_Umans_Multimod();

    return 0;
}