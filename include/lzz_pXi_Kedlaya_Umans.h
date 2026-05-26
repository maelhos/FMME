#pragma once

#include <NTL/lzz_pX.h>
#include <NTL/ZZ.h>
#include <NTL/tools.h>
#include <unordered_map>
#include <algorithm>
#include "lzz_pXi.h"
#include "lzz_pX_CRT.h"
#include "lzz_pX_Multipoint_Field.h"
#include "prime_sieve.h"
#include "lzz_pXY_Kedlaya_Umans.h"

NTL_CLIENT
PML_CLIENT

#define get_time GetWallTime

class zz_pXi_Kedlaya_Umans
{
private:
    unordered_map<long, zz_pX_Multipoint_Field> evaluators;
    prime_sieve ps;
public:
    zz_pXi_Kedlaya_Umans() : evaluators{}, ps{}
        { }
    ~zz_pXi_Kedlaya_Umans() {}

    template <auto N>
    void evaluateLookupNaive(Vec<zz_p>& val, const zz_pXi<N>& pol, const zz_pX_Multipoint_Field& ev, const Vec<array<zz_p, N>>& pts)
    {
        // todo
    }
    
    template <auto N>
    long computePhi(const zz_pXi<N>& pol)
    {
        // todo
    }

    template <auto N>
    ZZ computeUpper(const zz_pXi<N>& pol)
    {
        // todo
    }

    void getEvaluator(zz_pX_Multipoint_Field& ev)
    {
        // todo
    }

    // asserts that degX,Y pol < p
    template <auto N>
    void evaluateLookupNaiveReduced(Vec<zz_p>& val, const zz_pXi<N>& pol, const zz_pX_Multipoint_Field& ev, const Vec<array<zz_p, N>>& pts)
    {
        // todo
    }

    // asserts that degX,Y pol < p
    template <auto N>
    void multimodEvaluate(Vec<ZZ>& valz, Vec<zz_p>& val, const zz_pXi<N>& pol, const Vec<array<zz_p, N>>& pts, long t)
    {
        // todo
    }
};

template <auto N>
ZZ eval_ZZ_zz_pXi(const zz_pXi<N>& f, const ZZ a)
{

}