#pragma once

#include <NTL/lzz_pX.h>
#include <NTL/ZZ.h>
#include <NTL/tools.h>
#include <unordered_map>
#include <algorithm>
#include "lzz_pX_CRT.h"
#include "lzz_pX_Multipoint_Field.h"
#include "kedlaya_prime_sieve.h"

NTL_CLIENT
PML_CLIENT

#define get_time GetWallTime

class zz_pXY_Kedlaya_Umans
{
private:
    unordered_map<long, zz_pX_Multipoint_Field> evaluators;
    kedlaya_prime_sieve ps;
public:
    zz_pXY_Kedlaya_Umans() : evaluators{}, ps{}
        { }
    ~zz_pXY_Kedlaya_Umans() {}

    //void evaluate(Vec<zz_p>& val, const zz_pXY& pol, const Vec<pair<zz_p, zz_p>>& pts) {} // stub
//private:
    void evaluateLookupNaive(Vec<zz_p>& val, const zz_pXY& pol, const zz_pX_Multipoint_Field& ev, const Vec<pair<zz_p, zz_p>>& pts);
    long computePhi(const zz_pXY& pol);
    ZZ computeUpper(const zz_pXY& pol);
    void getEvaluator(zz_pX_Multipoint_Field& ev);

    // asserts that degX,Y pol < p
    void evaluateLookupNaiveReduced(Vec<zz_p>& val, const zz_pXY& pol, const zz_pX_Multipoint_Field& ev, const Vec<pair<zz_p, zz_p>>& pts);

    // asserts that degX,Y pol < p
    void multimodEvaluate(Vec<ZZ>& valz, Vec<zz_p>& val, const zz_pXY& pol, const Vec<pair<zz_p, zz_p>>& pts, long t);
};

ZZ eval_ZZ_zz_pX(const zz_pX& f, const ZZ a);
ZZ eval_ZZ_zz_pXY(const zz_pXY& f, const ZZ a, const ZZ b);