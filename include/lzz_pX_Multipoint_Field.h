#pragma once

#include <NTL/lzz_pX.h>
#include "lzz_pX_CRT.h"
#include "lzz_pXY.h"
#include <array>

NTL_CLIENT
PML_CLIENT

class zz_pX_Multipoint_Field
{
private:
    long r;
    zz_p q;
    zz_p qs;
    zz_pX_Multipoint_Geometric ev_even;
    zz_pX_Multipoint_Geometric ev_odd;
public:
    zz_pX_Multipoint_Field();
    ~zz_pX_Multipoint_Field() {}

    void evaluate(Vec<zz_p>& val, const zz_pX& pol) const;
    void evaluateReduced(Vec<Vec<zz_p>>& val, const zz_pXY& pol) const; // assume degX < p, degY < p
    inline void evaluate(Vec<Vec<zz_p>>& val, const zz_pXY& pol) const
    {
        zz_pXY r_pol;
        reduceFermat(r_pol, pol);
        evaluateReduced(val, r_pol);
    }
    

    static void reduceFermat(zz_pX& a, const zz_pX& b);
    static void reduceFermat(zz_pXY& a, const zz_pXY& b);
};

void evaluateWholeFieldNaive(Vec<zz_p>& val, const zz_pX& pol); 