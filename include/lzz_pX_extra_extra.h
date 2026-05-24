#pragma once

#include <NTL/lzz_pX.h>
#include <sstream>
#include "lzz_pX_CRT.h"
#include <NTL/lzz_pXFactoring.h>
#include "factorial_store.h"
#include "lzz_pX_middle_product.h"

NTL_CLIENT
PML_CLIENT

/** Pretty print for testing purposes */
string pp(const zz_pX& a);

void CompNaive(zz_pX& c, const zz_pX& a, const zz_pX& b);

// this is very slightly faster than CompNaive(c, a, x + b) ...
void ShiftX(zz_pX& c, const zz_pX& a, zz_p b);
void CompInterpolate(zz_pX& c, const zz_pX& a, const zz_pX& b);

struct zz_pXComp{
    zz_p a0, b0;
    long db;
    zz_pX_Multipoint_Geometric evg;
    Vec<zz_p> val_b;
    
    zz_pXComp(const zz_pX& b, long max_deg_a);  
};

void CompInterpolate(zz_pX& c, const zz_pX& a, const zz_pXComp& b);
vec_zz_p zz_pX_roots(const zz_pX& P);
