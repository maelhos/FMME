#pragma once
#include <NTL/lzz_p.h>

NTL_CLIENT

extern Vec<zz_p> fac_store;
extern Vec<zz_p> fac_inv_instore;


void extend_fac(long i);
void init_fac_store();

inline zz_p fac(long i)
{
    if (i >= fac_store.length()) extend_fac(i);
    return fac_store[i];
}

inline zz_p inv_fac(long i)
{
    if (i >= fac_inv_instore.length()) extend_fac(i);
    return fac_inv_instore[i];
}