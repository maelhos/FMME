#include "factorial_store.h"

Vec<zz_p> fac_store;
Vec<zz_p> fac_inv_instore;

void init_fac_store()
{
    fac_store.SetLength(1); fac_store[0] = zz_p(1);
    fac_inv_instore.SetLength(1); fac_inv_instore[0] = zz_p(1);
}

void extend_fac(long i)
{
    long old_l = fac_store.length();

    long new_size = 2*(i + 1);
    fac_store.SetLength(new_size);
    fac_inv_instore.SetLength(new_size);

    for (long j = old_l; j < new_size; j++)
    {
        fac_store[j] = zz_p(j) * fac_store[j - 1];
        fac_inv_instore[j] = 1 / fac_store[j];
    }
    
}