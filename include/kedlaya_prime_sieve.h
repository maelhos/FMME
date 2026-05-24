#pragma once

#include <NTL/lzz_pX.h>
#include <NTL/ZZ.h>

NTL_CLIENT

class kedlaya_prime_sieve
{
private:
    long size;
    Vec<long> primeList;
    Vec<bool> isPrime;

public:
    kedlaya_prime_sieve(long n);
    kedlaya_prime_sieve() : kedlaya_prime_sieve(64) {}

    inline Vec<long> primes() const
    {
        return primeList;
    }

    inline long biggestPrime() const
    {
        return primeList[primeList.length() - 1];
    }

    inline long getPrime(long i)
    {
        while (i >= primeList.length()) extend(2*size);
        return primeList[i];
    }

    void extend(long n);
    void boundedMultiModularBasis(Vec<long>& ret, long beta, const ZZ& beta_bar);
};

