#include "kedlaya_prime_sieve.h"

kedlaya_prime_sieve::kedlaya_prime_sieve(long s) : size(2)
{
    isPrime.SetLength(size);

    isPrime[0] = false;
    isPrime[1] = false;

    extend(max(64, s));
}

void kedlaya_prime_sieve::extend(long n)
{
    if (n <= size) 
    {
        return; // Already have primes up to n
    }
    
    long oldSize = size;
    size = n;
    
    // Resize the isPrime vector
    isPrime.SetLength(n);
    for (long i = oldSize; i < n; i++)
        isPrime[i] = true;
    
    for (long prime : primeList) {
        long start = ((oldSize + prime - 1) / prime) * prime;
        if (start == prime)
            start += prime;
        
        for (long multiple = start; multiple < size; multiple += prime)
            isPrime[multiple] = false;
    }
    
    for (long i = oldSize; i < size; i++) 
    {
        if (isPrime[i]) 
        {
            primeList.append(i);
            
            for (long multiple = i * i; multiple < size; multiple += i)
                isPrime[multiple] = false;
        }
    }
}

void kedlaya_prime_sieve::boundedMultiModularBasis(Vec<long>& ret, long beta, const ZZ& beta_bar)
{
    //cout << "[+] Called mutlimodular basis" << endl;
    while (biggestPrime() <= beta) extend(2*size); 

    long a = 0;
    while (primeList[a] <= beta) a++;// TODO: dicho search, beta is usually small though...

    ZZ prod(1);
    ret.SetLength(0);
    
    for (long i = a; i < primeList.length(); i++)
    {
        ret.append(primeList[i]);
        prod *= primeList[i];
        if (prod >= beta_bar) return;
    }
    
    //cout << "[+] Start climbing from prod = " << prod << " to " << beta_bar << endl;
    while (prod < beta_bar)
    {
        long oldPs = primeList.length();
        extend(2*size);

        for (long i = oldPs; i < primeList.length(); i++)
        {
            ret.append(primeList[i]);
            prod *= primeList[i];
            if (prod >= beta_bar) return;
        }
    }
}