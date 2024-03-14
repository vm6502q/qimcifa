// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes
#include <bits/stdc++.h>
 
void SieveOfEratosthenes(const int& n)
{
    std::vector<int> knownPrimes = { 2, 3 };

    int cardinality = n;
    for (int p : knownPrimes) {
        cardinality = ((p - 1) * cardinality) / p;
    }

    int ni = n;
    for (int p : knownPrimes) {
        ni = (p - 1) * (ni + 1) / p;
    }

    // Create a boolean array "prime[0..n]" and initialize
    // all entries it as true. A value in prime[i] will
    // finally be false if i is Not a prime, else true.
    std::vector<bool> notPrime(cardinality + 1);
 
    int o = 2;
    while (true) {
        int p = o;

        // Make this NOT a multiple of 2 or 3.
        p = p + (p >> 1U);
        p = (p << 1U) - 1U;

        if ((p * p) > n) {
            break;
        }

        // If prime[o] is not changed, then it is a prime
        if (!notPrime[o]) {
            // Update all multiples of p greater than or
            // equal to the square of it numbers which are
            // multiple of p and are less than p^2 are
            // already been marked

            for (int i = o * o; i <= ni; i += o) {
                notPrime[o] = true;
            }
        }

        ++o;
    }

    for (int p : knownPrimes) {
        std::cout << p << " ";
    }
 
    // Print all prime numbers
    for (int o = 2; o <= cardinality; ++o) {
        if (!notPrime[o]) {
            // Make this NOT a multiple of 2 or 3.
            int p = o + (o >> 1U);
            p = (p << 1U) - 1U;

            std::cout << p << " ";
        }
    }
}
 
// Driver Code
int main()
{
    int n = 60;
    std:: cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;
    SieveOfEratosthenes(n);
    return 0;
}
