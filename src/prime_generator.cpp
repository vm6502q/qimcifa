// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes
#include <bits/stdc++.h>
 
void SieveOfEratosthenes(const int& n)
{
    std::vector<int> knownPrimes = { 2, 3, 5 };

    int cardinality = n;
    for (int p : knownPrimes) {
        cardinality = (p - 1) * (cardinality + 1) / p;
    }

    // Create a boolean array "prime[0..n]" and initialize
    // all entries it as true. A value in prime[i] will
    // finally be false if i is Not a prime, else true.
    bool prime[cardinality + 1];
    memset(prime, true, sizeof(prime));
 
    for (int p = 7; p * p <= n; p++) {
        int o = p;
        for (int j = knownPrimes.size() - 1; j >= 0; --j) {
            const int q = knownPrimes[j];
            o = ((q - 1) * (o + 1)) / q;
        }
        // If prime[o] is not changed, then it is a prime
        if (prime[o] == true) {
            // Update all multiples of p greater than or
            // equal to the square of it numbers which are
            // multiple of p and are less than p^2 are
            // already been marked.
            for (int i = p * p; i <= n; i += p) {
                int o = i;
                for (int j = knownPrimes.size() - 1; j >= 0; --j) {
                    const int q = knownPrimes[j];
                    o = ((q - 1) * (o + 1)) / q;
                }
                prime[o] = false;
            }
        }
    }

    for (int p : knownPrimes) {
        std::cout << p << " ";
    }
 
    // Print all prime numbers
    for (int o = 2; o <= cardinality; o+=10) {
        for (int batchItem = 1; batchItem < 7; ++batchItem) {
            int p = o + batchItem;
            if (prime[p]) {
                // Make this NOT a multiple of 2 or 3.
                p = p + (p >> 1U);
                p = (p << 1U) - 1U;
                std::cout << p << " ";
            }
        }

        for (int batchItem = 8; batchItem < 10; ++batchItem) {
            int p = o + batchItem;
            if (prime[p]) {
                // Make this NOT a multiple of 2 or 3.
                p = p + (p >> 1U);
                p = (p << 1U) - 1U;
                std::cout << p << " ";
            }
        }
    }
}
 
// Driver Code
int main()
{
    int n = 30;
    std:: cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;
    SieveOfEratosthenes(n);
    return 0;
}
