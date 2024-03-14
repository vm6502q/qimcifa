// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes
#include <bits/stdc++.h>

std::vector<int> knownPrimes = { 2, 3 };

int backward(int ni) {
    for (int p : knownPrimes) {
        ni = ((p - 1) * (ni + 1)) / p;
    }

    return ni;
}

int forward(int p) {
    // Make this NOT a multiple of 2 or 3.
    p += (p >> 1U);
    return (p << 1U) - 1U;
}
 
void SieveOfEratosthenes(const int& n)
{
    const int cardinality = backward(n);

    // Create a boolean array "prime[0..n]" and initialize
    // all entries it as true. A value in prime[i] will
    // finally be false if i is Not a prime, else true.
    std::vector<bool> notPrime(cardinality + 1);
 
    int o = 2;
    while (true) {
        const int p = forward(o);

        if ((p * p) > n) {
            break;
        }

        // If prime[o] is not changed, then it is a prime
        if (!notPrime[o]) {
            // Update all multiples of p greater than or
            // equal to the square of it numbers which are
            // multiple of p and are less than p^2 are
            // already been marked

            for (int i = p * p; i <= n; i += p) {
                for (int j : knownPrimes) {
                    if ((p % j) == 0) {
                        continue;
                    }
                }
                notPrime[backward(i)] = true;
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
            std::cout << forward(o) << " ";
        }
    }
    std::cout << std::endl;
}
 
// Driver Code
int main()
{
    int n = 300;
    std:: cout << "Following are select (biased) prime numbers smaller than or equal to " << n << ":" << std::endl;
    SieveOfEratosthenes(n);
    return 0;
}
