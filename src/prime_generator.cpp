// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes
#include <bits/stdc++.h>

// 1/3 overall space complexity!
// 27/35 reduction in time complexity (and higher)!
std::vector<int> knownPrimes = { 2, 3, 5, 7 };

int backward(int ni) {
    ni = (ni + 1) >> 1;
    ni = ((ni + 1) << 1) / 3;
    return ni;
}

int forward(int p) {
    // Make this NOT a multiple of 2 or 3.
    p += (p >> 1U);
    return (p << 1U) - 1U;
}
 
void SieveOfEratosthenes(const int& n)
{
    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2 and 3, this reduces complexity by 2/3.
    const int cardinality = (n & ~1) / 3;

    // Create a boolean array "prime[0..cardinality]"
    // and initialize all entries it as true. Rather,
    // reverse the true/false meaning, so we can use
    // default initialization. A value in notPrime[i]
    // will finally be false only if i is a prime.
    std::vector<bool> notPrime(cardinality + 1);
 
    int o = 4;
    while (true) {
        const int p = forward(o);
        if ((p * p) > n) {
            break;
        }

        // This reduces time complexity by
        // skipping multiples of 5.
        if ((p % 5) == 0) {
            ++o;
            continue;
        }
        // This reduces time complexity by
        // skipping multiples of 7.
        if ((p % 7) == 0) {
            ++o;
            continue;
        }
        // (This process can continue, for other primes.)

        // If prime[o] is not changed, then it is a prime
        if (!notPrime[o]) {
            // Update all multiples of p greater than or
            // equal to the square of it numbers which are
            // multiple of p and are less than p^2 are
            // already been marked.

            for (int i = p * p; i <= n; i += p) {

                // If this is a multiple of one of the
                // filtered primes, then backwards(i)
                // will not return the correct number,
                // but this multiple has already been
                // struck from the set.
                bool isMultiple = false;
                for (int j : knownPrimes) {
                    if ((i % j) == 0) {
                        isMultiple = true;
                        break;
                    }
                }
                if (isMultiple) {
                    continue;
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
    for (int o = 3; o <= cardinality; ++o) {
        if (!notPrime[o]) {
            const int p = forward(o);

            // Skip multiples of 5.
            if ((p % 5) == 0) {
                continue;
            }
            // Skip multiples of 7.
            if ((p % 7) == 0) {
                continue;
            }
            // (This can continue, for other primes.)

            std::cout << p << " ";
        }
    }
    std::cout << std::endl;
}
 
// Driver Code
int main()
{
    int n = 100;
    std:: cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;
    SieveOfEratosthenes(n);
    return 0;
}
