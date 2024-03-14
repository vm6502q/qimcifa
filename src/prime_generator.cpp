// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes

#include <iostream>
#include <set>
#include <vector>

// 1/3 overall space complexity!
// log reduction in time complexity!
std::vector<size_t> knownPrimes = { 2, 3 };

size_t backward(size_t ni) {
    ni = (ni + 1) >> 1;
    ni = ((ni + 1) << 1) / 3;
    return ni;
}

size_t forward(size_t p) {
    // Make this NOT a multiple of 2 or 3.
    p += (p >> 1U);
    return (p << 1U) - 1U;
}

bool isTimeOrSpaceMultiple(size_t p) {
    for (size_t i : knownPrimes) {
        if ((p % i) == 0) {
            return true;
        }
    }
    return false;
}

bool isTimeMultiple(size_t p) {
    for (size_t i = 2U; i < knownPrimes.size(); ++i) {
        if ((p % knownPrimes[i]) == 0) {
            return true;
        }
    }
    return false;
}
 
void SieveOfEratosthenes(const size_t& n)
{
    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2 and 3, this reduces complexity by 2/3.
    const size_t cardinality = (n & ~1) / 3;

    // Create a boolean array "prime[0..cardinality]"
    // and initialize all entries it as true. Rather,
    // reverse the true/false meaning, so we can use
    // default initialization. A value in notPrime[i]
    // will finally be false only if i is a prime.
    std::vector<bool> notPrime(cardinality + 1);
 
    size_t o = 2;
    while (true) {
        const size_t p = forward(o);
        if ((p * p) > n) {
            break;
        }

        if (isTimeMultiple(p)) {
            ++o;
            continue;
        }

        // If prime[o] is not changed, then it is a prime
        if (!notPrime[o]) {
            // Update all multiples of p greater than or
            // equal to the square of it numbers which are
            // multiple of p and are less than p^2 are
            // already been marked.
            for (size_t i = p * p; i <= n; i += p) {
                // If this is a multiple of one of the
                // filtered primes, then backwards(i)
                // will not return the correct number,
                // but this multiple has already been
                // struck from the set.
                if (isTimeOrSpaceMultiple(i)) {
                    continue;
                }

                notPrime[backward(i)] = true;
            }
            
            knownPrimes.push_back(p);
        }

        ++o;
    }

    for (size_t p : knownPrimes) {
        std::cout << p << " ";
    }
    
    std::set<size_t> outputPrimes(knownPrimes.begin(), knownPrimes.end());
 
    // Print all prime numbers
    for (size_t o = 2; o <= cardinality; ++o) {
        if (!notPrime[o]) {
            const size_t p = forward(o);

            if (isTimeMultiple(p)) {
                continue;
            }

            outputPrimes.insert(p);
        }
    }

    for (size_t p : outputPrimes) {
        std::cout << p << " ";
    }
    std::cout << std::endl;
}
 
// Driver Code
int main()
{
    size_t n = 100;
    std:: cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;
    SieveOfEratosthenes(n);
    return 0;
}
