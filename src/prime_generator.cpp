// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes

#include <iostream>
#include <set>
#include <vector>

// log overall space complexity!
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
 
std::set<size_t> SieveOfEratosthenes(const size_t& n)
{
    if (n < 2) {
        return std::set<size_t>();
    }
    if (n < 3) {
        return std::set<size_t>(knownPrimes.begin(), knownPrimes.begin() + 1);
    }
    if (n < 5) {
        return std::set<size_t>(knownPrimes.begin(), knownPrimes.begin() + 2);
    }

    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2 and 3, this reduces complexity by 2/3.
    const size_t cardinality = (n & ~1) / 3;
 
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
        knownPrimes.push_back(p);

        ++o;
    }
    
    std::set<size_t> outputPrimes(knownPrimes.begin(), knownPrimes.end());
 
    // Get the remaining prime numbers.
    for (size_t o = 2; o <= cardinality; ++o) {
        const size_t p = forward(o);

        if (isTimeMultiple(p)) {
            continue;
        }

        outputPrimes.insert(p);
    }

    return outputPrimes;
}
 
// Driver Code
int main()
{
    size_t n = 100;

    std:: cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;

    const std::set<size_t> primes = SieveOfEratosthenes(n);

    for (size_t p : primes) {
        std::cout << p << " ";
    }
    std::cout << std::endl;

    return 0;
}
