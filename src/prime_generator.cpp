// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes

#include "config.h"

#include <iostream>
#include <vector>

#if USE_GMP
#include <boost/multiprecision/gmp.hpp>
#elif USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include "big_integer.hpp"
#endif

#if BIG_INT_BITS < 32
typedef uint32_t BigInteger
#elif BIG_INT_BITS < 64
typedef uint64_t BigInteger
#else
#if USE_GMP
typedef boost::multiprecision::mpz_int BigInteger;
#elif USE_BOOST
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<BIG_INT_BITS, BIG_INT_BITS,
    boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
    BigInteger;
#else
typedef BigInteger BigInteger;
#endif
#endif

// Improvements by Dan Strano of Unitary Fund, 2024:
// log overall space complexity!
// log reduction in time complexity!
std::vector<BigInteger> knownPrimes = { 2, 3 };

BigInteger backward(BigInteger ni) {
    ni = (ni + 1) >> 1;
    ni = ((ni + 1) << 1) / 3;
    return ni;
}

BigInteger forward(BigInteger p) {
    // Make this NOT a multiple of 2 or 3.
    p += (p >> 1U);
    return (p << 1U) - 1U;
}

bool isTimeOrSpaceMultiple(BigInteger p) {
    for (BigInteger i : knownPrimes) {
        if ((p % i) == 0) {
            return true;
        }
    }
    return false;
}

bool isTimeMultiple(BigInteger p) {
    for (size_t i = 2U; i < knownPrimes.size(); ++i) {
        if ((p % knownPrimes[i]) == 0) {
            return true;
        }
    }
    return false;
}

std::vector<BigInteger> TrialDivision(const BigInteger& n)
{
    if (n < 2) {
        return std::vector<BigInteger>();
    }
    if (n < 3) {
        return std::vector<BigInteger>(knownPrimes.begin(), knownPrimes.begin() + 1);
    }
    if (n < 5) {
        return std::vector<BigInteger>(knownPrimes.begin(), knownPrimes.begin() + 2);
    }

    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2 and 3, this reduces complexity by 2/3.
    const BigInteger cardinality = (~((~n) | 1)) / 3;

    // Get the remaining prime numbers.
    for (BigInteger o = 2; o <= cardinality; ++o) {
        const BigInteger p = forward(o);

        if (isTimeMultiple(p)) {
            continue;
        }

        knownPrimes.push_back(p);
    }

    return knownPrimes;
}
 
// Driver Code
int main()
{
    BigInteger n = 100;

    std:: cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;

    const std::vector<BigInteger> primes = TrialDivision(n);

    for (BigInteger p : primes) {
        std::cout << p << " ";
    }
    std::cout << std::endl;

    return 0;
}
