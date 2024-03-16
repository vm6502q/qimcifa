// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes

// Improved by Dan Strano of Unitary Fund, 2024.
// We can think of trial division as exact inverse of
// Sieve of Eratosthenes, with log space and log time.
// The modular division part is a costly atomic operation.
// It need only be carried out up the square root of the
// number under trial. Multiples of 2, 3, 5, 7, and 11 can
// be entirely skipped in loop enumeration.

#include "config.h"

#include <iostream>
#include <vector>

#include <pybind11/pybind11.h>

#if USE_GMP
#include <boost/multiprecision/gmp.hpp>
#elif USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include "big_integer.hpp"
#endif

#if BIG_INT_BITS < 33
typedef uint32_t BigInteger;
#elif BIG_INT_BITS < 65
typedef uint64_t BigInteger;
#else
#if USE_GMP
typedef boost::multiprecision::mpz_int BigInteger;
#elif USE_BOOST
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<BIG_INT_BITS, BIG_INT_BITS,
    boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
    BigInteger;
#endif
#endif

inline BigInteger gcd(BigInteger n1, BigInteger n2)
{
#if USE_GMP || USE_BOOST
    while (n2) {
#else
    if (bi_compare_0(n2) != 0) {
#endif
        const BigInteger t = n1;
        n1 = n2;
        n2 = t % n2;
    }

    return n1;
}

inline BigInteger sqrt(const BigInteger& toTest)
{
    // Otherwise, find b = sqrt(b^2).
    BigInteger start = 1U, end = toTest >> 1U, ans = 0U;
    do {
        const BigInteger mid = (start + end) >> 1U;

        // If toTest is a perfect square
        const BigInteger sqr = mid * mid;
        if (sqr == toTest) {
            ans = mid;
            break;
        }

        if (sqr < toTest) {
            // Since we need floor, we update answer when mid*mid is smaller than p, and move closer to sqrt(p).
            start = mid + 1U;
            ans = mid;
        } else {
            // If mid*mid is greater than p
            end = mid - 1U;
        }
    } while (start <= end);

    return ans;
}

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

#if 0
bool isTimeOrSpaceMultiple(BigInteger p, const std::vector<BigInteger>& knownPrimes) {
    const BigInteger sqrtP = sqrt(p);
    if ((sqrtP * sqrtP) == p) {
        return true;
    }
    for (BigInteger kp : knownPrimes) {
        if (kp >= sqrtP) {
            return false;
        }
        if ((p % kp) == 0) {
            return true;
        }
    }
    return false;
}
#endif

bool isTimeMultiple(BigInteger p, const std::vector<BigInteger>& knownPrimes) {
    const BigInteger sqrtP = sqrt(p);
    if ((sqrtP * sqrtP) == p) {
        return true;
    }
    for (size_t i = 5U; i < knownPrimes.size(); ++i) {
        const BigInteger kp = knownPrimes[i];
        if (kp >= sqrtP) {
            return false;
        }
        if ((p % kp) == 0) {
            return true;
        }
    }
    return false;
}

std::vector<BigInteger> TrialDivision(const BigInteger& n)
{
    std::vector<BigInteger> knownPrimes =
        { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 77, 79, 83, 91, 97, 103, 107, 109, 119, 127, 133, 137, 139, 149, 151, 161, 163, 167 };

    if (n < 2) {
        return std::vector<BigInteger>();
    }

    if (n < 170) {
        for (size_t i = 1U; i < knownPrimes.size(); ++i) {
            if (n < knownPrimes[i]) {
               return std::vector<BigInteger>(knownPrimes.begin(), knownPrimes.begin() + (i - 1U));
            }
        }
    }

    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2 and 3, this reduces complexity by 2/3.
    // const BigInteger cardinality = (~((~n) | 1)) / 3;

    // Get the remaining prime numbers.
    bool isWorking = true;
    int lcv7 = -11;
    int lcv11 = -16;
    BigInteger o = 3;
    BigInteger wheel = 17;
    while (isWorking) {
        for (int i = 0; i < 6; ++i) {
            if (lcv7 == 11) {
                lcv7 = 1;
                continue;
            }
            if (lcv7 == 7) {
                lcv7 = 8;
                continue;
            }
            ++lcv7;

            if (lcv11 == 17) {
                lcv11 = 1;
                continue;
            }
            if (lcv11 == 7) {
                lcv11 = 8;
                continue;
            }
            ++lcv11;

            BigInteger p = forward(o + i);

            if (p > n) {
                isWorking = false;
                break;
            }

            // **Hear me out**: We've "solved" up to multiples of 11.
            // It's trivial to know much higher primes than this.
            // At any such boundary of our knowledge, we can assume
            // that the highest prime necessary to know, to skip the
            // beginning work of the algorithm, would be the square
            // of the highest "inside-out" Wheel Factorization prime.
            //
            // Grant me only one step further, that the least expensive
            // way to remove 13 from here might be n % 13. For the edge
            // case, < 170 (13*13+1=169+1) is skipped, if we can know
            // that many primes (or obviously higher, hard storage).
            if (p < 170) {
                continue;
            }

            if (!(p % 13)) {
                // Skip
                continue;
            }

            if (gcd(p, wheel) != 1) {
                // Skip
                continue;
            }

            wheel *= p;
            knownPrimes.push_back(p);
        }

        for (int i = 7; i < 9; ++i) {
            if (lcv7 == 11) {
                lcv7 = 1;
                continue;
            }
            if (lcv7 == 7) {
                lcv7 = 8;
                continue;
            }
            ++lcv7;

            if (lcv11 == 17) {
                lcv11 = 1;
                continue;
            }
            if (lcv11 == 7) {
                lcv11 = 8;
                continue;
            }
            ++lcv11;

            BigInteger p = forward(o + i);

            if (p > n) {
                isWorking = false;
                break;
            }

            // **SEE LONG NOTE ABOVE**
            if (p < 170) {
                continue;
            }

            if (!(p % 13)) {
                // Skip
                continue;
            }

            if (gcd(p, wheel) != 1) {
                // Skip
                continue;
            }

            wheel *= p;
            knownPrimes.push_back(p);
        }

        o = o + 10;
    }

    return knownPrimes;
}

PYBIND11_MODULE(prime_gen, m) {
    m.doc() = "pybind11 plugin to generate prime numbers";
    m.def("prime_gen", &TrialDivision, "A function that returns all primes up to the value of its argument");
}
