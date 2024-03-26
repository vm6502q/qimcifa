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

#include <boost/dynamic_bitset.hpp>
#if USE_GMP
#include <boost/multiprecision/gmp.hpp>
#elif USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include "big_integer.hpp"
#endif

#include "dispatchqueue.hpp"

namespace qimcifa {
const size_t BATCH_SIZE = 1 << 10;
const size_t cpuCount = std::thread::hardware_concurrency();
DispatchQueue dispatch(cpuCount);

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
#else
typedef BigInteger BigInteger;
#endif
#endif

inline BigInteger sqrt(const BigInteger& toTest)
{
    // Otherwise, find b = sqrt(b^2).
    BigInteger start = 1U, end = toTest >> 1U, ans = 0U;
    do {
        const BigInteger mid = (start + end) >> 1U;

        // If toTest is a perfect square
        const BigInteger sqr = mid * mid;
        if (sqr == toTest) {
            return mid;
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

inline BigInteger backward(const BigInteger& n) {
    return ((~((~n) | 1U)) / 3U) + 1U;
}

inline BigInteger forward(const BigInteger& p) {
    // Make this NOT a multiple of 2 or 3.
    return (p << 1U) + (~(~p | 1U)) - 1U;
}

inline BigInteger backward5(BigInteger n) {
    n = ((n + 1U) << 2U) / 5U;
    n = ((n + 1U) << 1U) / 3U;
    return (n + 1U) >> 1U;
}

bool isMultipleParallel(const BigInteger& p, const size_t& nextPrimeIndex, const size_t& highestIndex,
    const std::vector<BigInteger>& knownPrimes) {
    const size_t _BATCH_SIZE = BATCH_SIZE;
    const size_t maxLcv = (highestIndex - nextPrimeIndex) / BATCH_SIZE;
    dispatch.resetResult();
    for (size_t i = 0; i < maxLcv; ++i) {
        size_t j = i * BATCH_SIZE + nextPrimeIndex;
        dispatch.dispatch([&knownPrimes, &p, _BATCH_SIZE, j]() {
            for (size_t k = 0; k < _BATCH_SIZE; ++k) {
                if ((p % knownPrimes[j + k]) == 0) {
                    return true;
                }
            }
            return false;
        });
    }

    return dispatch.finish();
}

bool isMultiple(const BigInteger& p, size_t nextIndex, const std::vector<BigInteger>& knownPrimes) {
    const BigInteger sqrtP = sqrt(p);
    const size_t highestIndex = std::distance(knownPrimes.begin(), std::upper_bound(knownPrimes.begin(), knownPrimes.end(), sqrtP));

    const size_t diff = highestIndex - nextIndex;
    if ((highestIndex > nextIndex) && ((diff >> 1U) > BATCH_SIZE)) {
        if (isMultipleParallel(p, nextIndex, highestIndex, knownPrimes)) {
            return true;
        }
    }
    nextIndex = diff % BATCH_SIZE;

    for (size_t i = nextIndex; i < highestIndex; ++i) {
        if ((p % knownPrimes[i]) == 0) {
            return true;
        }
    }
    return false;
}

template <typename BigInteger>
bool isMultiple(const BigInteger& p, const std::vector<BigInteger>& knownPrimes) {
    for (const BigInteger& prime : knownPrimes) {
        if ((p % prime) == 0) {
            return true;
        }
    }
    return false;
}

template <typename BigInteger>
boost::dynamic_bitset<size_t> wheel_inc(std::vector<BigInteger> primes, BigInteger limit) {
    BigInteger radius = 1U;
    for (const BigInteger& i : primes) {
        radius *= i;
    }
    if (limit < radius) {
        radius = limit;
    }
    const BigInteger prime = primes.back();
    primes.pop_back();
    std::vector<bool> o;
    for (BigInteger i = 1U; i <= radius; ++i) {
        if (!isMultiple(i, primes)) {
            o.push_back((i % prime) == 0);
        }
    }

    boost::dynamic_bitset<size_t> output(o.size());
    for (size_t i = 0U; i < o.size(); ++i) {
        output[i] = o[i];
    }
    output >>= 1U;

    return output;
}

template <typename BigInteger>
std::vector<boost::dynamic_bitset<size_t>> wheel_gen(const std::vector<BigInteger>& primes, BigInteger limit) {
    std::vector<boost::dynamic_bitset<size_t>> output;
    std::vector<BigInteger> wheelPrimes;
    for (const BigInteger& p : primes) {
        wheelPrimes.push_back(p);
        if (wheelPrimes.back() > 3U) {
            output.push_back(wheel_inc(wheelPrimes, limit));
        }
    }
    return output;
}

inline size_t GetWheelIncrement(std::vector<boost::dynamic_bitset<size_t>>& inc_seqs) {
    size_t wheelIncrement = 0U;
    bool is_wheel_multiple = false;
    do {
        for (size_t i = 0; i < inc_seqs.size(); ++i) {
            boost::dynamic_bitset<size_t>& wheel = inc_seqs[i];
            is_wheel_multiple = wheel.test(0U);
            wheel >>= 1U;
            if (is_wheel_multiple) {
                wheel[wheel.size() - 1U] = true;
                break;
            }
        }
        wheelIncrement++;
    } while (is_wheel_multiple);

    return wheelIncrement;
}

inline size_t GetWheel5Increment(uint32_t& wheel5) {
    size_t wheelIncrement = 0U;
    bool is_wheel_multiple = false;
    do {
        is_wheel_multiple = (bool)(wheel5 & 1U);
        wheel5 >>= 1U;
        if (is_wheel_multiple) {
            wheel5 |= 1U << 9U;
        }
        wheelIncrement++;
    } while (is_wheel_multiple);

    return wheelIncrement;
}

inline BigInteger makeNotMultiple(BigInteger n) {
    n |= 1U;
    if ((n % 3U) == 0U) {
        n -= 2U;
    }

    return n;
}

std::vector<BigInteger> TrialDivision(const BigInteger& n);
std::vector<BigInteger> SieveOfEratosthenes(const BigInteger& n);
std::vector<BigInteger> SegmentedSieveOfEratosthenes(const BigInteger& n);
} // namespace qimcifa
