////////////////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2024. All rights reserved.
//
// "A quantum-inspired Monte Carlo integer factoring algorithm"
//
// This example originally demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer
// factoring. It has since been developed into a general factoring algorithm and tool.
//
// The only potentially "original" part of this factoring algorithm is the "reverse trial
// division," as far as I can tell. The idea is, instead of performing typical trial division,
// we collect a short list of the first primes and remove all of their multiples from a
// "brute-force" guessing range by mapping a dense contiguous integer set, to a set without these
// multiples, by successively applying `guess = guess + guess / (p[i] - 1U) + 1U` for prime "`p`"
// in ascending (or any) order. Each prime applied this way effectively multiplies the
// brute-force guessing cardinality by a fraction (p-1)/p. Whatever "level" of primes we use, the
// cost per "guess" becomes higher.
//
// Then, we have a tuner that empirically estimates the cost per guess, and we multiply this by
// the (known) total cardinality of potential guesses. Whichever reverse trial division level has
// the lowest product of average cost per guess times guessing set cardinality should have the
// best performance, and the best level increases with the scale of the problem.
//
// Beyond this, we gain a functional advantage of a square-root over a more naive approach, by
// setting the brute force guessing range only between the highest prime in reverse trial
// division and the (modular) square root of the number to factor: if the number is semiprime,
// there is exactly one correct answer in this range, but including both factors in the range to
// search would cost us the square root advantage.
//
// Beyond that, we observed that many simple and well-known factoring techniques just don't pay
// dividends, for semiprime factoring. There's basically no point in checking either congruence
// of squares or even for a greatest common divisor, as these techniques require some dynamically
// variable overhead, and it tends to be faster (for semiprimes) just to check if a guess is an
// exact factor, on the smallest range we can identify that contains at least one correct answer.
//
// So, this is actually quite rudimentary and just "brute force," except for "reverse trial
// division" and the upper bound on the guessing range. It just better work entirely in CPU cache,
// then, but it only requires de minimis maximum memory footprint. (There are congruence of squares
// and greatest common divisor checks available for numbers besides semiprimes.)
//
// Theoretically, this algorithm might return to its original "quantum-inspired" design with the
// providence of a high-quality, high-throughput generator of uniform random bit strings. If we were
// to use the algorithm as-is, except guessing according to a uniform random distribution instead of
// systematically ascending through every possible "guess," then the average time to solution can be
// realized in any case, unlike the deterministic version of the algorithm. Then, no work towards
// the solution can ever be lost in event of interruption of the program, because every single guess
// (even the first) has the same probability (in the ideal) of leading to successful factoring.
//
// (This file was heavily adapted from
// https://github.com/ProjectQ-Framework/ProjectQ/blob/develop/examples/shor.py,
// with thanks to ProjectQ!)
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "config.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <float.h>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <stdlib.h>
#include <string>
#include <time.h>

#include <boost/dynamic_bitset.hpp>

#if USE_GMP
#include <boost/multiprecision/gmp.hpp>
#elif USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include "big_integer.hpp"
#endif

namespace Qimcifa {

// Make this a multiple of 2, 3, 5, 7, 11, 13, and 17.
constexpr int BIGGEST_WHEEL = 510510;
constexpr int MIN_RTD_LEVEL = 2;

#if USE_GMP
typedef boost::multiprecision::mpz_int BigIntegerInput;
#elif USE_BOOST
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<8192, 8192,
    boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
    BigIntegerInput;
#else
typedef BigInteger BigIntegerInput;
#endif

BigIntegerInput batchNumber;
BigIntegerInput batchBound;
BigIntegerInput batchCount;
std::mutex batchMutex;

inline void finish() {
    std::lock_guard<std::mutex> lock(batchMutex);
    batchNumber = batchBound;
}

inline BigIntegerInput getNextBatch() {
    std::lock_guard<std::mutex> lock(batchMutex);

#if IS_SQUARES_CONGRUENCE_CHECK
    BigIntegerInput result = batchNumber;

    if (batchNumber < batchBound) {
        ++batchNumber;
    }

    return result;
#else
    BigIntegerInput result = batchCount - (batchNumber + 1U);

    if (batchNumber == batchBound) {
        return batchBound;
    }

    ++batchNumber;

    return result;
#endif
}

// See https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
template <typename BigInteger> BigInteger ipow(BigInteger base, unsigned exp)
{
    BigInteger result = 1U;
    for (;;)
    {
        if (exp & 1U) {
            result *= base;
        }
        exp >>= 1U;
        if (!exp) {
            break;
        }
        base *= base;
    }

    return result;
}

template <typename BigInteger> inline BigInteger sqrt(const BigInteger& toTest)
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

template <typename BigInteger> inline uint64_t log2(BigInteger n) {
#if USE_GMP || USE_BOOST
    uint64_t pow = 0U;
    BigInteger p = n >> 1U;
    while (p) {
        p >>= 1U;
        ++pow;
    }
    return pow;
#else
    return bi_log2(n);
#endif
}

template <typename BigInteger> inline bool isPowerOfTwo(const BigInteger& x)
{
    // Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
#if USE_GMP || USE_BOOST
    return (x && !(x & (x - 1ULL)));
#else
    BigInteger y = x;
    bi_decrement(&y, 1U);
    bi_and_ip(&y, x);
    return (bi_compare_0(x) != 0) && (bi_compare_0(y) == 0);
#endif
}

template <typename BigInteger> inline BigInteger gcd(BigInteger n1, BigInteger n2)
{
    while (n2 != 0) {
        const BigInteger t = n1;
        n1 = n2;
        n2 = t % n2;
    }

    return n1;
}

template <typename BigInteger>
void printSuccess(const BigInteger& f1, const BigInteger& f2, const BigInteger& toFactor, const std::string& message,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock)
{
    std::cout << message << f1 << " * " << f2 << " = " << toFactor << std::endl;
    auto tClock =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
    // Report in seconds
    std::cout << "(Time elapsed: " << (tClock.count() / 1000000.0) << " seconds)" << std::endl;
    std::cout << "(Waiting to join other threads...)" << std::endl;
}

template <typename BigInteger> inline BigInteger backward(BigInteger n) {
    return ((~((~n) | 1U)) / 3U) + 1U;
}

template <typename BigInteger> inline BigInteger forward(BigInteger p) {
    // Make this NOT a multiple of 2 or 3.
    return (p << 1U) + (~(~p | 1U)) - 1U;
}

template <typename BigInteger> bool isMultiple(const BigInteger& p, const std::vector<BigInteger>& knownPrimes) {
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
        output.push_back(wheel_inc(wheelPrimes, limit));
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

#if IS_SQUARES_CONGRUENCE_CHECK
template <typename BigInteger>
inline bool checkCongruenceOfSquares(const BigInteger& toFactor, const BigInteger& toTest,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock)
{
    // The basic idea is "congruence of squares":
    // a^2 = b^2 mod N
    // If we're lucky enough that the above is true, for a^2 = toTest and (b^2 mod N) = remainder,
    // then we can immediately find a factor.

    // Consider a to be equal to "toTest."
    const BigInteger bSqr = (toTest * toTest) % toFactor;
    const BigInteger b = sqrt(bSqr);
    if ((b * b) != bSqr) {
        return false;
    }

    BigInteger f1 = gcd(toTest + b, toFactor);
    BigInteger f2 = gcd(toTest - b, toFactor);
    BigInteger fmul = f1 * f2;
    while ((fmul > 1U) && (fmul != toFactor) && ((toFactor % fmul) == 0)) {
        fmul = f1;
        f1 = f1 * f2;
        f2 = toFactor / (fmul * f2);
        fmul = f1 * f2;
    }
    if ((fmul == toFactor) && (f1 > 1U) && (f2 > 1U)) {
        // Inform the other threads on this node that we've succeeded and are done:
        printSuccess<BigInteger>(f1, f2, toFactor, "Congruence of squares: Found ", iterClock);
        return true;
    }

    return false;
}
#endif

template <typename BigInteger>
inline bool getSmoothNumbersIteration(const BigInteger& toFactor, const BigInteger& base,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock) {
#if IS_RSA_SEMIPRIME
    if ((toFactor % base) == 0U) {
        printSuccess<BigInteger>(base, toFactor / base, toFactor, "Exact factor: Found ", iterClock);
        return true;
    }
#else
    BigInteger n = gcd(base, toFactor);
    if (n != 1U) {
        printSuccess<BigInteger>(n, toFactor / n, toFactor, "Has common factor: Found ", iterClock);
        return true;
    }
#endif

#if IS_SQUARES_CONGRUENCE_CHECK
    return checkCongruenceOfSquares<BigInteger>(toFactor, base, iterClock);
#else
    return false;
#endif
}

template <typename BigInteger>
bool getSmoothNumbers(const BigInteger& toFactor, std::vector<boost::dynamic_bitset<uint64_t>>& inc_seqs, const BigInteger& offset,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock)
{
    for (BigInteger batchNum = (BigInteger)getNextBatch(); batchNum < batchBound; batchNum = (BigInteger)getNextBatch()) {
        const BigInteger batchStart = batchNum * BIGGEST_WHEEL + offset;
        const BigInteger batchEnd = (batchNum + 1U) * BIGGEST_WHEEL + offset;
        for (size_t p = batchStart; p < batchEnd;) {
            p += GetWheelIncrement(inc_seqs);
            if (getSmoothNumbersIteration<BigInteger>(toFactor, forward(p), iterClock)) {
                return true;
            }
        }
    }

    return false;
}
} // namespace Qimcifa
