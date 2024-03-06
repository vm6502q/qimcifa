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

#include <chrono>
#include <cmath>
#include <float.h>
#include <fstream>
#include <iomanip> // For setw
#include <iostream> // For cout
#include <random>
#include <stdlib.h>
#include <string>
#include <time.h>

#include <algorithm>
#include <future>
#include <map>
#include <mutex>

#if USE_GMP
#include <boost/multiprecision/gmp.hpp>
#elif USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include "big_integer.hpp"
#endif

namespace Qimcifa {

constexpr int BASE_TRIALS = 1U << 16U;
constexpr int MIN_RTD_LEVEL = 1;
constexpr int MIN_RTD_INDEX = 0;

#if !(USE_GMP || USE_BOOST)
typedef BigInteger bitCapIntInput;
typedef BigInteger bitCapInt;
const bitCapInt ZERO_BCI = 0U;

std::ostream& operator<<(std::ostream& os, bitCapInt b)
{
    if (bi_compare_0(b) == 0) {
        os << "0";
        return os;
    }

    // Calculate the base-10 digits, from lowest to highest.
    std::vector<std::string> digits;
    while (bi_compare_0(b) != 0) {
        bitCapInt quo;
        BIG_INTEGER_HALF_WORD rem;
        bi_div_mod_small(b, 10U, &quo, &rem);
        digits.push_back(std::to_string((unsigned char)rem));
        b = quo;
    }

    // Reversing order, print the digits from highest to lowest.
    for (size_t i = digits.size() - 1U; i > 0; --i) {
        os << digits[i];
    }
    // Avoid the need for a signed comparison.
    os << digits[0];

    return os;
}

std::istream& operator>>(std::istream& is, bitCapInt& b)
{
    // Get the whole input string at once.
    std::string input;
    is >> input;

    // Start the output address value at 0.
    b = ZERO_BCI;
    for (size_t i = 0; i < input.size(); ++i) {
        // Left shift by 1 base-10 digit.
        b = b * 10;
        // Add the next lowest base-10 digit.
        bi_increment(&b, (input[i] - 48U));
    }

    return is;
}

inline bool isPowerOfTwo(const bitCapInt& x)
{
    bitCapInt y = x;
    bi_decrement(&y, 1U);
    bi_and_ip(&y, x);
    return (bi_compare_0(x) != 0) && (bi_compare_0(y) == 0);
}

template <typename bitCapInt> inline bitCapInt sqrt(const bitCapInt& toTest)
{
    // Otherwise, find b = sqrt(b^2).
    bitCapInt start = 1U, end = toTest >> 1U, ans = 0U;
    do {
        const bitCapInt mid = (start + end) >> 1U;

        // If toTest is a perfect square
        const bitCapInt sqr = mid * mid;
        if (bi_compare(sqr, toTest) == 0) {
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
    } while (bi_compare(start, end) <= 0);

    return ans;
}
#endif

template <typename bitCapInt> inline bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
#if USE_GMP || USE_BOOST
    while (n2) {
#else
    while (bi_compare_0(n2) != 0) {
#endif
        const bitCapInt t = n1;
        n1 = n2;
        n2 = t % n2;
    }

    return n1;
}

template <typename bitCapInt> inline size_t pickTrialDivisionLevel(const int64_t& tdLevel, const size_t& nodeCount)
{
    if (tdLevel >= 0) {
        return tdLevel;
    }

    std::ifstream settingsFile ("qimcifa_calibration.ssv");
    std::string header;
    std::getline(settingsFile, header);
    size_t bestLevel = MIN_RTD_LEVEL;
    double bestCost = DBL_MAX;
    while (settingsFile.peek() != EOF)
    {
        size_t level;
        bitCapInt cardinality;
        double batchTime, cost;
        settingsFile >> level;
        settingsFile >> cardinality;
        settingsFile >> batchTime;
        settingsFile >> cost;

        if ((level >= MIN_RTD_LEVEL) && (cost < bestCost)) {
            bestLevel = level;
            bestCost = cost;
        }
    }
    settingsFile.close();

    std::cout << "Calibrated reverse trial division level: " << bestLevel << std::endl;
    std::cout << "Average sample time, times worst-case guess count, divided by thread count: " << (bestCost / (std::thread::hardware_concurrency() * nodeCount)) << " seconds" << std::endl;

    return bestLevel;
}

template <typename bitCapInt>
void printSuccess(const bitCapInt& f1, const bitCapInt& f2, const bitCapInt& toFactor, const std::string& message,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock)
{
    std::cout << message << f1 << " * " << f2 << " = " << toFactor << std::endl;
    auto tClock =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
    // Report in seconds
    std::cout << "(Time elapsed: " << (tClock.count() / 1000000.0) << " seconds)" << std::endl;
    std::cout << "(Waiting to join other threads...)" << std::endl;
}

template <typename bitCapInt> inline void finish(bitCapInt& batchNumber, std::mutex& batchMutex) {
    std::lock_guard<std::mutex> lock(batchMutex);
    batchNumber = (bitCapInt)-1;
}

template <typename bitCapInt> inline bitCapInt getNextBatch(bitCapInt& batchNumber, std::mutex& batchMutex) {
    std::lock_guard<std::mutex> lock(batchMutex);
    const bitCapInt result = batchNumber;
#if USE_GMP || USE_BOOST
    if (batchNumber != (bitCapInt)-1) {
        ++batchNumber;
    }
#else
    if (bi_compare(batchNumber, (bitCapInt)-1) != 0) {
        bi_increment(&batchNumber, 1U);
    }
#endif

    return result;
}

#if IS_SQUARES_CONGRUENCE_CHECK
template <typename bitCapInt>
inline bool checkCongruenceOfSquares(const bitCapInt& toFactor, const bitCapInt& toTest,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock)
{
    // The basic idea is "congruence of squares":
    // a^2 = b^2 mod N
    // If we're lucky enough that the above is true, for a^2 = toTest and (b^2 mod N) = remainder,
    // then we can immediately find a factor.

    bitCapInt remainder = (toTest * toTest) % toFactor;

    // The sqrt() algorithm is adapted from Gaurav Ahirwar's suggestion on
    // https://www.geeksforgeeks.org/square-root-of-an-integer/
    // It's a binary search for floor(sqrt(toTest)).

    // If a^2 = 1 mod N, then b = 1.
#if USE_GMP || USE_BOOST
    if (remainder > 1U) {
#else
    if (bi_compare_1(remainder) > 0) {
#endif
        // Otherwise, find b = sqrt(b^2).
        bitCapInt start = 1U, ans = 0U;
#if USE_GMP || USE_BOOST
        remainder >>= 1U;
#else
        bi_rshift_ip(&remainder, 1U);
#endif
        do {
            const bitCapInt mid = (start + remainder) >> 1U;

            // If toTest is a perfect square
            const bitCapInt sqr = mid * mid;
#if USE_GMP || USE_BOOST
            if (sqr == toTest) {
#else
            if (bi_compare(sqr, toTest) == 0) {
#endif
                ans = mid;
                break;
            }

            if (sqr < toTest) {
                // Since we need floor, we update answer when mid*mid is smaller than p, and move closer to sqrt(p).
                start = mid + 1U;
                ans = mid;
            } else {
                // If mid*mid is greater than p
                remainder = mid - 1U;
            }
#if USE_GMP || USE_BOOST
        } while (start <= remainder);
        if (start > remainder) {
#else
        } while (bi_compare(start, remainder) <= 0);
        if (bi_compare(start, remainder) > 0) {
#endif
            // Must be a perfect square.
            return false;
        }

        remainder = ans;
    }

    bitCapInt f1 = gcd(toTest + remainder, toFactor);
    bitCapInt f2 = gcd(toTest - remainder, toFactor);
    bitCapInt fmul = f1 * f2;
#if USE_GMP || USE_BOOST
    while ((fmul > 1U) && (fmul != toFactor) && ((toFactor % fmul) == 0)) {
#else
    while ((bi_compare_1(fmul) > 0) && (bi_compare(fmul, toFactor) != 0) && (bi_compare_0(toFactor % fmul) == 0)) {
#endif
        fmul = f1;
        f1 = f1 * f2;
        f2 = toFactor / (fmul * f2);
        fmul = f1 * f2;
    }
#if USE_GMP || USE_BOOST
    if ((fmul == toFactor) && (f1 > 1U) && (f2 > 1U)) {
#else
    if ((bi_compare(fmul, toFactor) == 0) && (bi_compare_1(f1) > 0) && (bi_compare_1(f2) > 0)) {
#endif
        // Inform the other threads on this node that we've succeeded and are done:
        printSuccess<bitCapInt>(f1, f2, toFactor, "Congruence of squares: Found ", iterClock);
        return true;
    }

    return false;
}
#endif

template <typename bitCapInt>
bool singleWordLoop(const bitCapInt& toFactor, const bitCapInt& range, const bitCapInt& threadMin, const bitCapInt& fullMinBase,
    const size_t& primeIndex, const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock,
    const std::vector<unsigned>& trialDivisionPrimes, bitCapInt& batchNumber, std::mutex& batchMutex)
{
    // Batching reduces mutex-waiting overhead, on the std::atomic broadcast.
    bitCapInt lcv = BASE_TRIALS * getNextBatch(batchNumber, batchMutex);

#if USE_GMP || USE_BOOST
    while ((lcv != (bitCapInt)-1) && (lcv < range)) {
#else
    while ((bi_compare(lcv, (bitCapInt)-1) != 0) && (bi_compare(lcv, range) < 0)) {
#endif
        for (int batchItem = 0U; batchItem < BASE_TRIALS; ++batchItem) {
            // Choose a base at random, >1 and <toFactor.
            bitCapInt base = threadMin + lcv + batchItem;

            for (size_t i = primeIndex; i > MIN_RTD_INDEX; --i) {
                // Make this NOT a multiple of prime "p", by adding it to itself divided by (p - 1), + 1.
                base = base + base / (trialDivisionPrimes[i] - 1U) + 1U;
            }

            // Make this not a multiple of 5.
            // base = base + (base >> 2U) + 1U;

            // Make this not a multiple of 3.
            // base = base + (base >> 1U) + 1U;

            // Make this odd, then shift the range.
            base = ((base << 1U) | 1U) + fullMinBase;

#if IS_RSA_SEMIPRIME
#if USE_GMP || USE_BOOST
            if ((toFactor % base) == 0U) {
#else
            if (bi_compare_0(toFactor % base) == 0U) {
#endif
                finish(batchNumber, batchMutex);
                printSuccess<bitCapInt>(base, toFactor / base, toFactor, "Exact factor: Found ", iterClock);
                return true;
            }
#else
            bitCapInt n = gcd(base, toFactor);
            if (n != 1U) {
                finish(batchNumber, batchMutex);
                printSuccess<bitCapInt>(n, toFactor / n, toFactor, "Has common factor: Found ", iterClock);
                return true;
            }
#endif

#if IS_SQUARES_CONGRUENCE_CHECK
            if (checkCongruenceOfSquares<bitCapInt>(toFactor, base, iterClock)) {
                finish(batchNumber, batchMutex);
                return true;
            }
#endif
        }

        // Check if finished, between batches.
        lcv = getNextBatch(batchNumber, batchMutex);
    }

    return true;
}

template <typename bitCapInt>
int mainBody(const bitCapInt& toFactor, const size_t& qubitCount, const size_t& nodeCount, const size_t& nodeId, const int64_t& tdLevel,
    const std::vector<unsigned>& trialDivisionPrimes)
{
    auto iterClock = std::chrono::high_resolution_clock::now();
    const bitCapInt fullMaxBase = sqrt(toFactor);
#if USE_GMP || USE_BOOST
    if (fullMaxBase * fullMaxBase == toFactor) {
#else
    if (bi_compare(fullMaxBase * fullMaxBase, toFactor) == 0) {
#endif
        std::cout << "Number to factor is a perfect square: " << fullMaxBase << " * " << fullMaxBase << " = " << toFactor;
        return 0;
    }

    int primeIndex = 0;
    while (primeIndex <= tdLevel) {
        const unsigned currentPrime = trialDivisionPrimes[primeIndex];
#if USE_GMP || USE_BOOST
        if ((toFactor % currentPrime) == 0) {
#else
        if (bi_compare_0(toFactor % currentPrime) == 0) {
#endif
            std::cout << "Factors: " << currentPrime << " * " << (toFactor / currentPrime) << " = " << toFactor
                      << std::endl;
            return 0;
        }
        ++primeIndex;
    }

    // We include potential factors as low as the next odd number after the highest trial division prime.
    bitCapInt fullMinBase = trialDivisionPrimes[tdLevel] + 2U;

    primeIndex = tdLevel - 1;
    while (primeIndex >= 0) {
        // The truncation here is a conservative bound, but it's exact if we
        // happen to be aligned to a perfect factor of all trial division.
        const unsigned currentPrime = trialDivisionPrimes[primeIndex];
        fullMinBase = (fullMinBase / currentPrime) * currentPrime;
        --primeIndex;
    }

    bitCapInt fullRange = fullMaxBase + 1U - fullMinBase;
    primeIndex = tdLevel - 1;
    while (primeIndex >= 0) {
        // The truncation here is a conservative bound, but it's exact if we
        // happen to be aligned to a perfect factor of all trial division.
        const unsigned currentPrime = trialDivisionPrimes[primeIndex];
        fullRange = ((fullRange + 1U) * (currentPrime - 1U)) / currentPrime;
        --primeIndex;
    }
    primeIndex = tdLevel - 1;

    bitCapInt batchNumber = 0U;
    std::mutex batchMutex;

    const auto workerFn = [toFactor, iterClock, primeIndex, qubitCount, fullMinBase, &trialDivisionPrimes, &batchNumber, &batchMutex]
        (bitCapInt threadMin, bitCapInt threadMax) {
        singleWordLoop<bitCapInt>(toFactor, threadMax - (threadMin + 1U), threadMin, fullMinBase, primeIndex, iterClock,
            trialDivisionPrimes, batchNumber, batchMutex);
    };

    const unsigned cpuCount = std::thread::hardware_concurrency();
    const bitCapInt nodeRange = (fullRange + nodeCount - 1U) / nodeCount;
    const bitCapInt nodeMin = fullMinBase + nodeRange * nodeId;
    const bitCapInt threadRange = (nodeRange + cpuCount - 1U) / cpuCount;
    std::vector<std::future<void>> futures(cpuCount);
    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        const bitCapInt threadMin = nodeMin + threadRange * cpu;
        const bitCapInt threadMax = threadMin + nodeRange;
        futures[cpu] = std::async(std::launch::async, workerFn, threadMin, threadMax);
    }

    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        futures[cpu].get();
    }

    return 0;
}
} // namespace Qimcifa

using namespace Qimcifa;

int main()
{
#if USE_GMP
    typedef boost::multiprecision::mpz_int bitCapIntInput;
#elif USE_BOOST
    typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<4096, 4096,
        boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
        bitCapIntInput;
#endif

    // First 1000 primes
    // (Only 100 included in program)
    // Source: https://gist.github.com/cblanc/46ebbba6f42f61e60666#file-gistfile1-txt
    const std::vector<unsigned> trialDivisionPrimes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
        61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
        181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
        311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439,
        443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541 }; //, 547, 557, 563, 569, 571, 577, 587,
#if 0
        593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727,
        733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877,
        881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021,
        1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
        1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283,
        1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429,
        1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549,
        1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667,
        1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
        1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973,
        1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089,
        2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243,
        2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377,
        2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
        2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677,
        2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791,
        2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939,
        2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089,
        3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
        3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391,
        3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541,
        3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677,
        3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833,
        3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
        4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129,
        4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271,
        4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447,
        4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597,
        4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
        4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931,
        4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051,
        5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227,
        5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399,
        5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
        5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683,
        5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839,
        5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007,
        6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151,
        6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
        6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451,
        6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637,
        6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791,
        6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949,
        6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
        7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253,
        7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477,
        7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589,
        7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741,
        7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 };
#endif

    // Print primes table by index:
    // for (size_t i = 0; i < trialDivisionPrimes.size(); ++i) {
    //     std::cout << i << ": " << trialDivisionPrimes[i] << ", ";
    // }

    bitCapIntInput toFactor;

    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    uint32_t qubitCount = 0;
    bitCapIntInput p = toFactor >> 1U;
#if USE_GMP || USE_BOOST
    while (p) {
        p >>= 1U;
#else
    while (bi_compare_0(p)) {
        bi_rshift_ip(&p, 1U);
#endif
        ++qubitCount;
    }
    // Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
#if USE_GMP || USE_BOOST
    if (!(toFactor && !(toFactor & (toFactor - 1ULL)))) {
#else
    if (!isPowerOfTwo(toFactor)) {
#endif
        qubitCount++;
    }
    std::cout << "Bits to factor: " << (int)qubitCount << std::endl;

    size_t nodeCount = 1U;
    size_t nodeId = 0U;
#if IS_DISTRIBUTED
    std::cout << "You can split this work across nodes, without networking!" << std::endl;
    do {
        std::cout << "Number of nodes (>=1): ";
        std::cin >> nodeCount;
        if (!nodeCount) {
            std::cout << "Invalid node count choice!" << std::endl;
        }
    } while (!nodeCount);
    if (nodeCount > 1U) {
        do {
            std::cout << "Which node is this? (0-" << (nodeCount - 1U) << "): ";
            std::cin >> nodeId;
            if (nodeId >= nodeCount) {
                std::cout << "Invalid node ID choice!" << std::endl;
            }
        } while (nodeId >= nodeCount);
    }
#endif

    int64_t tdLevel = -1;
    std::cout << "Reverse trial division level (minimum of " << MIN_RTD_LEVEL << ", or -1 for calibration file): ";
    std::cin >> tdLevel;
    if ((tdLevel > -1) && (tdLevel < MIN_RTD_LEVEL)) {
        tdLevel = MIN_RTD_LEVEL;
    }
    if (tdLevel >= (int64_t)trialDivisionPrimes.size()) {
        tdLevel = trialDivisionPrimes.size() - 1U;
    }
    tdLevel = pickTrialDivisionLevel<bitCapIntInput>(tdLevel, nodeCount);

    const unsigned highestPrime = trialDivisionPrimes[tdLevel];
    size_t primeFactorBits = 1U;
    p = highestPrime >> 1U;
#if USE_GMP || USE_BOOST
    while (p) {
        p >>= 1U;
#else
    while (bi_compare_0(p)) {
        bi_rshift_ip(&p, 1U);
#endif
        ++primeFactorBits;
    }
#if !(USE_GMP || USE_BOOST)
    typedef BigInteger bitCapInt;
    return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
#else
    const size_t QBCAPBITS = primeFactorBits + (((qubitCount >> 5U) + 1U) << 5U);
    if (QBCAPBITS < 64) {
        typedef uint64_t bitCapInt;
        return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
#if USE_GMP
    } else {
        return mainBody<bitCapIntInput>(toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
    }
#elif USE_BOOST
    } else if (QBCAPBITS < 128) {
        typedef boost::multiprecision::uint128_t bitCapInt;
        return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
    } else if (QBCAPBITS < 192) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<192, 192,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            bitCapInt;
        return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
    } else if (QBCAPBITS < 256) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<256, 256,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            bitCapInt;
        return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
    } else if (QBCAPBITS < 512) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<512, 512,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            bitCapInt;
        return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
    } else if (QBCAPBITS < 1024) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<1024, 1024,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            bitCapInt;
        return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
    } else if (QBCAPBITS < 2048) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<2048, 2048,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            bitCapInt;
        return mainBody<bitCapInt>((bitCapInt)toFactor, qubitCount, nodeCount, nodeId, tdLevel, trialDivisionPrimes);
    }
#endif
#endif
}
