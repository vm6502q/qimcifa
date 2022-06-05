////////////////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2022. All rights reserved.
//
// "A quantum-inspired Monte Carlo integer factoring algorithm"
//
// This example demonstrates a (Shor's-like) "quantum-inspired" algorithm for integer factoring.
// This approach is similar to Shor's algorithm, except with a uniformly random output from the
// quantum period-finding subroutine. Therefore, we don't need quantum computer simulation for
// this algorithm at all!
//
// (This file was heavily adapted from
// https://github.com/ProjectQ-Framework/ProjectQ/blob/develop/examples/shor.py,
// with thanks to ProjectQ!)
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include <chrono>
#include <cmath>
#include <iomanip> // For setw
#include <iostream> // For cout
#include <random>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <atomic>
#include <future>
#include <map>
#include <mutex>

// Turn this off, if you're not factoring a semi-prime number with equal-bit-width factors.
#define IS_RSA_SEMIPRIME 1
// Turn this off, if you don't want to coordinate across multiple (quasi-independent) nodes.
#define IS_DISTRIBUTED 1
// The maximum number of bits in Boost big integers is 2^QBCAPPOW.
// (2^7, only, needs custom std::cout << operator implementation.)
#define QBCAPPOW 7U

#if QBCAPPOW < 32U
#define bitLenInt uint32_t
#else
#define bitLenInt uint64_t
#endif

#if QBCAPPOW < 6U
#define bitCapInt uint32_t
#define ONE_BCI 1UL
#elif QBCAPPOW < 7U
#define bitCapInt uint64_t
#define ONE_BCI 1ULL
#elif QBCAPPOW < 8U
#include <boost/multiprecision/cpp_int.hpp>
#define bitCapInt boost::multiprecision::uint128_t
#define ONE_BCI 1ULL
#else
#include <boost/multiprecision/cpp_int.hpp>
#define bitCapInt                                                                                                      \
    boost::multiprecision::number<boost::multiprecision::cpp_int_backend<1ULL << QBCAPPOW, 1ULL << QBCAPPOW,           \
        boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
#define ONE_BCI 1ULL
#endif

#define WORD uint64_t
#define WORD_SIZE 64U

namespace Qimcifa {

#if QBCAPPOW == 7U
std::ostream& operator<<(std::ostream& os, bitCapInt b)
{
    // Calculate the base-10 digits, from lowest to highest.
    std::vector<std::string> digits;
    while (b) {
        digits.push_back(std::to_string((unsigned char)(b % 10U)));
        b /= 10U;
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
    b = 0;
    for (size_t i = 0; i < input.size(); ++i) {
        // Left shift by 1 base-10 digit.
        b *= 10;
        // Add the next lowest base-10 digit.
        b += (input[i] - 48U);
    }

    return is;
}
#endif

// Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
inline bool isPowerOfTwo(const bitCapInt& x) { return (x && !(x & (x - ONE_BCI))); }

inline bitLenInt log2(const bitCapInt& n)
{
#if __GNUC__ && QBCAPPOW < 7
// Source: https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers#answer-11376759
#if QBCAPPOW < 6
    return (bitLenInt)((sizeof(unsigned int) << 3U) - __builtin_clz((unsigned int)n) - 1U);
#else
    return (bitLenInt)((sizeof(unsigned long long) << 3U) - __builtin_clzll((unsigned long long)n) - 1U);
#endif
#else
    bitLenInt pow = 0U;
    bitCapInt p = n >> 1U;
    while (p) {
        p >>= 1U;
        ++pow;
    }
    return pow;
#endif
}

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
bitCapInt uipow(bitCapInt base, bitCapInt exp)
{
    bitCapInt result = 1U;
    for (;;) {
        if (base & 1U) {
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

bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
    while (n2) {
        const bitCapInt t = n1;
        n1 = n2;
        n2 = t % n2;
    }

    return n1;
}

typedef std::uniform_int_distribution<WORD> rand_dist;

std::vector<rand_dist> rangeRange(bitCapInt range) {
    std::vector<rand_dist> distToReturn;
    while (range) {
        distToReturn.push_back(rand_dist(0U, (WORD)range));
        range >>= WORD_SIZE;
    }
    std::reverse(distToReturn.begin(), distToReturn.end());

    return distToReturn;
}

} // namespace Qimcifa

using namespace Qimcifa;

int main()
{
    bitCapInt toFactor;
    bitCapInt nodeCount = 1U;
    bitCapInt nodeId = 0U;

    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    const bitLenInt qubitCount = log2(toFactor) + (isPowerOfTwo(toFactor) ? 0U : 1U);
    std::cout << "Bits to factor: " << (int)qubitCount << std::endl;

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
            std::cout << "Which node is this? (0-" << (int)(nodeCount - 1U) << "):";
            std::cin >> nodeId;
            if (nodeId >= nodeCount) {
                std::cout << "Invalid node ID choice!" << std::endl;
            }
        } while (nodeId >= nodeCount);
    }
#endif

    auto iterClock = std::chrono::high_resolution_clock::now();

    std::random_device rand_dev;
    std::mt19937 rand_gen(rand_dev());

    const unsigned cpuCount = std::thread::hardware_concurrency();
    std::atomic<bool> isFinished;
    isFinished = false;

#if IS_RSA_SEMIPRIME
    // Minimum prime number in bit width catalog
    std::map<bitLenInt, bitCapInt> primeDict = {
        { 16U, 32771U },
        { 28U, 134217757U },
        { 32U, 2147483659U },
        { 64U, 9223372036854775837U }
    };

    const bitLenInt primeBits = (qubitCount + 1U) >> 1U;
    const bitCapInt minPrime = primeDict[primeBits] ? primeDict[primeBits] : ((ONE_BCI << (primeBits - 1U)) + 1);
    const bitCapInt fullMinBase = minPrime;
    const bitCapInt fullMaxBase = toFactor / minPrime;
#else
    const bitCapInt fullMinBase = 2U;
    const bitCapInt fullMaxBase = toFactor / 2;
#endif
    const bitCapInt nodeRange = (fullMaxBase + 1U - fullMinBase) / nodeCount;
    const bitCapInt nodeMin = fullMinBase + nodeRange * nodeId;
    const bitCapInt nodeMax = ((nodeId + 1U) == nodeCount) ? fullMaxBase : (fullMinBase + nodeRange * (nodeId + 1U) - 1U);

    auto workerFn = [&toFactor, &nodeMin, &nodeMax, &iterClock, &rand_gen, &isFinished](int cpu) {
        // These constants are semi-redundant, but they're only defined once per thread,
        // and compilers differ on lambda expression capture of constants.

        // Batching reduces mutex-waiting overhead, on the std::atomic broadcast.
        // Batch size is BASE_TRIALS * PERIOD_TRIALS.

        // Number of times to reuse a random base:
        const size_t BASE_TRIALS = 1U << 16U;

        const double clockFactor = 1.0 / 1000.0; // Report in ms
        const unsigned threads = std::thread::hardware_concurrency();

        const bitCapInt threadRange = (nodeMax + 1U - nodeMin) / threads;
        const bitCapInt rMin = nodeMin + threadRange * cpu;
        const bitCapInt rMax = ((cpu + 1U) == threads) ? nodeMax : (nodeMin + threadRange * (cpu + 1U) - 1U);
        std::vector<rand_dist> rDist(rangeRange(rMax - rMin));

        for (;;) {
            for (size_t batchItem = 0U; batchItem < BASE_TRIALS; ++batchItem) {
                // Choose a base at random, >1 and <toFactor.
                bitCapInt base = rDist[0U](rand_gen);
                for (size_t i = 1U; i < rDist.size(); ++i) {
                    base <<= WORD_SIZE;
                    base |= rDist[i](rand_gen);
                }
                base += rMin;

#define PRINT_SUCCESS(f1, f2, toFactor, message)                                                                       \
    std::cout << message << (f1) << " * " << (f2) << " = " << (toFactor) << std::endl;                                 \
    auto tClock =                                                                                                      \
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);  \
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;                            \
    std::cout << "(Waiting to join other threads...)" << std::endl

#define TEST_GCD(testFactor, toFactor, message)                                                                        \
    if ((testFactor) != 1U) {                                                                                          \
        isFinished = true;                                                                                             \
        PRINT_SUCCESS(testFactor, toFactor / testFactor, toFactor, message);                                           \
        return;                                                                                                        \
    }

#define TEST_DIVIDE(testFactor, toFactor, message)                                                                     \
    if (((toFactor) % (testFactor)) == 0U) {                                                                           \
        isFinished = true;                                                                                             \
        PRINT_SUCCESS(testFactor, toFactor / testFactor, toFactor, message);                                           \
        return;                                                                                                        \
    }

                bitCapInt testFactor = gcd(toFactor, base);
                TEST_GCD(testFactor, toFactor, "Base has common factor: Found ");
            }

            // Check if finished, between batches.
            if (isFinished) {
                return;
            }
        }
    };

    std::vector<std::future<void>> futures(cpuCount);
    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        futures[cpu] = std::async(std::launch::async, workerFn, cpu);
    };

    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        futures[cpu].get();
    }

    return 0;
}
