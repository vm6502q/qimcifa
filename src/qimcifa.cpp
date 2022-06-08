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
// For now, this turns on/off elimination of factors of 5, for testing.
#define IS_HI_TRIAL_DIVISION 1
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

#if (QBCAPPOW < 6U) || (IS_RSA_SEMIPRIME && (QBCAPPOW < 7U))
#define WORD uint32_t
#define WORD_SIZE 32U
#else
#define WORD uint64_t
#define WORD_SIZE 64U
#endif

namespace Qimcifa {

#if QBCAPPOW == 7U
std::ostream& operator<<(std::ostream& os, bitCapInt b)
{
    if (b == 0) {
        os << "0";
        return os;
    }

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
inline bitCapInt uipow(bitCapInt base, bitCapInt exp)
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

inline bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
    while (n2) {
        const bitCapInt t = n1;
        n1 = n2;
        n2 = t % n2;
    }

    return n1;
}

typedef std::uniform_int_distribution<WORD> rand_dist;

std::vector<rand_dist> randRange(bitCapInt range)
{
    std::vector<rand_dist> distToReturn;
    while (range) {
        distToReturn.push_back(rand_dist(0U, (WORD)range));
        range >>= WORD_SIZE;
    }
    std::reverse(distToReturn.begin(), distToReturn.end());

    return distToReturn;
}

void printSuccess(bitCapInt f1, bitCapInt f2, bitCapInt toFactor, std::string message,
    std::chrono::time_point<std::chrono::high_resolution_clock> iterClock)
{
    const double clockFactor = 1.0 / 1000.0; // Report in ms

    std::cout << message << f1 << " * " << f2 << " = " << toFactor << std::endl;
    auto tClock =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;
    std::cout << "(Waiting to join other threads...)" << std::endl;
}

} // namespace Qimcifa

using namespace Qimcifa;

int main()
{
    bitCapInt toFactor;
    size_t nodeCount = 1U;
    size_t nodeId = 0U;

    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    const bitLenInt qubitCount = log2(toFactor) + (isPowerOfTwo(toFactor) ? 0U : 1U);
    std::cout << "Bits to factor: " << (int)qubitCount << std::endl;

    // We've eliminated all factors of 2 and 3 from the RNG.
    if ((toFactor & 1) == 0) {
        std::cout << "Factors: 2 * " << (toFactor >> 1) << " = " << toFactor << std::endl;
        return 0;
    }
    if ((toFactor % 3) == 0) {
        std::cout << "Factors: 3 * " << (toFactor / 3) << " = " << toFactor << std::endl;
        return 0;
    }
    if ((toFactor % 5) == 0) {
        std::cout << "Factors: 5 * " << (toFactor / 5) << " = " << toFactor << std::endl;
        return 0;
    }
    if ((toFactor % 7) == 0) {
        std::cout << "Factors: 7 * " << (toFactor / 7) << " = " << toFactor << std::endl;
        return 0;
    }

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
            std::cout << "Which node is this? (0-" << (nodeCount - 1U) << "):";
            std::cin >> nodeId;
            if (nodeId >= nodeCount) {
                std::cout << "Invalid node ID choice!" << std::endl;
            }
        } while (nodeId >= nodeCount);
    }
#endif

    auto iterClock = std::chrono::high_resolution_clock::now();

#if IS_RSA_SEMIPRIME
    std::map<bitLenInt, const std::vector<bitCapInt>> primeDict = { { 16U, { 16411U, 65521U } },
        { 28U, { 67108879U, 536870909U } }, { 32U, { 1073741827U, 8589934583U } } };

    const bitLenInt primeBits = (qubitCount + 1U) >> 1U;
    const bitCapInt minPrime =
        primeDict[primeBits].size() ? primeDict[primeBits][0] : ((ONE_BCI << (primeBits - 2U)) | 1U);
    const bitCapInt maxPrime =
        primeDict[primeBits].size() ? primeDict[primeBits][1] : ((ONE_BCI << (primeBits + 1U)) - 1U);
    const bitCapInt fullMinBase = ((toFactor / maxPrime) < minPrime) ? minPrime : ((toFactor / maxPrime) | 1U);
    const bitCapInt fullMaxBase = ((toFactor / minPrime) > maxPrime) ? maxPrime : ((toFactor / minPrime) | 1U);
#else
    // We include potential factors as low as 11.
    const bitCapInt fullMinBase = 11U;
    // We include potential factors as high as toFactor / 11.
    const bitCapInt fullMaxBase = toFactor / 11U;
#endif
    const bitCapInt nodeRange = ((nodeCount - 1U) + (fullMaxBase - fullMinBase)) / nodeCount;
    const bitCapInt nodeMin = fullMinBase + nodeRange * nodeId;
    const bitCapInt nodeMax = nodeMin + nodeRange;

    std::random_device rand_dev;
    std::mt19937 rand_gen(rand_dev());

    const unsigned cpuCount = std::thread::hardware_concurrency();
    std::atomic<bool> isFinished;
    isFinished = false;

    const auto workerFn = [toFactor, nodeMin, nodeMax, iterClock, &rand_gen, &isFinished](int cpu, unsigned cpuCount) {
        // These constants are semi-redundant, but they're only defined once per thread,
        // and compilers differ on lambda expression capture of constants.

        // Batching reduces mutex-waiting overhead, on the std::atomic broadcast.
        // Batch size is BASE_TRIALS * PERIOD_TRIALS.

        // Number of times to reuse a random base:
        const int BASE_TRIALS = 1U << 16U;

        const bitCapInt threadRange = (cpuCount + nodeMax - nodeMin) / cpuCount;
        // Make sure this is even multiple of 3, minus 1:
        const bitCapInt threadMin = ((nodeMin + threadRange * cpu + 41U) / 42U) * 42U - 1U;
        const bitCapInt threadMax = threadMin + threadRange + 1U;

#if IS_HI_TRIAL_DIVISION
        // First, we uniformly randomly guess offset from any multiple of 5.
        // Then, we're uniformly composing only numbers that are not multiples of 2 or 3.
        std::vector<rand_dist> baseDist(randRange((6U * (threadMax - threadMin) + 104U) / 105U));

        WORD randBitCache = 0;
        rand_dist fiveDist;
        const int maxBatch = (BASE_TRIALS << 1U) / WORD_SIZE;
        const int subBatchSize = WORD_SIZE >> 1U;
#else
        // We're choosing only even multiples of 2 or 3.
        std::vector<rand_dist> baseDist(randRange((threadMax - threadMin + 2U) / 3U));
        const int maxBatch = BASE_TRIALS;
#endif

        for (;;) {
            for (int batchItem = 0U; batchItem < maxBatch; ++batchItem) {
#if IS_HI_TRIAL_DIVISION
                randBitCache = fiveDist(rand_gen);
                for (int subBatchItem = 0U; subBatchItem < subBatchSize; ++subBatchItem) {
#endif
                    // Choose a base at random, >1 and <toFactor.
                    bitCapInt base = baseDist[0U](rand_gen);
#if (QBCAPPOW > 6U) && (!IS_RSA_SEMIPRIME || (QBCAPPOW > 7U))
                    for (size_t i = 1U; i < baseDist.size(); ++i) {
                        base <<= WORD_SIZE;
                        base |= baseDist[i](rand_gen);
                    }
#endif

#if IS_HI_TRIAL_DIVISION
                    // Make this NOT multiple of 7, by adding it to itself divided by six, + 1.
                    base += base / 6U + 1U;

                    // Make this a multiple of 5 (or 0), then randomly make it NOT one,
                    // by adding 1 and uniformly randomly guessing and adding 2 bits.
                    base += (base << 2U) + 1U + (randBitCache & 3U);
                    randBitCache >>= 2U;
#endif
                    // From this, we're composing numbers that are only NOT multiples of 2 or 3.
                    base += threadMin + (base << 1U) - (base & 1U);

#if IS_RSA_SEMIPRIME
                    if ((toFactor % base) == 0U) {
                        isFinished = true;
                        printSuccess(base, toFactor / base, toFactor, "Base has common factor: Found ", iterClock);
                        return;
                    }
#else
                    bitCapInt testFactor = gcd(toFactor, base);
                    if (testFactor != 1U) {
                        isFinished = true;
                        printSuccess(
                            testFactor, toFactor / testFactor, toFactor, "Base has common factor: Found ", iterClock);
                        return;
                    }
#endif
#if IS_HI_TRIAL_DIVISION
                    // We have an extra level of for loop nesting in this macro branch, for IS_HI_TRIAL_DIVISION.
                }
#endif
            }

            // Check if finished, between batches.
            if (isFinished) {
                return;
            }
        }
    };

    std::vector<std::future<void>> futures(cpuCount);
    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        futures[cpu] = std::async(std::launch::async, workerFn, cpu, cpuCount);
    }

    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        futures[cpu].get();
    }

    return 0;
}
