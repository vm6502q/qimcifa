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

#if (QBCAPPOW < 6U) || (IS_RSA_SEMIPRIME && (QBCAPPOW < 7U))
#define WORD uint32_t
#define WORD_SIZE 32U
#else
#define WORD uint64_t
#define WORD_SIZE 64U
#endif


#if QBCAPPOW < 6U
#define bitsInCap 32
#define bitCapInt uint32_t
#elif QBCAPPOW < 7U
#define bitsInCap 64
#define bitCapInt uint64_t
#else
#define bitsInCap (8U * (((bitLenInt)1U) << QBCAPPOW))
#include "big_integer.h"
#define bitCapInt BigInteger
#endif

#if QBCAPPOW < 7U
#define bci_create(a) (a)
#define bci_copy(a, o) *(o) = a
#define bci_set_0(a) *(a) = 0
#define bci_increment(l, r) *(l) += r
#define bci_add_i(l, r) *(l) += (r)
#define bci_add(l, r, o) *(o) = ((l) + (r))
#define bci_decrement(l, r) *(l) -= r
#define bci_sub_i(l, r) *(l) -= (r)
#define bci_sub(l, r, o) *(o) = ((l) - (r))
#define bci_mul(l, r, o) *(o) = ((l) * (r))
#define bci_div(l, r, o) *(o) = ((l) / (r))
#define bci_mod(l, r, o) *(o) = ((l) % (r))
#define bci_lshift(l, r, o) *(o) = ((l) << (r))
#define bci_lshift_ip(l, r) *(l) <<= (r)
#define bci_rshift(l, r, o) *(o) = ((l) >> (r))
#define bci_rshift_ip(l, r) *(l) >>= (r)
#define bci_and(l, r, o) *(o) = ((l) & (r))
#define bci_eq(l, r) ((l) == (r))
#define bci_neq(l, r) ((l) != (r))
#define bci_lt(l, r) ((l) < (r))
#define bci_gt(l, r) ((l) > (r))
#define bci_geq(l, r) ((l) >= (r))
#define bci_leq(l, r) ((l) <= (r))
#define bci_and_1(l) ((l) & 1U)
#define bci_compare(a, b) ((a > b) ? 1 : ((a < b) ? -1 : 0))
#define bci_eq_0(a) ((a) == 0)
#define bci_neq_0(a) ((a) != 0)
#define bci_eq_1(a) ((a) == 1)
#define bci_neq_1(a) ((a) != 1)
#define bci_gt_1(a) ((a) > 1)
#define bci_low64(a) (a)
#else
#define bci_create(a) bi_create(a)
#define bci_copy(a, o) bi_copy(&(a), o)
#define bci_set_0(a) bi_set_0(a)
#define bci_increment(l, r) bi_increment(l, r)
#define bci_add_ip(l, r) bi_add_ip(l, &(r))
#define bci_add(l, r, o) bi_add(&(l), &(r), o)
#define bci_decrement(l, r) bi_decrement(l, r)
#define bci_sub(l, r, o) bi_sub(&(l), &(r), o)
#define bci_sub_ip(l, r) bi_sub_ip(l, &(r))
#define bci_mul(l, r, o) bi_mul(&(l), &(r), o)
#define bci_div(l, r, o) bi_div_mod(&(l), &(r), o, NULL)
#define bci_mod(l, r, o) bi_div_mod(&(l), &(r), NULL, o)
#define bci_lshift(l, r, o) bi_lshift(&(l), r, o)
#define bci_lshift_ip(l, r) bi_lshift_ip(l, r)
#define bci_rshift(l, r, o) bi_rshift(&(l), r, o)
#define bci_rshift_ip(l, r) bi_rshift_ip(l, r)
#define bci_and(l, r, o) bi_and(&(l), &(r), o)
#define bci_eq(l, r) (bi_compare(&(l), &(r)) == 0)
#define bci_neq(l, r) (bi_compare(&(l), &(r)) != 0)
#define bci_lt(l, r) (bi_compare(&(l), &(r)) < 0)
#define bci_gt(l, r) (bi_compare(&(l), &(r)) > 0)
#define bci_geq(l, r) (bi_compare(&(l), &(r)) >= 0)
#define bci_leq(l, r) (bi_compare(&(l), &(r)) <= 0)
#define bci_and_1(l) bi_and_1(&(l))
#define bci_compare(a, b) (bi_compare(&(a), &(b)))
#define bci_eq_0(a) (bi_compare_0(&(a)) == 0)
#define bci_neq_0(a) (bi_compare_0(&(a)) != 0)
#define bci_eq_1(a) (bi_compare_1(&(a)) == 0)
#define bci_neq_1(a) (bi_compare_1(&(a)) != 0)
#define bci_gt_1(a) (bi_compare_1(&(a)) > 0)
#define bci_low64(a) ((a).bits[0])
#endif

namespace Qimcifa {

const bitCapInt ZERO_BCI = bci_create(0U);
const bitCapInt ONE_BCI = bci_create(1U);
const bitCapInt BIG_INT_10 = bci_create(10U);

#if QBCAPPOW > 6U
std::ostream& operator<<(std::ostream& os, bitCapInt b)
{
    bitCapInt t;

    // Calculate the base-10 digits, from lowest to highest.
    std::vector<std::string> digits;
    while (bci_neq_0(b)) {
        digits.push_back(std::to_string((unsigned char)(bci_low64(b) % 10U)));
        bci_copy(b, &t);
        bci_div(t, BIG_INT_10, &b);
    }

    if (digits.size() == 0) {
        return os;
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
    bitCapInt t;

    // Get the whole input string at once.
    std::string input;
    is >> input;

    // Start the output address value at 0.
    bi_set_0(&b);
    for (size_t i = 0; i < input.size(); ++i) {
        // Left shift by 1 base-10 digit.
        bci_copy(b, &t);
        bci_mul(t, BIG_INT_10, &b);
        // Add the next lowest base-10 digit.
        bci_increment(&b, input[i] - 48U);
    }

    return is;
}
#endif

// Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
bool isPowerOfTwo(const bitCapInt& x) {
    bitCapInt t1;
    bci_sub(x, ONE_BCI, &t1);
    bitCapInt t2;
    bci_and(x, t1, &t2);

    return bci_neq_0(x) && bci_eq_0(t2);
}

bitLenInt log2(const bitCapInt& n)
{
#if __GNUC__ && QBCAPPOW < 7
// Source: https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers#answer-11376759
#if QBCAPPOW < 6
    return (bitLenInt)(bitsInByte * sizeof(unsigned int) - __builtin_clz((unsigned int)n) - 1U);
#else
    return (bitLenInt)(bitsInByte * sizeof(unsigned long long) - __builtin_clzll((unsigned long long)n) - 1U);
#endif
#else
    return bi_log2(&n);
#endif
}

void gcd(bitCapInt n1, bitCapInt n2, bitCapInt* result)
{
    bitCapInt t1, t2;
    while (bci_neq_0(n2)) {
        bci_copy(n1, &t1);
        bci_copy(n2, &t2);
        bci_copy(n2, &n1);
        bci_mod(t1, t2, &n2);
    }
    bci_copy(n1, result);
}

typedef std::uniform_int_distribution<WORD> rand_dist;

std::vector<rand_dist> rangeRange(bitCapInt range) {
    std::vector<rand_dist> distToReturn;
    while (bci_neq_0(range)) {
        distToReturn.push_back(rand_dist(0U, (WORD)bci_low64(range)));
        bci_rshift_ip(&range, WORD_SIZE);
    }
    std::reverse(distToReturn.begin(), distToReturn.end());

    return distToReturn;
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

#if IS_RSA_SEMIPRIME
    std::map<bitLenInt, const std::vector<bitCapInt>> primeDict = {
        { 16U, { bci_create(16411U), bci_create(65521U) } },
        { 28U, { bci_create(67108879U), bci_create(536870909U) } },
        { 32U, { bci_create(1073741827U), bci_create(8589934583U) } }
    };

    const bitLenInt primeBits = (qubitCount + 1U) >> 1U;
    bitCapInt minPrime = bci_create(1), maxPrime = bci_create(1);
    if (primeDict[primeBits].size()) {
        bci_copy(primeDict[primeBits][0], &minPrime);
        bci_copy(primeDict[primeBits][1], &maxPrime);
    } else {
        bci_lshift_ip(&minPrime, primeBits - 2U);
        bci_increment(&minPrime, 1U);

        bci_lshift_ip(&maxPrime, primeBits + 1U);
        bci_decrement(&maxPrime, 1U);
    }

    // We're only picking odd numbers, so divide the boundaries both by 2.
    bitCapInt fullMinBase;
    bci_div(toFactor, maxPrime, &fullMinBase);
    if (bci_compare(fullMinBase, minPrime) == -1) {
        bci_copy(minPrime, &fullMinBase);
    }

    bitCapInt fullMaxBase;
    bci_div(toFactor, minPrime, &fullMaxBase);
    if (bci_compare(fullMinBase, maxPrime) == 1) {
        bci_copy(maxPrime, &fullMinBase);
    }
#else
    bitCapInt fullMinBase = bci_create(5U);
    bitCapInt fullMaxBase;
    bci_div(toFactor, fullMinBase, &fullMaxBase);
#endif
    const bitCapInt nodeCountBci = bci_create(nodeCount);
    const bitCapInt nodeIdBci = bci_create(nodeId);

    bitCapInt nodeRange, t;
    bci_sub(fullMaxBase, fullMinBase, &t);
    bci_increment(&t, nodeCount - 1U);
    bci_div(t, nodeCountBci, &nodeRange);

    bitCapInt nodeMin;
    bci_mul(nodeRange, nodeIdBci, &nodeMin);
    bci_add_ip(&nodeMin, fullMinBase);

    bitCapInt nodeMax;
    bci_add(nodeMin, nodeRange, &nodeMax);

    const double clockFactor = 1.0 / 1000.0; // Report in ms
    auto iterClock = std::chrono::high_resolution_clock::now();

    std::random_device rand_dev;
    std::mt19937 rand_gen(rand_dev());

    const unsigned cpuCount = std::thread::hardware_concurrency();
    std::atomic<bool> isFinished;
    isFinished = false;

    auto workerFn = [&toFactor, &nodeMin, &nodeMax, &iterClock, &rand_gen, &isFinished](int cpu, unsigned cpuCount) {
        // These constants are semi-redundant, but they're only defined once per thread,
        // and compilers differ on lambda expression capture of constants.

        // Batching reduces mutex-waiting overhead, on the std::atomic broadcast.
        // Batch size is BASE_TRIALS * PERIOD_TRIALS.

        // Number of times to reuse a random base:
        const size_t BASE_TRIALS = 1U << 8U;

        const bitCapInt cpuCountBci = bci_create(cpuCount);
        const bitCapInt cpuBci = bci_create(cpu);
        const bitCapInt BIG_INT_6 = bci_create(6U);

        bitCapInt threadRange, t;
        bci_sub(nodeMax, nodeMin, &t);
        bci_increment(&t, cpuCount);
        bci_div(t, cpuCountBci, &threadRange);

        bitCapInt threadMin;
        bci_mul(threadRange, cpuBci, &threadMin);
        bci_add_ip(&threadMin, nodeMin);
        bci_increment(&threadMin, 5U);
        bci_div(threadMin, BIG_INT_6, &t);
        bci_mul(t, BIG_INT_6, &threadMin);
        bci_decrement(&threadMin, 1U);

        bitCapInt threadMax;
        bci_add(threadMin, threadRange, &threadMax);
        bci_increment(&threadMax, 1U);

        std::vector<rand_dist> rDist(rangeRange(threadRange));

        bitCapInt base;

        for (;;) {
            for (size_t batchItem = 0U; batchItem < BASE_TRIALS; ++batchItem) {
                // Choose a base at random, >1 and <toFactor.
                bci_set_0(&base);
                for (size_t i = 0U; i < rDist.size(); i++) {
                    bci_lshift_ip(&base, WORD_SIZE);
                    bci_low64(base) = rDist[i](rand_gen);
                }
                bci_copy(base, &t);
                bci_lshift_ip(&base, 1U);
                bci_add_ip(&base, t);
                bci_add_ip(&base, threadMin);
                if (bci_low64(t) & 1U) {
                    bci_decrement(&base, 1U);
                }

                gcd(toFactor, base, &t);
                if (bci_neq_1(t)) {
                    isFinished = true;
                    return t;
                }
            }

            // Check if finished, between batches.
            if (isFinished) {
                return ZERO_BCI;
            }
        }
    };

    std::vector<std::future<bitCapInt>> futures(cpuCount);
    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        futures[cpu] = std::async(std::launch::async, workerFn, cpu, cpuCount);
    }

    bitCapInt f1, f2;
    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        f1 = futures[cpu].get();
        if (bci_neq_0(f1)) {
            bci_div(toFactor, f1, &f2);
            std::cout << "Base has common factor: Found ";
            // TODO: Serious problem with the output stream memory usage and execution time:
            // std::cout << f1 << " * " << f2 << " = " << toFactor << std::endl;
            auto tClock = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
            std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;
            std::cout << "(Waiting to join other threads...)" << std::endl;
        }
    }

    return 0;
}
