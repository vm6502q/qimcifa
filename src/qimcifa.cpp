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
#define bitsInByte 8U

#if QBCAPPOW < 8U
#define bitLenInt uint8_t
#elif QBCAPPOW < 16U
#define bitLenInt uint16_t
#elif QBCAPPOW < 32U
#define bitLenInt uint32_t
#else
#define bitLenInt uint64_t
#endif

#if QBCAPPOW < 6U
#define bitsInCap 32
#define bitCapInt uint32_t
#elif QBCAPPOW < 7U
#define bitsInCap 64
#define bitCapInt uint64_t
#else
#define bitsInCap (8U * (((bitLenInt)1U) << QBCAPPOW))
#include "big_integer.hpp"
#define bitCapInt BigInteger
#endif

#if QBCAPPOW < 7U
#define bci_copy(a) (a)
#define bci_create(a) (a)
#define bci_add(l, r) ((l) + (r))
#define bci_sub(l, r) ((l) - (r))
#define bci_mul(l, r) ((l) * (r))
#define bci_div(l, r) ((l) / (r))
#define bci_mod(l, r) ((l) % (r))
#define bci_lshift(l, r) ((l) << (r))
#define bci_rshift(l, r) ((l) >> (r))
#define bci_or(l, r) ((l) | (r))
#define bci_and(l, r) ((l) & (r))
#define bci_eq(l, r) ((l) == (r))
#define bci_neq(l, r) ((l) != (r))
#define bci_lt(l, r) ((l) < (r))
#define bci_gt(l, r) ((l) > (r))
#define bci_geq(l, r) ((l) >= (r))
#define bci_leq(l, r) ((l) <= (r))
#define bci_and_1(l) ((l) & 1U)
#define bci_compare(a, b) ((a > b) ? 1 : ((a < b) ? -1 : 0))
#define bci_compare_0(a) ((a > 0) ? 1 : ((a < 0) ? -1 : 0))
#define bci_eq_0(a) ((a) == 0)
#define bci_neq_0(a) ((a) != 0)
#define bci_eq_1(a) ((a) == 1)
#define bci_gt_1(a) ((a) > 1)
#define bci_first_word(a) (a)
#define bci_low64(a) (a)
#else
#define bci_copy(a) BigInteger(a)
#define bci_create(a) bi_create(a)
#define bci_add(l, r) bi_add(l, r)
#define bci_sub(l, r) bi_sub(l, r)
#define bci_mul(l, r) bi_mul(l, r)
#define bci_div(l, r) bi_div(l, r)
#define bci_mod(l, r) bi_mod(l, r)
#define bci_lshift(l, r) bi_lshift(l, r)
#define bci_rshift(l, r) bi_rshift(l, r)
#define bci_or(l, r) bi_or(l, r)
#define bci_and(l, r) bi_and(l, r)
#define bci_eq(l, r) (bi_compare(l, r) == 0)
#define bci_neq(l, r) (bi_compare(l, r) != 0)
#define bci_lt(l, r) (bi_compare(l, r) < 0)
#define bci_gt(l, r) (bi_compare(l, r) > 0)
#define bci_geq(l, r) (bi_compare(l, r) >= 0)
#define bci_leq(l, r) (bi_compare(l, r) <= 0)
#define bci_and_1(l) bi_and_1(l)
#define bci_compare(a, b) (bi_compare(a, b))
#define bci_compare_0(a) (bi_compare_0(a))
#define bci_eq_0(a) (bi_compare(a, ZERO_BCI) == 0)
#define bci_neq_0(a) (bi_compare(a, ZERO_BCI) != 0)
#define bci_eq_1(a) (bi_compare(a, ONE_BCI) == 0)
#define bci_gt_1(a) (bi_compare(a, ONE_BCI) > 0)
#define bci_first_word(a) (a.bits[0])
#define bci_low64(a) (a.bits[0])
#endif
namespace Qimcifa {

const bitCapInt ZERO_BCI = bci_create(0U);
const bitCapInt ONE_BCI = bci_create(1U);

#if QBCAPPOW > 6U
std::ostream& operator<<(std::ostream& os, const bitCapInt& bci) {
    const bitCapInt bci_10 = bci_create(10);

    bitCapInt b = bci_copy(bci);
    std::vector<std::string> digits;
    while (bci_neq_0(b)) {
        digits.push_back(std::to_string(b.bits[0] % 10));
        b = bci_div(b, bci_10);
    }
    for (int i = digits.size() - 1; i >= 0; i--) {
        os << digits[i];
    }

    return os;
}

std::istream &operator>>(std::istream &is, bitCapInt& b)
{
    const bitCapInt bci_10 = bci_create(10);

    std::string input;
    is >> input;

    b = ZERO_BCI;
    for (unsigned i = 0; i < input.size(); i++) {
        b = bci_mul(b, bci_10);
        b = bci_add(b, bci_create(input[i] - 48));
    }

    return is;
}
#endif

// Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
inline bool isPowerOfTwo(const bitCapInt& x) { return (bci_neq_0(x) && bci_eq_0(bci_and(x, bci_sub(x, ONE_BCI)))); }

inline bitLenInt log2(const bitCapInt& n)
{
#if __GNUC__ && QBCAPPOW < 7
// Source: https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers#answer-11376759
#if QBCAPPOW < 6
    return (bitLenInt)(bitsInByte * sizeof(unsigned int) - __builtin_clz((unsigned int)n) - 1U);
#else
    return (bitLenInt)(bitsInByte * sizeof(unsigned long long) - __builtin_clzll((unsigned long long)n) - 1U);
#endif
#else
    bitLenInt pow = 0U;
    bitCapInt p = bci_rshift(n, 1U);
    while (bci_neq_0(p)) {
        p = bci_rshift(p, 1U);
        pow++;
    }
    return pow;
#endif
}

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
inline bitCapInt uipow(const bitCapInt& base, const bitCapInt& exp)
{
    bitCapInt result = ONE_BCI;
    bitCapInt b = base;
    bitCapInt e = exp;
    for (;;) {
        if (bci_and_1(b)) {
            result = bci_mul(result, b);
        }
        e = bci_rshift(e, 1U);
        if (bci_eq_0(e)) {
            break;
        }
        b = bci_mul(b, b);
    }

    return result;
}

// It's fine if this is not exact for the whole bitCapInt domain, so long as it is <= the exact result.
inline bitLenInt intLog(const bitCapInt& base, const bitCapInt& arg)
{
    bitLenInt result = 0U;
    for (bitCapInt x = bci_copy(arg); bci_geq(x, base); x = bci_div(x, base)) {
        result++;
    }
    return result;
}

// Adapted from Gaurav Ahirwar's suggestion on https://www.geeksforgeeks.org/square-root-of-an-integer/
bitCapInt floorSqrt(const bitCapInt& x)
{
    // Base cases
    if (bci_eq_0(x) || bci_eq_1(x)) {
        return x;
    }

    // Binary search for floor(sqrt(x))
    bitCapInt start = ONE_BCI, end = bci_rshift(x, 1U), ans = ZERO_BCI;
    while (bci_leq(start, end)) {
        bitCapInt mid = bci_rshift(bci_add(start, end), 1U);

        // If x is a perfect square
        bitCapInt sqr = bci_mul(mid, mid);
        if (bci_eq(sqr, x)) {
            return mid;
        }

        if (bci_lt(sqr, x)) {
            // Since we need floor, we update answer when mid*mid is smaller than x, and move closer to sqrt(x).
            start = bci_add(mid, ONE_BCI);
            ans = bci_copy(mid);
        } else {
            // If mid*mid is greater than x
            end = bci_sub(mid, ONE_BCI);
        }
    }
    return ans;
}

bitCapInt gcd(const bitCapInt& n1, const bitCapInt& n2)
{
    if (bci_neq_0(n2)) {
        return gcd(n2, bci_mod(n1, n2));
    }
    return n1;
}

} // namespace Qimcifa

using namespace Qimcifa;

int main()
{
    typedef std::uniform_int_distribution<uint64_t> rand_dist;

    bitCapInt toFactor;
    uint64_t nodeCount = 1U;
    uint64_t nodeId = 0U;

    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    auto iterClock = std::chrono::high_resolution_clock::now();

    const bitLenInt qubitCount = log2(toFactor) + (isPowerOfTwo(toFactor) ? 0U : 1U);
    // const bitCapInt qubitPower = ONE_BCI << qubitCount;
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

    std::random_device rand_dev;
    std::mt19937 rand_gen(rand_dev());

    std::atomic<bool> isFinished;
    isFinished = false;

#if IS_RSA_SEMIPRIME
    std::map<bitLenInt, const std::vector<bitCapInt>> primeDict = { { 16U, { bci_create(32771U), bci_create(65521U) } },
        { 28U, { bci_create(134217757U), bci_create(268435399U) } }, { 32U, { bci_create(2147483659U), bci_create(4294967291U) } },
        { 64U, { bci_create(9223372036854775837U), bci_create(1844674407370955143U) } } };

    // If n is semiprime, \phi(n) = (p - 1) * (q - 1), where "p" and "q" are prime.
    // The minimum value of this formula, for our input, without consideration of actual
    // primes in the interval, is as follows:
    // (See https://www.mobilefish.com/services/rsa_key_generation/rsa_key_generation.php)
    const bitLenInt primeBits = (qubitCount + 1U) >> 1U;
    const bitCapInt fullMin = bci_lshift(ONE_BCI, primeBits - 1U);
    const bitCapInt fullMax = bci_sub(bci_lshift(ONE_BCI, primeBits), ONE_BCI);
    const bitCapInt minPrime = primeDict[primeBits].size() ? primeDict[primeBits][0] : bci_add(fullMin, ONE_BCI);
    const bitCapInt maxPrime = primeDict[primeBits].size() ? primeDict[primeBits][1] : fullMax;
    const bitCapInt fullMinR = bci_mul(bci_sub(minPrime, ONE_BCI), bci_sub(bci_div(toFactor, minPrime), ONE_BCI));
    const bitCapInt fullMaxR = bci_mul(bci_sub(maxPrime, ONE_BCI), bci_sub(bci_div(toFactor, maxPrime), ONE_BCI));
#else
    // \phi(n) is Euler's totient for n. A loose lower bound is \phi(n) >= sqrt(n/2).
    const bitCapInt fullMinR = floorSqrt(bci_rshift(toFactor, 1U));
    // A better bound is \phi(n) >= pow(n / 2, log(2)/log(3))
    // const bitCapInt fullMinR = pow(toFactor / 2, PHI_EXPONENT);

    // It can be shown that the period of this modular exponentiation can be no higher than 1
    // less than the modulus, as in https://www2.math.upenn.edu/~mlazar/math170/notes06-3.pdf.
    // Further, an upper bound on Euler's totient for composite numbers is n - sqrt(n). (See
    // https://math.stackexchange.com/questions/896920/upper-bound-for-eulers-totient-function-on-composite-numbers)
    const bitCapInt fullMaxR = bci_sub(toFactor, floorSqrt(toFactor));
#endif

    const unsigned cpuCount = std::thread::hardware_concurrency();
    std::vector<std::future<void>> futures(cpuCount);
    for (unsigned cpu = 0U; cpu < cpuCount; cpu++) {
        futures[cpu] = std::async(
            std::launch::async, [cpu, nodeId, nodeCount, toFactor, fullMinR, fullMaxR, &iterClock, &rand_gen, &isFinished] {
                // These constants are semi-redundant, but they're only defined once per thread,
                // and compilers differ on lambda expression capture of constants.

                // Batching reduces mutex-waiting overhead, on the std::atomic broadcast.
                // Batch size is BASE_TRIALS * PERIOD_TRIALS.

                // Number of times to reuse a random base:
                const size_t BASE_TRIALS = 1U << 9U;
                // Number of random period guesses per random base:
                const size_t PERIOD_TRIALS = 1U;

                const double clockFactor = 1.0 / 1000.0; // Report in ms
                const unsigned threads = std::thread::hardware_concurrency();

                const BigInteger BIG_INT_1 = bci_create(1);
                const BigInteger BIG_INT_2 = bci_create(2);

                const bitCapInt fullRange = bci_sub(bci_add(fullMaxR, BIG_INT_1), fullMinR);
                const bitCapInt nodeRange = bci_div(fullRange, bci_create(nodeCount));
                const bitCapInt nodeMin = bci_add(fullMinR, bci_mul(nodeRange, bci_create(nodeId)));
                const bitCapInt nodeMax = ((nodeId + 1U) == nodeCount)
                    ? fullMaxR
                    : bci_sub(bci_add(fullMinR, bci_mul(nodeRange, bci_create(nodeId + 1U))), BIG_INT_1);
                const bitCapInt threadRange = bci_div(bci_sub(bci_add(nodeMax, BIG_INT_1), nodeMin), bci_create(threads));
                const bitCapInt rMin = bci_add(nodeMin, bci_mul(threadRange, bci_create(cpu)));
                const bitCapInt rMax = ((cpu + 1U) == threads)
                    ? nodeMax
                    : bci_sub(bci_add(nodeMin, bci_mul(threadRange, bci_create(cpu + 1U))), BIG_INT_1);

                std::vector<rand_dist> baseDist;
                std::vector<rand_dist> rDist;
#if QBCAPPOW < 6U
                baseDist.push_back(rand_dist(2U, bci_sub(toFactor, BIG_INT_1));
                rDist.push_back(rand_dist(rMin, rMax));
#else
                const bitLenInt wordSize = 64U;
                const uint64_t wordMask = 0xFFFFFFFFFFFFFFFF;
                bitCapInt distPart = bci_sub(toFactor, bci_create(3));
                while (bci_compare_0(distPart) != 0) {
                    baseDist.push_back(rand_dist(0U, bci_low64(distPart)));
                    bci_rshift(distPart, wordSize);
                }
                std::reverse(rDist.begin(), rDist.end());

                distPart = bci_sub(rMax, rMin);
                while (bci_compare_0(distPart) != 0) {
                    rDist.push_back(rand_dist(0U, (uint64_t)(distPart.bits[0] & wordMask)));
                    distPart = bci_rshift(distPart, wordSize);
                }
                std::reverse(rDist.begin(), rDist.end());
#endif

                for (;;) {
                    for (size_t batchItem = 0U; batchItem < BASE_TRIALS; batchItem++) {
                        // Choose a base at random, >1 and <toFactor.
                        bitCapInt base = bci_create(baseDist[0U](rand_gen));
#if QBCAPPOW > 5U
                        for (size_t i = 1U; i < baseDist.size(); i++) {
                            base = bci_lshift(base, wordSize);
                            base.bits[0] = baseDist[i](rand_gen);
                        }
                        base = bci_add(base, BIG_INT_2);
#endif

                        const bitCapInt testFactor = gcd(toFactor, base);
                        if (bci_compare(testFactor, BIG_INT_1) != 0) {
                            // Inform the other threads on this node that we've succeeded and are done:
                            isFinished = true;

                            std::cout << "Chose non-relative prime: " << testFactor << " * " << bci_div(toFactor, testFactor)
                                      << std::endl;
                            auto tClock = std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now() - iterClock);
                            std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;
                            std::cout << "(Waiting to join other threads...)" << std::endl;
                            return;
                        }

                        // This would be where we perform the quantum period finding algorithm.
                        // However, we don't have a quantum computer!
                        // Instead, we "throw dice" for a guess to the output of the quantum subroutine.
                        // This guess will usually be wrong, at least for semi-prime inputs.
                        // If we try many times, though, this can be a practically valuable factoring method.

                        // y is meant to be close to some number c * qubitPower / r, where r is the period.
                        // c is a positive integer or 0, and we don't want the 0 case.
                        // y is truncated by the number of qubits in the register, at most.
                        // The maximum value of c before truncation is no higher than r.

                        // The period of ((base ^ x) MOD toFactor) can't be smaller than log_base(toFactor).
                        // (Also, toFactor is definitely NOT an exact multiple of base.)
                        // const bitCapInt logBaseToFactor = (bitCapInt)intLog(base, toFactor) + 1U;
                        // Euler's Theorem tells us, if gcd(a, n) = 1, then a^\phi(n) = 1 MOD n,
                        // where \phi(n) is Euler's totient for n.
                        // const bitCapInt fullMinR = (minPhi < logBaseToFactor) ? logBaseToFactor : minPhi;

                        // c is basically a harmonic degeneracy factor, and there might be no value in testing
                        // any case except c = 1, without loss of generality.

                        // This sets a nonuniform distribution on our y values to test.
                        // y values are close to qubitPower / rGuess, and we midpoint round.

                        // However, results are better with uniformity over r, rather than y.

                        // So, we guess r, between fullMinR and fullMaxR.
                        for (size_t rTrial = 0U; rTrial < PERIOD_TRIALS; rTrial++) {
                            // Choose a base at random, >1 and <toFactor.
                            bitCapInt r = bci_create(rDist[0U](rand_gen));
#if QBCAPPOW > 5U
                            for (size_t i = 1U; i < rDist.size(); i++) {
                                r = bci_lshift(r, wordSize);
                                r.bits[0] = rDist[i](rand_gen);
                            }
                            r = bci_add(r, rMin);
#endif
                            // Since our output is r rather than y, we can skip the continued fractions step.
                            const bitCapInt p = bci_and_1(r) ? r : bci_rshift(r, 1U);

#define PRINT_SUCCESS(f1, f2, toFactor, message)                                                                       \
    std::cout << message << (f1) << " * " << (f2) << " = " << (toFactor) << std::endl;                                 \
    auto tClock =                                                                                                      \
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);  \
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;                            \
    std::cout << "(Waiting to join other threads...)" << std::endl;

#if IS_RSA_SEMIPRIME
#define RGUESS p
#else
#define RGUESS r
#endif

                            // As a "classical" optimization, since \phi(toFactor) and factor bounds overlap,
                            // we first check if our guess for r is already a factor.
                            if ((bci_compare(RGUESS, BIG_INT_1) > 0) && (bci_compare(toFactor, bci_mul(bci_div(toFactor, RGUESS), RGUESS)) == 0)) {
                                // Inform the other threads on this node that we've succeeded and are done:
                                isFinished = true;

                                PRINT_SUCCESS(
                                    RGUESS, bci_div(toFactor, RGUESS), toFactor, "Success (on r trial division): Found ");
                                return;
                            }

                            const bitCapInt apowrhalf = bci_mod(uipow(base, p), toFactor);
                            bitCapInt f1 = gcd(bci_add(apowrhalf, BIG_INT_1), toFactor);
                            bitCapInt f2 = gcd(bci_sub(apowrhalf, BIG_INT_1), toFactor);
                            bitCapInt fmul = bci_mul(f1, f2);
                            while ((bci_compare(fmul, BIG_INT_1) > 0) && (bci_compare(fmul, toFactor) != 0) && (bci_compare(toFactor, bci_mul(bci_div(toFactor, fmul), fmul)) == 0)) {
                                fmul = f1;
                                f1 = bci_mul(fmul, f2);
                                f2 = bci_div(toFactor, bci_mul(fmul, f2));
                                fmul = bci_mul(f1, f2);
                            }
                            if ((bci_compare(fmul, BIG_INT_1) > 0) && (bci_compare(fmul, toFactor) == 0) && (bci_compare(f1, BIG_INT_1) > 0) && (bci_compare(f2, BIG_INT_1) > 0)) {
                                // Inform the other threads on this node that we've succeeded and are done:
                                isFinished = true;

                                PRINT_SUCCESS(f1, f2, toFactor,
                                    "Success (on r difference of squares): Found ");
                                return;
                            }
                        }
                    }

                    // Check if finished, between batches.
                    if (isFinished) {
                        return;
                    }
                }
            });
    };

    for (unsigned cpu = 0U; cpu < cpuCount; cpu++) {
        futures[cpu].get();
    }

    return 0;
}
