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
#define bci_create(a) (a)
#define bci_copy(a, o) *(o) = a
#define bci_add(l, r, o) *(o) = ((l) + (r))
#define bci_sub(l, r, o) *(o) = ((l) - (r))
#define bci_mul(l, r, o) *(o) = ((l) * (r))
#define bci_div(l, r, o) *(o) = ((l) / (r))
#define bci_mod(l, r, o) *(o) = ((l) % (r))
#define bci_lshift(l, r, o) *(o) = ((l) << (r))
#define bci_rshift(l, r, o) *(o) = ((l) >> (r))
#define bci_or(l, r, o) *(o) = ((l) | (r))
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
#define bci_first_word(a) (a)
#define bci_low64(a) (a)
#else
#define bci_create(a) BigInteger(a)
#define bci_copy(a, o) bi_copy(a, o)
#define bci_add(l, r, o) bi_add(l, r, o)
#define bci_sub(l, r, o) bi_sub(l, r, o)
#define bci_mul(l, r, o) bi_mul(l, r, o)
#define bci_div(l, r, o) bi_div_mod(l, r, o, NULL)
#define bci_mod(l, r, o) bi_div_mod(l, r, NULL, o)
#define bci_lshift(l, r, o) bi_lshift(l, r, o)
#define bci_rshift(l, r, o) bi_rshift(l, r, o)
#define bci_or(l, r, o) bi_or(l, r, o)
#define bci_and(l, r, o) bi_and(l, r, o)
#define bci_eq(l, r) (bi_compare(l, r) == 0)
#define bci_neq(l, r) (bi_compare(l, r) != 0)
#define bci_lt(l, r) (bi_compare(l, r) < 0)
#define bci_gt(l, r) (bi_compare(l, r) > 0)
#define bci_geq(l, r) (bi_compare(l, r) >= 0)
#define bci_leq(l, r) (bi_compare(l, r) <= 0)
#define bci_and_1(l) bi_and_1(l)
#define bci_compare(a, b) (bi_compare(a, b))
#define bci_eq_0(a) (bi_compare_0(a) == 0)
#define bci_neq_0(a) (bi_compare(a, ZERO_BCI) != 0)
#define bci_eq_1(a) (bi_compare(a, ONE_BCI) == 0)
#define bci_neq_1(a) (bi_compare(a, ONE_BCI) != 0)
#define bci_gt_1(a) (bi_compare(a, ONE_BCI) > 0)
#define bci_first_word(a) (a.bits[0])
#define bci_low64(a) (a.bits[0])
#endif
namespace Qimcifa {

const bitCapInt ZERO_BCI(0U);
const bitCapInt ONE_BCI(1U);

#if QBCAPPOW > 6U
std::ostream& operator<<(std::ostream& os, bitCapInt b) {
    const bitCapInt bci_10(10);

    std::vector<std::string> digits;
    while (bci_neq_0(b)) {
        digits.push_back(std::to_string(b.bits[0] % 10));
        bitCapInt t(b);
        bci_div(t, bci_10, &b);
    }
    for (int i = digits.size() - 1; i >= 0; i--) {
        os << digits[i];
    }

    return os;
}

std::istream &operator>>(std::istream &is, bitCapInt& b)
{
    const bitCapInt bci_10(10);

    std::string input;
    is >> input;

    b = ZERO_BCI;
    for (unsigned i = 0; i < input.size(); i++) {
        bitCapInt t(b);
        bci_mul(t, bci_10, &b);
        bci_copy(b, &t);
        bci_add(t, bci_create(input[i] - 48), &b);
    }

    return is;
}
#endif

// Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
inline bool isPowerOfTwo(const bitCapInt& x) {
    bitCapInt t1;
    bci_sub(x, ONE_BCI, &t1);
    bitCapInt t2;
    bci_and(x, t1, &t2);

    return bci_neq_0(x) && bci_eq_0(t2);
}

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
    return bi_log2(n);
#endif
}

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
void uipow(const bitCapInt& base, const bitCapInt& exp, bitCapInt* result)
{
    bi_copy(ONE_BCI, result);
    bitCapInt b(base);
    bitCapInt e(exp);
    for (;;) {
        if (bci_and_1(b)) {
            bitCapInt t(*result);
            bci_mul(t, b, result);
        }
        bitCapInt t(e);
        bci_rshift(t, 1U, &e);
        if (bci_eq_0(e)) {
            break;
        }
        bci_copy(b, &t);
        bci_mul(t, t, &b);
    }
}

// It's fine if this is not exact for the whole bitCapInt domain, so long as it is <= the exact result.
inline bitLenInt intLog(const bitCapInt& base, const bitCapInt& arg)
{
    bitLenInt result = 0U;
    bitCapInt t;
    for (bitCapInt x(arg); bci_geq(x, base); bci_div(t, base, &x)) {
        result++;
        bci_copy(x, &t);
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
    bitCapInt ans = ZERO_BCI, start = ONE_BCI;
    bitCapInt end;
    bci_rshift(x, 1U, &end);
    while (bci_leq(start, end)) {
        bitCapInt t, mid;
        bci_add(start, end, &t);
        bci_rshift(t, 1U, &mid);

        // If x is a perfect square
        bci_mul(mid, mid, &t);
        if (bci_eq(t, x)) {
            return mid;
        }

        if (bci_lt(t, x)) {
            // Since we need floor, we update answer when mid*mid is smaller than x, and move closer to sqrt(x).
            bci_add(mid, ONE_BCI, &start);
            bci_copy(mid, &ans);
        } else {
            // If mid*mid is greater than x
            bci_sub(mid, ONE_BCI, &end);
        }
    }
    return ans;
}

void gcd(bitCapInt n1, bitCapInt n2, bitCapInt* result)
{
    while (bci_neq_0(n2)) {
        bitCapInt t1(n1), t2(n2);
        bci_copy(n2, &n1);
        bci_mod(t1, t2, &n2);
    }
    bci_copy(n1, result);
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

    bitCapInt fullMin;
    bci_lshift(ONE_BCI, primeBits - 1U, &fullMin);

    bitCapInt fullMax, t1;
    bci_lshift(ONE_BCI, primeBits, &t1);
    bci_sub(t1, ONE_BCI, &fullMax);

    bitCapInt minPrime;
    if (primeDict[primeBits].size()) {
        bci_copy(primeDict[primeBits][0], &minPrime);
    } else {
        bci_add(fullMin, ONE_BCI, &minPrime);
    }

    bitCapInt maxPrime;
    if (primeDict[primeBits].size()) {
        bci_copy(primeDict[primeBits][1], &maxPrime);
    } else {
        bci_add(fullMin, ONE_BCI, &maxPrime);
    }

    bitCapInt fullMinR, t2, t3;
    bci_sub(minPrime, ONE_BCI, &t1);
    bci_div(toFactor, minPrime, &t2);
    bci_sub(t2, ONE_BCI, &t3);
    bci_mul(t1, t3, &fullMinR);

    bitCapInt fullMaxR;
    bci_sub(maxPrime, ONE_BCI, &t1);
    bci_div(toFactor, maxPrime, &t2);
    bci_sub(t2, ONE_BCI, &t3);
    bci_mul(t1, t3, &fullMaxR);
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

                const bitCapInt BIG_INT_1(1);

                bitCapInt t1, t2, t3;

                bitCapInt fullRange;
                bci_add(fullMaxR, BIG_INT_1, &t1);
                bci_sub(t1, fullMinR, &fullRange);

                bitCapInt nodeRange, nodeCountBci(nodeCount);
                bci_div(fullRange, nodeCountBci, &nodeRange);

                bitCapInt nodeMin, nodeIdBci(nodeId);
                bci_mul(nodeRange, nodeIdBci, &t1);
                bci_add(fullMinR, t1, &nodeMin);

                bitCapInt nodeMax;
                if ((nodeId + 1U) == nodeCount) {
                    bci_copy(fullMaxR, &nodeMax);
                } else {
                    bitCapInt nodeIdPlus1Bci(nodeId + 1U);
                    bci_mul(nodeRange, nodeIdPlus1Bci, &t1);
                    bci_add(fullMinR, t1, &t2);
                    bci_sub(t2, BIG_INT_1, &nodeMax);
                }

                bitCapInt threadRange, threadsBci(threads);
                bci_add(nodeMax, BIG_INT_1, &t1);
                bci_sub(t1, nodeMin, &t2);
                bci_div(t2, threadsBci, &threadRange);

                bitCapInt rMin, cpuBci(cpu);
                bci_mul(threadRange, cpuBci, &t1);
                bci_add(nodeMin, t1, &rMin);

                bitCapInt rMax;
                if ((cpu + 1U) == threads) {
                    bci_copy(nodeMax, &rMax);
                } else {
                    bitCapInt cpuPlus1Bci(cpu + 1U);
                    bci_mul(threadRange, cpuPlus1Bci, &t1);
                    bci_add(nodeMin, t1, &t2);
                    bci_sub(t2, BIG_INT_1, &rMax);
                }

                std::vector<rand_dist> baseDist;
                std::vector<rand_dist> rDist;
#if QBCAPPOW > 6U
                const bitLenInt wordSize = 64U;
                const uint64_t wordMask = 0xFFFFFFFFFFFFFFFF;
                const bitCapInt BIG_INT_3(3);
                bitCapInt distPart;
                bci_sub(toFactor, BIG_INT_3, &distPart);
                while (bci_neq_0(distPart)) {
                    baseDist.push_back(rand_dist(0U, bci_low64(distPart)));
                    bci_copy(distPart, &t1);
                    bci_rshift(t1, wordSize, &distPart);
                }
                std::reverse(rDist.begin(), rDist.end());

                bci_sub(rMax, rMin, &distPart);
                while (bci_neq_0(distPart)) {
                    rDist.push_back(rand_dist(0U, (uint64_t)(distPart.bits[0] & wordMask)));
                    bci_copy(distPart, &t1);
                    bci_rshift(t1, wordSize, &distPart);
                }
                std::reverse(rDist.begin(), rDist.end());
#else
                baseDist.push_back(rand_dist(2U, toFactor.bits[0]));
                rDist.push_back(rand_dist(rMin.bits[0], rMax.bits[0]));
#endif

                for (;;) {
                    for (size_t batchItem = 0U; batchItem < BASE_TRIALS; batchItem++) {
                        // Choose a base at random, >1 and <toFactor.
                        bitCapInt base(baseDist[0U](rand_gen));
#if QBCAPPOW > 6U
                        for (size_t i = 1U; i < baseDist.size(); i++) {
                            bi_copy(base, &t1);
                            bci_lshift(t1, wordSize, &base);
                            base.bits[0] = baseDist[i](rand_gen);
                        }
                        bi_copy(base, &t1);
                        const bitCapInt BIG_INT_2(2);
                        bci_add(t1, BIG_INT_2, &base);
#endif

                        gcd(toFactor, base, &t1);
                        if (bci_neq_1(t1)) {
                            // Inform the other threads on this node that we've succeeded and are done:
                            isFinished = true;

                            bci_div(toFactor, t1, &t2);

                            std::cout << "Chose non-relative prime: " << t1 << " * " << t2
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
                            bitCapInt r(rDist[0U](rand_gen));
#if QBCAPPOW > 6U
                            for (size_t i = 1U; i < rDist.size(); i++) {
                                bi_copy(r, &t1);
                                bci_lshift(t1, wordSize, &r);
                                r.bits[0] = rDist[i](rand_gen);
                            }
                            bi_copy(r, &t1);
                            bci_add(t1, rMin, &r);
#endif
                            // Since our output is r rather than y, we can skip the continued fractions step.
                            if (bci_and_1(r)) {
                                bci_copy(r, &t1);
                            } else {
                                bci_rshift(r, 1U, &t1);
                            }

#define PRINT_SUCCESS(f1, f2, toFactor, message)                                                                       \
    std::cout << message << (f1) << " * " << (f2) << " = " << (toFactor) << std::endl;                                 \
    auto tClock =                                                                                                      \
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);  \
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;                            \
    std::cout << "(Waiting to join other threads...)" << std::endl;

#if IS_RSA_SEMIPRIME
#define RGUESS t1
#else
#define RGUESS r
#endif

                            // As a "classical" optimization, since \phi(toFactor) and factor bounds overlap,
                            // we first check if our guess for r is already a factor.
                            bci_mod(toFactor, RGUESS, &t2);
                            if (bci_gt_1(RGUESS) && bci_eq_0(t2)) {
                                // Inform the other threads on this node that we've succeeded and are done:
                                isFinished = true;

                                bci_div(toFactor, RGUESS, &t2);

                                PRINT_SUCCESS(
                                    RGUESS, t2, toFactor, "Success (on r trial division): Found ");
                                return;
                            }

                            bitCapInt f1, f2, fmul;
                            uipow(base, t1, &t2);
                            bci_mod(t2, toFactor, &t1);

                            bci_add(t1, BIG_INT_1, &t2);
                            //gcd(t2, toFactor, &f1);

                            bci_sub(t1, BIG_INT_1, &t2);
                            //gcd(t2, toFactor, &f2);
#if 0
                            bci_mul(f1, f2, &fmul);

                            bci_mod(toFactor, fmul, &t1);

                            while (bci_gt_1(fmul) && bci_neq(fmul, toFactor) && bci_eq_0(t1)) {
                                bci_copy(f1, &fmul);
                                bci_mul(fmul, f2, &f1);

                                bci_mul(fmul, f2, &t1);
                                bci_div(toFactor, t1, &f2);

                                bci_mul(f1, f2, &fmul);

                                bci_mod(toFactor, fmul, &t1);
                            }

                            if (bci_gt_1(fmul) && bci_eq(fmul, toFactor) && bci_gt_1(f1) && bci_gt_1(f2)) {
                                // Inform the other threads on this node that we've succeeded and are done:
                                isFinished = true;

                                PRINT_SUCCESS(f1, f2, toFactor,
                                    "Success (on r difference of squares): Found ");
                                return;
                            }
#endif
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
