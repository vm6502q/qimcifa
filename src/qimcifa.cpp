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

#include "qimcifa.hpp"

namespace Qimcifa {

int mainBody(const BigInteger& toFactor)
{
    // First 9 primes
    std::vector<unsigned> trialDivisionPrimes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

    const unsigned cpuCount = std::thread::hardware_concurrency();

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

    int64_t tdLevel = 7;
    std::cout << "Wheel factorization level (minimum of " << MIN_RTD_LEVEL << ", max of " << trialDivisionPrimes.size() << ", or -1 for calibration file): ";
    std::cin >> tdLevel;
    if ((tdLevel > -1) && (tdLevel < MIN_RTD_LEVEL)) {
        tdLevel = MIN_RTD_LEVEL;
    }
    if (tdLevel > 9) {
        tdLevel = 9;
    }
    if (tdLevel < 0) {
        std::ifstream settingsFile ("qimcifa_calibration.ssv");
        std::string header;
        std::getline(settingsFile, header);
        double bestCost = DBL_MAX;
        while (settingsFile.peek() != EOF)
        {
            size_t level;
            BigInteger cardinality;
            double batchTime, cost;
            settingsFile >> level;
            settingsFile >> cardinality;
            settingsFile >> batchTime;
            settingsFile >> cost;

            if ((level >= MIN_RTD_LEVEL) && (cost < bestCost)) {
                tdLevel = level;
                bestCost = cost;
            }
        }
        settingsFile.close();

        std::cout << "Calibrated wheel factorization level: " << tdLevel << std::endl;
        std::cout << "Estimated average time to exit: " << (bestCost / (2 * cpuCount * nodeCount)) << " seconds" << std::endl;
    }

    // Starting clock right after user input is finished
    auto iterClock = std::chrono::high_resolution_clock::now();

    const BigInteger fullMaxBase = sqrt(toFactor);
    if (fullMaxBase * fullMaxBase == toFactor) {
        std::cout << "Number to factor is a perfect square: " << fullMaxBase << " * " << fullMaxBase << " = " << toFactor;
        return 0;
    }

    for (int64_t primeIndex = 0; primeIndex < tdLevel; ++primeIndex) {
        const unsigned currentPrime = trialDivisionPrimes[primeIndex];
        if ((toFactor % currentPrime) == 0) {
            std::cout << "Factors: " << currentPrime << " * " << (toFactor / currentPrime) << " = " << toFactor
                      << std::endl;
            return 0;
        }
        ++primeIndex;
    }

#if IS_SQUARES_CONGRUENCE_CHECK
    const BigInteger offset = (fullMaxBase / BIGGEST_WHEEL) * BIGGEST_WHEEL + 2U;
    const BigInteger fullRange = backward(1U + toFactor - offset);
#else
    const BigInteger offset = 1U;
    const BigInteger fullRange = backward(fullMaxBase);
#endif

    std::vector<boost::dynamic_bitset<uint64_t>> inc_seqs = wheel_gen(std::vector<BigInteger>(trialDivisionPrimes.begin(), trialDivisionPrimes.begin() + tdLevel), toFactor);
    inc_seqs.erase(inc_seqs.begin(), inc_seqs.begin() + MIN_RTD_LEVEL);

    const BigInteger nodeRange = (((fullRange + nodeCount - 1U) / nodeCount) + BIGGEST_WHEEL - 1U) / BIGGEST_WHEEL;
    batchNumber = nodeId * nodeRange;
    batchBound = (nodeId + 1) * nodeRange;
    batchCount = nodeCount * nodeRange;

    const auto workerFn = [toFactor, &inc_seqs, &offset, &iterClock] {
        std::vector<boost::dynamic_bitset<uint64_t>> inc_seqs_clone;
        inc_seqs_clone.reserve(inc_seqs.size());
        for (const auto& b : inc_seqs) {
            inc_seqs_clone.emplace_back(b);
        }
        getSmoothNumbers(toFactor, inc_seqs_clone, offset, iterClock);
    };

    std::vector<std::future<void>> futures;
    futures.reserve(cpuCount);

    for (unsigned cpu = 0U; cpu < cpuCount; ++cpu) {
        futures.push_back(std::async(std::launch::async, workerFn));
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
    BigInteger toFactor;

    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    uint32_t qubitCount = 0;
    BigInteger p = toFactor >> 1U;
    while (p != 0) {
        p >>= 1U;
        ++qubitCount;
    }
    if (!isPowerOfTwo(toFactor)) {
        qubitCount++;
    } else {
        std::cout << "Number to factor is a power of 2: 2 * " << (toFactor >> 1U) << " = " << toFactor;
    }
    std::cout << "Bits to factor: " << (int)qubitCount << std::endl;

    return mainBody(toFactor);
}
