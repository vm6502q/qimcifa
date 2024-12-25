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

template <typename BigInteger>
double mainBody(const BigInteger& toFactor, const uint64_t& tdLevel, const std::vector<unsigned>& trialDivisionPrimes)
{
    // When we factor this number, we split it into two factors (which themselves may be composite).
    // Those two numbers are either equal to the square root, or in a pair where one is higher and one lower than the square root.

    std::vector<boost::dynamic_bitset<uint64_t>> inc_seqs = wheel_gen(std::vector<BigInteger>(trialDivisionPrimes.begin(), trialDivisionPrimes.begin() + tdLevel), toFactor);
    inc_seqs.erase(inc_seqs.begin(), inc_seqs.begin() + 2U);

    const BigInteger fullMaxBase = backward(sqrt<BigInteger>(toFactor));
    const BigInteger offset = (fullMaxBase / BIGGEST_WHEEL) * BIGGEST_WHEEL + 1U;

    std::vector<BigInteger> smoothNumbers;
    auto iterClock = std::chrono::high_resolution_clock::now();
    batchNumber = 0U;
    batchBound = 1U;
    getSmoothNumbers(toFactor, inc_seqs, offset, iterClock);

    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - iterClock).count() * 1e-10;
}
} // namespace Qimcifa

using namespace Qimcifa;

double mainCase(BigIntegerInput toFactor, int tdLevel)
{
    uint32_t qubitCount = 0;
    BigIntegerInput p = toFactor;
    while (p) {
        p >>= 1U;
        ++qubitCount;
    }
    if (!isPowerOfTwo(toFactor)) {
        qubitCount++;
    }

    // First 9 primes
    std::vector<unsigned> trialDivisionPrimes = { 2, 3, 5, 7, 11, 13, 17, 19, 23 };

    if (qubitCount < 64) {
        typedef uint64_t BigInteger;
        return mainBody<BigInteger>((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
#if USE_GMP
    } else {
        return mainBody<BigIntegerInput>(toFactor, tdLevel, trialDivisionPrimes);
    }
#else
    } else if (qubitCount < 128) {
        typedef boost::multiprecision::uint128_t BigInteger;
        return mainBody<BigInteger>((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    } else if (qubitCount < 192) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<192, 192,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            BigInteger;
        return mainBody<BigInteger>((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    } else if (qubitCount < 256) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<256, 256,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            BigInteger;
        return mainBody<BigInteger>((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    } else if (qubitCount < 512) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<512, 512,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            BigInteger;
        return mainBody<BigInteger>((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    } else if (qubitCount < 1024) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<1024, 1024,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            BigInteger;
        return mainBody<BigInteger>((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    } else if (qubitCount < 2048) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<2048, 2048,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            BigInteger;
        return mainBody((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    } else if (qubitCount < 4096) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<4096, 4096,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            BigInteger;
        return mainBody((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    } else if (qubitCount < 8192) {
        typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<8192, 8192,
            boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
            BigInteger;
        return mainBody((BigInteger)toFactor, tdLevel, trialDivisionPrimes);
    }

    if (qubitCount >= 8192) {
        throw std::runtime_error("Number to factor exceeds 8192 templated max bit width!");
    }
#endif

    return -999.0;
}

int main() {
    BigIntegerInput toFactor;

    std::cout << "The qimcifa_tuner number need not be the exact number that will be factored with qimcifa, but closer is better." << std::endl;
    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    uint32_t qubitCount = 0;
    BigIntegerInput p = toFactor >> 1U;
    while (p) {
        p >>= 1U;
        ++qubitCount;
    }
    if (!isPowerOfTwo(toFactor)) {
        qubitCount++;
    }
    std::cout << "Bits to factor: " << (int)qubitCount << std::endl;

    size_t threadCount = 1;
    std::cout << "Total thread count (across all nodes): ";
    std::cin >> threadCount;

    // First 9 primes
    std::vector<unsigned> trialDivisionPrimes = { 2, 3, 5, 7, 11, 13, 17, 19, 23 };

    const BigIntegerInput range = backward(sqrt(toFactor));
    std::ofstream oSettingsFile ("qimcifa_calibration.ssv");
    oSettingsFile << "level, cardinality, batch time (ns), cost (s)" << std::endl;
    int wheelSize = 30;
    for (size_t i = MIN_RTD_LEVEL; i < 10U; ++i) {
        // Test
        const double time = mainCase(toFactor, i);
        wheelSize *= trialDivisionPrimes[i];
#if BIG_INTEGER_BITS > 64 && !USE_BOOST && !USE_GMP
        oSettingsFile << i << " " << range << " " << time << " " << (bi_to_double(range) * (time / wheelSize)) << std::endl;
#else
        oSettingsFile << i << " " << range << " " << time << " " << (range.convert_to<double>() * (time / wheelSize)) << std::endl;
#endif
    }
    oSettingsFile.close();

    std::ifstream iSettingsFile ("qimcifa_calibration.ssv");
    std::string header;
    std::getline(iSettingsFile, header);
    size_t bestLevel = -1;
    double bestCost = DBL_MAX;
    while (iSettingsFile.peek() != EOF)
    {
        size_t level;
        BigIntegerInput cardinality;
        double batchTime, cost;
        iSettingsFile >> level;
        iSettingsFile >> cardinality;
        iSettingsFile >> batchTime;
        iSettingsFile >> cost;

        if ((level >= MIN_RTD_LEVEL) && (cost < bestCost)) {
            bestLevel = level;
            bestCost = cost;
        }
    }
    iSettingsFile.close();

    std::cout << "Calibrated reverse trial division level: " << bestLevel << std::endl;
    std::cout << "Estimated average time to exit: " << (bestCost / (2 * threadCount)) << " seconds" << std::endl;

    return 0;
}
