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

// Make this a multiple of 2, 3, 5, 7, 11, 13, 17, 19, and 23.
constexpr int BIGGEST_WHEEL = 223092870;
// Make this a multiple of 2, 3, 5, 7, or 11.
constexpr int SMALLEST_WHEEL = 2310;
constexpr int MIN_RTD_LEVEL = 5;

#if USE_GMP
typedef boost::multiprecision::mpz_int BigInteger;
#elif USE_BOOST
typedef boost::multiprecision::cpp_int BigInteger;
#endif

BigInteger batchNumber;
BigInteger batchBound;
BigInteger batchCount;
std::mutex batchMutex;

inline void finish() {
    std::lock_guard<std::mutex> lock(batchMutex);
    batchNumber = batchBound;
}

inline BigInteger getNextBatch() {
    std::lock_guard<std::mutex> lock(batchMutex);

#if IS_SQUARES_CONGRUENCE_CHECK
    BigInteger result = batchNumber;

    if (batchNumber < batchBound) {
        ++batchNumber;
    }

    return result;
#else
    BigInteger result = batchCount - (batchNumber + 1U);

    if (batchNumber >= batchBound) {
        return batchBound;
    }

    ++batchNumber;

    return result;
#endif
}

// See https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
inline BigInteger ipow(BigInteger base, unsigned exp)
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

inline uint64_t log2(BigInteger n) {
#if USE_GMP || USE_BOOST
    uint64_t pow = 0U;
    while (n >>= 1U) {
        ++pow;
    }
    return pow;
#else
    return bi_log2(n);
#endif
}

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

inline bool isPowerOfTwo(const BigInteger& x)
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

inline BigInteger gcd(BigInteger n1, BigInteger n2)
{
    while (n2 != 0) {
        const BigInteger t = n1;
        n1 = n2;
        n2 = t % n2;
    }

    return n1;
}

inline void printSuccess(const BigInteger& f1, const BigInteger& f2, const BigInteger& toFactor, const std::string& message,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock)
{
    finish();
    std::cout << message << f1 << " * " << f2 << " = " << toFactor << std::endl;
    auto tClock =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
    // Report in seconds
    std::cout << "(Time elapsed: " << (tClock.count() / 1000000.0) << " seconds)" << std::endl;
    std::cout << "(Waiting to join other threads...)" << std::endl;
}

constexpr unsigned short wheel11[480U] = {
    1U, 13U, 17U, 19U, 23U, 29U, 31U, 37U, 41U, 43U, 47U, 53U, 59U, 61U, 67U, 71U, 73U, 79U, 83U, 89U, 97U, 101U,
    103U, 107U, 109U, 113U, 127U, 131U, 137U, 139U, 149U, 151U, 157U, 163U, 167U, 169U, 173U, 179U, 181U, 191U,
    193U, 197U, 199U, 211U, 221U, 223U, 227U, 229U, 233U, 239U, 241U, 247U, 251U, 257U, 263U, 269U, 271U, 277U,
    281U, 283U, 289U, 293U, 299U, 307U, 311U, 313U, 317U, 323U, 331U, 337U, 347U, 349U, 353U, 359U, 361U, 367U,
    373U, 377U, 379U, 383U, 389U, 391U, 397U, 401U, 403U, 409U, 419U, 421U, 431U, 433U, 437U, 439U, 443U, 449U,
    457U, 461U, 463U, 467U, 479U, 481U, 487U, 491U, 493U, 499U, 503U, 509U, 521U, 523U, 527U, 529U, 533U, 541U,
    547U, 551U, 557U, 559U, 563U, 569U, 571U, 577U, 587U, 589U, 593U, 599U, 601U, 607U, 611U, 613U, 617U, 619U,
    629U, 631U, 641U, 643U, 647U, 653U, 659U, 661U, 667U, 673U, 677U, 683U, 689U, 691U, 697U, 701U, 703U, 709U,
    713U, 719U, 727U, 731U, 733U, 739U, 743U, 751U, 757U, 761U, 767U, 769U, 773U, 779U, 787U, 793U, 797U, 799U,
    809U, 811U, 817U, 821U, 823U, 827U, 829U, 839U, 841U, 851U, 853U, 857U, 859U, 863U, 871U, 877U, 881U, 883U,
    887U, 893U, 899U, 901U, 907U, 911U, 919U, 923U, 929U, 937U, 941U, 943U, 947U, 949U, 953U, 961U, 967U, 971U,
    977U, 983U, 989U, 991U, 997U, 1003U, 1007U, 1009U, 1013U, 1019U, 1021U, 1027U, 1031U, 1033U, 1037U, 1039U,
    1049U, 1051U, 1061U, 1063U, 1069U, 1073U, 1079U, 1081U, 1087U, 1091U, 1093U, 1097U, 1103U, 1109U, 1117U, 1121U,
    1123U, 1129U, 1139U, 1147U, 1151U, 1153U, 1157U, 1159U, 1163U, 1171U, 1181U, 1187U, 1189U, 1193U, 1201U, 1207U,
    1213U, 1217U, 1219U, 1223U, 1229U, 1231U, 1237U, 1241U, 1247U, 1249U, 1259U, 1261U, 1271U, 1273U, 1277U, 1279U,
    1283U, 1289U, 1291U, 1297U, 1301U, 1303U, 1307U, 1313U, 1319U, 1321U, 1327U, 1333U, 1339U, 1343U, 1349U, 1357U,
    1361U, 1363U, 1367U, 1369U, 1373U, 1381U, 1387U, 1391U, 1399U, 1403U, 1409U, 1411U, 1417U, 1423U, 1427U, 1429U,
    1433U, 1439U, 1447U, 1451U, 1453U, 1457U, 1459U, 1469U, 1471U, 1481U, 1483U, 1487U, 1489U, 1493U, 1499U, 1501U,
    1511U, 1513U, 1517U, 1523U, 1531U, 1537U, 1541U, 1543U, 1549U, 1553U, 1559U, 1567U, 1571U, 1577U, 1579U, 1583U,
    1591U, 1597U, 1601U, 1607U, 1609U, 1613U, 1619U, 1621U, 1627U, 1633U, 1637U, 1643U, 1649U, 1651U, 1657U, 1663U,
    1667U, 1669U, 1679U, 1681U, 1691U, 1693U, 1697U, 1699U, 1703U, 1709U, 1711U, 1717U, 1721U, 1723U, 1733U, 1739U,
    1741U, 1747U, 1751U, 1753U, 1759U, 1763U, 1769U, 1777U, 1781U, 1783U, 1787U, 1789U, 1801U, 1807U, 1811U, 1817U,
    1819U, 1823U, 1829U, 1831U, 1843U, 1847U, 1849U, 1853U, 1861U, 1867U, 1871U, 1873U, 1877U, 1879U, 1889U, 1891U,
    1901U, 1907U, 1909U, 1913U, 1919U, 1921U, 1927U, 1931U, 1933U, 1937U, 1943U, 1949U, 1951U, 1957U, 1961U, 1963U,
    1973U, 1979U, 1987U, 1993U, 1997U, 1999U, 2003U, 2011U, 2017U, 2021U, 2027U, 2029U, 2033U, 2039U, 2041U, 2047U,
    2053U, 2059U, 2063U, 2069U, 2071U, 2077U, 2081U, 2083U, 2087U, 2089U, 2099U, 2111U, 2113U, 2117U, 2119U, 2129U,
    2131U, 2137U, 2141U, 2143U, 2147U, 2153U, 2159U, 2161U, 2171U, 2173U, 2179U, 2183U, 2197U, 2201U, 2203U, 2207U,
    2209U, 2213U, 2221U, 2227U, 2231U, 2237U, 2239U, 2243U, 2249U, 2251U, 2257U, 2263U, 2267U, 2269U, 2273U, 2279U,
    2281U, 2287U, 2291U, 2293U, 2297U, 2309U
};

inline BigInteger forward(const BigInteger& p) {
    // Make this NOT a multiple of 2, 3, 5, 7, or 11.
    return wheel11[(size_t)(p % 480U)] + (p / 480U) * 2310U;
}

inline BigInteger backward(const BigInteger& n) {
    return std::distance(wheel11, std::lower_bound(wheel11, wheel11 + 480U, size_t(n % 2310U))) + 480U * (size_t)(n / 2310U) + 1U;
}

inline bool isMultiple(const BigInteger& p, const std::vector<BigInteger>& knownPrimes) {
    for (const BigInteger& prime : knownPrimes) {
        if ((p % prime) == 0) {
            return true;
        }
    }
    return false;
}

inline boost::dynamic_bitset<size_t> wheel_inc(std::vector<BigInteger> primes, BigInteger limit) {
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

inline std::vector<boost::dynamic_bitset<size_t>> wheel_gen(const std::vector<BigInteger>& primes, BigInteger limit) {
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
        ++wheelIncrement;
    } while (is_wheel_multiple);

    return wheelIncrement;
}

#if IS_SQUARES_CONGRUENCE_CHECK
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

inline bool getSmoothNumbersIteration(const BigInteger& toFactor, const BigInteger& base,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock) {
#if IS_SQUARES_CONGRUENCE_CHECK
    return checkCongruenceOfSquares<BigInteger>(toFactor, base, iterClock);
#elif IS_RSA_SEMIPRIME
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

    return false;
}

inline bool getSmoothNumbers(const BigInteger& toFactor, std::vector<boost::dynamic_bitset<uint64_t>>& inc_seqs, const BigInteger& offset,
    const std::chrono::time_point<std::chrono::high_resolution_clock>& iterClock)
{
    for (BigInteger batchNum = (BigInteger)getNextBatch(); batchNum < batchBound; batchNum = (BigInteger)getNextBatch()) {
        const BigInteger batchStart = batchNum * BIGGEST_WHEEL + offset;
        const BigInteger batchEnd = (batchNum + 1U) * BIGGEST_WHEEL + offset;
        for (BigInteger p = batchStart; p < batchEnd;) {
            p += GetWheelIncrement(inc_seqs);
            if (getSmoothNumbersIteration<BigInteger>(toFactor, forward(p), iterClock)) {
                return true;
            }
        }
    }

    return false;
}
} // namespace Qimcifa
