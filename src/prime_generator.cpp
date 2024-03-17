// Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
// C++ program to print all primes smaller than or equal to
// n using Sieve of Eratosthenes

// Improved by Dan Strano of Unitary Fund, 2024.
// We can think of trial division as exact inverse of
// Sieve of Eratosthenes, with log space and log time.
// The modular division part is a costly atomic operation.
// It need only be carried out up the square root of the
// number under trial. Multiples of 2, 3, 5, 7, and 11 can
// be entirely skipped in loop enumeration.

#include "config.h"

#include <future>
#include <iostream>
#include <vector>

#if USE_GMP
#include <boost/multiprecision/gmp.hpp>
#elif USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include "big_integer.hpp"
#endif

#if BIG_INT_BITS < 33
typedef uint32_t BigInteger;
#elif BIG_INT_BITS < 65
typedef uint64_t BigInteger;
#else
#if USE_GMP
typedef boost::multiprecision::mpz_int BigInteger;
#elif USE_BOOST
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<BIG_INT_BITS, BIG_INT_BITS,
    boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>
    BigInteger;
#else
typedef BigInteger BigInteger;
#endif
#endif

#if 0
inline size_t log2(BigInteger n) {
    size_t log2n = 0;
    BigInteger _n = n;
    while (_n >>= 1) {
        ++log2n;
    }

    return log2n;
}

inline BigInteger gcd(BigInteger n1, BigInteger n2)
{
#if USE_GMP || USE_BOOST
    while (n2) {
#else
    if (bi_compare_0(n2) != 0) {
#endif
        const BigInteger t = n1;
        n1 = n2;
        n2 = t % n2;
    }

    return n1;
}
#endif

inline BigInteger sqrt(const BigInteger& toTest)
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

BigInteger backward(BigInteger ni) {
    ni = (ni + 1) >> 1;
    ni = ((ni + 1) << 1) / 3;
    return ni;
}

BigInteger forward(BigInteger p) {
    // Make this NOT a multiple of 2 or 3.
    p += (p >> 1U);
    return (p << 1U) - 1U;
}

const size_t BATCH_SIZE = 1 << 12;

bool isMultipleParallel(const BigInteger& p, const size_t& nextPrimeIndex, const size_t& highestPrimeIndex,
    const std::vector<BigInteger>& knownPrimes) {
    const size_t _BATCH_SIZE = BATCH_SIZE;
    const unsigned cpuCount = std::thread::hardware_concurrency();
    const size_t batchSize = cpuCount * _BATCH_SIZE;
    const size_t maxLcv = (highestPrimeIndex - nextPrimeIndex) - batchSize;
    std::vector<std::future<bool>> futures(cpuCount);
    for (size_t i = 0;;) {
        for (unsigned cpu = 0; cpu < cpuCount; ++cpu) {
            i += BATCH_SIZE;
            futures[cpu] = std::async(std::launch::async,
                [&knownPrimes, _BATCH_SIZE, cpu, i, nextPrimeIndex](const BigInteger& p) {
                    for (size_t j = 0; j < _BATCH_SIZE; ++j) {
                        if ((p % knownPrimes[nextPrimeIndex + i + (cpu * _BATCH_SIZE) + j]) != 1) {
                            return true;
                        }
                    }
                    return false;
                }, p);
        }
        for (unsigned cpu = 0; cpu < cpuCount; ++cpu) {
            if (futures[cpu].get()) {
                return true;
            }
            if (i <= maxLcv) {
                i += BATCH_SIZE;
                futures[cpu] = std::async(std::launch::async,
                    [&knownPrimes, _BATCH_SIZE, cpu, i, nextPrimeIndex](const BigInteger& p) {
                        for (size_t j = 0; j < _BATCH_SIZE; ++j) {
                            if ((p % knownPrimes[nextPrimeIndex + i + (cpu * _BATCH_SIZE) + j]) != 1) {
                                return true;
                            }
                        }
                        return false;
                    }, p);
            }
        }
    }

    return false;
}

bool isMultiple(const BigInteger& p, size_t nextIndex, const std::vector<BigInteger>& knownPrimes) {
    const BigInteger sqrtP = sqrt(p);
    if ((sqrtP * sqrtP) == p) {
        return true;
    }
    size_t highestIndex = 0U;
    for (size_t m = (knownPrimes.size() + 1) >> 1U; m > 1; m = (m + 1) >> 1) {
        if (knownPrimes[highestIndex + m] <= sqrtP) {
            highestIndex += m;
        }
    }

    /*const size_t diff = highestIndex - nextIndex;
    const unsigned cpuCount = std::thread::hardware_concurrency();
    if ((diff / cpuCount) > BATCH_SIZE) {
        if (isMultipleParallel(p, nextIndex, highestIndex, knownPrimes)) {
            return true;
        }
    }
    nextIndex = diff % (BATCH_SIZE * cpuCount);*/

    for (size_t i = nextIndex; i <= highestIndex; ++i) {
        if ((p % knownPrimes[i]) == 0) {
            return true;
        }
    }
    return false;
}

std::list<bool> wheel_inc(std::vector<BigInteger> primes) {
    BigInteger prime = primes.back();
    primes.pop_back();
    BigInteger radius = 1U;
    for (BigInteger i : primes) {
        radius *= i;
    }
    std::list<bool> output;
    for (BigInteger i = 1; i < radius; ++i) {
        if (!isMultiple(i, 2, primes)) {
            output.push_back((i % prime) == 0);
        }
    }

    return output;
}

#if 0
def wheel_gen(primes):
    output = []
    for i in range(3, len(primes) + 1):
        output.append(wheel_inc(primes[:i]))
        output[-1] = output[-1][1:] + output[-1][:1]
    return output
#endif

std::vector<BigInteger> TrialDivision(const BigInteger& n)
{
    // First 2 primes
    std::vector<BigInteger> knownPrimes= { 2, 3 };

    if (n < 2) {
        return std::vector<BigInteger>();
    }

    if (n < (knownPrimes.back() + 2)) {
        const BigInteger sqrtN = sqrt(n);
        size_t highestIndex = 0U;
        for (size_t m = (knownPrimes.size() + 1) >> 1U; m > 1; m = (m + 1) >> 1) {
            if (knownPrimes[highestIndex + m] <= sqrtN) {
                highestIndex += m;
            }
        }
        return std::vector<BigInteger>(knownPrimes.begin(), knownPrimes.begin() + highestIndex);
    }

    std::vector<BigInteger> wheelPrimes= { 2, 3 };

    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2 and 3, this reduces complexity by 2/3.
    // const BigInteger cardinality = (~((~n) | 1)) / 3;

    // Get the remaining prime numbers.
    std::vector<std::list<bool>> inc_seqs;
    BigInteger o = 1U;
    size_t wheel_limit = 11;
    while (true) {
        ++o;
        bool is_wheel_multiple = false;
        for (size_t i = 0; i < inc_seqs.size(); ++i) {
            std::list<bool>& wheel = inc_seqs[i];
            is_wheel_multiple = wheel.front();
            std::rotate(wheel.begin(), next(wheel.begin()), wheel.end());
            if (is_wheel_multiple) {
                break;
            }
        }

        BigInteger p = forward(o);

        if (is_wheel_multiple) {
            continue;
        }

        if (p > n) {
            break;
        }
        if (isMultiple(p, wheelPrimes.size() - 1, knownPrimes)) {
            // Skip
            continue;
        }

        knownPrimes.push_back(p);
        if (p <= wheel_limit) {
            wheelPrimes.push_back(p);
            inc_seqs.push_back(wheel_inc(knownPrimes));
            std::list<bool>& wheel = inc_seqs.back();
            std::rotate(wheel.begin(), next(next(wheel.begin())), wheel.end());
        }
    }

    return knownPrimes;
}
 
// Driver Code
int main()
{
    BigInteger n; // = 1000000;

    std::cout << "Primes up to number: ";
    std::cin >> n;

    std:: cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;

    const std::vector<BigInteger> primes = TrialDivision(n);

    for (BigInteger p : primes) {
        std::cout << p << " ";
    }
    std::cout << std::endl;

    return 0;
}
