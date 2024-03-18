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

#include <boost/dynamic_bitset.hpp>
#if USE_GMP
#include <boost/multiprecision/gmp.hpp>
#elif USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include "big_integer.hpp"
#endif

#include "dispatchqueue.hpp"

DispatchQueue dispatch(std::thread::hardware_concurrency());

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

inline BigInteger backward(BigInteger ni) {
    ni = (ni + 1) >> 1;
    ni = ((ni + 1) << 1) / 3;
    return ni;
}

inline BigInteger forward(BigInteger p) {
    // Make this NOT a multiple of 2 or 3.
    p += (p >> 1U);
    return (p << 1U) - 1U;
}

const size_t BATCH_SIZE = 1 << 12;
bool isMultipleParallel(const BigInteger& p, const size_t& nextPrimeIndex, const size_t& highestIndex,
    const std::vector<BigInteger>& knownPrimes) {
    const size_t _BATCH_SIZE = BATCH_SIZE;
    const size_t maxLcv = (highestIndex - nextPrimeIndex) / BATCH_SIZE;
    dispatch.resetResult();
    for (size_t i = 0; i < maxLcv; ++i) {
        size_t j = i * BATCH_SIZE + nextPrimeIndex;
        dispatch.dispatch([&knownPrimes, &p, _BATCH_SIZE, j]() {
            for (size_t k = 0; k < _BATCH_SIZE; ++k) {
                if ((p % knownPrimes[j + k]) == 0) {
                    return true;
                }
            }
            return false;
        });
    }

    return dispatch.finish();
}

bool isMultiple(const BigInteger& p, size_t nextIndex, const std::vector<BigInteger>& knownPrimes) {
    const BigInteger sqrtP = sqrt(p);
    const size_t highestIndex = std::distance(knownPrimes.begin(), std::upper_bound(knownPrimes.begin(), knownPrimes.end(), sqrtP));

    const size_t diff = highestIndex - nextIndex;
    if ((highestIndex > nextIndex) && ((diff / std::thread::hardware_concurrency()) > BATCH_SIZE)) {
        if (isMultipleParallel(p, nextIndex, highestIndex, knownPrimes)) {
            return true;
        }
    }
    nextIndex = diff % BATCH_SIZE;

    for (size_t i = nextIndex; i < highestIndex; ++i) {
        if ((p % knownPrimes[i]) == 0) {
            return true;
        }
    }
    return false;
}

bool isMultiple(const BigInteger& p, const std::vector<BigInteger>& knownPrimes) {
    for (const BigInteger& prime : knownPrimes) {
        if ((p % prime) == 0) {
            return true;
        }
    }
    return false;
}

boost::dynamic_bitset<uint32_t> wheel_inc(std::vector<BigInteger> primes) {
    BigInteger radius = 1U;
    for (const BigInteger& i : primes) {
        radius *= i;
    }
    const BigInteger prime = primes.back();
    primes.pop_back();
    std::vector<bool> o;
    size_t count = 0;
    for (size_t i = 1U; i < radius; ++i) {
        if (!isMultiple(i, primes)) {
            o.push_back((i % prime) == 0);
            ++count;
        }
    }

    boost::dynamic_bitset<uint32_t> output(count);
    for (size_t i = 0U; i < count; ++i) {
        output[i] = o[i];
    }
    output >>= 1U;

    return output;
}

std::vector<boost::dynamic_bitset<uint32_t>> wheel_gen(const std::vector<BigInteger>& primes) {
    std::vector<boost::dynamic_bitset<uint32_t>> output;
    std::vector<BigInteger> wheelPrimes;
    for (const BigInteger p : primes) {
        wheelPrimes.push_back(p);
        if (wheelPrimes.back() > 3) {
            output.push_back(wheel_inc(wheelPrimes));
        }
    }
    return output;
}


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
    std::vector<boost::dynamic_bitset<uint32_t>> inc_seqs;
    BigInteger o = 1U;
    size_t wheel_limit = 13U;
    while (true) {
        ++o;
        bool is_wheel_multiple = false;
        for (size_t i = 0; i < inc_seqs.size(); ++i) {
            boost::dynamic_bitset<uint32_t>& wheel = inc_seqs[i];
            is_wheel_multiple = wheel[0U];
            wheel >>= 1U;
            if (is_wheel_multiple) {
                wheel[wheel.size() - 1U] = true;
                break;
            }
        }

        if (is_wheel_multiple) {
            continue;
        }

        BigInteger p = forward(o);
        if (p > n) {
            break;
        }
        if (isMultiple(p, wheelPrimes.size(), knownPrimes)) {
            // Skip
            continue;
        }

        knownPrimes.push_back(p);
        if (p <= wheel_limit) {
            wheelPrimes.push_back(p);
            inc_seqs.push_back(wheel_inc(knownPrimes));
            boost::dynamic_bitset<uint32_t>& wheel = inc_seqs.back();
            wheel >>= 1U;
            wheel[wheel.size() - 1U] = true;
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
