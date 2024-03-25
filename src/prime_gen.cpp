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

#include "prime_generator.hpp"

#include <pybind11/pybind11.h>

namespace qimcifa {
std::vector<BigInteger> TrialDivision(const BigInteger& n)
{
    // First 3 primes
    std::vector<BigInteger> knownPrimes = { 2, 3 };

    if (n < 2) {
        return std::vector<BigInteger>();
    }

    if (n < (knownPrimes.back() + 2)) {
        const auto highestPrimeIt = std::upper_bound(knownPrimes.begin(), knownPrimes.end(), n);
        return std::vector<BigInteger>(knownPrimes.begin(), highestPrimeIt);
    }

    std::vector<size_t> wheelPrimes= { 2, 3 };

    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2 and 3, this reduces complexity by 2/3.
    // const BigInteger cardinality = (~((~n) | 1)) / 3;

    // Get the remaining prime numbers.
    std::vector<boost::dynamic_bitset<size_t>> inc_seqs;
    const size_t wheel_limit = 17U;
    for (BigInteger o = 1U; forward(o) < n;) {
        o += GetWheelIncrement(inc_seqs);

        const BigInteger p = forward(o);
        if (p > n) {
            break;
        }
        if (isMultiple(p, wheelPrimes.size(), knownPrimes)) {
            // Skip
            continue;
        }

        knownPrimes.push_back(p);
        if (p <= wheel_limit) {
            wheelPrimes.push_back((size_t)p);
            inc_seqs.push_back(wheel_inc(knownPrimes, n));
            boost::dynamic_bitset<size_t>& wheel = inc_seqs.back();
            wheel >>= 1U;
            wheel[wheel.size() - 1U] = true;
        }
    }

    return knownPrimes;
}

std::vector<BigInteger> SieveOfEratosthenes(const BigInteger& n)
{
    std::vector<BigInteger> knownPrimes = { 2, 3, 5 };
    if (n < 2) {
        return std::vector<BigInteger>();
    }

    if (n < (knownPrimes.back() + 2)) {
        const auto highestPrimeIt = std::upper_bound(knownPrimes.begin(), knownPrimes.end(), n);
        return std::vector<BigInteger>(knownPrimes.begin(), highestPrimeIt);
    }

    BigInteger threadLimit = 26U;

    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2, 3, and 5 this reduces complexity to 4/15.
    const size_t cardinality = (size_t)backward5(n);

    // Create a boolean array "prime[0..cardinality]"
    // and initialize all entries it as true. Rather,
    // reverse the true/false meaning, so we can use
    // default initialization. A value in notPrime[i]
    // will finally be false only if i is a prime.
    std::vector<bool> notPrime(cardinality + 1);

    // Get the remaining prime numbers.
    dispatch.resetResult();
    std::vector<boost::dynamic_bitset<size_t>> inc_seqs = wheel_gen(knownPrimes, n);
    size_t o = 1U;
    for (;;) {
        o += GetWheelIncrement(inc_seqs);

        const BigInteger p = forward(o);
        if ((p * p) > n) {
            break;
        }

        if (threadLimit < p) {
            dispatch.finish();
            threadLimit *= threadLimit;
        }


        const size_t q = (size_t)backward5(p);
        if (notPrime[q] == true) {
            continue;
        }

        knownPrimes.push_back(p);

        dispatch.dispatch([&n, p, &notPrime]() {
            // We are skipping multiples of 2, 3, and 5
            // for space complexity, for 4/15 the bits.
            // More are skipped by the wheel for time.
            const BigInteger p2 = p << 1U;
            const BigInteger p4 = p << 2U;
            BigInteger i = p * p;

            // "p" already definitely not a multiple of 3.
            // Its remainder when divided by 3 can be 1 or 2.
            // If it is 2, we can do a "half iteration" of the
            // loop that would handle remainder of 1, and then
            // we can proceed with the 1 remainder loop.
            // This saves 2/3 of updates (or modulo).
            if ((p % 3U) == 2U) {
                notPrime[(size_t)backward5(i)] = true;
                i += p2;
                if (i > n) {
                    return false;
                }
            }

            std::vector<bool> wheel30;
            wheel30.reserve(30);
            for (int j = 0; j < 15; ++j) {
                wheel30.push_back(i % 5);
                if (wheel30.back()) {
                    notPrime[(size_t)backward5(i)] = true;
                }
                i += p4;
                if (i > n) {
                    return false;
                }

                wheel30.push_back(i % 5);
                if (wheel30.back()) {
                    notPrime[(size_t)backward5(i)] = true;
                }
                i += p2;
                if (i > n) {
                    return false;
                }
            }

            for (;;) {
                for (int j = 0; j < 30; j+=2) {
                    if (wheel30[j]) {
                        notPrime[(size_t)backward5(i)] = true;
                    }
                    i += p4;
                    if (i > n) {
                        return false;
                    }

                    if (wheel30[j + 1]) {
                        notPrime[(size_t)backward5(i)] = true;
                    }
                    i += p2;
                    if (i > n) {
                        return false;
                    }
                }
            }

            return false;
        });
    }
    dispatch.finish();

    inc_seqs = wheel_gen(knownPrimes, n);
    o = 1U;
    for (;;) {
        o += GetWheelIncrement(inc_seqs);

        const BigInteger p = forward(o);
        if (p > n) {
            break;
        }

        if (notPrime[(size_t)backward5(p)] == true) {
            continue;
        }

        knownPrimes.push_back(p);
    }

    return knownPrimes;
}

std::vector<BigInteger> SegmentedSieveOfEratosthenes(const BigInteger& n)
{
    // TODO: This should scale to the system.
    // It's 16 GB in bytes.
    const BigInteger limit = BigInteger(137438953472ULL);

    // `backward(n)` counts assuming that multiples
    // of 2 and 3 have been removed.
    if (backward(n) <= limit) {
        return SieveOfEratosthenes(n);
    }

    // Process segments of length `limit` at a time.
    BigInteger low = limit | 1U;
    if ((low % 3U) == 0U) {
        low -= 2U;
    }
    BigInteger high = (limit << 1U) | 1U;
    if ((high % 3U) == 0U) {
        high -= 2U;
    }

    // Compute all primes smaller than or equal to limit using simple sieve
    std::vector<BigInteger> knownPrimes = SieveOfEratosthenes(limit);
    dispatch.resetResult();

    // Process one segment at a time until we pass n
    while (low < n) {
        if (high >= n) {
            high = n | 1U;
            if ((high % 3U) == 0U) {
                high -= 2U;
            }
        }

        // Cardinality with multiples of 2 and 3 removed is 1/3 of total.
        const BigInteger bLow = backward(low);
        std::vector<bool> notPrime((size_t)(backward(high) - bLow) + 1U);

        // Use the found primes by simpleSieve() to find
        // primes in current range
        for (size_t k = 2U; k < knownPrimes.size(); k++) {
            // Find the minimum number in [low..high] that is
            // a multiple of prime[i] (divisible by prime[i])
            // For example, if low is 31 and prime[i] is 3,
            // we start with 33.
            const BigInteger& p = knownPrimes[k];

            dispatch.dispatch([&bLow, &high, &low, p, &notPrime]() {
                // We are skipping multiples of 2, 3, and 5
                // for space complexity, for 4/15 the bits.
                // More are skipped by the wheel for time.
                const BigInteger p2 = p << 1U;
                const BigInteger p4 = p << 2U;
                BigInteger i = ((low + p - 1U) / p) * p;
                while (((i & 1U) == 0U) || ((i % 3U) == 0U)) {
                    i += p;
                }

                // "p" already definitely not a multiple of 3.
                // Its remainder when divided by 3 can be 1 or 2.
                // If it is 2, we can do a "half iteration" of the
                // loop that would handle remainder of 1, and then
                // we can proceed with the 1 remainder loop.
                // This saves 2/3 of updates (or modulo).
                if ((i % 3U) == 2U) {
                    const size_t q = (size_t)(backward(i) - bLow);
                    if (q >= notPrime.size()) {
                        return false;
                    }
                    notPrime[q] = true;
                    i += p2;
                }

                for (;;) {
                    size_t q = (size_t)(backward(i) - bLow);
                    if (q >= notPrime.size()) {
                        return false;
                    }
                    notPrime[q] = true;
                    i += p4;

                    q = (size_t)(backward(i) - bLow);
                    if (q >= notPrime.size()) {
                        return false;
                    }
                    notPrime[q] = true;
                    i += p2;
                }

                return false;
            });
        }
        dispatch.finish();

        // Numbers which are not marked as false are prime
        for (size_t i = 0; i < notPrime.size(); ++i) {
            if (notPrime[i] == false) {
                knownPrimes.push_back(forward(i + bLow));
            }
        }

        // Update low and high for next segment
        low = (low + limit) | 1U;
        if ((low % 3U) == 0U) {
            low -= 2U;
        }
        high = (high + limit) | 1U;
        if ((high % 3U) == 0U) {
            high -= 2U;
        }
    }

    return knownPrimes;
}
} // namespace qimcifa

using namespace qimcifa;

PYBIND11_MODULE(prime_gen, m) {
    m.doc() = "pybind11 plugin to generate prime numbers";
    m.def("prime_gen", &SieveOfEratosthenes, "A function that returns all primes up to the value of its argument");
}
