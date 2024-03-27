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

#include <cmath>

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

    knownPrimes.reserve(std::expint(log(n)) - std::expint(log(2)));

    // We are excluding multiples of the first few
    // small primes from outset. For multiples of
    // 2, 3, and 5 this reduces complexity to 4/15.
    const size_t cardinality = (size_t)backward5(n);

    // Create a boolean array "prime[0..cardinality]"
    // and initialize all entries it as true. Rather,
    // reverse the true/false meaning, so we can use
    // default initialization. A value in notPrime[i]
    // will finally be false only if i is a prime.
    size_t bitsPerWord = 64U;
    const size_t arrayWidth = (cardinality + bitsPerWord) / bitsPerWord;
    std::unique_ptr<uint64_t> uNotPrime(new uint64_t[arrayWidth]());
    std::unique_ptr<std::mutex> uNotPrimeMutex(new std::mutex[arrayWidth]());
    uint64_t* notPrime = uNotPrime.get();
    std::mutex* notPrimeMutex = uNotPrimeMutex.get();

    // We dispatch multiple marking asynchronously.
    // If we've already marked all primes up to x,
    // we're free to continue to up to x * x,
    // then we synchronize.
    BigInteger threadBoundary = 36U;

    // Get the remaining prime numbers.
    uint32_t wheel5 = (1U << 7U) | 1U;
    size_t o = 1U;
    for (;;) {
        o += GetWheel5Increment(wheel5);

        const BigInteger p = forward(o);
        if ((p * p) > n) {
            break;
        }

        if (threadBoundary < p) {
            dispatch.finish();
            threadBoundary *= threadBoundary;
        }

        const size_t q = (size_t)backward5(p);
        if (notPrime[q / bitsPerWord] == 1ULL << (q % bitsPerWord)) {
            continue;
        }

        knownPrimes.push_back(p);

        dispatch.dispatch([&n, p, &bitsPerWord, &notPrime, &notPrimeMutex]() {
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
                const size_t q = (size_t)backward5(i);
                const size_t w = q / bitsPerWord;
                const uint64_t b = 1ULL << (q % bitsPerWord);
                if (true) {
                    std::lock_guard<std::mutex> lock(notPrimeMutex[w]);
                    notPrime[w] |= b;
                }
                i += p2;
                if (i > n) {
                    return false;
                }
            }

            for (;;) {
                if (i % 5U) {
                    const size_t q = (size_t)backward5(i);
                    const size_t w = q / bitsPerWord;
                    const uint64_t b = 1ULL << (q % bitsPerWord);
                    std::lock_guard<std::mutex> lock(notPrimeMutex[w]);
                    notPrime[w] |= b;
                }
                i += p4;
                if (i > n) {
                    return false;
                }

                if (i % 5U) {
                    const size_t q = (size_t)backward5(i);
                    const size_t w = q / bitsPerWord;
                    const uint64_t b = 1ULL << (q % bitsPerWord);
                    std::lock_guard<std::mutex> lock(notPrimeMutex[w]);
                    notPrime[w] |= b;
                }
                i += p2;
                if (i > n) {
                    return false;
                }
            }

            return false;
        });
    }

    dispatch.finish();
    notPrimeMutex = NULL;
    uNotPrimeMutex.reset();

    for (;;) {
        const BigInteger p = forward(o);
        if (p > n) {
            break;
        }

        o += GetWheel5Increment(wheel5);

        const size_t q = (size_t)backward5(p);
        if (notPrime[q / bitsPerWord] == 1ULL << (q % bitsPerWord)) {
            continue;
        }

        knownPrimes.push_back(p);
    }

    return knownPrimes;
}

std::vector<BigInteger> SegmentedSieveOfEratosthenes(const BigInteger& n, const size_t& limit)
{
    // `backward(n)` counts assuming that multiples
    // of 2 and 3 have been removed.
    if ((n < 7U) || (backward5(n) <= limit)) {
        return SieveOfEratosthenes(n);
    }

    // Process segments of length `limit` at a time.
    BigInteger low = makeNotMultiple(limit);
    BigInteger high = makeNotMultiple(low + limit);

    // Compute all primes smaller than or equal to limit using simple sieve
    std::vector<BigInteger> knownPrimes = SieveOfEratosthenes(limit);
    knownPrimes.reserve(std::expint(log(n)) - std::expint(log(2)));
    dispatch.resetResult();

    // Process one segment at a time until we pass n
    while (low < n) {
        if (high >= n) {
            high = makeNotMultiple(n);
        }

        // Cardinality with multiples of 2 and 3 removed is 1/3 of total.
        const BigInteger bLow = backward5(low);
        const size_t cardinality = (size_t)(backward5(high) - bLow);
        size_t bitsPerWord = 64U;
        const size_t arrayWidth = (cardinality + bitsPerWord) / bitsPerWord;
        std::unique_ptr<uint64_t> uNotPrime(new uint64_t[arrayWidth]);
        std::unique_ptr<std::mutex> uNotPrimeMutex(new std::mutex[arrayWidth]);
        uint64_t* notPrime = uNotPrime.get();
        std::mutex* notPrimeMutex = uNotPrimeMutex.get();

        // Use the found primes by simpleSieve() to find
        // primes in current range
        for (size_t k = 3U; k < knownPrimes.size(); k++) {
            // Find the minimum number in [low..high] that is
            // a multiple of prime[i] (divisible by prime[i])
            // For example, if low is 31 and prime[i] is 3,
            // we start with 33.
            const BigInteger& p = knownPrimes[k];

            dispatch.dispatch([&bLow, &high, &low, &cardinality, &bitsPerWord, p, &notPrime, &notPrimeMutex]() {
                // We are skipping multiples of 2, 3, and 5
                // for space complexity, for 4/15 the bits.
                // More are skipped by the wheel for time.
                const BigInteger p2 = p << 1U;
                const BigInteger p4 = p << 2U;
                BigInteger i = ((low + p - 1U) / p) * p;
                if ((i & 1U) == 0U) {
                    i += p;
                }
                while (((i % 3U) == 0U) || ((i % 5U) == 0U)) {
                    i += p2;
                }

                // "p" already definitely not a multiple of 3.
                // Its remainder when divided by 3 can be 1 or 2.
                // If it is 2, we can do a "half iteration" of the
                // loop that would handle remainder of 1, and then
                // we can proceed with the 1 remainder loop.
                // This saves 2/3 of updates (or modulo).
                if ((i % 3U) == 2U) {
                    const size_t q = (size_t)(backward5(i) - bLow);
                    if (q > cardinality) {
                        return false;
                    }
                    const size_t w = q / bitsPerWord;
                    const uint64_t b = 1ULL << (q % bitsPerWord);
                    if (true) {
                        std::lock_guard<std::mutex> lock(notPrimeMutex[w]);
                        notPrime[w] |= b;
                    }
                    i += p2;
                }

                for (;;) {
                    if (i % 5U) {
                        size_t q = (size_t)(backward5(i) - bLow);
                        if (q > cardinality) {
                            return false;
                        }
                        const size_t w = q / bitsPerWord;
                        const uint64_t b = 1ULL << (q % bitsPerWord);
                        if (true) {
                            std::lock_guard<std::mutex> lock(notPrimeMutex[w]);
                            notPrime[w] |= b;
                        }
                    }
                    i += p4;

                    if (i % 5U) {
                        size_t q = (size_t)(backward5(i) - bLow);
                        if (q > cardinality) {
                            return false;
                        }
                        const size_t w = q / bitsPerWord;
                        const uint64_t b = 1ULL << (q % bitsPerWord);
                        if (true) {
                            std::lock_guard<std::mutex> lock(notPrimeMutex[w]);
                            notPrime[w] |= b;
                        }
                    }
                    i += p2;
                }

                return false;
            });
        }
        dispatch.finish();

        // Numbers which are not marked as false are prime
        size_t q = 0U;
        for (size_t o = 0U; ; ++o) {
            const size_t p = forward(o + bLow);
            if (p > n) {
                break;
            }
            if ((p % 5U) == 0U) {
                continue;
            }
            if (notPrime[q]) {
                knownPrimes.push_back(p);
            }
            ++q;
        }

        // Update low and high for next segment
        low = makeNotMultiple(low + limit);
        high = makeNotMultiple(low + limit);
    }

    return knownPrimes;
}
} // namespace qimcifa

using namespace qimcifa;

// Driver Code
int main()
{
    BigInteger n = 1000000000U; // 1e9

    std::cout << "Primes up to number: ";
    std::cin >> n;

    std::cout << "Following are the prime numbers smaller than or equal to " << n << ":" << std::endl;

    // const std::vector<BigInteger> primes = TrialDivision(n);
    const std::vector<BigInteger> primes = SieveOfEratosthenes(n);

    for (BigInteger p : primes) {
        std::cout << p << " ";
    }
    std::cout << std::endl;

    return 0;
}
