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
    // Primes up to 36
    std::vector<BigInteger> knownPrimes = { 2, 3, 5, 11, 13, 17, 19, 23, 29, 31 };
    if (n < 2) {
        return std::vector<BigInteger>();
    }

    if (n < (knownPrimes.back() + 2)) {
        const auto highestPrimeIt = std::upper_bound(knownPrimes.begin(), knownPrimes.end(), n);
        return std::vector<BigInteger>(knownPrimes.begin(), highestPrimeIt);
    }

    // Wheels up to 5
    knownPrimes = { 2, 3, 5 };

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
    inc_seqs.erase(inc_seqs.begin(), inc_seqs.begin() + 2U);
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
} // namespace qimcifa

using namespace qimcifa;

// Driver Code
int main()
{
    BigInteger n = 10000000000U; // 1e10

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
