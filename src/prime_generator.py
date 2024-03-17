# Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
# C++ program to print all primes smaller than or equal to
# n using Sieve of Eratosthenes

# Improved by Dan Strano of Unitary Fund, 2024.
# We can think of trial division as exact inverse of
# Sieve of Eratosthenes, with log space and log time.
# The modular division part is a costly atomic operation.
# It need only be carried out up the square root of the
# number under trial. Multiples of 2, 3, 5, 7, and 11 can
# be entirely skipped in loop enumeration.

import math

def backward(ni):
    ni = (ni + 1) >> 1
    ni = int(((ni + 1) << 1) / 3)
    return ni;

def forward(p):
    # Make this NOT a multiple of 2 or 3.
    p = p + (p >> 1)
    return (p << 1) - 1

def isTrialDivisionMultiple(p, nextIndex, knownPrimes):
    sqrtP = math.isqrt(p);
    if (sqrtP * sqrtP) == p:
        return True

    highestIndex = 0
    m = (len(knownPrimes) + 1) >> 1
    while m > 1:
        if knownPrimes[highestIndex + m] <= sqrtP:
            highestIndex = highestIndex + m
        m = (m + 1) >> 1

    for i in knownPrimes[nextIndex:(highestIndex + 1)]:
        if (p % i) == 0:
            return True

    return False
 
def wheel_inc(primes):
    wheelPrimes = primes[:-1]
    prime = primes[-1]
    radius = 1
    for i in primes:
        radius *= i
    output = []
    counter = 1
    for i in range(1, radius):
        if not isTrialDivisionMultiple(i, 2, wheelPrimes):
            output.append((i % prime) == 0)
            counter = counter + 1

    return output

def wheel_gen(primes):
    output = []
    for i in range(3, len(primes) + 1):
        output.append(wheel_inc(primes[:i]))
        output[-1] = output[-1][1:] + output[-1][:1]
    return output
 
def TrialDivision(n):
    knownPrimes = [ 2, 3, 5 ]
    wheelPrimes = [ 2, 3 ]

    if n < (knownPrimes[-1] + 2):
        return [p for p in knownPrimes if p <= n]

    # We are excluding multiples of the first few
    # small primes from outset. For multiples of
    # 2 and 3, this reduces complexity by 2/3.
    # cardinality = int((~((~n) | 1)) / 3)

    # From here, for each new prime we find, if it is
    # less than or equal to wheel_limit, we build a
    # new "inside-out" wheel.
    inc_seqs = wheel_gen(wheelPrimes) if len(wheelPrimes) > 2 else []
    wheel_limit = 11

    # Get the remaining prime numbers.
    o = 1
    while True:
        o = o + 1
        is_wheel_multiple = False
        for i in range(len(inc_seqs)):
            is_wheel_multiple = inc_seqs[i][0]
            inc_seqs[i] = inc_seqs[i][1:] + inc_seqs[i][:1]
            if is_wheel_multiple:
                break

        p = forward(o)

        if is_wheel_multiple:
            continue

        if p > n:
            break
        if isTrialDivisionMultiple(p, len(wheelPrimes) - 1, knownPrimes):
            # Skip
            continue

        knownPrimes.append(p)
        if p <= wheel_limit:
            wheelPrimes.append(p)
            inc_seqs.append(wheel_inc(knownPrimes))
            inc_seqs[-1] = inc_seqs[-1][2:] + inc_seqs[-1][:2]

    return knownPrimes

# Driver Code
if __name__ == '__main__':
    n = 1000000

    print("Following are the prime numbers smaller than or equal to " + str(n) + ":")
    print(TrialDivision(n))
