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

def isTimeOrSpaceMultiple(p, knownPrimes):
    sqrt_p = math.isqrt(p)
    if (sqrt_p * sqrt_p) == p:
        return True
    for kp in knownPrimes:
        if kp >= sqrt_p:
            return False
        if (p % kp) == 0:
            return True
    return False

def isTimeMultiple(p, knownPrimes):
    sqrt_p = math.isqrt(p)
    if (sqrt_p * sqrt_p) == p:
        return True
    for kp in knownPrimes[5:]:
        if kp >= sqrt_p:
            return False
        if (p % kp) == 0:
            return True
    return False
 
def TrialDivision(n):
    knownPrimes = [ 2, 3, 5 ]

    if n < 7:
        return [p for p in knownPrimes if p <= n]

    # We are excluding multiples of the first few
    # small primes from outset. For multiples of
    # 2 and 3, this reduces complexity by 2/3.
    # cardinality = int((~((~n) | 1)) / 3)
 
    # Get the remaining prime numbers.
    o = 3
    lcv7 = -11
    lcv11 = -16
    # wheel = 1
    isWorking = True
    while isWorking:
        for i in range(0, 6):
            if lcv7 == 11:
                lcv7 = 1
                continue
            if lcv7 == 7:
                lcv7 = 8
                continue
            lcv7 = lcv7 + 1

            if lcv11 == 17:
                lcv11 = 1
                continue
            if lcv11 == 7:
                lcv11 = 8
                continue
            lcv11 = lcv11 + 1

            p = forward(o + i)

            if p > n:
                isWorking = False
                break

            if isTimeMultiple(p, knownPrimes):
                # Skip
                continue

            knownPrimes.append(p)

        for i in range(7, 9):
            if lcv7 == 11:
                lcv7 = 1
                continue
            if lcv7 == 7:
                lcv7 = 8
                continue
            lcv7 = lcv7 + 1

            if lcv11 == 17:
                lcv11 = 1
                continue
            if lcv11 == 7:
                lcv11 = 8
                continue
            lcv11 = lcv11 + 1

            p = forward(o + i)

            if p > n:
                isWorking = False
                break

            if isTimeMultiple(p, knownPrimes):
                # Skip
                continue

            knownPrimes.append(p)

        o = o + 10

    return knownPrimes

# Driver Code
if __name__ == '__main__':
    n = 100

    print("Following are the prime numbers smaller than or equal to " + str(n) + ":")
    print(TrialDivision(n))
