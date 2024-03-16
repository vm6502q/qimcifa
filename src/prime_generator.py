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
    knownPrimes = [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 77, 79, 83, 91, 97, 103, 107, 109, 119, 127, 133, 137, 139, 149, 151, 161, 163, 167 ]

    if n < 170:
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

            # **Hear me out**: We've "solved" up to multiples of 11.
            # It's trivial to know much higher primes than this.
            # At any such boundary of our knowledge, we can assume
            # that the highest prime necessary to know, to skip the
            # beginning work of the algorithm, would be the square
            # of the highest "inside-out" Wheel Factorization prime.
            #
            # Grant me only one step further, that the least expensive
            # way to remove 13 from here might be n % 13. For the edge
            # case, < 170 (13*13+1=169+1) is skipped, if we can know
            # that many primes (or obviously higher, hard storage).
            if p < 170:
                # Skip
                continue

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

            # **SEE LONG NOTE ABOVE**
            if p < 170:
                # Skip
                continue

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
