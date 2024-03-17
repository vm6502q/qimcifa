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

def isTrialDivisionMultiple(p, nextPrimeIndex, knownPrimes):
    sqrtP = math.isqrt(p);
    if (sqrtP * sqrtP) == p:
        return True

    for i in knownPrimes[nextPrimeIndex:]:
        if (i >= sqrtP):
            return False
        if (p % i) == 0:
            return True

    return False

def isMultiple(p, knownPrimes):
    for kp in knownPrimes:
        if (p % kp) == 0:
            return True
    return False
 
def wheel_inc(primes):
    wheelPrimes = primes[:-1]
    radius = 1
    for i in primes:
        radius *= i
    output = []
    counter = 1
    for i in range(1, radius):
        if not isMultiple(i, wheelPrimes):
            output.append(isMultiple(i, primes))
            counter = counter + 1

    output = output[1:] + output[:1]

    return output

def wheel_gen(primes):
    output = []
    for i in range(2, len(primes)):
        output.append(wheel_inc(primes[:i+1]))
    return output
 
def TrialDivision(n):
    wheelPrimes = [ 2, 3, 5, 7, 11 ]
    knownPrimes = [ 2, 3, 5, 7, 11, 13, 17, 19 ]

    if n < (knownPrimes[-1] + 2):
        return [p for p in knownPrimes if p <= n]

    # We are excluding multiples of the first few
    # small primes from outset. For multiples of
    # 2 and 3, this reduces complexity by 2/3.
    # cardinality = int((~((~n) | 1)) / 3)
 
    # Get the remaining prime numbers.
    o = 2
    inc_seqs = wheel_gen(wheelPrimes)
    isWorking = True
    while isWorking:
        for i in range(len(inc_seqs)):
            o = o + (2 if inc_seqs[i][0] else 1)
            inc_seqs[i] = inc_seqs[i][1:] + inc_seqs[i][:1]
    
        p = forward(o)

        if p > n:
            isWorking = False
            break

        if isTrialDivisionMultiple(p, 2, knownPrimes):
            # Skip
            continue

        knownPrimes.append(p)

    return knownPrimes

# Driver Code
if __name__ == '__main__':
    n = 100

    print("Following are the prime numbers smaller than or equal to " + str(n) + ":")
    print(TrialDivision(n))
