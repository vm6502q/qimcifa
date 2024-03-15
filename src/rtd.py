# Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
# C++ program to print all primes smaller than or equal to
# n using Sieve of Eratosthenes

# Improved by Dan Strano of Unitary Fund, 2024.
# We can think of trial division as exact inverse of
# Sieve of Eratosthenes, with log space and log time.
# The modular division part is a costly atomic operation.
# It need only be carried out up the square root of the
# number under trial. Multiples of 2, 3, 5, and 7 can be
# entirely skipped in loop enumeration.

def isTimeOrSpaceMultiple(p, knownPrimes):
    for kp in knownPrimes:
        if (p % kp) == 0:
            return True
    return False
 
def rtd(n):
    knownPrimes = [ 2, 3, 5, 7 ]
    output = []
    counter = 1
    for i in range(1, n):
        if not isTimeOrSpaceMultiple(i, knownPrimes):
            output.append((counter, i))
            counter = counter + 1
    print(output)

# Driver Code
if __name__ == '__main__':
    n = 1000
    print(rtd(n))
