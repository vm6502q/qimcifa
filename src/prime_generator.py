# Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
# C++ program to print all primes smaller than or equal to
# n using Sieve of Eratosthenes

# 1/3 overall space complexity!
# log reduction in time complexity!
knownPrimes = [ 2, 3 ]

def backward(ni):
    ni = (ni + 1) >> 1
    ni = int(((ni + 1) << 1) / 3)
    return ni;

def forward(p):
    # Make this NOT a multiple of 2 or 3.
    p = p + (p >> 1)
    return (p << 1) - 1

def isTimeOrSpaceMultiple(p):
    for i in knownPrimes:
        if (p % i) == 0:
            return True
    return False

def isTimeMultiple(p):
    for i in knownPrimes[2:]:
        if (p % i) == 0:
            return True
    return False
 
def SieveOfEratosthenes(n):
    if n < 2:
        return set([])
    if n < 3:
        return set([2])
    if n < 5:
        return set([2, 3])

    # We are excluding multiples of the first few
    # small primes from outset. For multiples of
    # 2 and 3, this reduces complexity by 2/3.
    cardinality = int((n & ~1) / 3)

    # Create a boolean array "prime[0..cardinality]"
    # and initialize all entries it as true. Rather,
    # reverse the true/false meaning, so we can use
    # default initialization. A value in notPrime[i]
    # will finally be false only if i is a prime.
    notPrime= [False] * (cardinality + 1)
 
    o = 2;
    while True:
        p = forward(o)
        if (p * p) > n:
            break

        if isTimeMultiple(p):
            # Skip
            o = o + 1
            continue

        # If prime[o] is not changed, then it is a prime
        if not notPrime[o]:            
            knownPrimes.append(p)

        # Increment "while" loop.
        o = o + 1;
    
    outputPrimes = set(knownPrimes)
 
    # Get the remaining prime numbers.
    for o in range(2, cardinality + 1):
        if not notPrime[o]:
            p = forward(o)

            if isTimeMultiple(p):
                continue

            outputPrimes.add(p);

    return outputPrimes
 
# Driver Code
if __name__ == '__main__':
    n = 100

    print("Following are the prime numbers smaller than or equal to ", n, ":")
    print(SieveOfEratosthenes(n))
