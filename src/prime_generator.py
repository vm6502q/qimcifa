# Source: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
# C++ program to print all primes smaller than or equal to
# n using Sieve of Eratosthenes

# Improvements by Dan Strano of Unitary Fund, 2024:
# log space complexity!
# log reduction in time complexity!
knownPrimes = [ 2, 3, 5 ]

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
 
def TrialDivision(n):
    if n < 2:
        return []
    if n < 3:
        return [2]
    if n < 5:
        return [2, 3]
    if n < 7:
        return [2, 3, 5]

    # We are excluding multiples of the first few
    # small primes from outset. For multiples of
    # 2 and 3, this reduces complexity by 2/3.
    cardinality = int((~((~n) | 1)) / 3)
 
    # Get the remaining prime numbers.
    o = 2
    isWorking = True
    while isWorking:
        for i in range(1, 7):
            p = forward(o + i)

            if p > n:
                isWorking = False
                break

            if isTimeMultiple(p):
                # Skip
                continue

            knownPrimes.append(p)

        for i in range(8, 10):
            p = forward(o + i)

            if p > n:
                isWorking = False
                break

            if isTimeMultiple(p):
                # Skip
                continue

            knownPrimes.append(p)

        o = o + 10

    return knownPrimes

# Driver Code
if __name__ == '__main__':
    n = 100

    print("Following are the prime numbers smaller than or equal to ", n, ":")
    print(TrialDivision(n))
