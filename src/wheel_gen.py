# This script demonstrates the sequence that is being
# enumerated by "reverse trial division" in the prime
# genertor script.

def isMultiple(p, knownPrimes):
    for kp in knownPrimes:
        if (p % kp) == 0:
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
        if not isMultiple(i, wheelPrimes):
            output.append((i % prime) == 0)
            print((counter, i, output[-1]))
            counter = counter + 1

    return output

def wheel_gen(primes):
    output = []
    for i in range(1, len(primes) + 1):
        output.append(wheel_inc(primes[:i]))
        print()
    return output

# Driver Code
if __name__ == '__main__':
    print(wheel_gen([2, 3, 5]))
