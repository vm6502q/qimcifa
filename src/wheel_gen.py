# This script demonstrates the sequence that is being
# enumerated by "reverse trial division" in the prime
# genertor script.

def isMultiple(p, knownPrimes):
    for kp in knownPrimes:
        if (p % kp) == 0:
            return True
    return False

def wheel_inc(primes):
    radius = 1
    for i in primes:
        radius *= i
    prime = primes[-1]
    wheelPrimes = primes[:-1]
    output = []
    counter = 1
    for i in range(1, radius + 1):
        if not isMultiple(i, wheelPrimes):
            isMult = (i % prime) == 0
            output.append(isMult)
            print((counter, i, isMult))
            counter = counter + 1

    output = output[1:] + output[:1]

    return output

def wheel_gen(primes):
    output = []
    for i in range(2, len(primes)):
        output.append(wheel_inc(primes[:i+1]))
        print()
    return output

# Driver Code
if __name__ == '__main__':
    print(wheel_gen([2, 3, 5, 7]))
