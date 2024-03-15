# This script demonstrates the sequence that is being
# enumerated by "reverse trial division" in the prime
# genertor script.

def isTimeOrSpaceMultiple(p, knownPrimes):
    for kp in knownPrimes:
        if (p % kp) == 0:
            return True
    return False
 
def rtd(n):
    knownPrimes = [ 2, 3, 5, 7, 11 ]
    output = []
    counter = 1
    for i in range(1, n):
        if not isTimeOrSpaceMultiple(i, knownPrimes):
            output.append((counter, i))
            counter = counter + 1
    print(output)

# Driver Code
if __name__ == '__main__':
    print(rtd(1000))
