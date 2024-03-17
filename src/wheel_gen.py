# This script demonstrates the sequence that is being
# enumerated by "reverse trial division" in the prime
# genertor script.

def isTimeOrSpaceMultiple(p, knownPrimes):
    for kp in knownPrimes:
        if (p % kp) == 0:
            return True
    return False
 
def wheel_gen():
    qimcifaPrimes = [ 2, 3, 5, 7, 11 ]
    wheelPrimes = [ 2, 3, 5, 7, 11, 13 ]
    radius = 1
    for i in wheelPrimes:
        radius *= i
    output = []
    counter = 1
    for i in range(1, radius):
        if not isTimeOrSpaceMultiple(i, qimcifaPrimes):
            output.append((counter, i, not isTimeOrSpaceMultiple(i, wheelPrimes)))
            counter = counter + 1
    print(output)

# Driver Code
if __name__ == '__main__':
    print(wheel_gen())
