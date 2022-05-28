////////////////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2022. All rights reserved.
//
// "A quantum-inspired Monte Carlo integer factoring algorithm"
//
// This example demonstrates a (Shor's-like) "quantum-inspired" algorithm for integer factoring.
// This approach is similar to Shor's algorithm, except with a uniformly random output from the
// quantum period-finding subroutine. Therefore, we don't need quantum computer simulation for
// this algorithm at all!
//
// (This file was heavily adapted from
// https://github.com/ProjectQ-Framework/ProjectQ/blob/develop/examples/shor.py,
// with thanks to ProjectQ!)
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "common/oclengine.hpp"

#include <chrono>
#include <cmath>
#include <iomanip> // For setw
#include <iostream> // For cout
#include <random>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <atomic>
#include <future>
#include <map>
#include <mutex>

// Turn this off, if you're not factoring a semi-prime number with equal-bit-width factors.
#define IS_RSA_SEMIPRIME 1
// Turn this off, if you don't want to coordinate across multiple (quasi-independent) nodes.
#define IS_DISTRIBUTED 1
// The maximum number of bits in Boost big integers is 2^QBCAPPOW.
// (2^7, only, needs custom std::cout << operator implementation.)
#define QBCAPPOW 8U
#define bitsInByte 8U

#if QBCAPPOW < 8U
#define bitLenInt uint8_t
#elif QBCAPPOW < 16U
#define bitLenInt uint16_t
#elif QBCAPPOW < 32U
#define bitLenInt uint32_t
#else
#define bitLenInt uint64_t
#endif

#if QBCAPPOW < 6U
#define bitsInCap 32
#define bitCapInt uint32_t
#elif QBCAPPOW < 7U
#define bitsInCap 64
#define bitCapInt uint64_t
#else
#define bitsInCap (8U * (((bitLenInt)1U) << QBCAPPOW))
#include "big_integer.hpp"
#define bitCapInt BigInteger
#endif

#if QBCAPPOW < 7U
#define bci_copy(a) (a)
#define bci_create(a) (a)
#define bci_add(l, r) (l + r)
#define bci_sub(l, r) (l - r)
#define bci_mul(l, r) (l * r)
#define bci_div(l, r) (l / r)
#define bci_mod(l, r) (l % r)
#define bci_lshift(l, r) (l << r)
#define bci_rshift(l, r) (l >> r)
#define bci_or(l, r) (l | r)
#define bci_and(l, r) (l & r)
#define bci_and_1(l) (l & 1U)
#define bci_eq_0(a) (a == 0)
#define bci_neq_0(a) (a != 0)
#define bci_eq_1(a) (a == 1)
#define bci_eq(l, r) (l == r)
#define bci_lt(l, r) (l < r)
#define bci_gt(l, r) (l > r)
#define bci_geq(l, r) (l >= r)
#define bci_leq(l, r) (l <= r)
#else
#define bci_copy(a) BigInteger(a)
#define bci_create(a) bi_create(a)
#define bci_add(l, r) bi_add(l, r)
#define bci_sub(l, r) bi_sub(l, r)
#define bci_mul(l, r) bi_mul(l, r)
#define bci_div(l, r) bi_div(l, r)
#define bci_mod(l, r) bi_mod(l, r)
#define bci_lshift(l, r) bi_lshift(l, r)
#define bci_rshift(l, r) bi_rshift(l, r)
#define bci_or(l, r) bi_or(l, r)
#define bci_and(l, r) bi_and(l, r)
#define bci_and_1(l) bi_and_1(l)
#define bci_eq_0(a) (bi_compare(a, bi_empty()) == 0)
#define bci_neq_0(a) (bi_compare(a, bi_empty()) != 0)
#define bci_eq_1(a) (bi_compare(a, bi_create(1)))
#define bci_eq(l, r) (bi_compare(l, r) == 0)
#define bci_lt(l, r) (bi_compare(l, r) < 0)
#define bci_gt(l, r) (bi_compare(l, r) > 0)
#define bci_geq(l, r) (bi_compare(l, r) >= 0)
#define bci_leq(l, r) (bi_compare(l, r) <= 0)
#endif

namespace Qimcifa {

const bitCapInt ONE_BCI = bci_create(1U);

std::ostream& operator<<(std::ostream& os, const bitCapInt& bci) {
    const bitCapInt bci_10 = bci_create(10);

    bitCapInt b = bci_copy(bci);
    std::vector<std::string> digits;
    while (bci_neq_0(b)) {
        digits.push_back(std::to_string(b.bits[0] % 10));
        b = bci_div(b, bci_10);
    }
    for (int i = digits.size() - 1; i >= 0; i--) {
        os << digits[i];
    }

    return os;
}

std::istream &operator>>(std::istream &is, bitCapInt& b)
{
    const bitCapInt bci_10 = bci_create(10);

    std::string input;
    is >> input;

    b = bci_create(0);
    for (int i = input.size() - 1; i >= 0; i--) {
        b = bci_mul(b, bci_10);
        b = bci_add(b, bci_create(input[i] - 48));
    }

    return is;
}

// Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
inline bool isPowerOfTwo(const bitCapInt& x) { return (bci_neq_0(x) && bci_eq_0(bci_and(x, bci_sub(x, ONE_BCI)))); }

inline bitLenInt log2(const bitCapInt& n)
{
#if __GNUC__ && QBCAPPOW < 7
// Source: https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers#answer-11376759
#if QBCAPPOW < 6
    return (bitLenInt)(bitsInByte * sizeof(unsigned int) - __builtin_clz((unsigned int)n) - 1U);
#else
    return (bitLenInt)(bitsInByte * sizeof(unsigned long long) - __builtin_clzll((unsigned long long)n) - 1U);
#endif
#else
    bitLenInt pow = 0U;
    bitCapInt p = bci_rshift(n, 1U);
    while (bci_neq_0(p)) {
        p = bci_rshift(p, 1U);
        pow++;
    }
    return pow;
#endif
}

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
inline bitCapInt uipow(const bitCapInt& base, const bitCapInt& exp)
{
    bitCapInt result = ONE_BCI;
    bitCapInt b = base;
    bitCapInt e = exp;
    for (;;) {
        if (bci_and_1(b)) {
            result = bci_mul(result, b);
        }
        e = bci_rshift(e, 1U);
        if (bci_eq_0(e)) {
            break;
        }
        b = bci_mul(b, b);
    }

    return result;
}

// It's fine if this is not exact for the whole bitCapInt domain, so long as it is <= the exact result.
inline bitLenInt intLog(const bitCapInt& base, const bitCapInt& arg)
{
    bitLenInt result = 0U;
    for (bitCapInt x = bci_copy(arg); bci_geq(x, base); x = bci_div(x, base)) {
        result++;
    }
    return result;
}

// Adapted from Gaurav Ahirwar's suggestion on https://www.geeksforgeeks.org/square-root-of-an-integer/
bitCapInt floorSqrt(const bitCapInt& x)
{
    // Base cases
    if (bci_eq_0(x) || bci_eq_1(x)) {
        return x;
    }

    // Binary search for floor(sqrt(x))
    bitCapInt start = ONE_BCI, end = bci_rshift(x, 1U), ans = bci_create(0U);
    while (bci_leq(start, end)) {
        bitCapInt mid = bci_rshift(bci_add(start, end), 1U);

        // If x is a perfect square
        bitCapInt sqr = bci_mul(mid, mid);
        if (bci_eq(sqr, x)) {
            return mid;
        }

        if (bci_lt(sqr, x)) {
            // Since we need floor, we update answer when mid*mid is smaller than x, and move closer to sqrt(x).
            start = bci_add(mid, ONE_BCI);
            ans = bci_copy(mid);
        } else {
            // If mid*mid is greater than x
            end = bci_sub(mid, ONE_BCI);
        }
    }
    return ans;
}

bitCapInt gcd(const bitCapInt& n1, const bitCapInt& n2)
{
    if (bci_neq_0(n2)) {
        return gcd(n2, bci_mod(n1, n2));
    }
    return n1;
}

class bad_alloc : public std::bad_alloc {
private:
    std::string m;

public:
    bad_alloc(std::string message)
        : m(message)
    {
        // Intentionally left blank.
    }

    const char* what() const noexcept { return m.c_str(); }
};

typedef std::shared_ptr<cl::Buffer> BufferPtr;

BufferPtr MakeBuffer(const cl::Context& context, cl_mem_flags flags, size_t size, void* host_ptr = NULL)
{
    cl_int error;
    BufferPtr toRet = std::make_shared<cl::Buffer>(context, flags, size, host_ptr, &error);
    if (error != CL_SUCCESS) {
        if (error == CL_MEM_OBJECT_ALLOCATION_FAILURE) {
            throw bad_alloc("CL_MEM_OBJECT_ALLOCATION_FAILURE in QEngineOCL::MakeBuffer()");
        }
        if (error == CL_OUT_OF_HOST_MEMORY) {
            throw bad_alloc("CL_OUT_OF_HOST_MEMORY in QEngineOCL::MakeBuffer()");
        }
        if (error == CL_INVALID_BUFFER_SIZE) {
            throw bad_alloc("CL_INVALID_BUFFER_SIZE in QEngineOCL::MakeBuffer()");
        }
        throw std::runtime_error("OpenCL error code on buffer allocation attempt: " + std::to_string(error));
    }

    return toRet;
}

} // namespace Qimcifa

using namespace Qimcifa;

int main()
{
    if (!(OCLEngine::Instance().GetDeviceCount())) {
        throw std::runtime_error("Tried to initialize QEngineOCL, but no available OpenCL devices.");
    }

    typedef std::uniform_int_distribution<uint64_t> rand_dist;

    bitCapInt toFactor;
    uint64_t nodeCount = 1U;
    uint64_t nodeId = 0U;

    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    auto iterClock = std::chrono::high_resolution_clock::now();

    const bitLenInt qubitCount = log2(toFactor) + (isPowerOfTwo(toFactor) ? 0U : 1U);
    // const bitCapInt qubitPower = ONE_BCI << qubitCount;
    std::cout << "Bits to factor: " << (int)qubitCount << std::endl;

#if IS_DISTRIBUTED
    std::cout << "You can split this work across nodes, without networking!" << std::endl;
    do {
        std::cout << "Number of nodes (>=1): ";
        std::cin >> nodeCount;
        if (!nodeCount) {
            std::cout << "Invalid node count choice!" << std::endl;
        }
    } while (!nodeCount);
    if (nodeCount > 1U) {
        do {
            std::cout << "Which node is this? (0-" << (int)(nodeCount - 1U) << "):";
            std::cin >> nodeId;
            if (nodeId >= nodeCount) {
                std::cout << "Invalid node ID choice!" << std::endl;
            }
        } while (nodeId >= nodeCount);
    }
#endif

    std::random_device rand_dev;
    std::mt19937 rand_gen(rand_dev());

    std::atomic<bool> isFinished;
    isFinished = false;

#if IS_RSA_SEMIPRIME
    std::map<bitLenInt, const std::vector<bitCapInt>> primeDict;// = { { 16U, { 32771U, 65521U } },
    //    { 28U, { 134217757U, 268435399U } }, { 32U, { 2147483659U, 4294967291U } },
    //    { 64U, { 9223372036854775837U, 1844674407370955143U } } };

    // If n is semiprime, \phi(n) = (p - 1) * (q - 1), where "p" and "q" are prime.
    // The minimum value of this formula, for our input, without consideration of actual
    // primes in the interval, is as follows:
    // (See https://www.mobilefish.com/services/rsa_key_generation/rsa_key_generation.php)
    const bitLenInt primeBits = (qubitCount + 1U) >> 1U;
    const bitCapInt fullMin = bci_lshift(ONE_BCI, primeBits - 1U);
    const bitCapInt fullMax = bci_sub(bci_lshift(ONE_BCI, primeBits), ONE_BCI);
    const bitCapInt minPrime = primeDict[primeBits].size() ? primeDict[primeBits][0] : bci_add(fullMin, ONE_BCI);
    const bitCapInt maxPrime = primeDict[primeBits].size() ? primeDict[primeBits][1] : fullMax;
    const bitCapInt fullMinR = bci_mul(bci_sub(minPrime, ONE_BCI), bci_sub(bci_div(toFactor, minPrime), ONE_BCI));
    const bitCapInt fullMaxR = bci_mul(bci_sub(maxPrime, ONE_BCI), bci_sub(bci_div(toFactor, maxPrime), ONE_BCI));
#else
    // \phi(n) is Euler's totient for n. A loose lower bound is \phi(n) >= sqrt(n/2).
    const bitCapInt fullMinR = floorSqrt(bci_rshift(toFactor, 1U));
    // A better bound is \phi(n) >= pow(n / 2, log(2)/log(3))
    // const bitCapInt fullMinR = pow(toFactor / 2, PHI_EXPONENT);

    // It can be shown that the period of this modular exponentiation can be no higher than 1
    // less than the modulus, as in https://www2.math.upenn.edu/~mlazar/math170/notes06-3.pdf.
    // Further, an upper bound on Euler's totient for composite numbers is n - sqrt(n). (See
    // https://math.stackexchange.com/questions/896920/upper-bound-for-eulers-totient-function-on-composite-numbers)
    const bitCapInt fullMaxR = bci_sub(toFactor, floorSqrt(toFactor));
#endif

    const DeviceContextPtr deviceContext = OCLEngine::Instance().GetDeviceContextPtr(-1);
    const cl::Context context = deviceContext->context;
    const cl::CommandQueue queue = deviceContext->queue;

    const size_t groupSize = deviceContext->GetPreferredSizeMultiple();
    const size_t groupCount = deviceContext->GetPreferredConcurrency();
    const size_t itemCount = groupSize * groupCount;

    const size_t batchSize = 1U << 9U;
    BufferPtr ulongArgsBufferPtr = MakeBuffer(context, CL_MEM_READ_ONLY, sizeof(uint64_t) * 4);
    uint64_t ulongArgs[4] = { batchSize, nodeId, nodeCount, (uint64_t)IS_RSA_SEMIPRIME };
    cl_int error = queue.enqueueWriteBuffer(*ulongArgsBufferPtr, CL_TRUE, 0U, sizeof(uint64_t) * 4, ulongArgs, NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }

    BufferPtr bciArgsBufferPtr = MakeBuffer(context, CL_MEM_READ_ONLY, sizeof(bitCapInt) * 3);
    bitCapInt bciArgs[3] = { toFactor, fullMinR, fullMaxR };
    error = queue.enqueueWriteBuffer(*bciArgsBufferPtr, CL_TRUE, 0U, sizeof(bitCapInt) * 3, bciArgs, NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }
    
    BufferPtr rngSeedBufferPtr = MakeBuffer(context, CL_MEM_READ_WRITE, sizeof(uint64_t) * itemCount * 4U);
    std::unique_ptr<uint64_t[]> rngSeeds(new uint64_t[itemCount * 4]);
    rand_dist seedDist(0, 0xFFFFFFFFFFFFFFFF);
    rand_dist cSeedDist(0, 0x3FFFFFFFFFFFFFF);
    for (size_t i = 0U; i < itemCount; i++) {
        rngSeeds.get()[i * 4 + 0] = seedDist(rand_gen);
        rngSeeds.get()[i * 4 + 1] = cSeedDist(rand_gen);
        rngSeeds.get()[i * 4 + 2] = seedDist(rand_gen);
        rngSeeds.get()[i * 4 + 3] = seedDist(rand_gen);
    }
    error = queue.enqueueWriteBuffer(*rngSeedBufferPtr, CL_TRUE, 0U, sizeof(uint64_t) * itemCount * 4U, rngSeeds.get(), NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }
    
    std::unique_ptr<bitCapInt[]> outputArray(new bitCapInt[itemCount]());
    BufferPtr outputBufferPtr = MakeBuffer(context, CL_MEM_WRITE_ONLY, sizeof(bitCapInt) * itemCount);
    error = queue.enqueueWriteBuffer(*outputBufferPtr, CL_TRUE, 0U, sizeof(bitCapInt) * itemCount, outputArray.get(), NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }

    // We have to reserve the kernel, because its argument hooks are unique. The same kernel therefore can't be used by
    // other QEngineOCL instances, until we're done queueing it.
    OCLDeviceCall ocl = deviceContext->Reserve(OCL_API_QIMCIFA_BATCH);

    ocl.call.setArg(0, *ulongArgsBufferPtr);
    ocl.call.setArg(1, *bciArgsBufferPtr);
    ocl.call.setArg(2, *rngSeedBufferPtr);
    ocl.call.setArg(3, *outputBufferPtr);

    bitCapInt testFactor = ONE_BCI;
    while (bci_leq(testFactor, ONE_BCI)) {
        cl::Event kernelEvent;
        error = queue.enqueueNDRangeKernel(ocl.call, cl::NullRange, // kernel, offset
            cl::NDRange(groupCount), // global number of work items
            cl::NDRange(groupSize), // local number (per group)
            NULL, // vector of events to wait for
            &kernelEvent); // handle to wait for the kernel

        if (error != CL_SUCCESS) {
            throw std::runtime_error("Failed to enqueue kernel, error code: " + std::to_string(error));
        }

        kernelEvent.wait();
        queue.enqueueReadBuffer(*outputBufferPtr, CL_TRUE, 0U, sizeof(bitCapInt) * itemCount, outputArray.get(), NULL);
        for (size_t i = 0; i < itemCount; i++) {
            testFactor = outputArray.get()[i];
            if (bci_geq(testFactor, ONE_BCI) && bci_eq(bci_mul(bci_div(toFactor, testFactor), testFactor), toFactor)) {
                break;
            }
        }
    }

    std::cout << "Success: " << testFactor << " * " << bci_div(toFactor, testFactor) << std::endl;
    const double clockFactor = 1.0 / 1000.0; // Report in ms
    auto tClock = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;

    return 0;
}
