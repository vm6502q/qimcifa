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
#define QBCAPPOW 7U

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
#include "big_integer.h"
#define bitCapInt BigInteger
#endif

#if QBCAPPOW < 7U
#define bci_create(a) (a)
#define bci_copy(a, o) *(o) = a
#define bci_add(l, r, o) *(o) = ((l) + (r))
#define bci_sub(l, r, o) *(o) = ((l) - (r))
#define bci_mul(l, r, o) *(o) = ((l) * (r))
#define bci_div(l, r, o) *(o) = ((l) / (r))
#define bci_mod(l, r, o) *(o) = ((l) % (r))
#define bci_lshift(l, r, o) *(o) = ((l) << (r))
#define bci_rshift(l, r, o) *(o) = ((l) >> (r))
#define bci_or(l, r, o) *(o) = ((l) | (r))
#define bci_and(l, r, o) *(o) = ((l) & (r))
#define bci_eq(l, r) ((l) == (r))
#define bci_neq(l, r) ((l) != (r))
#define bci_lt(l, r) ((l) < (r))
#define bci_gt(l, r) ((l) > (r))
#define bci_geq(l, r) ((l) >= (r))
#define bci_leq(l, r) ((l) <= (r))
#define bci_and_1(l) ((l) & 1U)
#define bci_compare(a, b) ((a > b) ? 1 : ((a < b) ? -1 : 0))
#define bci_eq_0(a) ((a) == 0)
#define bci_neq_0(a) ((a) != 0)
#define bci_eq_1(a) ((a) == 1)
#define bci_neq_1(a) ((a) != 1)
#define bci_gt_1(a) ((a) > 1)
#else
#define bci_create(a) bi_create(a)
#define bci_copy(a, o) bi_copy(&(a), o)
#define bci_add(l, r, o) bi_add(&(l), &(r), o)
#define bci_sub(l, r, o) bi_sub(&(l), &(r), o)
#define bci_mul(l, r, o) bi_mul(&(l), &(r), o)
#define bci_div(l, r, o) bi_div_mod(&(l), &(r), o, NULL)
#define bci_mod(l, r, o) bi_div_mod(&(l), &(r), NULL, o)
#define bci_lshift(l, r, o) bi_lshift(&(l), r, o)
#define bci_rshift(l, r, o) bi_rshift(&(l), r, o)
#define bci_or(l, r, o) bi_or(&(l), &(r), o)
#define bci_and(l, r, o) bi_and(&(l), &(r), o)
#define bci_eq(l, r) (bi_compare(&(l), &(r)) == 0)
#define bci_neq(l, r) (bi_compare(&(l), &(r)) != 0)
#define bci_lt(l, r) (bi_compare(&(l), &(r)) < 0)
#define bci_gt(l, r) (bi_compare(&(l), &(r)) > 0)
#define bci_geq(l, r) (bi_compare(&(l), &(r)) >= 0)
#define bci_leq(l, r) (bi_compare(&(l), &(r)) <= 0)
#define bci_and_1(l) bi_and_1(&(l))
#define bci_compare(a, b) (bi_compare(&(a), &(b)))
#define bci_eq_0(a) (bi_compare_0(&(a)) == 0)
#define bci_neq_0(a) (bi_compare(&(a), &ZERO_BCI) != 0)
#define bci_eq_1(a) (bi_compare(&(a), &ONE_BCI) == 0)
#define bci_neq_1(a) (bi_compare(&(a), &ONE_BCI) != 0)
#define bci_gt_1(a) (bi_compare(&(a), &ONE_BCI) > 0)
#endif

namespace Qimcifa {

const bitCapInt ZERO_BCI = bci_create(0U);
const bitCapInt ONE_BCI = bci_create(1U);

#if QBCAPPOW > 6U
std::ostream& operator<<(std::ostream& os, bitCapInt b) {
    const bitCapInt bci_10 = bci_create(10);

    std::vector<std::string> digits;
    while (bci_neq_0(b)) {
        digits.push_back(std::to_string(b.bits[0] % 10));
        bitCapInt t(b);
        bci_div(t, bci_10, &b);
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

    b = ZERO_BCI;
    for (unsigned i = 0; i < input.size(); i++) {
        bitCapInt t(b);
        bci_mul(t, bci_10, &b);
        bci_copy(b, &t);
        bitCapInt c = bci_create(input[i] - 48);
        bci_add(t, c, &b);
    }

    return is;
}
#endif

// Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
bool isPowerOfTwo(const bitCapInt& x) {
    bitCapInt t1;
    bci_sub(x, ONE_BCI, &t1);
    bitCapInt t2;
    bci_and(x, t1, &t2);

    return bci_neq_0(x) && bci_eq_0(t2);
}

bitLenInt log2(const bitCapInt& n)
{
#if __GNUC__ && QBCAPPOW < 7
// Source: https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers#answer-11376759
#if QBCAPPOW < 6
    return (bitLenInt)(bitsInByte * sizeof(unsigned int) - __builtin_clz((unsigned int)n) - 1U);
#else
    return (bitLenInt)(bitsInByte * sizeof(unsigned long long) - __builtin_clzll((unsigned long long)n) - 1U);
#endif
#else
    return bi_log2(&n);
#endif
}

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
void uipow(bitCapInt b, bitCapInt e, bitCapInt* result)
{
    bci_copy(ONE_BCI, result);
    for (;;) {
        if (bci_and_1(b)) {
            bitCapInt t(*result);
            bci_mul(t, b, result);
        }
        bitCapInt t;
        bci_copy(e, &t);
        bci_rshift(t, 1U, &e);
        if (bci_eq_0(e)) {
            break;
        }
        bci_copy(b, &t);
        bci_mul(t, t, &b);
    }
}

// It's fine if this is not exact for the whole bitCapInt domain, so long as it is <= the exact result.
bitLenInt intLog(const bitCapInt& base, const bitCapInt& arg)
{
    bitLenInt result = 0U;
    bitCapInt t;
    for (bitCapInt x(arg); bci_geq(x, base); bci_div(t, base, &x)) {
        result++;
        bci_copy(x, &t);
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
    bitCapInt ans = ZERO_BCI, start = ONE_BCI;
    bitCapInt end;
    bci_rshift(x, 1U, &end);
    while (bci_leq(start, end)) {
        bitCapInt t, mid;
        bci_add(start, end, &t);
        bci_rshift(t, 1U, &mid);

        // If x is a perfect square
        bci_mul(mid, mid, &t);
        if (bci_eq(t, x)) {
            return mid;
        }

        if (bci_lt(t, x)) {
            // Since we need floor, we update answer when mid*mid is smaller than x, and move closer to sqrt(x).
            bci_add(mid, ONE_BCI, &start);
            bci_copy(mid, &ans);
        } else {
            // If mid*mid is greater than x
            bci_sub(mid, ONE_BCI, &end);
        }
    }
    return ans;
}

void gcd(bitCapInt n1, bitCapInt n2, bitCapInt* result)
{
    while (bci_neq_0(n2)) {
        bitCapInt t1(n1), t2(n2);
        bci_copy(n2, &n1);
        bci_mod(t1, t2, &n2);
    }
    bci_copy(n1, result);
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
    std::map<bitLenInt, const std::vector<bitCapInt>> primeDict = { { 16U, { bci_create(32771U), bci_create(65521U) } },
        { 28U, { bci_create(134217757U), bci_create(268435399U) } }, { 32U, { bci_create(2147483659U), bci_create(4294967291U) } },
        { 64U, { bci_create(9223372036854775837U), bci_create(1844674407370955143U) } } };

    // If n is semiprime, \phi(n) = (p - 1) * (q - 1), where "p" and "q" are prime.
    // The minimum value of this formula, for our input, without consideration of actual
    // primes in the interval, is as follows:
    // (See https://www.mobilefish.com/services/rsa_key_generation/rsa_key_generation.php)
    const bitLenInt primeBits = (qubitCount + 1U) >> 1U;

    bitCapInt fullMin;
    bci_lshift(ONE_BCI, primeBits - 1U, &fullMin);

    bitCapInt fullMax, t1;
    bci_lshift(ONE_BCI, primeBits, &t1);
    bci_sub(t1, ONE_BCI, &fullMax);

    bitCapInt minPrime;
    if (primeDict[primeBits].size()) {
        bci_copy(primeDict[primeBits][0], &minPrime);
    } else {
        bci_add(fullMin, ONE_BCI, &minPrime);
    }

    bitCapInt maxPrime;
    if (primeDict[primeBits].size()) {
        bci_copy(primeDict[primeBits][1], &maxPrime);
    } else {
        bci_add(fullMin, ONE_BCI, &maxPrime);
    }

    bitCapInt fullMinR, t2, t3;
    bci_sub(minPrime, ONE_BCI, &t1);
    bci_div(toFactor, minPrime, &t2);
    bci_sub(t2, ONE_BCI, &t3);
    bci_mul(t1, t3, &fullMinR);

    bitCapInt fullMaxR;
    bci_sub(maxPrime, ONE_BCI, &t1);
    bci_div(toFactor, maxPrime, &t2);
    bci_sub(t2, ONE_BCI, &t3);
    bci_mul(t1, t3, &fullMaxR);
#else
    // \phi(n) is Euler's totient for n. A loose lower bound is \phi(n) >= sqrt(n/2).
    bitCapInt fullMinR;
    bci_rshift(toFactor, 1U, &t1);
    floorSqrt(t1, &fullMinR);
    // A better bound is \phi(n) >= pow(n / 2, log(2)/log(3))
    // const bitCapInt fullMinR = pow(toFactor / 2, PHI_EXPONENT);

    // It can be shown that the period of this modular exponentiation can be no higher than 1
    // less than the modulus, as in https://www2.math.upenn.edu/~mlazar/math170/notes06-3.pdf.
    // Further, an upper bound on Euler's totient for composite numbers is n - sqrt(n). (See
    // https://math.stackexchange.com/questions/896920/upper-bound-for-eulers-totient-function-on-composite-numbers)
    bitCapInt fullMaxR;
    floorSqrt(toFactor, &t1);
    bci_sub(toFactor, t1, &fullMaxR);
#endif

    const DeviceContextPtr deviceContext = OCLEngine::Instance().GetDeviceContextPtr(-1);
    const cl::Context context = deviceContext->context;
    const cl::CommandQueue queue = deviceContext->queue;

    const size_t groupSize = deviceContext->GetPreferredSizeMultiple();
    const size_t groupCount = deviceContext->GetPreferredConcurrency();
    const size_t itemCount = groupSize * groupCount;
    const uint64_t batchSize = 1U << 9U;
    
    const size_t ulongArgsCount = 5U;
    BufferPtr ulongArgsBufferPtr = MakeBuffer(context, CL_MEM_READ_ONLY, sizeof(bitCapInt) * ulongArgsCount);
    uint64_t ulongArgs[5] = { batchSize, nodeId, nodeCount, (uint64_t)(log2(toFactor) >> 6), (uint64_t)IS_RSA_SEMIPRIME };
    cl_int error = queue.enqueueWriteBuffer(*ulongArgsBufferPtr, CL_TRUE, 0U, sizeof(bitCapInt) * ulongArgsCount, ulongArgs, NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }

    const size_t bciArgsCount = 3U;
    BufferPtr bciArgsBufferPtr = MakeBuffer(context, CL_MEM_READ_ONLY, sizeof(bitCapInt) * bciArgsCount);
    bitCapInt bciArgs[3] = { toFactor, fullMinR, fullMaxR };
    error = queue.enqueueWriteBuffer(*bciArgsBufferPtr, CL_TRUE, 0U, sizeof(bitCapInt) * bciArgsCount, bciArgs, NULL);
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

    bitCapInt testFactor = bci_create(1U);
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
        // std::cout << "Finished batch." << std::endl;

        queue.enqueueReadBuffer(*outputBufferPtr, CL_TRUE, 0U, sizeof(bitCapInt) * itemCount, outputArray.get(), NULL);
        for (size_t i = 0; i < itemCount; i++) {
            testFactor = outputArray.get()[i];
            BigInteger remainder;
            bci_mod(toFactor, testFactor, &remainder);
            if (bci_gt_1(testFactor) && bci_eq_0(remainder)) {
                break;
            }
        }
    }

    bci_div(toFactor, testFactor, &t1);
    std::cout << "Success: " << testFactor << " * " << t1 << std::endl;
    const double clockFactor = 1.0 / 1000.0; // Report in ms
    auto tClock = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;

    return 0;
}
