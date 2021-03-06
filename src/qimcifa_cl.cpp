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

#include "config.h"
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

#include <boost/random.hpp>

#include "big_integer.h"
#define bitCapInt BigInteger
#define bitLenInt uint32_t
#define WORD uint64_t
#define bitsInByte 8U

#define bci_create(a) bi_create(a)
#define bci_copy_ip(a, o) bi_copy_ip(&(a), o)
#define bci_set_0(a) bi_set_0(a)
#define bci_increment(l, r) bi_increment(l, r)
#define bci_add_ip(l, r) bi_add_ip(l, &(r))
#define bci_add(l, r) bi_add(&(l), &(r))
#define bci_decrement(l, r) bi_decrement(l, r)
#define bci_sub(l, r) bi_sub(&(l), &(r))
#define bci_sub_ip(l, r) bi_sub_ip(l, &(r))
#define bci_mul(l, r) bi_mul(&(l), &(r))
#define bci_mul_small(l, r) bi_mul_small(&(l), r)
#define bci_div(l, r, o) bi_div_mod(&(l), &(r), o, NULL)
#define bci_div_small(l, r, o) bi_div_mod_small(&(l), r, o, NULL)
#define bci_mod(l, r, o) bi_div_mod(&(l), &(r), NULL, o)
#define bci_mod_small(l, r, o) bi_div_mod_small(&(l), r, NULL, o)
#define bci_lshift(l, r) bi_lshift(&(l), r)
#define bci_lshift_ip(l, r) bi_lshift_ip(l, r)
#define bci_rshift(l, r) bi_rshift(&(l), r)
#define bci_rshift_ip(l, r) bi_rshift_ip(l, r)
#define bci_and(l, r) bi_and(&(l), &(r))
#define bci_eq(l, r) (bi_compare(&(l), &(r)) == 0)
#define bci_neq(l, r) (bi_compare(&(l), &(r)) != 0)
#define bci_lt(l, r) (bi_compare(&(l), &(r)) < 0)
#define bci_gt(l, r) (bi_compare(&(l), &(r)) > 0)
#define bci_geq(l, r) (bi_compare(&(l), &(r)) >= 0)
#define bci_leq(l, r) (bi_compare(&(l), &(r)) <= 0)
#define bci_and_1(l) bi_and_1(&(l))
#define bci_compare(a, b) (bi_compare(&(a), &(b)))
#define bci_eq_0(a) (bi_compare_0(&(a)) == 0)
#define bci_neq_0(a) (bi_compare_0(&(a)) != 0)
#define bci_eq_1(a) (bi_compare_1(&(a)) == 0)
#define bci_neq_1(a) (bi_compare_1(&(a)) != 0)
#define bci_gt_1(a) (bi_compare_1(&(a)) > 0)
#define bci_low64(a) ((a).bits[0])

namespace Qimcifa {

const bitCapInt ONE_BCI = bci_create(1U);

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

std::ostream& operator<<(std::ostream& os, bitCapInt b)
{
    bitCapInt t;

    // Calculate the base-10 digits, from lowest to highest.
    std::vector<std::string> digits;
    while (bci_neq_0(b)) {
        digits.push_back(std::to_string((unsigned char)(bci_low64(b) % 10U)));
        bci_copy_ip(b, &t);
        bci_div_small(t, 10, &b);
    }

    if (digits.size() == 0) {
        os << "0";
        return os;
    }

    // Reversing order, print the digits from highest to lowest.
    for (size_t i = digits.size() - 1U; i > 0; --i) {
        os << digits[i];
    }
    // Avoid the need for a signed comparison.
    os << digits[0];

    return os;
}

std::istream& operator>>(std::istream& is, bitCapInt& b)
{
    // Get the whole input string at once.
    std::string input;
    is >> input;

    // Start the output address value at 0.
    bi_set_0(&b);
    for (size_t i = 0; i < input.size(); ++i) {
        // Left shift by 1 base-10 digit.
        b = bci_mul_small(b, 10);
        // Add the next lowest base-10 digit.
        bci_increment(&b, input[i] - 48U);
    }

    return is;
}

// Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
bool isPowerOfTwo(const bitCapInt& x)
{
    bitCapInt t1 = bci_sub(x, ONE_BCI);
    bitCapInt t2 = bci_and(x, t1);

    return bci_neq_0(x) && bci_eq_0(t2);
}

bitLenInt log2(const bitCapInt& n) { return bi_log2(&n); }

// Count of distinct primes increases logarithmically, over the integers increasing from 0. The cost of additional
// trial division factors is linear. The complexity asymptote of "multiples elimination" (complement to trial
// division by primes) is O(log), with the grant of a primes table that scales linearly in query count for cost.
// However, the O(log) asymptote is FAR practically slower, (at least for now). The empirical level follows:
inline size_t pickTrialDivisionLevel(size_t qubitCount)
{
#if defined(TRIAL_DIVISION_LEVEL_OVERRIDE)
    if (TRIAL_DIVISION_LEVEL_OVERRIDE >= 0) {
        return TRIAL_DIVISION_LEVEL_OVERRIDE;
    }
#endif

#if IS_RSA_SEMIPRIME
    if (qubitCount < 56) {
        return 0;
    }

    return (qubitCount >> 3U) - 6U;
#else
    return (qubitCount >> 4U) + 34U;
#endif
}

int mainBody(bitCapInt toFactor, size_t qubitCount, size_t nodeCount, size_t nodeId,
    const std::vector<unsigned>& trialDivisionPrimes)
{
    bitCapInt t;
    auto iterClock = std::chrono::high_resolution_clock::now();
    const int TRIAL_DIVISION_LEVEL = pickTrialDivisionLevel(qubitCount);
#if IS_RSA_SEMIPRIME
    int primeIndex = TRIAL_DIVISION_LEVEL;
    unsigned currentPrime = trialDivisionPrimes[primeIndex];

    const bitLenInt primeBits = (qubitCount + 1U) >> 1U;
    bitCapInt fullMinBase = bci_create(1);
    bci_lshift_ip(&fullMinBase, primeBits - 2U);
    bci_increment(&fullMinBase, 1U);

    bitCapInt fullMaxBase = bci_create(1);
    bci_lshift_ip(&fullMaxBase, primeBits + 1U);
    bci_decrement(&fullMaxBase, 1U);
#else
    int primeIndex = 0;
    unsigned currentPrime = 2U;
    while (primeIndex <= TRIAL_DIVISION_LEVEL) {
        currentPrime = trialDivisionPrimes[primeIndex];
        unsigned t2;
        bci_mod_small(toFactor, currentPrime, &t2);
        if (!t2) {
            bci_div_small(toFactor, currentPrime, &t);
            std::cout << "Factors: " << currentPrime << " * " << t << " = " << toFactor << std::endl;
            return 0;
        }
        ++primeIndex;
    }

    // We include potential factors as low as the next odd number after the highest trial division prime.
    currentPrime += 2U;
    bitCapInt fullMinBase = bci_create(currentPrime);
    // We include potential factors as high as toFactor / nextPrime.
    bitCapInt fullMaxBase;
    bci_div_small(toFactor, currentPrime, &fullMaxBase);
#endif

    primeIndex = TRIAL_DIVISION_LEVEL;
    while (primeIndex >= 0) {
        // The truncation here is a conservative bound, but it's exact if we
        // happen to be aligned to a perfect factor of all trial division.
        currentPrime = trialDivisionPrimes[primeIndex];
        bci_div_small(fullMinBase, currentPrime, &t);
        fullMinBase = bci_mul_small(t, currentPrime);
        --primeIndex;
    }

    bitCapInt fullRange = bci_sub(fullMaxBase, fullMinBase);
    bci_increment(&fullRange, 1U);
    primeIndex = TRIAL_DIVISION_LEVEL;
    while (primeIndex >= 0) {
        // The truncation here is a conservative bound, but it's exact if we
        // happen to be aligned to a perfect factor of all trial division.
        currentPrime = trialDivisionPrimes[primeIndex];
        t = bci_mul_small(fullRange, currentPrime - 1U);
        bci_div_small(t, currentPrime, &fullRange);
        --primeIndex;
    }
    primeIndex = TRIAL_DIVISION_LEVEL;

    bitCapInt nodeRange;
    bci_copy_ip(fullRange, &t);
    bci_increment(&t, nodeCount - 1U);
    bci_div_small(t, nodeCount, &nodeRange);

    bitCapInt nodeMin = bci_mul_small(nodeRange, nodeId);
    bci_add_ip(&nodeMin, fullMinBase);

    const DeviceContextPtr deviceContext = OCLEngine::Instance().GetDeviceContextPtr(-1);
    const cl::Context context = deviceContext->context;
    const cl::CommandQueue queue = deviceContext->queue;

    const size_t groupSize = deviceContext->GetPreferredSizeMultiple();
    const size_t groupCount = deviceContext->GetPreferredConcurrency();
    const size_t itemCount = groupSize * groupCount;

    BufferPtr primesBufferPtr = MakeBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned) * primeIndex);
    cl_int error = queue.enqueueWriteBuffer(*primesBufferPtr, CL_TRUE, 0U, sizeof(unsigned) * primeIndex, &(trialDivisionPrimes[0]), NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }

    const size_t uArgsCount = 4U;
    const unsigned batchSize = 1U << 9U;
    BufferPtr uArgsBufferPtr = MakeBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned) * uArgsCount);
    unsigned uArgs[4] = { batchSize, (unsigned)nodeId, (unsigned)nodeCount, (unsigned)primeIndex };
    error = queue.enqueueWriteBuffer(*uArgsBufferPtr, CL_TRUE, 0U, sizeof(unsigned) * uArgsCount, uArgs, NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }

    const size_t bciArgsCount = 4U;
    BufferPtr bciArgsBufferPtr = MakeBuffer(context, CL_MEM_READ_ONLY, sizeof(bitCapInt) * bciArgsCount);
    bitCapInt bciArgs[4] = { toFactor, nodeMin, nodeRange, fullMinBase };
    error = queue.enqueueWriteBuffer(*bciArgsBufferPtr, CL_TRUE, 0U, sizeof(bitCapInt) * bciArgsCount, bciArgs, NULL);
    if (error != CL_SUCCESS) {
        throw std::runtime_error("Failed to enqueue buffer write, error code: " + std::to_string(error));
    }

    typedef std::uniform_int_distribution<uint64_t> rand_dist;
    std::random_device rand_dev;
    boost::taus88 rand_gen(rand_dev());

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
#if IS_RSA_SEMIPRIME
    OCLDeviceCall ocl = deviceContext->Reserve(OCL_API_QIMCIFA_RSA_BATCH);
#else
    OCLDeviceCall ocl = deviceContext->Reserve(OCL_API_QIMCIFA_BATCH);
#endif

    ocl.call.setArg(0, *rngSeedBufferPtr);
    ocl.call.setArg(1, *primesBufferPtr);
    ocl.call.setArg(2, *uArgsBufferPtr);
    ocl.call.setArg(3, *bciArgsBufferPtr);
    ocl.call.setArg(4, *outputBufferPtr);

    bitCapInt testFactor = bci_create(1U);
    while (!bci_gt_1(testFactor)) {
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
            if (bci_gt_1(testFactor)) {
                bitCapInt rmndr;
                bci_mod(toFactor, testFactor, &rmndr);
                if (bci_eq_0(rmndr)) {
                    break;
                }
            }
        }
    }

    bitCapInt f2;
    bci_div(toFactor, testFactor, &f2);
    std::cout << "Success: " << testFactor << " * " << f2 << " = " << toFactor << std::endl;
    const double clockFactor = 1.0 / 1000.0; // Report in ms
    auto tClock = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;

    return 0;
}
} // namespace Qimcifa

using namespace Qimcifa;

int main()
{
    // First 1000 primes
    // Source: https://gist.github.com/cblanc/46ebbba6f42f61e60666#file-gistfile1-txt
    const std::vector<unsigned> trialDivisionPrimes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
        61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
        181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
        311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439,
        443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587,
        593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727,
        733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877,
        881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021,
        1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
        1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283,
        1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429,
        1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549,
        1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667,
        1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
        1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973,
        1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089,
        2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243,
        2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377,
        2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
        2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677,
        2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791,
        2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939,
        2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089,
        3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
        3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391,
        3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541,
        3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677,
        3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833,
        3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
        4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129,
        4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271,
        4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447,
        4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597,
        4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
        4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931,
        4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051,
        5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227,
        5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399,
        5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
        5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683,
        5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839,
        5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007,
        6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151,
        6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
        6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451,
        6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637,
        6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791,
        6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949,
        6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
        7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253,
        7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477,
        7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589,
        7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741,
        7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 };

    // Print primes table by index:
    // for (size_t i = 0; i < trialDivisionPrimes.size(); ++i) {
    //     std::cout << i << ": " << trialDivisionPrimes[i] << ", ";
    // }

    bitCapInt toFactor;
    size_t nodeCount = 1U;
    size_t nodeId = 0U;

    std::cout << "Number to factor: ";
    std::cin >> toFactor;

    uint32_t qubitCount = log2(toFactor);
    // Source: https://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
    if (!isPowerOfTwo(toFactor)) {
        qubitCount++;
    }
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
            std::cout << "Which node is this? (0-" << (nodeCount - 1U) << "): ";
            std::cin >> nodeId;
            if (nodeId >= nodeCount) {
                std::cout << "Invalid node ID choice!" << std::endl;
            }
        } while (nodeId >= nodeCount);
    }
#endif

    return mainBody(toFactor, qubitCount, nodeCount, nodeId, trialDivisionPrimes);
}
