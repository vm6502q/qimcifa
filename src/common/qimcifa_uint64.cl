////////////////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qimcifa contributors 2017-2022. All rights reserved.
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

// TODO: Consider this pseudocode, copied from qimcifa.cpp.

// #include "kiss09.cl"
// #include "big_integer.cl"

#define bitCapInt BigInteger

bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
    bitCapInt t1, t2;
    while (bi_compare_0(&n2) != 0) {
        bi_copy(&n1, &t1);
        bi_copy(&n2, &t2);
        bi_copy(&n2, &n1);
        bi_div_mod(&t1, &t2, NULL, &n2);
    }
    return n1;
}

bool isPowerOfTwo(bitCapInt x)
{
    bitCapInt t;
    bi_copy(&x, &t);
    bi_decrement(&t, 1);
    bi_and_ip(&t, &x);

    return (bi_compare_0(&x) != 0) && (bi_compare_0(&t) == 0);
}

void kernel qimcifa_batch(global ulong* rngSeeds, global unsigned* trialDivisionPrimes, global unsigned* unsignedArgs, global bitCapInt* bitCapIntArgs, global bitCapInt* outputs)
{
    //const unsigned wordSize = 64;
    const unsigned thread = get_global_id(0);
    const unsigned threadCount = get_global_size(0);
    const unsigned batchSize = unsignedArgs[0];
    const unsigned nodeId = unsignedArgs[1];
    const unsigned nodeCount = unsignedArgs[2];
    const unsigned primesLength = unsignedArgs[3];
    const bitCapInt toFactor = bitCapIntArgs[0];
    const bitCapInt nodeMin = bitCapIntArgs[1];
    const bitCapInt nodeMax = bitCapIntArgs[2];
    bitCapInt t;

    const ulong4 rngLoad = vload4(thread, rngSeeds);
    kiss09_state rngState;
    rngState.x = rngLoad.x;
    rngState.c = rngLoad.y;
    rngState.y = rngLoad.z;
    rngState.z = rngLoad.w;

    bitCapInt threadRange;
    t = bi_sub(&nodeMax, &nodeMin);
    bi_increment(&t, threadCount - 1U);
    bi_div_mod_small(&t, threadCount, &threadRange, NULL);

    bitCapInt threadMin = bi_mul_small(&threadRange, thread);
    bi_add_ip(&threadMin, &nodeMin);
    threadMin.bits[0] |= 1U;
    bitCapInt threadMax = bi_add(&threadMin, &threadRange);

    const size_t rngWordCount = (bi_log2(&threadRange) >> 6) + (isPowerOfTwo(threadRange) ? 0 : 1);

    // Align the lower limit to a multiple of ALL trial division factors.
    unsigned privPrimes[64];
    for (unsigned i = 0; i < primesLength; ++i) {
        unsigned currentPrime = trialDivisionPrimes[i];
        privPrimes[i] = currentPrime;
        bi_div_mod_small(&threadMin, currentPrime, &t, NULL);
        threadMin = bi_mul_small(&t, currentPrime);
    }
    threadMin.bits[0] |= 1;
    bi_increment(&threadMin, 2);

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        for (size_t i = 0; i < rngWordCount; i++) {
            t.bits[i] = kiss09_ulong(rngState);
        }
        bitCapInt base;
        bi_div_mod(&t, &threadRange, NULL, &base);

        for (size_t i = primesLength - 1; i > 2; --i) {
            // Make this NOT a multiple of prime "p", by adding it to itself divided by (p - 1), + 1.
            bi_div_mod_small(&base, privPrimes[i] - 1, &t, NULL);
            bi_add_ip(&base, &t);
            bi_increment(&base, 1);
        }

        // Make this NOT a multiple of 5, by adding it to itself divided by 4, + 1.
        t = bi_rshift(&base, 2);
        bi_add_ip(&base, &t);
        bi_increment(&base, 1);

        // We combine the 2 and 3 multiple removal steps.
        // Make this NOT a multiple of 3, by adding it to itself divided by 2, + 1.
        // Then, make this odd, when added to the minimum.
        t = bi_lshift(&base, 1);
        base.bits[0] &= ~1;
        bi_add_ip(&base, &t);
        bi_add_ip(&base, &threadMin);

        bitCapInt n = gcd(base, toFactor);
        if (bi_compare_1(&n) != 0) {
            outputs[thread] = n;

            return;
        }
    }

    // We need to update the global seeds, but only if we didn't succeeed.
    const size_t threadOffset = thread * 4;
    rngSeeds[threadOffset + 0] = rngState.x;
    rngSeeds[threadOffset + 1] = rngState.c;
    rngSeeds[threadOffset + 2] = rngState.y;
    rngSeeds[threadOffset + 3] = rngState.z;
}

void kernel qimcifa_rsa_batch(global ulong* rngSeeds, global unsigned* trialDivisionPrimes, global unsigned* unsignedArgs, global bitCapInt* bitCapIntArgs, global bitCapInt* outputs)
{
    //const unsigned wordSize = 64;
    const unsigned thread = get_global_id(0);
    const unsigned threadCount = get_global_size(0);
    const unsigned batchSize = unsignedArgs[0];
    const unsigned nodeId = unsignedArgs[1];
    const unsigned nodeCount = unsignedArgs[2];
    const unsigned primesLength = unsignedArgs[3];
    const bitCapInt toFactor = bitCapIntArgs[0];
    const bitCapInt nodeMin = bitCapIntArgs[1];
    const bitCapInt nodeMax = bitCapIntArgs[2];
    bitCapInt t;

    const ulong4 rngLoad = vload4(thread, rngSeeds);
    kiss09_state rngState;
    rngState.x = rngLoad.x;
    rngState.c = rngLoad.y;
    rngState.y = rngLoad.z;
    rngState.z = rngLoad.w;

    bitCapInt threadRange;
    t = bi_sub(&nodeMax, &nodeMin);
    bi_increment(&t, threadCount - 1U);
    bi_div_mod_small(&t, threadCount, &threadRange, NULL);

    bitCapInt threadMin = bi_mul_small(&threadRange, thread);
    bi_add_ip(&threadMin, &nodeMin);
    threadMin.bits[0] |= 1U;
    bitCapInt threadMax = bi_add(&threadMin, &threadRange);

    const size_t rngWordCount = (bi_log2(&threadRange) >> 6) + (isPowerOfTwo(threadRange) ? 0 : 1);

    // Align the lower limit to a multiple of ALL trial division factors.
    unsigned privPrimes[64];
    for (unsigned i = 0; i < primesLength; ++i) {
        unsigned currentPrime = trialDivisionPrimes[i];
        privPrimes[i] = currentPrime;
        bi_div_mod_small(&threadMin, currentPrime, &t, NULL);
        threadMin = bi_mul_small(&t, currentPrime);
    }
    threadMin.bits[0] |= 1;
    bi_increment(&threadMin, 2);

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        for (size_t i = 0; i < rngWordCount; i++) {
            t.bits[i] = kiss09_ulong(rngState);
        }
        bitCapInt base;
        bi_div_mod(&t, &threadRange, NULL, &base);

        for (size_t i = primesLength - 1; i > 2; --i) {
            // Make this NOT a multiple of prime "p", by adding it to itself divided by (p - 1), + 1.
            bi_div_mod_small(&base, privPrimes[i] - 1, &t, NULL);
            bi_add_ip(&base, &t);
            bi_increment(&base, 1);
        }

        // Make this NOT a multiple of 5, by adding it to itself divided by 4, + 1.
        t = bi_rshift(&base, 2);
        bi_add_ip(&base, &t);
        bi_increment(&base, 1);

        // We combine the 2 and 3 multiple removal steps.
        // Make this NOT a multiple of 3, by adding it to itself divided by 2, + 1.
        // Then, make this odd, when added to the minimum.
        t = bi_lshift(&base, 1);
        base.bits[0] &= ~1;
        bi_add_ip(&base, &t);
        bi_add_ip(&base, &threadMin);

        bitCapInt n;
        bi_div_mod(&toFactor, &base, NULL, &n);
        if (bi_compare_0(&n) == 0) {
            outputs[thread] = base;

            return;
        }
    }

    // We need to update the global seeds, but only if we didn't succeeed.
    const size_t threadOffset = thread * 4;
    rngSeeds[threadOffset + 0] = rngState.x;
    rngSeeds[threadOffset + 1] = rngState.c;
    rngSeeds[threadOffset + 2] = rngState.y;
    rngSeeds[threadOffset + 3] = rngState.z;
}
