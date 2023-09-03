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

// #include "kiss09.cl"
// #include "big_integer.cl"

#define bitCapInt BigInteger

bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
    while (bi_compare_0(&n2) != 0) {
        BigInteger t1 = bi_copy(&n1);
        BigInteger t2 = bi_copy(&n2);
        bi_copy_ip(&n2, &n1);
        bi_div_mod(&t1, &t2, 0, &n2);
    }
    return n1;
}

bool isPowerOfTwo(bitCapInt x)
{
    bitCapInt t = bi_copy(&x);
    bi_decrement(&t, 1);
    bi_and_ip(&t, &x);

    return (bi_compare_0(&x) != 0) && (bi_compare_0(&t) == 0);
}

void kernel qimcifa_batch(global ulong* rngSeeds, constant unsigned* trialDivisionPrimes, constant unsigned* unsignedArgs, constant bitCapInt* bitCapIntArgs, global bitCapInt* outputs)
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
    const bitCapInt nodeRange = bitCapIntArgs[2];
    const bitCapInt fullMinBase = bitCapIntArgs[3];

    const ulong4 rngLoad = vload4(thread, rngSeeds);
    kiss09_state rngState;
    rngState.x = rngLoad.x;
    rngState.c = rngLoad.y;
    rngState.y = rngLoad.z;
    rngState.z = rngLoad.w;

    bitCapInt threadRange, t = bi_copy(&nodeRange);
    bi_increment(&t, threadCount - 1U);
    bi_div_mod_small(&t, threadCount, &threadRange, 0);

    bitCapInt threadMin = bi_mul_small(&threadRange, thread);
    bi_add_ip(&threadMin, &nodeMin);
    bitCapInt threadMax = bi_add(&threadMin, &threadRange);

    const size_t rngBitCount = bi_log2(&threadRange) + (isPowerOfTwo(threadRange) ? 0 : 1);
    const size_t rngWordCount = (rngBitCount < 64) ? 1 : ((rngBitCount >> 6) + 1);

    // Align the lower limit to a multiple of ALL trial division factors.
    unsigned privPrimes[64];
    for (int i = primesLength - 1; i >= 0; --i) {
        privPrimes[i] = trialDivisionPrimes[i] - 1;
    }

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        for (size_t i = 0; i < rngWordCount; i++) {
            t.bits[i] = kiss09_ulong(rngState);
        }
        for (size_t i = rngWordCount; i < BIG_INTEGER_WORD_SIZE; i++) {
            t.bits[i] = 0;
        }
        bitCapInt base;
        bi_div_mod(&t, &threadRange, 0, &base);

        for (size_t i = primesLength - 1; i > 0; --i) {
            // Make this NOT a multiple of prime "p", by adding it to itself divided by (p - 1), + 1.
            bi_div_mod_small(&base, privPrimes[i], &t, 0);
            bi_add_ip(&base, &t);
            bi_increment(&base, 1);
        }
        // Then, make this odd, when added to the minimum.
        bi_lshift_ip(&base, 1);
        base.bits[0] |= 1;
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

void kernel qimcifa_rsa_batch(global ulong* rngSeeds, constant unsigned* trialDivisionPrimes, constant unsigned* unsignedArgs, constant bitCapInt* bitCapIntArgs, global bitCapInt* outputs)
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
    const bitCapInt nodeRange = bitCapIntArgs[2];
    const bitCapInt fullMinBase = bitCapIntArgs[3];

    const ulong4 rngLoad = vload4(thread, rngSeeds);
    kiss09_state rngState;
    rngState.x = rngLoad.x;
    rngState.c = rngLoad.y;
    rngState.y = rngLoad.z;
    rngState.z = rngLoad.w;

    bitCapInt threadRange, t = bi_copy(&nodeRange);
    bi_increment(&t, threadCount - 1U);
    bi_div_mod_small(&t, threadCount, &threadRange, 0);

    bitCapInt threadMin = bi_mul_small(&threadRange, thread);
    bi_add_ip(&threadMin, &nodeMin);
    bitCapInt threadMax = bi_add(&threadMin, &threadRange);

    const size_t rngBitCount = bi_log2(&threadRange) + (isPowerOfTwo(threadRange) ? 0 : 1);
    const size_t rngWordCount = (rngBitCount < 64) ? 1 : ((rngBitCount >> 6) + 1);

    // Align the lower limit to a multiple of ALL trial division factors.
    unsigned privPrimes[64];
    for (int i = primesLength - 1; i >= 0; --i) {
        privPrimes[i] = trialDivisionPrimes[i] - 1;
    }

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        for (size_t i = 0; i < rngWordCount; i++) {
            t.bits[i] = kiss09_ulong(rngState);
        }
        for (size_t i = rngWordCount; i < BIG_INTEGER_WORD_SIZE; i++) {
            t.bits[i] = 0;
        }
        bitCapInt base;
        bitCapInt quotient;
        bi_div_mod(&t, &threadRange, &quotient, &base);

        for (size_t i = primesLength - 1; i > 0; --i) {
            // Make this NOT a multiple of prime "p", by adding it to itself divided by (p - 1), + 1.
            bi_div_mod_small(&base, privPrimes[i], &t, 0);
            bi_add_ip(&base, &t);
            bi_increment(&base, 1);
        }
        // Then, make this odd, when added to the minimum.
        bi_lshift_ip(&base, 1);
        base.bits[0] |= 1;
        bi_add_ip(&base, &threadMin);

        bitCapInt n;
        bi_div_mod(&toFactor, &base, &quotient, &n);
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
