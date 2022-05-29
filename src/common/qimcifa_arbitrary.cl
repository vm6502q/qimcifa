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

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
void uipow(BigInteger b, BigInteger e, BigInteger* result)
{
    bi_set_0(result);
    bi_increment(result, 1);
    for (;;) {
        if (bi_and_1(&b)) {
            BigInteger t;
            bi_copy(result, &t);
            bi_mul(&t, &b, result);
        }
        BigInteger t;
        bi_copy(&e, &t);
        bi_rshift(&t, 1U, &e);
        if (bi_compare_0(&e) == 0) {
            break;
        }
        bi_copy(&b, &t);
        bi_mul(&t, &t, &b);
    }
}


void gcd(BigInteger n1, BigInteger n2, BigInteger* result)
{
    while (bi_compare_0(&n2) != 0) {
        BigInteger t1, t2;
        bi_copy(&n1, &t1);
        bi_copy(&n2, &t2);
        bi_copy(&n2, &n1);
        bi_div_mod(&t1, &t2, 0, &n2);
    }
    bi_copy(&n1, result);
}

void kernel qimcifa_batch(global ulong* ulongArgs, global BigInteger* BigIntegerArgs, global ulong* rngSeeds, global BigInteger* outputs)
{
    const ulong thread = get_global_id(0);
    const ulong threadCount = get_global_size(0);

    const ulong batchSize = ulongArgs[0];
    const ulong nodeId = ulongArgs[1];
    const ulong nodeCount = ulongArgs[2];
    const ulong byteCount = ulongArgs[3];
    // This can become a bit flag register, in the future:
    const ulong isRsaSemiprime = ulongArgs[4];

    const BigInteger toFactor = bi_load((global ulong*)BigIntegerArgs);
    const BigInteger fullMinR = bi_load((global ulong*)(BigIntegerArgs + 1));
    const BigInteger fullMaxR = bi_load((global ulong*)(BigIntegerArgs + 2));

    const ulong4 rngLoad = vload4(thread, rngSeeds);
    kiss09_state rngState;
    rngState.x = rngLoad.x;
    rngState.c = rngLoad.y;
    rngState.y = rngLoad.z;
    rngState.z = rngLoad.w;

    const BigInteger BIG_INT_1 = bi_create(1);
    const BigInteger BIG_INT_2 = bi_create(2);
    const BigInteger BIG_INT_3 = bi_create(3);

    BigInteger t1, t2, t3;

    BigInteger fullRange;
    bi_add(&fullMaxR, &BIG_INT_1, &t1);
    bi_sub(&t1, &fullMinR, &fullRange);

    BigInteger nodeRange, nodeCountBci = bi_create(nodeCount);
    bi_div_mod(&fullRange, &nodeCountBci, &nodeRange, 0);

    BigInteger nodeMin, nodeIdBci = bi_create(nodeId);
    bi_mul(&nodeRange, &nodeIdBci, &t1);
    bi_add(&fullMinR, &t1, &nodeMin);

    BigInteger nodeMax;
    if ((nodeId + 1U) == nodeCount) {
        bi_copy(&fullMaxR, &nodeMax);
    } else {
        BigInteger nodeIdPlus1Bci = bi_create(nodeId + 1U);
        bi_mul(&nodeRange, &nodeIdPlus1Bci, &t1);
        bi_add(&fullMinR, &t1, &t2);
        bi_sub(&t2, &BIG_INT_1, &nodeMax);
    }

    BigInteger threadRange, threadCountBci = bi_create(nodeId);
    bi_add(&nodeMax, &BIG_INT_1, &t1);
    bi_sub(&t1, &nodeMin, &t2);
    bi_div_mod(&t2, &threadCountBci, &threadRange, 0);

    BigInteger rMin, threadBci = bi_create(thread);
    bi_mul(&threadRange, &threadBci, &t1);
    bi_add(&nodeMin, &t1, &rMin);

    BigInteger rMax;
    if ((thread + 1U) == threadCount) {
        bi_copy(&nodeMax, &rMax);
    } else {
        BigInteger threadPlus1Bci = bi_create(thread + 1U);
        bi_mul(&threadRange, &threadPlus1Bci, &t1);
        bi_add(&nodeMin, &t1, &t2);
        bi_sub(&t2, &BIG_INT_1, &rMax);
    }

    BigInteger rRange;
    bi_copy(&rMax, &t1);
    bi_increment(&t1, 1);
    bi_sub(&t1, &rMin, &rRange);

    BigInteger toFactorMin3;
    bi_sub(&toFactor, &BIG_INT_3, &toFactorMin3);

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        BigInteger base;
        for (size_t i = 0; i < byteCount; i++) {
            base.bits[i] = kiss09_ulong(rngState);
        }
        bi_div_mod(&base, &toFactorMin3, 0, &t1);
        bi_add(&t1, &BIG_INT_2, &base);

        gcd(toFactor, base, &t1);
        if (bi_compare(&t1, &BIG_INT_1) != 0) {
            for (int i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
                ((global ulong*)outputs)[BIG_INTEGER_WORD_SIZE * thread + i] = t1.bits[i];
            }

            return;
        }

        // This would be where we perform the quantum period finding algorithm.
        // However, we don't have a quantum computer!
        // Instead, we "throw dice" for a guess to the output of the quantum subroutine.
        // This guess will usually be wrong, at least for semi-prime inputs.
        // If we try many times, though, this can be a practically valuable factoring method.

        // y is meant to be close to some number c * qubitPower / r, where r is the period.
        // c is a positive integer or 0, and we don't want the 0 case.
        // y is truncated by the number of qubits in the register, at most.
        // The maximum value of c before truncation is no higher than r.

        // The period of ((base ^ x) MOD toFactor) can't be smaller than log_base(toFactor).
        // (Also, toFactor is definitely NOT an exact multiple of base.)
        // const BigInteger logBaseToFactor = (BigInteger)intLog(base, toFactor) + 1;
        // Euler's Theorem tells us, if gcd(a, n) = 1, then a^\phi(n) = 1 MOD n,
        // where \phi(n) is Euler's totient for n.
        // const BigInteger fullMinR = (minPhi < logBaseToFactor) ? logBaseToFactor : minPhi;

        // c is basically a harmonic degeneracy factor, and there might be no value in testing
        // any case except c = 1, without loss of generality.

        // This sets a nonuniform distribution on our y values to test.
        // y values are close to qubitPower / rGuess, and we midpoint round.

        // However, results are better with uniformity over r, rather than y.

        // So, we guess r, between fullMinR and fullMaxR.
        // Choose a base at random, >1 and <toFactor.
        BigInteger r;
        for (size_t i = 0; i < byteCount; i++) {
            r.bits[i] = kiss09_ulong(rngState);
        }
        bi_div_mod(&r, &rRange, 0, &t1);
        bi_add(&t1, &rMin, &r);

        // Since our output is r rather than y, we can skip the continued fractions step.
        if (bi_and_1(&r)) {
            bi_copy(&r, &t1);
        } else {
            bi_rshift(&r, 1U, &t1);
        }

        BigInteger rGuess;
        if (isRsaSemiprime) {
            bi_copy(&t1, &rGuess);
        } else {
            bi_copy(&r, &rGuess);
        }

        // As a "classical" optimization, since \phi(toFactor) and factor bounds overlap,
        // we first check if our guess for r is already a factor.
        bi_div_mod(&toFactor, &rGuess, 0, &t2);
        if ((bi_compare(&rGuess, &BIG_INT_1) > 0) && (bi_compare_0(&t2) == 0)) {
            for (int i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
                ((global ulong*)outputs)[BIG_INTEGER_WORD_SIZE * thread + i] = rGuess.bits[i];
            }

            return;
        }
        
        BigInteger f1, f2, fmul;
        uipow(base, t1, &t2);
        bi_div_mod(&t2, &toFactor, 0, &t1);

        bi_add(&t1, &BIG_INT_1, &t2);
        gcd(t2, toFactor, &f1);

        bi_sub(&t1, &BIG_INT_1, &t2);
        gcd(t2, toFactor, &f2);

        bi_mul(&f1, &f2, &fmul);

        bi_div_mod(&toFactor, &fmul, 0, &t1);

        while ((bi_compare(&fmul, &BIG_INT_1) > 0) && (bi_compare(&fmul, &toFactor) != 0) && (bi_compare_0(&t1) == 0)) {
            bi_copy(&f1, &fmul);
            bi_mul(&fmul, &f2, &f1);

            bi_mul(&fmul, &f2, &t1);
            bi_div_mod(&toFactor, &t1, &f2, 0);

            bi_mul(&f1, &f2, &fmul);

            bi_div_mod(&toFactor, &fmul, 0, &t1);
        }
        if ((bi_compare(&fmul, &BIG_INT_1) > 0) && (bi_compare(&fmul, &toFactor) == 0) && (bi_compare(&f1, &BIG_INT_1) > 0) && (bi_compare(&f2, &BIG_INT_1) > 0)) {
            for (int i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
                ((global ulong*)outputs)[BIG_INTEGER_WORD_SIZE * thread + i] = f1.bits[i];
            }

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
