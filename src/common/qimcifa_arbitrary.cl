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

#define bitLenInt uint
#define bitCapInt BigInteger
#define BIT_AND_1(x) (x.bits[0] & 1)

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
bitCapInt uipow(bitCapInt b, bitCapInt e)
{
    const BigInteger BIGINT_0 = bi_create(0);
    bitCapInt result = bi_create(1);
    for (;;) {
        if (BIT_AND_1(b)) {
            result = bi_mul(result, b);
        }
        bi_rshift(e, 1);
        if (bi_compare(e, BIGINT_0) == 0) {
            break;
        }
        b = bi_mul(b, b);
    }

    return result;
}

bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
    const BigInteger BIGINT_0 = bi_create(0);
    while (bi_compare(n2, BIGINT_0) != 0) {
        const bitCapInt temp = bi_copy(n2);
        n2 = bi_mod(n1, n2);
        n1 = temp;
    }
    return n1;
}

void kernel qimcifa_batch(global ulong* ulongArgs, global ulong* bitCapIntArgs, global ulong* rngSeeds, global ulong* outputs)
{
    const bitLenInt wordSize = 64;

    const ulong thread = get_global_id(0);
    const ulong threadCount = get_global_size(0);
    const ulong batchSize = ulongArgs[0];
    const ulong nodeId = bitCapIntArgs[1];
    const ulong nodeCount = bitCapIntArgs[2];
    // This can become a bit flag register, in the future:
    const ulong isRsaSemiprime = bitCapIntArgs[3];
    const bitCapInt toFactor = bi_load(bitCapIntArgs);
    const bitCapInt fullMinR = bi_load(bitCapIntArgs + 1);
    const bitCapInt fullMaxR = bi_load(bitCapIntArgs + 2);

    const ulong4 rngLoad = vload4(thread, rngSeeds);
    kiss09_state rngState;
    rngState.x = rngLoad.x;
    rngState.c = rngLoad.y;
    rngState.y = rngLoad.z;
    rngState.z = rngLoad.w;

    const bitCapInt BIGINT_1 = bi_create(1);
    const bitCapInt BIGINT_2 = bi_create(2);
    const bitCapInt BIGINT_3 = bi_create(3);

    const bitCapInt fullRange = bi_subtract(bi_add(fullMaxR, BIGINT_1), fullMinR);
    const bitCapInt nodeRange = bi_div(fullRange, bi_create(nodeCount));
    const bitCapInt nodeMin = bi_add(fullMinR, bi_mul(nodeRange, bi_create(nodeId)));
    const bitCapInt nodeMax = ((nodeId + 1) == nodeCount)
        ? fullMaxR
        : bi_add(fullMinR, bi_subtract(bi_mul(nodeRange, bi_create(nodeId + 1)), BIGINT_1));
    const bitCapInt threadRange = bi_div(bi_subtract(bi_add(nodeMax, BIGINT_1), nodeMin), bi_create(threadCount));
    const bitCapInt rMin = bi_add(nodeMin, bi_mul(threadRange, bi_create(thread)));
    const bitCapInt rMax = ((thread + 1) == threadCount)
        ? nodeMax
        : bi_add(nodeMin, bi_subtract(bi_mul(threadRange, bi_create(thread + 1)), BIGINT_1));
    const bitCapInt rRange = bi_subtract(bi_add(rMax, BIGINT_1), rMin);

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        bitCapInt base = bi_create(kiss09_ulong(rngState));
        for (size_t i = 1; i < BIG_INTEGER_WORD_SIZE; i++) {
            bi_lshift_word(base, 1);
            base.bits[0] = kiss09_ulong(rngState);
        }
        base = bi_add(bi_mod(base, bi_subtract(toFactor, BIGINT_3)), BIGINT_2);

        const bitCapInt testFactor = gcd(toFactor, base);
        if (bi_compare(testFactor, BIGINT_1) != 0) {
            const ulong threadOffset = thread * BIG_INTEGER_WORD_SIZE;
            for (int i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
                outputs[threadOffset + i] = testFactor.bits[i];
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
        // const bitCapInt logBaseToFactor = (bitCapInt)intLog(base, toFactor) + 1;
        // Euler's Theorem tells us, if gcd(a, n) = 1, then a^\phi(n) = 1 MOD n,
        // where \phi(n) is Euler's totient for n.
        // const bitCapInt fullMinR = (minPhi < logBaseToFactor) ? logBaseToFactor : minPhi;

        // c is basically a harmonic degeneracy factor, and there might be no value in testing
        // any case except c = 1, without loss of generality.

        // This sets a nonuniform distribution on our y values to test.
        // y values are close to qubitPower / rGuess, and we midpoint round.

        // However, results are better with uniformity over r, rather than y.

        // So, we guess r, between fullMinR and fullMaxR.
        // Choose a base at random, >1 and <toFactor.
        bitCapInt r = bi_create(kiss09_ulong(rngState));
        for (size_t i = 1; i < BIG_INTEGER_WORD_SIZE; i++) {
            bi_lshift_word(r, 1);
            r.bits[0] = kiss09_ulong(rngState);
        }
        r = bi_add(bi_mod(r, rRange), rMin);

        // Since our output is r rather than y, we can skip the continued fractions step.
        const bitCapInt p = BIT_AND_1(r) ? r : bi_rshift(r, 1);

        const bitCapInt rGuess = (isRsaSemiprime & 1) ? p : r;

        // As a "classical" optimization, since \phi(toFactor) and factor bounds overlap,
        // we first check if our guess for r is already a factor.
        if ((bi_compare(rGuess, BIGINT_1) != 0)
            && (bi_compare(bi_mul(bi_div(toFactor, rGuess), rGuess), toFactor) == 0)) {
            const ulong threadOffset = thread * BIG_INTEGER_WORD_SIZE;
            for (int i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
                outputs[threadOffset + i] = rGuess.bits[i];
            }

            return;
        }

        const bitCapInt apowrhalf = bi_mod(uipow(base, p), toFactor);
        bitCapInt f1 = gcd(bi_add(apowrhalf, BIGINT_1), toFactor);
        bitCapInt f2 = gcd(bi_subtract(apowrhalf, BIGINT_1), toFactor);
        bitCapInt fmul = bi_mul(f1, f2);
        while ((bi_compare(fmul, BIGINT_1) != 0)
            && (bi_compare(fmul, toFactor) != 0)
            && (bi_compare(bi_mul(bi_div(toFactor, fmul), fmul), toFactor) == 0)) {
            fmul = f1;
            f1 = bi_mul(fmul, f2);
            f2 = bi_div(toFactor, bi_mul(fmul, f2));
            fmul = bi_mul(f1, f2);
        }
        if ((bi_compare(fmul, BIGINT_1) != 0)
            && (bi_compare(fmul, toFactor) == 0)
            && (bi_compare(f1, BIGINT_1) != 0)
            && (bi_compare(f2, BIGINT_1) != 0)) {
            const ulong threadOffset = thread * BIG_INTEGER_WORD_SIZE;
            for (int i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
                outputs[threadOffset + i] = f1.bits[i];
            }

            return;
        }
    }

    // We need to update the global seeds, but only if we didn't succeeed.
    const ulong threadOffset = thread * 4;
    rngSeeds[threadOffset + 0] = rngState.x;
    rngSeeds[threadOffset + 1] = rngState.c;
    rngSeeds[threadOffset + 2] = rngState.y;
    rngSeeds[threadOffset + 3] = rngState.z;
}
