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
#define BIT_AND_1(x) (x.data.bits[0] & 1)

// Source:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int#answer-101613
bitCapInt uipow(bitCapInt b, bitCapInt e)
{
    const BigInteger BIGINT_0 = big_integer_create(0);
    bitCapInt result = big_integer_create(1);
    for (;;) {
        if (BIT_AND_1(b)) {
            result = big_integer_mul(result, b);
        }
        big_integer_rshift(e, 1);
        if (big_integer_compare(e, BIGINT_0) == 0) {
            break;
        }
        b = big_integer_mul(b, b);
    }

    return result;
}

bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
    const BigInteger BIGINT_0 = big_integer_create(0);
    while (big_integer_compare(n2, BIGINT_0) != 0) {
        const bitCapInt temp = big_integer_copy(n2);
        n2 = big_integer_mod(n1, n2);
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
    const bitCapInt isRsaSemiprime = bitCapIntArgs[3];
    const bitCapInt toFactor = bitCapIntArgs[0];
    const bitCapInt fullMinR = bitCapIntArgs[1];
    const bitCapInt fullMaxR = bitCapIntArgs[2];

    const ulong4 rngLoad = vload4(thread, rngSeeds);
    kiss09_state rngState;
    rngState.x = rngLoad.x;
    rngState.c = rngLoad.y;
    rngState.y = rngLoad.z;
    rngState.z = rngLoad.w;

    const bitCapInt BIGINT_1 = big_integer_create(1);
    const bitCapInt BIGINT_2 = big_integer_create(2);
    const bitCapInt BIGINT_3 = big_integer_create(3);

    const bitCapInt fullRange = big_integer_subtract(big_integer_add(fullMaxR, BIGINT_1), fullMinR);
    const bitCapInt nodeRange = big_integer_div(fullRange, big_integer_create(nodeCount));
    const bitCapInt nodeMin = big_integer_add(fullMinR, big_integer_mul(nodeRange, big_integer_create(nodeId)));
    const bitCapInt nodeMax = ((nodeId + 1) == nodeCount)
        ? fullMaxR
        : big_integer_add(fullMinR, big_integer_subtract(big_integer_mul(nodeRange, big_integer_create(nodeId + 1)), BIGINT_1));
    const bitCapInt threadRange = big_integer_div(big_integer_subtract(big_integer_add(nodeMax, BIGINT_1), nodeMin), big_integer_create(threadCount));
    const bitCapInt rMin = big_integer_add(nodeMin, big_integer_mul(threadRange, big_integer_create(thread)));
    const bitCapInt rMax = ((thread + 1) == threadCount)
        ? nodeMax
        : big_integer_add(nodeMin, big_integer_subtract(big_integer_mul(threadRange, big_integer_create(thread + 1)), BIGINT_1));
    const bitCapInt rRange = big_integer_subtract(big_integer_add(rMax, BIGINT_1), rMin);

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        bitCapInt base = big_integer_create(kiss09_ulong(rngState));
        for (size_t i = 1; i < toFactor.data.length; i++) {
            big_integer_lshift_64(base, 1);
            base.data.bits[0] = kiss09_ulong(rngState);
        }
        base = big_integer_add(big_integer_mod(base, big_integer_subtract(toFactor, BIGINT_3)), BIGINT_2);

        const bitCapInt testFactor = gcd(toFactor, base);
        if (big_integer_compare(testFactor, BIGINT_1) != 0) {
            outputs[thread] = testFactor;

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
        bitCapInt r =big_integer_create(kiss09_ulong(rngState));
        for (size_t i = 1; i < fullMaxR.data.length; i++) {
            big_integer_lshift_64(r, 1);
            r.data.bits[0] = kiss09_ulong(rngState);
        }
        r = big_integer_add(big_integer_mod(r, rRange), rMin);

        // Since our output is r rather than y, we can skip the continued fractions step.
        const bitCapInt p = BIT_AND_1(r) ? r : big_integer_rshift(r, 1);

        const bitCapInt rGuess = BIT_AND_1(isRsaSemiprime) ? p : r;

        // As a "classical" optimization, since \phi(toFactor) and factor bounds overlap,
        // we first check if our guess for r is already a factor.
        if ((big_integer_compare(rGuess, BIGINT_1) != 0)
            && (big_integer_compare(big_integer_mul(big_integer_div(toFactor, rGuess), rGuess), toFactor) == 0) {
            outputs[thread] = rGuess;

            return;
        }

        const bitCapInt apowrhalf = big_integer_mod(uipow(base, p), toFactor);
        bitCapInt f1 = gcd(big_integer_add(apowrhalf, BIGINT_1), toFactor);
        bitCapInt f2 = gcd(big_integer_subtract(apowrhalf, BIGINT_1), toFactor);
        bitCapInt fmul = big_integer_mul(f1, f2);
        while ((big_integer_compare(fmul, BIGINT_1) != 0)
            && (big_integer_compare(fmul, toFactor) != 0)
            && (big_integer_compare(big_integer_mul(big_integer_div(toFactor, fmul), fmul), toFactor) == 0) {
            fmul = f1;
            f1 = big_integer_mul(fmul, f2);
            f2 = big_integer_div(toFactor, big_integer_mul(fmul, f2));
            fmul = big_integer_mul(f1, f2);
        }
        if ((big_integer_compare(fmul, BIGINT_1) != 0)
            && (big_integer_compare(fmul, toFactor) == 0)
            && (big_integer_compare(f1, BIGINT_1) != 0)
            && (big_integer_compare(f2, BIGINT_1) != 0)) {
            outputs[thread] = f1;

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
