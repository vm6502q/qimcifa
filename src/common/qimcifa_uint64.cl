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

#include "kiss09.cl"

#define bitLenInt unsigned char
#define bitCapInt unsigned long

bitCapInt gcd(bitCapInt n1, bitCapInt n2)
{
    while (n2) {
        const bitCapInt temp = n2;
        n2 = n1 % n2;
        n1 = temp;
    }
    return n1;
}

void kernel qimcifa_batch(global bitCapInt* bitCapIntArgs, global bitCapInt* rngSeeds)
{
    const bitCapInt threads = get_global_size(0);
    const bitCapInt toFactor = bitCapIntArgs[0];
    const bitCapInt batchSize = bitCapIntArgs[1];
    const bitCapInt nodeId = bitCapIntArgs[2];
    const bitCapInt nodeCount = bitCapIntArgs[3];
    const bitCapInt fullMinR = bitCapIntArgs[4];
    const bitCapInt fullMaxR = bitCapIntArgs[5];
    // This can become a bit flag register, in the future:
    const bitCapInt isRsaSemiprime = bitCapIntArgs[6];

    const kiss09_state rngState = vload4(threads, rngSeeds);

    const bitLenInt byteCount = (log2(toFactor) >> 3) + 1;

    const bitCapInt fullRange = fullMaxR + 1 - fullMinR;
    const bitCapInt nodeRange = fullRange / nodeCount;
    const bitCapInt nodeMin = fullMinR + nodeRange * nodeId;
    const bitCapInt nodeMax = ((nodeId + 1) == nodeCount) ? fullMaxR : (fullMinR + nodeRange * (nodeId + 1) - 1);
    const bitCapInt threadRange = (nodeMax + 1 - nodeMin) / threads;
    const bitCapInt rMin = nodeMin + threadRange * cpu;
    const bitCapInt rMax = ((cpu + 1) == threads) ? nodeMax : (nodeMin + threadRange * (cpu + 1) - 1);

    // TODO: Replace the initialization of random base and r distributions, equivalently,
    // with functions for random number generation.

    for (size_t batchItem = 0; batchItem < batchSize; batchItem++) {
        // Choose a base at random, >1 and <toFactor.
        bitCapInt base = baseDist[0](rand_gen);
#if QBCAPPOW > 5U
        for (size_t i = 1; i < byteCount; i++) {
            base <<= wordSize;
            base |= baseDist[i](rand_gen);
        }
        base += 2U;
#endif

        const bitCapInt testFactor = gcd(toFactor, base);
        if (testFactor != 1) {
            // Inform the other threads on this node that we've succeeded and are done, if feasible.
            std::cout << "Chose non-relative prime: " << testFactor << " * " << (toFactor / testFactor) << std::endl;
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
        bitCapInt r = rDist[0](rand_gen);
#if QBCAPPOW > 5U
        for (size_t i = 1; i < byteCount; i++) {
            r <<= wordSize;
            r |= rDist[i](rand_gen);
        }
        r += rMin;
#endif
        // Since our output is r rather than y, we can skip the continued fractions step.
        const bitCapInt p = (r & 1) ? r : (r >> 1);

#define PRINT_SUCCESS(f1, f2, toFactor, message)                                                                       \
    std::cout << message << (f1) << " * " << (f2) << " = " << (toFactor) << std::endl;

        const bitCapInt rGuess = isRsaSemiprime ? p : r;

        // As a "classical" optimization, since \phi(toFactor) and factor bounds overlap,
        // we first check if our guess for r is already a factor.
        if ((rGuess > 1) && ((toFactor / rGuess) * rGuess) == toFactor) {
            // Inform the other threads on this node that we've succeeded and are done, if feasible.
            PRINT_SUCCESS(rGuess, toFactor / rGuess, toFactor, "Success (on r trial division): Found ");

            return;
        }

        const bitCapInt apowrhalf = uipow(base, p) % toFactor;
        bitCapInt f1 = (bitCapInt)gcd(apowrhalf + 1, toFactor);
        bitCapInt f2 = (bitCapInt)gcd(apowrhalf - 1, toFactor);
        bitCapInt fmul = f1 * f2;
        while ((fmul > 1) && (fmul != toFactor) && (((toFactor / fmul) * fmul) == toFactor)) {
            fmul = f1;
            f1 = fmul * f2;
            f2 = toFactor / (fmul * f2);
            fmul = f1 * f2;
        }
        if ((fmul > 1) && (fmul == toFactor) && (f1 > 1) && (f2 > 1)) {
            // Inform the other threads on this node that we've succeeded and are done, if feasible.
            PRINT_SUCCESS(rGuess, toFactor / rGuess, toFactor, "Success (on r difference of squares): Found ");

            return;
        }
    }
}
