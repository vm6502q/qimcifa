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

// TODO: Consider this pseudocode, copied directly from qimcifa.cpp.

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


    std::vector<std::future<void>> futures(cpuCount);
    for (unsigned cpu = 0U; cpu < cpuCount; cpu++) {
        futures[cpu] = std::async(std::launch::async,
            [cpu, nodeId, nodeCount, toFactor, fullMinR, fullMaxR, &iterClock, &rand_gen, &isFinished] {
                // These constants are semi-redundant, but they're only defined once per thread,
                // and compilers differ on lambda expression capture of constants.

                // Batching reduces mutex-waiting overhead, on the std::atomic broadcast.
                // Batch size is BASE_TRIALS * PERIOD_TRIALS.

                // Number of times to reuse a random base:
                const size_t BASE_TRIALS = 1U << 9U;
                // Number of random period guesses per random base:
                const size_t PERIOD_TRIALS = 1U;

                const double clockFactor = 1.0 / 1000.0; // Report in ms
                const unsigned threads = std::thread::hardware_concurrency();

                const bitCapInt fullRange = fullMaxR + 1U - fullMinR;
                const bitCapInt nodeRange = fullRange / nodeCount;
                const bitCapInt nodeMin = fullMinR + nodeRange * nodeId;
                const bitCapInt nodeMax =
                    ((nodeId + 1U) == nodeCount) ? fullMaxR : (fullMinR + nodeRange * (nodeId + 1U) - 1U);
                const bitCapInt threadRange = (nodeMax + 1U - nodeMin) / threads;
                const bitCapInt rMin = nodeMin + threadRange * cpu;
                const bitCapInt rMax = ((cpu + 1U) == threads) ? nodeMax : (nodeMin + threadRange * (cpu + 1U) - 1U);

                std::vector<rand_dist> baseDist;
                std::vector<rand_dist> rDist;
#if QBCAPPOW < 6U
                baseDist.push_back(rand_dist(2U, toFactor - 1U));
                rDist.push_back(rand_dist(rMin, rMax));
#else
                const bitLenInt wordSize = 32U;
                const bitCapInt wordMask = 0xFFFFFFFF;

                bitCapInt distPart = toFactor - 3U;
                while (distPart) {
                    baseDist.push_back(rand_dist(0U, (uint32_t)(distPart & wordMask)));
                    distPart >>= wordSize;
                }
                std::reverse(rDist.begin(), rDist.end());

                distPart = rMax - rMin;
                while (distPart) {
                    rDist.push_back(rand_dist(0U, (uint32_t)(distPart & wordMask)));
                    distPart >>= wordSize;
                }
                std::reverse(rDist.begin(), rDist.end());
#endif

                for (;;) {
                    for (size_t batchItem = 0U; batchItem < BASE_TRIALS; batchItem++) {
                        // Choose a base at random, >1 and <toFactor.
                        bitCapInt base = baseDist[0U](rand_gen);
#if QBCAPPOW > 5U
                        for (size_t i = 1U; i < baseDist.size(); i++) {
                            base <<= wordSize;
                            base |= baseDist[i](rand_gen);
                        }
                        base += 2U;
#endif

                        const bitCapInt testFactor = gcd(toFactor, base);
                        if (testFactor != 1U) {
                            // Inform the other threads on this node that we've succeeded and are done:
                            isFinished = true;

                            std::cout << "Chose non-relative prime: " << testFactor << " * " << (toFactor / testFactor)
                                      << std::endl;
                            auto tClock = std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now() - iterClock);
                            std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;
                            std::cout << "(Waiting to join other threads...)" << std::endl;
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
                        // const bitCapInt logBaseToFactor = (bitCapInt)intLog(base, toFactor) + 1U;
                        // Euler's Theorem tells us, if gcd(a, n) = 1, then a^\phi(n) = 1 MOD n,
                        // where \phi(n) is Euler's totient for n.
                        // const bitCapInt fullMinR = (minPhi < logBaseToFactor) ? logBaseToFactor : minPhi;

                        // c is basically a harmonic degeneracy factor, and there might be no value in testing
                        // any case except c = 1, without loss of generality.

                        // This sets a nonuniform distribution on our y values to test.
                        // y values are close to qubitPower / rGuess, and we midpoint round.

                        // However, results are better with uniformity over r, rather than y.

                        // So, we guess r, between fullMinR and fullMaxR.
                        for (size_t rTrial = 0U; rTrial < PERIOD_TRIALS; rTrial++) {
                            // Choose a base at random, >1 and <toFactor.
                            bitCapInt r = rDist[0U](rand_gen);
#if QBCAPPOW > 5U
                            for (size_t i = 1U; i < rDist.size(); i++) {
                                r <<= wordSize;
                                r |= rDist[i](rand_gen);
                            }
                            r += rMin;
#endif
                            // Since our output is r rather than y, we can skip the continued fractions step.
                            const bitCapInt p = (r & 1U) ? r : (r >> 1U);

#define PRINT_SUCCESS(f1, f2, toFactor, message)                                                                       \
    std::cout << message << (f1) << " * " << (f2) << " = " << (toFactor) << std::endl;                                 \
    auto tClock =                                                                                                      \
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - iterClock);  \
    std::cout << "(Time elapsed: " << (tClock.count() * clockFactor) << "ms)" << std::endl;                            \
    std::cout << "(Waiting to join other threads...)" << std::endl;

#if IS_RSA_SEMIPRIME
#define RGUESS p
#else
#define RGUESS r
#endif

                            // As a "classical" optimization, since \phi(toFactor) and factor bounds overlap,
                            // we first check if our guess for r is already a factor.
                            if ((RGUESS > 1U) && ((toFactor / RGUESS) * RGUESS) == toFactor) {
                                // Inform the other threads on this node that we've succeeded and are done:
                                isFinished = true;

                                PRINT_SUCCESS(
                                    RGUESS, toFactor / RGUESS, toFactor, "Success (on r trial division): Found ");
                                return;
                            }

                            const bitCapInt apowrhalf = uipow(base, p) % toFactor;
                            bitCapInt f1 = (bitCapInt)gcd(apowrhalf + 1U, toFactor);
                            bitCapInt f2 = (bitCapInt)gcd(apowrhalf - 1U, toFactor);
                            bitCapInt fmul = f1 * f2;
                            while ((fmul > 1U) && (fmul != toFactor) && (((toFactor / fmul) * fmul) == toFactor)) {
                                fmul = f1;
                                f1 = fmul * f2;
                                f2 = toFactor / (fmul * f2);
                                fmul = f1 * f2;
                            }
                            if ((fmul > 1U) && (fmul == toFactor) && (f1 > 1U) && (f2 > 1U)) {
                                // Inform the other threads on this node that we've succeeded and are done:
                                isFinished = true;

                                PRINT_SUCCESS(f1, f2, toFactor, "Success (on r difference of squares): Found ");
                                return;
                            }
                        }
                    }

                    // Check if finished, between batches.
                    if (isFinished) {
                        return;
                    }
                }
            });
    };
