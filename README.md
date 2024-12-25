# qimcifa
Quantum-inspired Monte Carlo integer factoring algorithm

## Copyright and license
(c) Daniel Strano and the Qrack contributors 2017-2024. All rights reserved.

Licensed under the GNU Lesser General Public License V3.
See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html for details.

## About 
This example originally demonstrated a (Shor's-like) "quantum-inspired" algorithm for integer factoring. It has since been developed into a general factoring algorithm and tool.

The only potentially "original" part of this factoring algorithm is the "reverse wheel factorization," as far as I can tell. The idea is, instead of performing typical wheel factorization or trial division, we collect a short list of the first primes and remove all of their multiples from a "brute-force" guessing range by mapping a dense contiguous integer set, to a set without these multiples, by successively applying `guess = guess + guess / (p[i] - 1U) + 1U` for prime "`p`" in ascending (or any) order. Each prime applied this way effectively multiplies the brute-force guessing cardinality by a fraction (p-1)/p. Whatever "level" of primes we use, the cost per "guess" becomes higher.

Then, we have a tuner that empirically estimates the cost per guess, and we multiply this by the (known) total cardinality of potential guesses. Whichever reverse wheel factorization level has the lowest product of average cost per guess times guessing set cardinality should have the best performance, and the best level increases with the scale of the problem.

Beyond this, we gain a functional advantage of a square-root over a more naive approach, by setting the brute force guessing range only between the highest prime in reverse wheel factorization and the (modular) square root of the number to factor: if the number is semiprime, there is exactly one correct answer in this range, but including both factors in the range to search would cost us the square root advantage.

Beyond that, we observed that many simple and well-known factoring techniques just don't pay dividends, for semiprime factoring. There's basically no point in checking either congruence of squares or even for a greatest common divisor, as these techniques require some dynamically-variable overhead, and it tends to be faster (for semiprimes) just to check if a guess is an exact factor, on the smallest range we can identify that contains at least one correct answer.

So, this is actually quite rudimentary and just "brute force," except for "reverse wheel factorization" and the upper bound on the guessing range. It just better work entirely in CPU cache, then, but it only requires de minimis maximum memory footprint. (There are congruence of squares and greatest common divisor checks available for numbers besides semiprimes.)

Theoretically, this algorithm might return to its original "quantum-inspired" design with the availability of a high-quality, high-throughput generator of uniform random bit strings. If we were to use the algorithm as-is, except guessing according to a uniform random distribution instead of systematically ascending through every possible "guess," then the average time to solution can be realized in any case, unlike the deterministic version of the algorithm. Then, no work towards the solution can ever be lost in event of interruption of the program, because every single guess (even the first) has the same probability (in the ideal) of leading to successful factoring.

A nearly "brute-force" technique like this has a surprising advantage: basically 0 network communication is needed to coordinate an arbitrarily high amount of parallelism to factor a single number. Each trial division instance is effectively 100% independent of all others (i.e. entirely "embarrassingly parallel"), so Qimcifa offers an interface that allows work to be split between an arbitrarily high number of nodes with absolutely no network communication at all. In terms of incentives of those running different, cooperating nodes in the context of this specific number of integer factoring, all one ultimately cares about is knowing the correct factorization answer _by any means._ For pratical applications, there is no point at all in factoring a number whose factors are already known. When a hypothetical answer is forwarded to the (0-communication) "network" of collaborating nodes, _it is trivial to check whether the answer is correct_ (such as by simply entering the multiplication and equality check with the original number into a Python shell console)! Hence, collaborating node operators only need to trust that all participants in the "network" are actually performing their alloted segment of guesses and would actually communicate the correct answer to the entire group of collaborating nodes if any specific invidual happened to find the answer, but any purported answer is still trivial to verify.
