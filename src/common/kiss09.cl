// BSD 3-Clause License
//
// Copyright (c) 2018, Tadej Ciglarič, Erik Štrumbelj, Rok Češnovar. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
// following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this list of conditions and the following
// disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
// disclaimer in the documentation and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
// products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
@file

Implements KISS (Keep It Simple, Stupid) generator, proposed in 2009.

G. Marsaglia, 64-bit kiss rngs, https://www.thecodingforums.com/threads/64-bit-kiss-rngs.673657.
*/

#define KISS09_FLOAT_MULTI 5.4210108624275221700372640e-20f
#define KISS09_DOUBLE_MULTI 5.4210108624275221700372640e-20

// https://www.thecodingforums.com/threads/64-bit-kiss-rngs.673657/

/**
State of kiss09 RNG.
*/
typedef struct {
    ulong x, c, y, z;
} kiss09_state;

/**
Generates a random 64-bit unsigned integer using kiss09 RNG.

@param state State of the RNG to use.
*/
#define kiss09_ulong(state)                                                                                            \
    (/*multiply with carry*/                                                                                           \
        state.c = state.x >> 6, state.x += (state.x << 58) + state.c,                                                  \
        state.c += state.x<(state.x << 58) + state.c, /*xorshift*/                                                     \
                       state.y ^= state.y << 13, state.y ^= state.y>> 17,                                              \
        state.y ^= state.y << 43, /*linear congruential*/                                                              \
        state.z = 6906969069UL * state.z + 1234567UL, state.x + state.y + state.z)

/**
Generates a random 64-bit unsigned integer using kiss09 RNG.

This is alternative implementation of kiss09 RNG as a function.

@param state State of the RNG to use.
*/
#define kiss09_func_ulong(state) _kiss09_func_ulong(&state)
ulong _kiss09_func_ulong(kiss09_state* state)
{
    // multiply with carry
    ulong t = (state->x << 58) + state->c;
    state->c = state->x >> 6;
    state->x += t;
    state->c += state->x < t;
    // xorshift
    state->y ^= state->y << 13;
    state->y ^= state->y >> 17;
    state->y ^= state->y << 43;
    // linear congruential
    state->z = 6906969069UL * state->z + 1234567UL;
    return state->x + state->y + state->z;
}

/**
Seeds kiss09 RNG.

@param state Variable, that holds state of the generator to be seeded.
@param seed Value used for seeding. Should be randomly generated for each instance of generator (thread).
*/
void kiss09_seed(kiss09_state* state, ulong j)
{
    state->x = 1234567890987654321UL ^ j;
    state->c = 123456123456123456UL ^ j;
    state->y = 362436362436362436UL ^ j;
    if (state->y == 0) {
        state->y = 1;
    }
    state->z = 1066149217761810UL ^ j;
}

/**
Generates a random 32-bit unsigned integer using kiss09 RNG.

@param state State of the RNG to use.
*/
#define kiss09_uint(state) ((uint)kiss09_ulong(state))

/**
Generates a random float using kiss09 RNG.

@param state State of the RNG to use.
*/
#define kiss09_float(state) (kiss09_ulong(state) * KISS09_FLOAT_MULTI)

/**
Generates a random double using kiss09 RNG.

@param state State of the RNG to use.
*/
#define kiss09_double(state) (kiss09_ulong(state) * KISS09_DOUBLE_MULTI)

/**
Generates a random double using kiss09 RNG. Since kiss09 returns 64-bit numbers this is equivalent to kiss09_double.

@param state State of the RNG to use.
*/
#define kiss09_double2(state) kiss09_double(state)
