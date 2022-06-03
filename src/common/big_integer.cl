//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qimcifa contributors, 2022. All rights reserved.
//
// This header has adapted for OpenCL and C, from big_integer.c by Andre Azevedo.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

/*
** big_integer.c
**     Description: "Arbitrary"-precision integer
**     Author: Andre Azevedo <http://github.com/andreazevedo>
**/

// The MIT License (MIT)
//
// Copyright (c) 2014 Andre Azevedo
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#define BIG_INTEGER_WORD_SIZE 2
#define BIG_INTEGER_MAX_WORD_INDEX 1
#define BIG_INTEGER_WORD_BITS 64
#define BIG_INTEGER_WORD_POWER 6
#define BIG_INTEGER_WORD unsigned long
#define BIG_INTEGER_BITS 128

typedef struct BigInteger {
    BIG_INTEGER_WORD bits[BIG_INTEGER_WORD_SIZE];
} BigInteger;

void bi_set_0(BigInteger* p)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        p->bits[i] = 0;
    }
}

void bi_copy(const BigInteger* in, BigInteger* out)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        out->bits[i] = in->bits[i];
    }
}

int bi_compare(const BigInteger* left, const BigInteger* right)
{
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
        if (left->bits[i] > right->bits[i]) {
            return 1;
        }
        if (left->bits[i] < right->bits[i]) {
            return -1;
        }
    }

    return 0;
}

int bi_compare_0(const BigInteger* left)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        if (left->bits[i]) {
            return 1;
        }
    }

    return 0;
}

void bi_add(const BigInteger* left, const BigInteger* right, BigInteger* result)
{
    bi_copy(left, result);

    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result->bits[i] += right->bits[i];
        if (result->bits[i] < left->bits[i]) {
            ++(result->bits[i + 1]);
        }
    }
    result->bits[BIG_INTEGER_MAX_WORD_INDEX] += right->bits[BIG_INTEGER_MAX_WORD_INDEX];
}

void bi_add_ip(BigInteger* left, const BigInteger* right)
{
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        BIG_INTEGER_WORD temp = left->bits[i];
        left->bits[i] += right->bits[i];
        if (left->bits[i] < temp) {
            ++(left->bits[i + 1]);
        }
    }
    left->bits[BIG_INTEGER_MAX_WORD_INDEX] += right->bits[BIG_INTEGER_MAX_WORD_INDEX];
}

void bi_sub(const BigInteger* left, const BigInteger* right, BigInteger* result)
{
    bi_copy(left, result);

    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result->bits[i] -= right->bits[i];
        if (result->bits[i] > left->bits[i]) {
            --(result->bits[i + 1]);
        }
    }
    result->bits[BIG_INTEGER_MAX_WORD_INDEX] -= right->bits[BIG_INTEGER_MAX_WORD_INDEX];
}

void bi_sub_ip(BigInteger* left, const BigInteger* right)
{
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        BIG_INTEGER_WORD temp = left->bits[i];
        left->bits[i] -= right->bits[i];
        if (left->bits[i] > temp) {
            --(left->bits[i + 1]);
        }
    }
    left->bits[BIG_INTEGER_MAX_WORD_INDEX] -= right->bits[BIG_INTEGER_MAX_WORD_INDEX];
}

void bi_increment(BigInteger* pBigInt, BIG_INTEGER_WORD value)
{
    BIG_INTEGER_WORD temp = pBigInt->bits[0];
    pBigInt->bits[0] += value;
    if (temp <= pBigInt->bits[0]) {
        return;
    }
    for (int i = 1; i < BIG_INTEGER_WORD_SIZE; i++) {
        temp = pBigInt->bits[i]++;
        if (temp <= pBigInt->bits[i]) {
            break;
        }
    }
}

void bi_decrement(BigInteger* pBigInt, BIG_INTEGER_WORD value)
{
    BIG_INTEGER_WORD temp = pBigInt->bits[0];
    pBigInt->bits[0] -= value;
    if (temp >= pBigInt->bits[0]) {
        return;
    }
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
        temp = pBigInt->bits[i]--;
        if (temp >= pBigInt->bits[i]) {
            break;
        }
    }
}

BigInteger bi_create(BIG_INTEGER_WORD val) {
    BigInteger result;
    result.bits[0] = val;
    for (int i = 1; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = 0;
    }

    return result;
}

BigInteger bi_load(global BIG_INTEGER_WORD* a)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = a[i];
    }

    return result;
}

void bi_lshift_word(const BigInteger* left, BIG_INTEGER_WORD rightMult, BigInteger* result)
{
    if (!rightMult) {
        bi_copy(left, result);
        return;
    }

    bi_set_0(result);
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i] = left->bits[i - rightMult];
    }
}

void bi_rshift_word(const BigInteger* left, BIG_INTEGER_WORD rightMult, BigInteger* result)
{
    if (!rightMult) {
        bi_copy(left, result);
        return;
    }

    bi_set_0(result);
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i - rightMult] = left->bits[i];
    }
}

void bi_lshift(const BigInteger* left, BIG_INTEGER_WORD right, BigInteger* result)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);
    bi_set_0(result);

    bi_lshift_word(left, rShift64, result);
    if (!rMod) {
        return;
    }

    BigInteger lWord;
    bi_copy(result, &lWord);
    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i] = carry | (lWord.bits[i] << rMod);
        carry = lWord.bits[i] >> rModComp;
    }
}

void bi_rshift(const BigInteger* left, BIG_INTEGER_WORD right, BigInteger* result)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);
    bi_set_0(result);

    bi_rshift_word(left, rShift64, result);
    if (!rMod) {
        return;
    }

    BigInteger lWord;
    bi_copy(result, &lWord);
    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
        result->bits[i] = carry | (lWord.bits[i] >> rMod);
        carry = lWord.bits[i] << rModComp;
    }
}

int bi_log2(const BigInteger* n)
{
    int pw = 0;
    BigInteger p;
    bi_rshift(n, 1U, &p);
    while (bi_compare_0(&p) != 0) {
        BigInteger temp;
        bi_copy(&p, &temp);
        bi_rshift(&temp, 1U, &p);
        pw++;
    }
    return pw;
}

int bi_and_1(const BigInteger* left) { return left->bits[0] & 1; }

void bi_and(const BigInteger* left, const BigInteger* right, BigInteger* result)
{
    bi_copy(left, result);
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i] &= right->bits[i];
    }
}

void bi_or(const BigInteger* left, const BigInteger* right, BigInteger* result)
{
    bi_copy(left, result);
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i] |= right->bits[i];
    }
}

void bi_xor(const BigInteger* left, const BigInteger* right, BigInteger* result)
{
    bi_copy(left, result);
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i] ^= right->bits[i];
    }
}

void bi_not(const BigInteger* left, BigInteger* result)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i] = ~(left->bits[i]);
    }
}

#if 0
// Complexity - O(x^2)

#define BIG_INTEGER_HALF_WORD_SIZE 1
#define BIG_INTEGER_HALF_WORD_BITS 32
#define BIG_INTEGER_HALF_WORD_MASK 0xFFFFFFFF
#define BIG_INT_LO_HALF_WORD(a) ((a)&BIG_INTEGER_HALF_WORD_MASK)
#define BIG_INT_HI_HALF_WORD(a) ((a) >> BIG_INTEGER_HALF_WORD_BITS)
BigInteger bi_mul(const BigInteger* left, const BigInteger* right)
{
    BigInteger result;
    BigInteger partResult;
    for (int i = 0; i < BIG_INTEGER_HALF_WORD_SIZE; ++i) {
        bi_set_0(partResult);
        const int maxJ = BIG_INTEGER_HALF_WORD_SIZE - i;
        for (int j = 0; j < maxJ; ++j) {
            const unsigned long lLoHalfWord = BIG_INT_LO_HALF_WORD(left->bits[j]);
            const unsigned long lHiHalfWord = BIG_INT_HI_HALF_WORD(left->bits[j]);
            const unsigned long rLoHalfWord = BIG_INT_LO_HALF_WORD(right->bits[i]);
            const unsigned long rHiHalfWord = BIG_INT_HI_HALF_WORD(right->bits[i]);

            partResult.bits[2 * j] += (lLoHalfWord * rLoHalfWord) |
                ((lHiHalfWord * rLoHalfWord + lLoHalfWord * rHiHalfWord) << BIG_INTEGER_HALF_WORD_BITS);
            partResult.bits[2 * j + 1] += lHiHalfWord * rHiHalfWord;
        }
        result = bi_add(result, bi_lshift_word(partResult, i));
    }

    return result;
}
#endif

// Adapted from Qrack! (The fundamental algorithm was discovered before.)
// Complexity - O(log)
void bi_mul(const BigInteger* left, const BigInteger* right, BigInteger* result)
{
    int rightLog2 = bi_log2(right);
    if (rightLog2 == 0) {
        // right == 1
        bi_copy(left, result);
        return;
    }
    int maxI = BIG_INTEGER_BITS - rightLog2;

    bi_set_0(result);
    for (int i = 0; i < maxI; ++i) {
        BigInteger partMul;
        bi_lshift(right, i, &partMul);
        if (bi_compare_0(&partMul) == 0) {
            break;
        }
        if (1 & (left->bits[i / BIG_INTEGER_WORD_BITS] >> (i % BIG_INTEGER_WORD_BITS))) {
            BigInteger temp;
            bi_copy(result, &temp);
            bi_add(&temp, &partMul, result);
        }
    }
}

// Adapted from Qrack! (The fundamental algorithm was discovered before.)
// Complexity - O(log)
void bi_div_mod(const BigInteger* left, const BigInteger* right, BigInteger* quotient, BigInteger* rmndr)
{
    const int lrCompare = bi_compare(left, right);

    if (lrCompare < 0) {
        // left < right
        if (quotient) {
            // quotient = 0
            bi_set_0(quotient);
        }
        if (rmndr) {
            // rmndr = left
            bi_copy(left, rmndr);
        }
        return;
    }

    BigInteger BIG_INT_1 = bi_create(1);
    if (lrCompare == 0) {
        // left == right
        if (quotient) {
            // quotient = 1
            bi_copy(&BIG_INT_1, quotient);
        }
        if (rmndr) {
            // rmndr = 0
            bi_set_0(rmndr);
        }
        return;
    }

    // Otherwise, past this point, left > right.

    int rightLog2 = bi_log2(right);
    if (rightLog2 == 0) {
        // right == 1
        if (quotient) {
            // quotient = left
            bi_copy(left, quotient);
        }
        if (rmndr) {
            // rmndr = 0
            bi_set_0(rmndr);
        }
        return;
    }

    // Past this point, left > right > 1.

    if (quotient) {
        bi_set_0(quotient);
    }
    BigInteger leftCopy;
    bi_copy(left, &leftCopy);
    for (int i = ((BIG_INTEGER_BITS - 1) - rightLog2); i >= 0; --i) {
        BigInteger partMul;
        bi_lshift(right, i, &partMul);

        const int c = bi_compare(&leftCopy, &partMul);
        if (c < 0) {
            continue;
        }

        // Past this point, leftCopy >= partMul.

        if (quotient) {
            BigInteger temp1, temp2;
            bi_copy(quotient, &temp2);
            bi_lshift(&BIG_INT_1, i, &temp1);
            bi_add(&temp1, &temp2, quotient);
        }

        if (c == 0) {
            // leftCopy == partMul
            if (rmndr) {
                bi_set_0(rmndr);
            }
            return;
        }

        // Otherwise, c > 1, meaning leftCopy > partMul.

        BigInteger temp1;
        bi_copy(&leftCopy, &temp1);
        bi_sub(&temp1, &partMul, &leftCopy);
    }

    if (rmndr) {
        bi_copy(&leftCopy, rmndr);
    }
}
