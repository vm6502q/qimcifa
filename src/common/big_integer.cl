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
    result->bits[0] = 0;
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result->bits[i] += left->bits[i] + right->bits[i];
        result->bits[i + 1] = (result->bits[i] < left->bits[i]) ? 1 : 0;
    }
    result->bits[BIG_INTEGER_MAX_WORD_INDEX] += right->bits[BIG_INTEGER_MAX_WORD_INDEX];
}

void bi_add_ip(BigInteger* left, const BigInteger* right)
{
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        BIG_INTEGER_WORD temp = left->bits[i];
        left->bits[i] += right->bits[i];
        int j = i;
        while (left->bits[j] < temp) {
            temp = left->bits[++j]++;
        }
    }
    left->bits[BIG_INTEGER_MAX_WORD_INDEX] += right->bits[BIG_INTEGER_MAX_WORD_INDEX];
}

void bi_sub(const BigInteger* left, const BigInteger* right, BigInteger* result)
{
    result->bits[0] = 0;
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result->bits[i] += left->bits[i] - right->bits[i];
        result->bits[i + 1] = (result->bits[i] > left->bits[i]) ? -1 : 0;
    }
    result->bits[BIG_INTEGER_MAX_WORD_INDEX] -= right->bits[BIG_INTEGER_MAX_WORD_INDEX];
}

void bi_sub_ip(BigInteger* left, const BigInteger* right)
{
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        BIG_INTEGER_WORD temp = left->bits[i];
        left->bits[i] -= right->bits[i];
        int j = i;
        while (left->bits[j] > temp) {
            temp = left->bits[++j]--;
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
    for (BIG_INTEGER_WORD i = 0; i < rightMult; ++i) {
        result->bits[i] = 0;
    }
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i] = left->bits[i - rightMult];
    }
}

void bi_lshift_word_ip(BigInteger* left, BIG_INTEGER_WORD rightMult)
{
    if (!rightMult) {
        return;
    }
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        left->bits[i] = left->bits[i - rightMult];
    }
    for (BIG_INTEGER_WORD i = 0; i < rightMult; ++i) {
        left->bits[i] = 0;
    }
}

void bi_rshift_word(const BigInteger* left, BIG_INTEGER_WORD rightMult, BigInteger* result)
{
    if (!rightMult) {
        bi_copy(left, result);
        return;
    }
    for (BIG_INTEGER_WORD i = 0; i < rightMult; ++i) {
        result->bits[BIG_INTEGER_WORD_SIZE - (i + 1)] = 0;
    }
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result->bits[i - rightMult] = left->bits[i];
    }
}

void bi_rshift_word_ip(BigInteger* left, BIG_INTEGER_WORD rightMult)
{
    if (!rightMult) {
        return;
    }
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        left->bits[i - rightMult] = left->bits[i];
    }
    for (BIG_INTEGER_WORD i = 0; i < rightMult; ++i) {
        left->bits[BIG_INTEGER_WORD_SIZE - (i + 1)] = 0;
    }
}

void bi_lshift(const BigInteger* left, BIG_INTEGER_WORD right, BigInteger* result)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

    bi_lshift_word(left, rShift64, result);
    if (!rMod) {
        return;
    }

    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        right = result->bits[i];
        result->bits[i] = carry | (right << rMod);
        carry = right >> rModComp;
    }
}

void bi_lshift_ip(BigInteger* left, BIG_INTEGER_WORD right)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

    bi_lshift_word_ip(left, rShift64);
    if (!rMod) {
        return;
    }

    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        right = left->bits[i];
        left->bits[i] = carry | (right << rMod);
        carry = right >> rModComp;
    }
}

void bi_rshift(const BigInteger* left, BIG_INTEGER_WORD right, BigInteger* result)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

    bi_rshift_word(left, rShift64, result);
    if (!rMod) {
        return;
    }

    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
        right = result->bits[i];
        result->bits[i] = carry | (right >> rMod);
        carry = right << rModComp;
    }
}

void bi_rshift_ip(BigInteger* left, BIG_INTEGER_WORD right)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

    bi_rshift_word_ip(left, rShift64);
    if (!rMod) {
        return;
    }

    BigInteger lWord;
    bi_copy(left, &lWord);
    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
        right = left->bits[i];
        left->bits[i] = carry | (right >> rMod);
        carry = right << rModComp;
    }
}

int bi_log2(const BigInteger* n)
{
    int pw = 0;
    BigInteger p;
    bi_rshift(n, 1U, &p);
    while (bi_compare_0(&p) != 0) {
        bi_rshift_ip(&p, 1U);
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
        const int iWord = i / BIG_INTEGER_WORD_BITS;
        if (1 & (left->bits[iWord] >> (i - (iWord * BIG_INTEGER_WORD_BITS)))) {
            for (int j = iWord; j < BIG_INTEGER_WORD_SIZE; j++) {
                BIG_INTEGER_WORD temp = result->bits[j];
                result->bits[j] += partMul.bits[j];
                int k = j;
                while ((k < BIG_INTEGER_WORD_SIZE) && (temp > result->bits[k])) {
                    temp = result->bits[++k]++;
                }
            }
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
            int iWord = i / BIG_INTEGER_WORD_BITS;
            const BIG_INTEGER_WORD shiftWord = 1 << (i - (iWord * BIG_INTEGER_WORD_BITS));
            BIG_INTEGER_WORD temp = quotient->bits[iWord];
            quotient->bits[iWord] += shiftWord;
            while ((iWord < BIG_INTEGER_WORD_SIZE) && (temp > quotient->bits[iWord])) {
                temp = quotient->bits[++iWord]++;
            }
        }

        if (c == 0) {
            // leftCopy == partMul
            if (rmndr) {
                bi_set_0(rmndr);
            }
            return;
        }

        // Otherwise, c > 1, meaning leftCopy > partMul.
        for (int j = i / BIG_INTEGER_WORD_BITS; j < BIG_INTEGER_MAX_WORD_INDEX; ++j) {
            BIG_INTEGER_WORD temp = leftCopy.bits[j];
            leftCopy.bits[j] -= partMul.bits[j];
            if (leftCopy.bits[j] > temp) {
                --(leftCopy.bits[j + 1]);
            }
        }
        leftCopy.bits[BIG_INTEGER_MAX_WORD_INDEX] -= partMul.bits[BIG_INTEGER_MAX_WORD_INDEX];
    }

    if (rmndr) {
        bi_copy(&leftCopy, rmndr);
    }
}
