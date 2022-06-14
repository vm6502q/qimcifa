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

inline void bi_set_0(BigInteger* p)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        p->bits[i] = 0;
    }
}

inline void bi_copy(const BigInteger* in, BigInteger* out)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        out->bits[i] = in->bits[i];
    }
}

inline int bi_compare(const BigInteger* left, const BigInteger* right)
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

inline int bi_compare_0(const BigInteger* left)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        if (left->bits[i]) {
            return 1;
        }
    }

    return 0;
}

inline int bi_compare_1(const BigInteger* left)
{
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i > 0; --i) {
        if (left->bits[i]) {
            return 1;
        }
    }
    if (left->bits[0] > 1) {
        return 1;
    }
    if (left->bits[0] < 1) {
        return -1;
    }

    return 0;
}

inline BigInteger bi_add(const BigInteger* left, const BigInteger* right)
{
    BigInteger result;
    result.bits[0] = 0;
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result.bits[i] += left->bits[i] + right->bits[i];
        result.bits[i + 1] = (result.bits[i] < left->bits[i]) ? 1 : 0;
    }
    result.bits[BIG_INTEGER_MAX_WORD_INDEX] += right->bits[BIG_INTEGER_MAX_WORD_INDEX];

    return result;
}

inline void bi_add_ip(BigInteger* left, const BigInteger* right)
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

inline BigInteger bi_sub(const BigInteger* left, const BigInteger* right)
{
    BigInteger result;
    result.bits[0] = 0;
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result.bits[i] += left->bits[i] - right->bits[i];
        result.bits[i + 1] = (result.bits[i] > left->bits[i]) ? -1 : 0;
    }
    result.bits[BIG_INTEGER_MAX_WORD_INDEX] -= right->bits[BIG_INTEGER_MAX_WORD_INDEX];

    return result;
}

inline void bi_sub_ip(BigInteger* left, const BigInteger* right)
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

inline void bi_increment(BigInteger* pBigInt, BIG_INTEGER_WORD value)
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

inline void bi_decrement(BigInteger* pBigInt, BIG_INTEGER_WORD value)
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

inline BigInteger bi_create(BIG_INTEGER_WORD val) {
    BigInteger result;
    result.bits[0] = val;
    for (int i = 1; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = 0;
    }

    return result;
}

inline BigInteger bi_load(BIG_INTEGER_WORD* a)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = a[i];
    }

    return result;
}

inline BigInteger bi_lshift_word(const BigInteger* left, BIG_INTEGER_WORD rightMult)
{
    if (!rightMult) {
        return *left;
    }

    BigInteger result;
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = left->bits[i - rightMult];
    }
    for (BIG_INTEGER_WORD i = 0; i < rightMult; ++i) {
        result.bits[i] = 0;
    }

    return result;
}

inline void bi_lshift_word_ip(BigInteger* left, BIG_INTEGER_WORD rightMult)
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

inline BigInteger bi_rshift_word(const BigInteger* left, BIG_INTEGER_WORD rightMult)
{
    if (!rightMult) {
        return *left;
    }

    BigInteger result;
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i - rightMult] = left->bits[i];
    }
    for (BIG_INTEGER_WORD i = 0; i < rightMult; ++i) {
        result.bits[BIG_INTEGER_WORD_SIZE - (i + 1)] = 0;
    }

    return result;
}

inline void bi_rshift_word_ip(BigInteger* left, BIG_INTEGER_WORD rightMult)
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

inline BigInteger bi_lshift(const BigInteger* left, BIG_INTEGER_WORD right)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

    BigInteger result = bi_lshift_word(left, rShift64);
    if (!rMod) {
        return result;
    }

    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        right = result.bits[i];
        result.bits[i] = carry | (right << rMod);
        carry = right >> rModComp;
    }

    return result;
}

inline void bi_lshift_ip(BigInteger* left, BIG_INTEGER_WORD right)
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

inline BigInteger bi_rshift(const BigInteger* left, BIG_INTEGER_WORD right)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

    BigInteger result = bi_rshift_word(left, rShift64);
    if (!rMod) {
        return result;
    }

    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
        right = result.bits[i];
        result.bits[i] = carry | (right >> rMod);
        carry = right << rModComp;
    }

    return result;
}

inline void bi_rshift_ip(BigInteger* left, BIG_INTEGER_WORD right)
{
    const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
    const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

    bi_rshift_word_ip(left, rShift64);
    if (!rMod) {
        return;
    }

    const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
    BIG_INTEGER_WORD carry = 0;
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
        right = left->bits[i];
        left->bits[i] = carry | (right >> rMod);
        carry = right << rModComp;
    }
}

inline int bi_log2(const BigInteger* n)
{
    int pw = 0;
    BigInteger p = bi_rshift(n, 1U);
    while (bi_compare_0(&p) != 0) {
        bi_rshift_ip(&p, 1U);
        ++pw;
    }
    return pw;
}

inline int bi_and_1(const BigInteger* left) { return left->bits[0] & 1; }

inline BigInteger bi_and(const BigInteger* left, const BigInteger* right)
{
    BigInteger result;
    bi_copy(left, &result);
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] &= right->bits[i];
    }

    return result;
}

inline void bi_and_ip(BigInteger* left, const BigInteger* right)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        left->bits[i] &= right->bits[i];
    }
}

inline BigInteger bi_or(const BigInteger* left, const BigInteger* right)
{
    BigInteger result;
    bi_copy(left, &result);
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] |= right->bits[i];
    }

    return result;
}

inline void bi_or_ip(BigInteger* left, const BigInteger* right)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        left->bits[i] |= right->bits[i];
    }
}

inline BigInteger bi_xor(const BigInteger* left, const BigInteger* right)
{
    BigInteger result;
    bi_copy(left, &result);
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] ^= right->bits[i];
    }

    return result;
}

inline void bi_xor_ip(BigInteger* left, const BigInteger* right)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        left->bits[i] ^= right->bits[i];
    }
}

inline BigInteger bi_not(const BigInteger* left)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = ~(left->bits[i]);
    }

    return result;
}

inline void bi_not_ip(BigInteger* left)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        left->bits[i] = ~(left->bits[i]);
    }
}


// "School book multiplication" (on half words)
// Complexity - O(x^2)
BigInteger bi_mul(const BigInteger* left, const BigInteger* right)
{
    int maxI = BIG_INTEGER_WORD_SIZE << 1;
    const BIG_INTEGER_WORD wordSize = BIG_INTEGER_WORD_BITS >> 1U;
    const BIG_INTEGER_WORD m = (1ULL << wordSize) - 1ULL;

    BigInteger result = bi_create(0);
    for (int i = 0; i < maxI; ++i) {
        BIG_INTEGER_WORD carry = 0;
        int i2 = i >> 1;
        for (int j = 0; j < (maxI - i); ++j) {
            int j2 = j >> 1;
            if (((i & 1) == 0) && ((j & 1) == 0)) {
                BIG_INTEGER_WORD temp = (right->bits[j2] & m) * (left->bits[i2] & m) + (result.bits[i2 + j2] & m) + carry;
                carry = temp >> wordSize;
                result.bits[i2 + j2] &= ~m;
                result.bits[i2 + j2] |= temp & m;
            } else if (((i & 1) == 0) && ((j & 1) == 1)) {
                BIG_INTEGER_WORD temp = (right->bits[j2] >> wordSize) * (left->bits[i2] & m) + (result.bits[i2 + j2] >> wordSize) + carry;
                carry = temp >> wordSize;
                result.bits[i2 + j2] &= ~(m << wordSize);
                result.bits[i2 + j2] |= (temp & m) << wordSize;
            } else if (((i & 1) == 1) && ((j & 1) == 0)) {
                BIG_INTEGER_WORD temp = (right->bits[j2] & m) * (left->bits[i2] >> wordSize) + (result.bits[i2 + j2] >> wordSize) + carry;
                carry = temp >> wordSize;
                result.bits[i2 + j2] &= ~(m << wordSize);
                result.bits[i2 + j2] |= (temp & m) << wordSize;
            } else {
                BIG_INTEGER_WORD temp = (right->bits[j2] >> wordSize) * (left->bits[i2] >> wordSize) + (result.bits[i2 + j2 + 1] & m) + carry;
                carry = temp >> wordSize;
                result.bits[i2 + j2 + 1] &= ~m;
                result.bits[i2 + j2 + 1] |= temp & m;
            }
        }
    }

    return result;
}

#if 0
// Adapted from Qrack! (The fundamental algorithm was discovered before.)
// Complexity - O(log)
BigInteger bi_mul(const BigInteger* left, const BigInteger* right)
{
    int rightLog2 = bi_log2(right);
    if (rightLog2 == 0) {
        // right == 1
        return *left;
    }
    int maxI = BIG_INTEGER_BITS - rightLog2;

    BigInteger result;
    bi_set_0(&result);
    for (int i = 0; i < maxI; ++i) {
        BigInteger partMul = bi_lshift(right, i);
        if (bi_compare_0(&partMul) == 0) {
            break;
        }
        const int iWord = i / BIG_INTEGER_WORD_BITS;
        if (1 & (left->bits[iWord] >> (i - (iWord * BIG_INTEGER_WORD_BITS)))) {
            for (int j = iWord; j < BIG_INTEGER_WORD_SIZE; j++) {
                BIG_INTEGER_WORD temp = result.bits[j];
                result.bits[j] += partMul.bits[j];
                int k = j;
                while ((k < BIG_INTEGER_WORD_SIZE) && (temp > result.bits[k])) {
                    temp = result.bits[++k]++;
                }
            }
        }
    }

    return result;
}
#endif

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

    if (lrCompare == 0) {
        // left == right
        if (quotient) {
            // quotient = 1
            bi_set_0(quotient);
            quotient->bits[0] = 1;
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
    if (!rmndr) {
        rmndr = &leftCopy;
    }
    bi_copy(left, rmndr);
    for (int i = ((BIG_INTEGER_BITS - 1) - rightLog2); i >= 0; --i) {
        BigInteger partMul = bi_lshift(right, i);

        const int c = bi_compare(rmndr, &partMul);
        if (c < 0) {
            continue;
        }

        // Past this point, rmndr >= partMul.

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
            // rmndr == partMul
            bi_set_0(rmndr);
            return;
        }

        // Otherwise, c > 1, meaning rmndr > partMul.
        for (int j = i / BIG_INTEGER_WORD_BITS; j < BIG_INTEGER_MAX_WORD_INDEX; ++j) {
            BIG_INTEGER_WORD temp = rmndr->bits[j];
            rmndr->bits[j] -= partMul.bits[j];
            if (rmndr->bits[j] > temp) {
                --(rmndr->bits[j + 1]);
            }
        }
        rmndr->bits[BIG_INTEGER_MAX_WORD_INDEX] -= partMul.bits[BIG_INTEGER_MAX_WORD_INDEX];
    }
}
