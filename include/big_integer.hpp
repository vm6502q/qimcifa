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

    BigInteger() {
        for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
            bits[i] = 0;
        }
    }

    BigInteger(const BigInteger& o) {
        for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
            bits[i] = o.bits[i];
        }
    }

    BigInteger& operator=(const BigInteger& o) {
        for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
            bits[i] = o.bits[i];
        }

        return *this;
    }
} BigInteger;

inline void bi_set_0(BigInteger* p)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        p->bits[i] = 0;
    }
}

int bi_compare(const BigInteger& left, const BigInteger& right)
{
    for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
        if (left.bits[i] > right.bits[i]) {
            return 1;
        }
        if (left.bits[i] < right.bits[i]) {
            return -1;
        }
    }

    return 0;
}

int bi_compare_0(const BigInteger& left)
{
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        if (left.bits[i]) {
            return 1;
        }
    }

    return 0;
}

BigInteger bi_add(const BigInteger& left, const BigInteger& right)
{
    BigInteger result(left);

    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result.bits[i] += right.bits[i];
        if (result.bits[i] < left.bits[i]) {
            ++result.bits[i + 1];
        }
    }
    result.bits[BIG_INTEGER_MAX_WORD_INDEX] += right.bits[BIG_INTEGER_MAX_WORD_INDEX];

    return result;
}

BigInteger bi_sub(const BigInteger& left, const BigInteger& right)
{
    BigInteger result(left);

    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result.bits[i] -= right.bits[i];
        if (result.bits[i] > left.bits[i]) {
            result.bits[i + 1]--;
        }
    }
    result.bits[BIG_INTEGER_MAX_WORD_INDEX] -= right.bits[BIG_INTEGER_MAX_WORD_INDEX];

    return result;
}

inline void bi_increment(BigInteger* pBigInt, const BIG_INTEGER_WORD& value)
{
    int i = 0;
    BIG_INTEGER_WORD temp = pBigInt->bits[i];
    pBigInt->bits[i] += value;
    while ((i < BIG_INTEGER_MAX_WORD_INDEX) && (temp > pBigInt->bits[i])) {
        ++i;
        temp = pBigInt->bits[i];
        ++(pBigInt->bits[i]);
    }
}

inline void bi_decrement(BigInteger* pBigInt, const BIG_INTEGER_WORD& value)
{
    int i = 0;
    BIG_INTEGER_WORD temp = pBigInt->bits[i];
    pBigInt->bits[i] -= value;
    while ((i < BIG_INTEGER_MAX_WORD_INDEX) && (temp < pBigInt->bits[i])) {
        ++i;
        temp = pBigInt->bits[i];
        pBigInt->bits[i]--;
    }
}

BigInteger bi_create(const BIG_INTEGER_WORD& value)
{
    BigInteger bigInt;
    bigInt.bits[0] = value;

    return bigInt;
}

BigInteger bi_load(BIG_INTEGER_WORD* a)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = a[i];
    }

    return result;
}

BigInteger bi_lshift_word(const BigInteger& left, const BIG_INTEGER_WORD& rightMult)
{
    if (!rightMult) {
        return left;
    }

    if (rightMult >= BIG_INTEGER_WORD_SIZE) {
        return BigInteger();
    }

    BigInteger result;
    for (BIG_INTEGER_WORD i = 0; i < rightMult; ++i) {
        result.bits[i] = 0;
    }
    for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i + rightMult] = left.bits[i];
    }

    return result;
}

BigInteger bi_rshift_word(const BigInteger& left, const BIG_INTEGER_WORD& rightMult)
{
    if (!rightMult) {
        return left;
    }

    if (rightMult >= BIG_INTEGER_WORD_SIZE) {
        return BigInteger();
    }

    const int maxLcv = BIG_INTEGER_WORD_SIZE - rightMult;
    BigInteger result;
    for (int i = 0; i < maxLcv; ++i) {
        result.bits[i] = left.bits[i + rightMult];
    }
    for (int i = maxLcv; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = 0;
    }

    return result;
}

BigInteger bi_lshift(const BigInteger& left, const BIG_INTEGER_WORD& right)
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
        result.bits[i] = carry | (left.bits[i] << rMod);
        carry = left.bits[i] >> rModComp;
    }

    return result;
}

BigInteger bi_rshift(const BigInteger& left, const BIG_INTEGER_WORD& right)
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
        result.bits[i] = carry | (left.bits[i] >> rMod);
        carry = left.bits[i] << rModComp;
    }

    return result;
}

int bi_log2(const BigInteger& n)
{
    int pow = 0;
    BigInteger p = bi_rshift(n, 1U);
    while (bi_compare_0(p) != 0) {
        p = bi_rshift(p, 1U);
        pow++;
    }
    return pow;
}

int bi_and_1(const BigInteger& left) { return left.bits[0] & 1; }

BigInteger bi_and(const BigInteger& left, const BigInteger& right)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = left.bits[i] & right.bits[i];
    }
    return result;
}

BigInteger bi_or(const BigInteger& left, const BigInteger& right)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = left.bits[i] | right.bits[i];
    }
    return result;
}

BigInteger bi_xor(const BigInteger& left, const BigInteger& right)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = left.bits[i] ^ right.bits[i];
    }
    return result;
}

BigInteger bi_not(const BigInteger& left)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = ~left.bits[i];
    }
    return result;
}

#if 0
// Complexity - O(x^2)

#define BIG_INTEGER_HALF_WORD_SIZE 1
#define BIG_INTEGER_HALF_WORD_BITS 32
#define BIG_INTEGER_HALF_WORD_MASK 0xFFFFFFFF
#define BIG_INT_LO_HALF_WORD(a) ((a)&BIG_INTEGER_HALF_WORD_MASK)
#define BIG_INT_HI_HALF_WORD(a) ((a) >> BIG_INTEGER_HALF_WORD_BITS)
BigInteger bi_mul(const BigInteger& left, const BigInteger& right)
{
    BigInteger result;
    BigInteger partResult;
    for (int i = 0; i < BIG_INTEGER_HALF_WORD_SIZE; ++i) {
        bi_set_0(partResult);
        const int maxJ = BIG_INTEGER_HALF_WORD_SIZE - i;
        for (int j = 0; j < maxJ; ++j) {
            const unsigned long lLoHalfWord = BIG_INT_LO_HALF_WORD(left.bits[j]);
            const unsigned long lHiHalfWord = BIG_INT_HI_HALF_WORD(left.bits[j]);
            const unsigned long rLoHalfWord = BIG_INT_LO_HALF_WORD(right.bits[i]);
            const unsigned long rHiHalfWord = BIG_INT_HI_HALF_WORD(right.bits[i]);

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
BigInteger bi_mul(const BigInteger& left, const BigInteger& right)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_BITS; ++i) {
        const BigInteger partMul = bi_lshift(right, i);
        if (bi_compare_0(partMul) == 0) {
            break;
        }
        if (1 & (left.bits[i / BIG_INTEGER_WORD_BITS] >> (i % BIG_INTEGER_WORD_BITS))) {
            result = bi_add(result, partMul);
        }
    }

    return result;
}

// Adapted from Qrack! (The fundamental algorithm was discovered before.)
// Complexity - O(log)
BigInteger bi_div(const BigInteger& left, const BigInteger& right)
{
    const BigInteger BIG_INT_1 = bi_create(1);
    BigInteger result;
    BigInteger leftCopy(left);
    for (int i = bi_log2(right); i >= 0; --i) {
        const BigInteger partMul = bi_lshift(right, i);
        const int c = bi_compare(leftCopy, partMul);
        if (c < 0) {
            continue;
        }
        result = bi_add(result, bi_lshift(BIG_INT_1, i));
        if (c == 0) {
            break;
        }
        leftCopy = bi_sub(leftCopy, partMul);
    }

    return result;
}

BigInteger bi_mod(const BigInteger& left, const BigInteger& right)
{
    BigInteger temp = bi_div(left, right);
    temp = bi_mul(temp, right);

    return bi_sub(left, temp);
}
