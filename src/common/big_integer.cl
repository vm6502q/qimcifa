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

#define BIG_INTEGER_DATA_MAX_SIZE 10

typedef struct BigIntegerData {
    unsigned int bits[BIG_INTEGER_DATA_MAX_SIZE];
    int length;
} BigIntegerData;

typedef struct BigInteger {
    char sign;
    BigIntegerData data;
} BigInteger;

#define UINT_NUM_BITS (sizeof(unsigned int) * 8)
// const int UINT_NUM_BITS =	(sizeof(unsigned int) * 8);

/* PRIVATE FUNCTIONS IMPLEMENTATION */
void big_integer_clear_trash_data(BigIntegerData* pBigIntData)
{
    int i;
    for (i = pBigIntData->length; i < BIG_INTEGER_DATA_MAX_SIZE; ++i)
        pBigIntData->bits[i] = 0;
};

BigIntegerData big_integer_empty_data()
{
    BigIntegerData bigIntData;
    bigIntData.length = 0;
    big_integer_clear_trash_data(&bigIntData);
    return bigIntData;
};

BigIntegerData big_integer_create_data(const unsigned int bits[], const int length)
{
    BigIntegerData bigIntData;
    if (bits && length > 0)
        for (int i = 0; i < length; i++)
            bigIntData.bits[i] = bits[i];
    bigIntData.length = length;

    big_integer_clear_trash_data(&bigIntData);

    return bigIntData;
};

BigInteger big_integer_create_internal(const char sign, const BigIntegerData data)
{
    BigInteger bigInt;
    bigInt.sign = sign;
    bigInt.data = data;

    return bigInt;
};

void big_integer_normalize_from(BigIntegerData* pBigIntData, const int from)
{
    int i;
    for (i = from; i >= 0; --i) {
        if (pBigIntData->bits[i] != 0) {
            pBigIntData->length = i + 1;
            break;
        }
    }
};

void big_integer_normalize(BigIntegerData* pBigIntData)
{
    big_integer_normalize_from(pBigIntData, BIG_INTEGER_DATA_MAX_SIZE - 1);
};

int big_integer_compare_data(const BigIntegerData* pLeft, const BigIntegerData* pRight)
{
    /* if the lengths are different */
    if (pLeft->length > pRight->length)
        return 1;
    if (pLeft->length < pRight->length)
        return -1;

    int length = pLeft->length;
    int i;
    for (i = (length - 1); i >= 0; --i) {
        if (pLeft->bits[i] > pRight->bits[i])
            return 1;
        if (pLeft->bits[i] < pRight->bits[i])
            return -1;
    }

    return 0;
};

int big_integer_compare_data_uint(const BigIntegerData* pBigIntData, unsigned int value)
{
    if (pBigIntData->length == 0)
        return -1;
    if (pBigIntData->length > 1)
        return 1;

    if (pBigIntData->bits[0] > value)
        return 1;
    else if (pBigIntData->bits[0] < value)
        return -1;

    return 0;
};

#define MAX(l, r) ((l < r) ? r : l)

BigIntegerData big_integer_add_data(const BigIntegerData left, const BigIntegerData right)
{
    int uIntNumBits = UINT_NUM_BITS;

    BigIntegerData result = big_integer_empty_data();

    int len = MAX(left.length, right.length);

    unsigned long sum = 0;
    int i;
    for (i = 0; i < len; ++i) {
        sum += (unsigned long)left.bits[i] + right.bits[i];
        result.bits[i] = (unsigned int)sum;
        sum >>= uIntNumBits;
    }

    if (sum > 0) {
        result.bits[i] = (unsigned int)sum;
        i++;
    }

    result.length = i;

    return result;
};

/* left > right always */
BigIntegerData big_integer_subtract_data(const BigIntegerData left, const BigIntegerData right)
{
    BigIntegerData result = big_integer_empty_data();

    int len = MAX(left.length, right.length);

    unsigned long borrow = 0;
    int i;
    for (i = 0; i < len; ++i) {
        /* what happens here is that, if left is less than right, borrow will become
           "negative" (not really because it is unsigned), and the bit pattern for that is
           the 1's complement (complementing it to get to 0), which is exactly the remainder
           of this term in the subtraction. */
        borrow = (unsigned long)left.bits[i] - right.bits[i] - borrow;

        result.bits[i] = (unsigned int)borrow;

        /* here we just want the first 1 after removing the lower order term */
        borrow = (borrow >> UINT_NUM_BITS) & 1;
    }

    big_integer_normalize_from(&result, i);

    return result;
};

void big_integer_increment_data(BigIntegerData* pBigIntData, const unsigned int value)
{
    unsigned long carry = value;
    int i = 0;
    while (carry > 0) {
        carry += (unsigned long)pBigIntData->bits[i];
        pBigIntData->bits[i] = (unsigned int)carry;
        carry >>= UINT_NUM_BITS;
        ++i;
    }

    if (i > pBigIntData->length)
        pBigIntData->length = i;
};

/* pBigIntData > value */
void big_integer_decrement_data(BigIntegerData* pBigIntData, const unsigned int value)
{
    unsigned long borrow = value;
    int i = 0;
    while (borrow > 0) {
        borrow = (unsigned long)pBigIntData->bits[i] - borrow;
        pBigIntData->bits[i] = (unsigned int)borrow;
        borrow = (borrow >> UINT_NUM_BITS) & 1;
        ++i;
    }

    big_integer_normalize_from(pBigIntData, i);
};

/* PUBLIC FUNCTIONS IMPLEMENTATION */
BigInteger big_integer_create(long value)
{
    BigInteger bigInt;
    int numBits = UINT_NUM_BITS;

    if (value == 0) {
        bigInt.sign = 0;
        bigInt.data.bits[0] = 0;
        bigInt.data.length = 1;
    } else {
        unsigned long uValue;
        if (value < 0) {
            bigInt.sign = -1;
            uValue = (unsigned long)-value;
        } else {
            bigInt.sign = 1;
            uValue = (unsigned long)value;
        }

        bigInt.data.length = 0;
        while (uValue > 0) {
            bigInt.data.bits[bigInt.data.length++] = (unsigned int)uValue;
            uValue >>= numBits;
        }
    }

    big_integer_clear_trash_data(&bigInt.data);

    return bigInt;
};

int big_integer_to_int(const BigInteger bigInt)
{
    if (bigInt.sign == 0)
        return 0;

    if (bigInt.sign == -1)
        return -(int)bigInt.data.bits[0];

    return (int)bigInt.data.bits[0];
};

long big_integer_to_long_long(const BigInteger bigInt)
{
    if (bigInt.sign == 0)
        return 0;

    int uIntNumBits = UINT_NUM_BITS;
    int maxLength = sizeof(long) / sizeof(unsigned int);

    unsigned long result = 0;
    int i = 0;
    for (i = 0; i < bigInt.data.length; ++i) {
        result |= ((unsigned long)bigInt.data.bits[i]) << (uIntNumBits * i);
    }

    if (bigInt.sign == -1)
        return -(long)result;

    return result;
};

int big_integer_compare(const BigInteger left, const BigInteger right)
{
    /* if one is positive and the other negative */
    if (left.sign > right.sign)
        return 1;
    if (left.sign < right.sign)
        return -1;

    /* if they have the same sign */
    char sign = left.sign;
    return sign * big_integer_compare_data(&left.data, &right.data);
};

BigInteger big_integer_add(const BigInteger left, const BigInteger right)
{
    if (left.sign == 0)
        return right;
    if (right.sign == 0)
        return left;

    if (left.sign == right.sign)
        return big_integer_create_internal(left.sign, big_integer_add_data(left.data, right.data));

    /* compare the MOD of the numbers */
    int compRes = big_integer_compare_data(&left.data, &right.data);

    if (compRes == 0)
        return big_integer_create(0);

    if (compRes > 0) /* left > right */
        return big_integer_create_internal(left.sign, big_integer_subtract_data(left.data, right.data));

    return big_integer_create_internal(right.sign, big_integer_subtract_data(right.data, left.data));
};

BigInteger big_integer_subtract(const BigInteger left, const BigInteger right)
{
    if (left.sign == 0)
        return big_integer_create_internal(-right.sign, right.data);
    if (right.sign == 0)
        return left;

    if (left.sign != right.sign)
        return big_integer_create_internal(left.sign, big_integer_add_data(left.data, right.data));

    /* compare the MOD of the numbers */
    int compRes = big_integer_compare_data(&left.data, &right.data);

    if (compRes == 0)
        return big_integer_create(0);

    if (compRes > 0) /* left > right */
        return big_integer_create_internal(left.sign, big_integer_subtract_data(left.data, right.data));

    return big_integer_create_internal(-right.sign, big_integer_subtract_data(right.data, left.data));
};

void big_integer_increment(BigInteger* bigInt, const unsigned int value)
{
    if (bigInt->sign >= 0) /* bigInt >= 0 */
    {
        if (bigInt->sign == 0 && value > 0)
            bigInt->sign = 1;
        big_integer_increment_data(&bigInt->data, value);
    } else /* bigInt < 0 */
    {
        int compRes = big_integer_compare_data_uint(&bigInt->data, value);

        if (compRes == 0) /* |bigInt| == |value| */
        {
            bigInt->sign = 0;
            bigInt->data.length = 0;
            big_integer_clear_trash_data(&bigInt->data);
        } else if (compRes > 0) /* |bigInt| > |value| */
        {
            big_integer_decrement_data(&bigInt->data, value);
        } else /* |bigInt| < |value| */
        {
            bigInt->sign = 1;
            bigInt->data.bits[0] = value - bigInt->data.bits[0];
        }
    }
};

void big_integer_decrement(BigInteger* bigInt, const unsigned int value)
{
    if (bigInt->sign <= 0) /* bigInt <= 0 */
    {
        if (bigInt->sign == 0 && value > 0)
            bigInt->sign = -1;
        big_integer_increment_data(&bigInt->data, value);
    } else /* bigInt > 0 */
    {
        int compRes = big_integer_compare_data_uint(&bigInt->data, value);

        if (compRes == 0) /* |bigInt| == |value| */
        {
            bigInt->sign = 0;
            bigInt->data.length = 0;
            big_integer_clear_trash_data(&bigInt->data);
        } else if (compRes > 0) /* |bigInt| > |value| */
        {
            big_integer_decrement_data(&bigInt->data, value);
        } else /* |bigInt| < |value| */
        {
            bigInt->sign = -1;
            bigInt->data.bits[0] = value - bigInt->data.bits[0];
        }
    }
};

/// Qimcifa API functions
BigInteger big_integer_copy(BigInteger orig) {
    BigInteger result;
    result.sign = orig.sign;
    result.data.length = orig.data.length;
    for (int i = 0; i < BIG_INTEGER_DATA_MAX_SIZE; i++) {
        result.data.bits[i] = orig.data.bits[i];
    }
    
    return result;
}

BigInteger big_integer_mul(const BigInteger left, const BigInteger right)
{
    BigInteger result = big_integer_create(0);
    for (BigInteger i = big_integer_create(0); big_integer_compare(i, right) == -1; big_integer_increment(&i, 1)) {
        result = big_integer_add(result, left);
    }

    return result;
}

BigInteger big_integer_div(const BigInteger left, const BigInteger right)
{
    BigInteger result = big_integer_create(0);
    for (BigInteger i = big_integer_copy(left); big_integer_compare(i, right) >= 0; i = big_integer_subtract(i, right)) {
       big_integer_increment(&result, 1);
    }

    return result;
}

BigInteger big_integer_mod(const BigInteger left, const BigInteger right)
{
    BigInteger temp = big_integer_div(left, right);
    temp = big_integer_mul(temp, right);

    return big_integer_subtract(left, temp);
}

BigInteger big_integer_lshift_64(const BigInteger left, const unsigned int rightMult)
{
    const int maxLcv = left.data.length - rightMult;

    BigInteger result = big_integer_copy(left);
    
    for (int i = 0; i < maxLcv; i++) {
        result.data.bits[i] = left.data.bits[i + rightMult];
    }
    for (int i = maxLcv; i < BIG_INTEGER_DATA_MAX_SIZE; i++) {
        result.data.bits[i] = 0;
    }
    result.sign = (char)(result.data.bits[left.data.length - 1] & 1);
    result.data.length += rightMult;

    return result;
}

BigInteger big_integer_rshift_64(const BigInteger left, const unsigned int rightMult)
{
    const int maxLcv = left.data.length - rightMult;

    BigInteger result = big_integer_copy(left);
    
    for (int i = 0; i < rightMult; i++) {
        result.data.bits[i] = 0;
    }
    for (int i = 0; i < maxLcv; i++) {
        result.data.bits[i + rightMult] = left.data.bits[i];
    }
    result.sign = (char)(result.data.bits[left.data.length - 1] & 1);
    result.data.length -= rightMult;

    return result;
}

BigInteger big_integer_lshift(const BigInteger left, const unsigned int right)
{
    const int maxLcv = left.data.length - 1;
    const unsigned int rShift64 = right / 64;
    const unsigned int rMod = right - (rShift64 * 64);
    const unsigned int rModComp = 64 - rMod;
    const unsigned long rMask = (1UL << rMod) - 1U;

    BigInteger result = big_integer_lshift_64(left, rShift64);

    unsigned long carry = 0;
    for (int i = 0; i < maxLcv; i++) {
        result.data.bits[i] = carry | (left.data.bits[i] << rMod);
        carry = left.data.bits[i] >> rModComp;
    }
    result.sign = (char)(result.data.bits[left.data.length - 1] & 1);
    if (carry) {
        result.data.length++;
    }

    return result;
}

BigInteger big_integer_rshift(const BigInteger left, const unsigned int right)
{
    const int maxLcv = left.data.length - 1;
    const unsigned int rShift64 = right / 64;
    const unsigned int rMod = right - (rShift64 * 64);
    const unsigned int rModComp = 64 - rMod;
    const unsigned long rMask = (1UL << rMod) - 1U;

    BigInteger result = big_integer_rshift_64(left, rShift64);

    unsigned long carry = 0;
    for (int i = maxLcv; i >= 0; i++) {
        result.data.bits[i] = carry | (left.data.bits[i] >> rMod);
        carry = left.data.bits[i] << rModComp;
    }
    result.sign = (char)(result.data.bits[left.data.length - 1] & 1);
    if (!carry) {
        result.data.length--;
    }

    return result;
}
