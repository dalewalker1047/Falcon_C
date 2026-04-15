#include <stdlib.h>
#include <stdint.h>

/*
 * flopt type is a representation of floating point number in a custom format.
 * the flopt number consists of the following components:
 * - 1 bit for the sign (s)
 * - 11 bits for the exponent (e)
 * - 52 bits for the mantissa (m)
 *
 */
typedef uint64_t flopt;


/*
 * flopt_create is a function that creates a flopt number from its components: sign (s), exponent (e), and mantissa (m).
 * The function takes three parameters:
 * - s: an integer representing the sign of the number (0 for positive, 1 for negative)
 * - e: an integer representing the exponent of the number
 * - m: the mantissa of the number, which must satisfy 2^54 <= m < 2^55
 * 
 * If m is 0, then a zero is returned with the same sign
 * If e < -1076, then zero is returned
 * If e >= 1076 and e != 0, m must be within the range (2^54, 2^55) 
 * 
*/
flopt flopt_create(int s, int e, uint64_t m) {
    flopt x;
    uint32_t t; 
    unsigned f;

    // If e < 1076, we need to change it to just be 0 
    e += 1076;
    t = (uint32_t)e >> 31; // t is 0 if e >= 0, and 0xFFFFFFFF if e < 0
    m &= (uint64_t)t - 1; // extract the mantissa if e >= 0, otherwise set m to 0

    // if m is 0, we want e to be 0 as well, but keep s the same. 
    t = (uint32_t)(m >> 54); 
    e &= -(int)t;

    // Now we can construct the flopt number by combining the sign, exponent, and mantissa.
    // The mantissa always has a leading 1 unlesss m is 0, in which case the leading 1 is not present.
    // If the mantissa is not zero, the leading 1 will increment the exponent
    // The final flopt has the sign bit in the MSB, followed by the exponent in the next 11 bits, 
    // and the mantissa in the remaining 52 bits.
    x = (((uint64_t)s << 63) | (m >> 2)) + ((uint64_t)(uint32_t)e << 52);

    // Use the last 3 bits of the mantissa for rounding
    // if the last three bits are 011, 110, or 111, we need to increment 
    f = (unsigned)m & 7U; // f is the last 3 bits of m, which are used for rounding
    x += (0xC8U >> f) & 1; // 0xC8U is 11001000, so if f is 3, 6, or 7, we add 1 to x for rounding

    return x;
}

/*
 * Time-constant right shift of an unsigned 64-bit number
 * 
 * shift must be in the range [0, 63]
*/
uint64_t unsigned_shift_right(uint64_t x, int shift) {
    x ^= (x ^ (x >> 32)) & -(uint64_t)(shift >> 5);
    return x >> (shift & 63);
}

/*
 * Time-constant right shift of a signed 64-bit number
 * 
 * shift must be in the range [0, 63]
*/
int64_t signed_shift_right(int64_t x, int shift) {
    x ^= (x ^ (x >> 32)) & -(int64_t)(shift >> 5);
    return x >> (shift & 63);
}

/*
 * Time-constant left shift of an unsigned 64-bit number
 * 
 * shift must be in the range [0, 63]
 */
uint64_t unsigned_shift_left(uint64_t x, int shift) {
    x ^= (x ^ (x << 32)) & -(uint64_t)(shift >> 5);
    return x << (shift & 63);
}

/*
 * Time-constant left shift of a signed 64-bit number
 * 
 * shift must be in the range [0, 63]
 */
int64_t signed_shift_left(int64_t x, int shift) {
    x ^= (x ^ (x << 32)) & -(int64_t)(shift >> 5);
    return x << (shift & 63);
}


/*
 * flopt_norm is a function that normalizes a flopt number x, 
 * which means it shifts the mantissa of x left until the leading 1 is in the correct position, 
 * while adjusting the exponent accordingly (so if we shift the value by n bits, then we need to subtract n from the exponent).
 * 
 * The function takes three parameters:
 * - s: an integer representing the sign of the number (0 for positive, 1 for negative)
 * - m: the mantissa of the number, which must satisfy 2^54 <= m < 2^55
 * - e: an integer representing the exponent of the number
 * 
 * The function returns a flopt number with the normalized mantissa and adjusted exponent, and the specified sign.
*/
flopt flopt_norm(int s, uint64_t m, int e) {
    // We're going to do a binary-search-like method to find the position of the leading 1 in the mantissa, 
    // and shift the mantissa left until the leading 1 is in the correct position, while adjusting the exponent accordingly

    uint32_t flag;

    e -= 63; 
    flag = (uint32_t)(m >> 32); // extract the top 32 bits of the mantissa
    flag = (flag | -flag) >> 31;  // flag is 0 if the top 32 bits of the mantissa are 0, and 1 if they are not 0
    m ^= (m ^ (m << 32)) & ((uint64_t)flag - 1);  // ;left shift m by s if flag is 0
    e += (int)(flag << 5); // add 32 to the exponent if flag is 1

    flag = (uint32_t)(m >> 48); // extract the top 16 bits of the mantissa
    flag = (flag | -flag) >> 31;  // flag is 0 if the top 16 bits of the mantissa are 0, and 1 if they are not 0
    m ^= (m ^ (m << 16)) & ((uint64_t)flag - 1); // left shift m by s if flag is 0
    e += (int)(flag << 4); // add 16 to the exponent if flag is 1

    flag = (uint32_t)(m >> 56); // extract the top 8 bits of the mantissa
    flag = (flag | -flag) >> 31;  // flag is 0 if the top 8 bits of the mantissa are 0, and 1 if they are not 0
    m ^= (m ^ (m << 8)) & ((uint64_t)flag - 1); // left shift m by s if flag is 0
    e += (int)(flag << 3); // add 8 to the exponent if flag is 1

    flag = (uint32_t)(m >> 60); // extract the top 4 bits of the mantissa
    flag = (flag | -flag) >> 31;  // flag is 0 if the top 4 bits of the mantissa are 0, and 1 if they are not 0
    m ^= (m ^ (m << 4)) & ((uint64_t)flag - 1); // left shift m by s if flag is 0
    e += (int)(flag << 2); // add 4 to the exponent if flag is 1

    flag = (uint32_t)(m >> 62); // extract the top 2 bits of the mantissa
    flag = (flag | -flag) >> 31;  // flag is 0 if the top 2 bits of the mantissa are 0, and 1 if they are not 0
    m ^= (m ^ (m << 2)) & ((uint64_t)flag - 1); // left shift m by s if flag is 0
    e += (int)(flag << 1); // add 2 to the exponent if flag is 1

    flag = (uint32_t)(m >> 63); // extract the top bit of the mantissa
    flag = (flag | -flag) >> 31;  // flag is 0 if the top bit of the mantissa is 0, and 1 if it is not 0
    m ^= (m ^ (m << 1)) & ((uint64_t)flag - 1); // left shift m by s if flag is 0
    e += (int)flag; // add 1 to the exponent if flag is 1

    return flopt_create(s, e, m); // create a flopt number with the normalized mantissa and adjusted exponent, and the specified sign

}


/*
 * flopt_add is a function that adds two flopt numbers x and y, and returns the result as a flopt number.
 * The function first ensures that abs(x) >= abs(y) by potentially swapping x and y
 * Then it extracts the sign, exponent, and mantissa of x and y, and performs the addition or subtraction of the mantissas depending on the signs of x and y.
 * Finally, it normalizes the result and rounds it to fit into the flopt format before returning it.
*/
flopt flopt_add(flopt x, flopt y) {

    /*
     * First we need to make sure that abs(x) >= abs(y)
     * Unless abs(x) = abs(y), the result will have the same sign as x
     * If abs(x) = abs(y) and they have different signs, the result will be zero with a positive sign
    */

    int sign_x, sign_y; // sign bits of x and y
    int exp_x, exp_y; // exponent bits of x and y
    uint64_t mant_x, mant_y; // mantissa bits of x and y

    uint64_t mask = ((uint64_t)1 << 63) - 1;
    uint64_t abs_diff; 
    abs_diff = (x & mask) - (y & mask);

    // flag to decide whether we need to swap x and y to ensure abs(x) >= abs(y)
    uint32_t switch_flag; 
    switch_flag = (uint32_t) (abs_diff >> 63) | ((1U - (uint32_t)(-abs_diff >> 63)) & (uint32_t)(x >> 63));

    mask = (x ^ y) & -(uint64_t)switch_flag; // mask is 0 if we don't need to switch, and 0xFFFFFFFFFFFFFFFF if we do need to switch

    // switch x and y if necessary to ensure abs(x) >= abs(y)
    x ^= mask;
    y ^= mask;

    // ------------------------------------------------------
    // extract the sign, exponent, and mantissa of x and y

    exp_x = (int)(x >> 52); 
    sign_x = exp_x >> 11; // extract the sign of x from the exponent
    exp_x &= 0x7FF; // extract the exponent of x

    // we need to add the implicit leading 1 to the mantissa if the exponent is not zero
    mask = (uint64_t)(uint32_t)((exp_x + 0x7FF) >> 11) << 52;

    // extract the mantissa of x and add the implicit leading 1 if necessary, then shift left by 3 to make room for rounding bits
    mant_x = ((x & (((uint64_t)1 << 52) - 1)) | mask) << 3; 
    exp_x -= 1078; // adjust the exponent of x to be unbiased
    exp_y = (int)(y >> 52);
    sign_y = exp_y >> 11; // extract the sign of y from the exponent
    exp_y &= 0x7FF; // extract the exponent of y

    mask = (uint64_t)(uint32_t)((exp_y + 0x7FF) >> 11) << 52; // mask for the implicit leading 1 of y
    mant_y = ((y & (((uint64_t)1 << 52) - 1)) | mask) << 3;
    exp_y -= 1078; // adjust the exponent of y to be unbiased


    int exp_diff = exp_x - exp_y; // difference between exponents
    mant_y &= -(uint64_t)((uint32_t)(exp_diff - 60) >> 31); // if the shift count is larger than 59 bits, then we set mant_y to 0
    exp_diff &= 63;


    // we need to shift mant_y to align with mant_x, but we also need to add the bits that will be shifted out of mant_y to the mantissa for rounding
    mask = unsigned_shift_left(1, exp_diff) - 1; // mask for the bits that will be shifted out of mant_y
    mant_y |= (mant_y & mask) + mask; // add the bits that will be shifted out of mant_y to the mantissa for rounding
    mant_y = unsigned_shift_right(mant_y, exp_diff); // shift mant_y to align with mant_x


    // if x and y have the same sign, we can just add the mantissas,
    // otherwise,, we need to subtract the smaller mantissa from the larger mantissa

    mant_x += mant_y - ((mant_y << 1) & -(sign_x ^ sign_y)); // if sign_x and sign_y are the same, we add mant_y to mant_x, otherwise we subtract mant_y from mant_x

    // We need to normalize before returning
    flopt result = flopt_norm((int)sign_x, mant_x, exp_x); 


    int normed_s = result >> 63; // extract the sign of the result
    int normed_e = (result >> 52) & 0x7FF; // extract the exponent of the result
    uint64_t normed_m = result & (((uint64_t)1 << 52) - 1); // extract the mantissa of the result

    // We need to scale the value down to be in (2^54..2^55-1)
    // the last bit is sticky
    normed_m |= ((uint32_t)normed_m & 0x1FF) + 0x1FF;
    normed_m >>= 9; 
    normed_e += 9; // we need to add 9 to the exponent to account for the right shift of the mantissa

    return flopt_create(normed_s, normed_e, normed_m); // create a flopt number with the normalized mantissa and adjusted exponent, and the specified sign

}

/*
 * flopt_sub is a function that subtracts two flopt numbers x and y, and returns the result as a flopt number.
 * The function simply negates the sign of y and adds it to x using the flopt_add function, since subtraction can be performed as
 * 
 * This function has inputs:
 * - x: a flopt number
 * - y: a flopt number
 * and returns a flopt number that is the difference of x and y.
*/
flopt flopt_sub(flopt x, flopt y) {
    // To subtract y from x, we can just negate the sign of y and add it to x
    y ^= (uint64_t)1 << 63; // negate the sign of y by flipping the MSB
    return flopt_add(x, y); // add x and the negated y using the flopt_add function
}


/*
 * flopt_mul is a function that multiplies two flopt numbers x and y, and returns the result as a flopt number.
 * The function first extracts the sign, exponent, and mantissa of x and y, and
 * then it multiplies the mantissas of x and y using a method that avoids overflow by splitting the mantissas into two parts to calculate the product.
 * Finally, it normalizes the result and rounds it to fit into the flopt format before returning it.
 * 
 * This function has inputs:
 * - x: a flopt number
 * - y: a flopt number
 * and returns a flopt number that is the product of x and y.
*/
flopt flopt_mul(flopt x, flopt y) {

    int sign_x, sign_y;         // sign bits of x and y
    int exp_x, exp_y;           // exponent bits of x and y
    uint64_t mant_x, mant_y;    // mantissa bits of x and y
    int sign_z;                 // sign of the result
    int exp_z;                  // exponent of the result
    uint64_t mant_z;            // mantissa of the result
    
    uint32_t x_top, y_top;          // top halves of the mantissas of x and y
    uint32_t x_bot, y_bot;          // bottom halves of the mantissas of x and y
    uint32_t z_bot, z_mid, z_top;   // parts of the product of the mantissas of x and y

    // Extract the mantissas as scaled uints
    mant_x = (x & (((uint64_t)1 << 52) -1)) | ((uint64_t)1 << 52); // extract the mantissa of x and add the implicit leading 1
    mant_y = (y & (((uint64_t)1 << 52) -1)) | ((uint64_t)1 << 52); // extract the mantissa of y and add the implicit leading 1

    // x and y are 53-bit ints, so to multiply, we need to split them into parts to avoid overflow
    x_bot = (uint32_t)mant_x & 0x01FFFFFF; // bottom 25 bits of mant_x
    x_top = (uint32_t)(mant_x >> 25); // top 28 bits of mant_x
    y_bot = (uint32_t)mant_y & 0x01FFFFFF; // bottom 25 bits of mant_y
    y_top = (uint32_t)(mant_y >> 25); // top 28 bits of mant_y


    // this is just FOIL for the product of mant_x and mant_y,
    // since we need to split the mantissas into two parts to multiply

    uint64_t shift_flag; // flag to determine if we need to shift the result of the multiplication
    shift_flag = (uint64_t)x_bot * (uint64_t)y_bot;

    // Multiplication of the last terms of the FOIL expansion 
    z_bot = (uint32_t)shift_flag & 0x01FFFFFF; // bottom 25 bits of the product of the bottom halves
    z_mid = (uint32_t)(shift_flag >> 25); // top 7 bits of the product of the bottom halves, which will be added to the middle part of the result

    // Multiplication of the outer terms of the FOIL expansion
    shift_flag = (uint64_t)x_bot * (uint64_t)y_top;
    z_mid += (uint32_t)shift_flag & 0x01FFFFFF; // add the product of the bottom half of x and the top half of y to z_mid
    z_top = (uint32_t)(shift_flag >> 25); // top 7 bits of the product of the bottom half of x and the top half of y

    // Multiplication of the inner terms of the FOIL expansion
    shift_flag = (uint64_t)x_top * (uint64_t)y_bot;
    z_mid += (uint32_t)shift_flag & 0x01FFFFFF; // add the product of the top half of x and the bottom half of y to z_mid
    z_top += (uint32_t)(shift_flag >> 25); // add the top 7 bits of the product of the top half of x and the bottom half of y to z_top

    // Multiplication of the first terms of the FOIL expansion
    mant_z = (uint64_t)x_top * (uint64_t)y_top; // product of the top halves of x and y
    z_top += (z_mid >> 25); // we need to add the carry from z_mid to z_top, since z_mid is 25 bits, the carry can only be 1, so we just need to add the top bit of z_mid to z_top
    z_mid &= 0x01FFFFFF; // we need to add the carry from z_mid to z_top, and then mask z_mid to be 25 bits
    mant_z += z_top; // add z_top to the product of the top halves (at this point, mant_z is the final product without z_bot)

    // mant_x and mant_y are in the range [2^52, 2^53), 
    // so their product is in the range [2^104, 2^106).
    // we want mant_z to be in the range [2^54, 2^55), 
    // so we need to reassemble and round the product of the mantissas accordingly

    // we need to add the bits that will be shifted out of z_mid and z_bot to the mantissa for rounding
    // z_bot and z_mid are treated as the sticky bits 
    mant_z |= ((z_bot | z_mid) + 0x01FFFFFF) >> 25; 

    // Now we need to normalize mant_z to be in the range [2^54, 2^55)
    uint64_t z_temp = (mant_z >> 1) | (mant_z & 1); // we need to check if we need to shift mant_z right by 1, which is the case if the top bit of mant_z is 1, in which case we also need to add the sticky bit to mant_z for rounding
    shift_flag = mant_z >> 55; // shift_flag is 1 if the top bit of mant_z is 1, and 0 otherwise
    mant_z ^= (mant_z ^ z_temp) & -shift_flag; // if shift_flag is 1, we shift mant_z right by 1 and add the sticky bit to mant_z for rounding

    /*
     * Now we can determine the exponent and sign of the result.
     * 
     * We first get the scaling factor: 
     * - the exponents are bised by 1023
     * - the mantissas are scalaed by 2^52, so add 52 bias for each exponent
     * - we right-shifted z by 50 bits to get mant_z, so we need to add 50 to the exponent bias
     * - and then the shift flag deterines whether we need to add another 1 or not
     * 
     * So the total bias = 2(1023+52) - 50 + shift_flag
    */

    exp_x = (int)(x >> 52) & 0x7FF; // extract the exponent of x
    exp_y = (int)(y >> 52) & 0x7FF; // extract the exponent of y
    exp_z = exp_x + exp_y - 2100 + (int)shift_flag; 

    sign_z = (int)((x^y) >> 63); // the sign of the result is the XOR of the signs of x and y

    // We've a problem if one of the inputs is 0, so we need a separate case to return 0
    int zero_flag = ((exp_x + 0x7FF) & (exp_y + 0x7FF)) >> 11; // zero_flag is 0 if neither x nor y is zero, and 1 if at least one of them is zero
    mant_z &= -(uint64_t)zero_flag; // if zero_flag is 1, we set mant_z to 0

    return flopt_create(sign_z, exp_z, mant_z); // create a flopt number with the calculated sign, exponent, and mantissa (it will deal with the 0 exponent)

}

/*
 * flopt_div is a function that divides two flopt numbers x and y, and returns the result as a flopt number.
 * The function first extracts the sign, exponent, and mantissa of x and y, and then it divides the mantissa of x by the mantissa of y using a method that calculates the quotient bit by bit to avoid overflow.
 * Finally, it normalizes the result and rounds it to fit into the flopt format before returning it.
 * 
 * This function has inputs:
 * - x: a flopt number (the dividend)
 * - y: a flopt number (the divisor)
 * 
 * and returns a flopt number that is the quotient of x and y.
*/
flopt flopt_div(flopt x, flopt y) {

    uint64_t mant_x, mant_y, mant_z;    // mantissa bits of x, y, and the result
    uint64_t quotient;                  // quotient of the division of the mantissas of x and y
    int exp_x, exp_y, exp_z;            // exponent bits of x, y, and the result
    int sign;                           // sign of the result

    int loop_counter;                   // counter for the loop to calculate the quotient
    uint64_t shift_flag;

    // Extract the mantissas as scaled uints
    mant_x = (x & (((uint64_t)1 << 52) -1)) | ((uint64_t)1 << 52); 
    mant_y = (y & (((uint64_t)1 << 52) -1)) | ((uint64_t)1 << 52); 


    // We need to do bit-by-bit division for 55 bits
    quotient = 0;
    for (loop_counter = 0; loop_counter < 55; loop_counter++) {
        int flag; 
        flag = ((mant_x - mant_y) >> 63) - 1; // flag is 0xFFFFFFFFFFFFFFFF if mant_x >= mant_y, and 0 otherwise
        mant_x -= mant_y & flag; // if mant_x >= mant_y, we subtract mant_y from mant_x
        quotient |= flag & 1; // if mant_x >= mant_y, we set the current bit of the quotient to 1
        mant_x <<= 1; // we need to shift mant_x left by 1 to get the next bit of the quotient
        quotient <<= 1; // we need to shift the quotient left by 1 to make room for the next bit
    }

    // After the for-loop, quotient is 55 bits, followed by an extra 0
    // We want that last bit to be sticky to represent if the remainder is non-zero for rounding purposes, so we need to add it to the quotient
    quotient |= (mant_x | -mant_x) >> 63; // we need to add the bits that will be shifted out of mant_x to the quotient for rounding

    uint64_t temp_q;
    temp_q = (quotient >> 1) | (quotient & 1); // we need to check if we need to shift the quotient right by 1, which is the case if the top bit of the quotient is 1, in which case we also need to add the sticky bit to the quotient for rounding

    shift_flag = quotient >> 55; // flag is 1 if the top bit of the quotient is 1, and 0 otherwise
    quotient ^= (quotient ^ temp_q) & -shift_flag; // if flag is 1, we shift the quotient right by 1 and add the sticky bit to the quotient for rounding


    // Now we can determine the exponent and sign of the result.
    exp_x = (int)(x >> 52) & 0x7FF; // extract the exponent of x
    exp_y = (int)(y >> 52) & 0x7FF; // extract the exponent of y
    exp_z = exp_x - exp_y - 55 + (int)shift_flag; // we need to subtract 55 from the exponent to account for the fact that we want the quotient to be in the range [2^54, 2^55)

    sign = (int)((x^y) >> 63); // the sign of the result is the XOR of the signs of x and y

    // We've a problem if the divisor is 0, so we need a separate case if x is 0
    int zero_flag = (exp_x + 0x7FF) >> 11; // zero_flag is 0 if x is not zero, and 1 if x is zero
    sign &= zero_flag; // if zero_flag is 1, we set the sign to 0 (since the result should be 0)
    exp_z &= -zero_flag; // if zero_flag is 1, we set the exponent to 0
    quotient &= -(uint64_t)zero_flag; // if zero_flag is 1, we set the quotient to 0

    return flopt_create(sign, exp_z, quotient); // create a flopt number with the calculated sign, exponent, and mantissa (it will deal with the 0 exponent)

}

/*
 * flopt_neg is a function that negates a flopt number x by flipping the sign bit of x, and returns the result as a flopt number.
 * This function has input:
 * - x: a flopt number
 * and returns a flopt number that is the negation of x.
*/
flopt flopt_neg(flopt x){
    x ^= (uint64_t)1 << 63; // negate the sign of x by flipping the MSB
    return x; // return the negated x
}

/*
 * flopt_double is a function that doubles a flopt number x by adding 1 to the exponent of x, and returns the result as a flopt number.
 * This function has input:
 * - x: a flopt number
 * and returns a flopt number that is the double of x.
*/
flopt flopt_double(flopt x) {
    int exp = (int)(x >> 52) & 0x7FF; // extract the exponent of x
    exp += 1; // we need to add 1 to the exponent to double the value
    return (x & ~(((uint64_t)0x7FF) << 52)) | ((uint64_t)exp << 52);  // we need to combine the new exponent with the original sign and mantissa of x to get the final result
}

/*
 * flopt_half is a function that halves a flopt number x by subtracting 1 from the exponent of x, and returns the result as a flopt number.
 * This function has input:
 * - x: a flopt number
 * and returns a flopt number that is the half of x.
*/
flopt flopt_half(flopt x) {
    int exp = (int)(x >> 52) & 0x7FF; // extract the exponent of x
    exp -= 1; // we need to subtract 1 from the exponent to half the value
    return (x & ~(((uint64_t)0x7FF) << 52)) | ((uint64_t)exp << 52);  // we need to combine the new exponent with the original sign and mantissa of x to get the final result
}

/*
 * flopt_inv is a function that calculates the multiplicative inverse of a flopt number x by dividing 1 by x using the flopt_div function, and returns the result as a flopt number.
 * This function has input:
 * - x: a flopt number
 * and returns a flopt number that is the multiplicative inverse of x.
*/
flopt flopt_inv(flopt x) {
    // To get the inverse of x, we can just divide 1 by x using the flopt_div function
    flopt one = flopt_create(0, 1076, (uint64_t)1 << 54); // create a flopt number that represents 1 (sign = 0, exponent = 1076 to account for the bias and scaling, and mantissa = 2^54 to account for the scaling)
    return flopt_div(one, x); // divide 1 by x using the flopt_div function to get the inverse of x
}