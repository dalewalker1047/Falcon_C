#ifndef restrict
#define restrict   __restrict
#endif

#ifndef FFT_H
#define FFT_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "flopt.c"

static const flopt flopt_zero = 0;

// Math for complex numbers that splits complex numbers into real and imaginary parts using our custom flopt type

//addition of two complex number(a+bi)+(c+di)=(a+c)+(b+d)i
void fft_add(flopt* out_re, flopt* out_im, flopt a_re, flopt a_im, flopt b_re, flopt b_im);

//subtraction of two complex number (a+bi)?(c+di)=(a?c)+(b?d)i
void fft_sub(flopt* out_re, flopt* out_im, flopt a_re, flopt a_im, flopt b_re, flopt b_im);

//complex number multiplication (a+bi)(c+di)=(ac?bd)+(ad+bc)i
void fft_mul(flopt* out_re, flopt* out_im, flopt a_re, flopt a_im, flopt b_re, flopt b_im);

//complex number squaring (a+bi)^2=(a^2?b^2)+2abi
void fft_sqr(flopt* out_re, flopt* out_im, flopt a_re, flopt a_im);

//complex number inversion 1/(a+bi) = (a-bi)/(a^2 + b^2)
void fft_inv(flopt* out_re, flopt* out_im, flopt a_re, flopt a_im);

//division of two complex number a/b = a (1/b)
void fft_div(flopt* out_re, flopt* out_im, flopt a_re, flopt a_im, flopt b_re, flopt b_im);

/*
 * FFT (falcon-fft.c).
 *
 * A real polynomial is represented as an array of N 'flopt' elements.
 * The FFT representation of a real polynomial contains N/2 complex
 * elements; each is stored as two real numbers, for the real and
 * imaginary parts, respectively. See falcon-fft.c for details on the
 * internal representation.
 * 
 * Compute FFT in-place: the source array should contain a real
 * polynomial (N coefficients); its storage area is reused to store
 * the FFT representation of that polynomial (N/2 complex numbers).
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */

void FFT(flopt* f, unsigned logn);

/*
 * Compute the inverse FFT in-place: the source array should contain the
 * FFT representation of a real polynomial (N/2 elements); the resulting
 * real polynomial (N coefficients of type 'flopt') is written over the
 * array.
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */
void iFFT(flopt* f, unsigned logn);

/*
 * Add polynomial b to polynomial a. a and b MUST NOT overlap. This
 * function works in both normal and FFT representations.
 */
void poly_add(flopt* restrict a, const flopt* restrict b, unsigned logn);

/*
 * Subtract polynomial b from polynomial a. a and b MUST NOT overlap. This
 * function works in both normal and FFT representations.
 */
void poly_sub(flopt* restrict a, const flopt* restrict b, unsigned logn);

/*
 * Negate polynomial a. This function works in both normal and FFT
 * representations.
 */
void poly_neg(flopt* a, unsigned logn);

/*
 * Compute adjoint of polynomial a. This function works only in FFT
 * representation.
 */
void poly_adj_fft(flopt* a, unsigned logn);

/*
 * Multiply polynomial a with polynomial b. a and b MUST NOT overlap.
 * This function works only in FFT representation.
 */
void poly_mul_fft(flopt* restrict a, const flopt* restrict b, unsigned logn);

/*
 * Multiply polynomial a with the adjoint of polynomial b. a and b MUST NOT
 * overlap. This function works only in FFT representation.
 */
void poly_muladj_fft(flopt* restrict a, const flopt* restrict b, unsigned logn);

/*
 * Multiply polynomial with its own adjoint. This function works only in FFT
 * representation.
 */
void poly_mulselfadj_fft(flopt* a, unsigned logn);

/*
 * Multiply polynomial with a real constant. This function works in both
 * normal and FFT representations.
 */
void poly_mulconst(flopt* a, flopt x, unsigned logn);

/*
 * Divide polynomial a by polynomial b, modulo X^N+1 (FFT representation).
 * a and b MUST NOT overlap.
 */
void poly_div_fft(flopt* restrict a, const flopt* restrict b, unsigned logn);

/*
 * Given f and g (in FFT representation), compute 1/(f*adj(f)+g*adj(g))
 * (also in FFT representation). Since the result is auto-adjoint, all its
 * coordinates in FFT representation are real; as such, only the first N/2
 * values of d[] are filled (the imaginary parts are skipped).
 *
 * Array d MUST NOT overlap with either a or b.
 */
void poly_invnorm2_fft(flopt* restrict d,
const flopt* restrict a, const flopt* restrict b, unsigned logn);

/*
 * Given F, G, f and g (in FFT representation), compute F*adj(f)+G*adj(g)
 * (also in FFT representation). Destination d MUST NOT overlap with
 * any of the source arrays.
 */
void poly_add_muladj_fft(flopt* restrict d,
	const flopt* restrict F, const flopt* restrict G,
	const flopt* restrict f, const flopt* restrict g, unsigned logn);

/*
 * Multiply polynomial a by polynomial b, where b is autoadjoint. Both
 * a and b are in FFT representation. Since b is autoadjoint, all its
 * FFT coefficients are real, and the array b contains only N/2 elements.
 * a and b MUST NOT overlap.
 */
void poly_mul_autoadj_fft(flopt* restrict a,
	const flopt* restrict b, unsigned logn);

/*
 * Divide polynomial a by polynomial b, where b is autoadjoint. Both
 * a and b are in FFT representation. Since b is autoadjoint, all its
 * FFT coefficients are real, and the array b contains only N/2 elements.
 * a and b MUST NOT overlap.
 */
void poly_div_autoadj_fft(flopt* restrict a,
	const flopt* restrict b, unsigned logn);

/*
 * Perform an LDL decomposition of an auto-adjoint matrix G, in FFT
 * representation. On input, g00, g01 and g11 are provided (where the
 * matrix G = [[g00, g01], [adj(g01), g11]]). On output, the d00, l10
 * and d11 values are written in g00, g01 and g11, respectively
 * (with D = [[d00, 0], [0, d11]] and L = [[1, 0], [l10, 1]]).
 * (In fact, d00 = g00, so the g00 operand is left unmodified.)
 */
void poly_LDL_fft(const flopt* restrict g00,
	flopt* restrict g01, flopt* restrict g11, unsigned logn);

/*
 * Perform an LDL decomposition of an auto-adjoint matrix G, in FFT
 * representation. This is identical to poly_LDL_fft() except that
 * g00, g01 and g11 are unmodified; the outputs d11 and l10 are written
 * in two other separate buffers provided as extra parameters.
 */
void poly_LDLmv_fft(flopt* restrict d11, flopt* restrict l10,
	const flopt* restrict g00, const flopt* restrict g01,
	const flopt* restrict g11, unsigned logn);

/*
 * Apply "split" operation on a polynomial in FFT representation:
 * f = f0(x^2) + x*f1(x^2), for half-size polynomials f0 and f1
 * (polynomials modulo X^(N/2)+1). f0, f1 and f MUST NOT overlap.
 */
void poly_split_fft(flopt* restrict f0, flopt* restrict f1,
	const flopt* restrict f, unsigned logn);

/*
 * Apply "merge" operation on two polynomials in FFT representation:
 * given f0 and f1, polynomials moduo X^(N/2)+1, this function computes
 * f = f0(x^2) + x*f1(x^2), in FFT representation modulo X^N+1.
 * f MUST NOT overlap with either f0 or f1.
 */
void poly_merge_fft(flopt* restrict f,
	const flopt* restrict f0, const flopt* restrict f1, unsigned logn);

#endif