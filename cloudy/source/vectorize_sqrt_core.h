/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_SQRT_CORE_H
#define VECTORIZE_SQRT_CORE_H

#include "vectorize_math.h"

//
// Written by Peter A.M. van Hoof, Royal Observatory of Belgium, Brussels
//
// this file contains vectorized versions of the single and double variants of the sqrt()
// and hypot() functions. They are vectorized using AVX instructions, but also make use of
// AVX2, FMA, and AVX512 instructions when available. The basic algorithms for calculating
// the sqrt() functions were derived from the algorithm for calculating rsqrt() described
// here: http://en.wikipedia.org/wiki/Fast_inverse_square_root
//
// Alternatively one can also use the sqrt hardware instruction, but on some hardware the
// the software implementation is faster... The hardware instruction is chosen as the
// default implementation below.
//

#ifdef __AVX__

VECLL_CONST(sqrt_mask1,0x7fffffffffffffff);

#ifdef __AVX2__
VECLL_CONST(sqrt_magic,0x5fe6eb50c7b537a9);
#else
VECL_CONST(sqrt_magic,0x5fe6eb50c7b537a9);
#endif

#ifdef __AVX512F__

inline v8df v1sqrtd_core(v8df x)
{
#if 0
	v8di ir = _mm512_castpd_si512(x);
	ir = _mm512_srli_epi64(ir, 1);
	ir = _mm512_sub_epi64(sqrt_magic, ir);
	v8df r = _mm512_castsi512_pd(ir);
	__mmask8 zmask = _mm512_cmp_pd_mask(x, dbl_min, _CMP_LT_OQ);
	r = _mm512_mask_load_pd(r, zmask, &zero);
	// do not precompute x/2. to avoid underflow to denormalized numbers
	// this may be flushed to zero, but is also very slow....
	v8df rx = _mm512_mul_pd(r, x);
	v8df rx2 = _mm512_mul_pd(rx, mhalf);
	v8df y = _mm512_fmadd_pd(rx2, r, c1p5);
	r = _mm512_mul_pd(r, y);
	rx = _mm512_mul_pd(r, x);
	rx2 = _mm512_mul_pd(rx, mhalf);
	y = _mm512_fmadd_pd(rx2, r, c1p5);
	r = _mm512_mul_pd(r, y);
	rx = _mm512_mul_pd(r, x);
	rx2 = _mm512_mul_pd(rx, mhalf);
	y = _mm512_fmadd_pd(rx2, r, c1p5);
	r = _mm512_mul_pd(r, y);
	rx = _mm512_mul_pd(r, x);
	rx2 = _mm512_mul_pd(rx, mhalf);
	y = _mm512_fmadd_pd(rx2, r, c1p5);
	r = _mm512_mul_pd(r, y);
	return _mm512_mul_pd(x, r);
#else
	return _mm512_sqrt_pd(x);
#endif
}

inline v8df v1hypotd_core(v8df x, v8df y)
{
	v8df xp = _mm512_max_pd(x, y);
	v8df yp = _mm512_min_pd(x, y);
	__mmask8 zmask = _mm512_cmp_pd_mask(xp, zero, _CMP_NEQ_OQ);
	v8df arg = _mm512_mask_div_pd(zero, zmask, yp, xp);
	arg = _mm512_fmadd_pd(arg, arg, one);
	v8df s = v1sqrtd_core(arg);
	return _mm512_mul_pd(xp, s);
}

#else

inline v4df v1sqrtd_core(v4df x)
{
#if 0
#ifdef __AVX2__
	v4di ir = _mm256_castpd_si256(x);
	ir = _mm256_srli_epi64(ir, 1);
	ir = _mm256_sub_epi64(sqrt_magic, ir);
#else
	v2df xl = _mm256_extractf128_pd(x, 0);
	v2df xh = _mm256_extractf128_pd(x, 1);
	v2di ixl = _mm_castpd_si128(xl);
	v2di ixh = _mm_castpd_si128(xh);
	ixl = _mm_srli_epi64(ixl, 1);
	ixh = _mm_srli_epi64(ixh, 1);
	ixl = _mm_sub_epi64(sqrt_magic, ixl);
	ixh = _mm_sub_epi64(sqrt_magic, ixh);
	v4di ir = _mm256_setzero_si256();
	ir = _mm256_insertf128_si256(ir, ixl, 0);
	ir = _mm256_insertf128_si256(ir, ixh, 1);
#endif
	v4df r = _mm256_castsi256_pd(ir);
	v4df zmask = _mm256_cmp_pd(x, dbl_min, _CMP_LT_OQ);
	r = _mm256_blendv_pd(r, zero, zmask);
	// do not precompute x/2. to avoid underflow to denormalized numbers
	// this may be flushed to zero, but is also very slow....
	v4df rx = _mm256_mul_pd(r, x);
	v4df rx2 = _mm256_mul_pd(rx, mhalf);
#ifdef __FMA__
	v4df y = _mm256_fmadd_pd(rx2, r, c1p5);
#else
	v4df y = _mm256_mul_pd(rx2, r);
	y = _mm256_add_pd(y, c1p5);
#endif
	r = _mm256_mul_pd(r, y);
	rx = _mm256_mul_pd(r, x);
	rx2 = _mm256_mul_pd(rx, mhalf);
#ifdef __FMA__
	y = _mm256_fmadd_pd(rx2, r, c1p5);
#else
	y = _mm256_mul_pd(rx2, r);
	y = _mm256_add_pd(y, c1p5);
#endif
	r = _mm256_mul_pd(r, y);
	rx = _mm256_mul_pd(r, x);
	rx2 = _mm256_mul_pd(rx, mhalf);
#ifdef __FMA__
	y = _mm256_fmadd_pd(rx2, r, c1p5);
#else
	y = _mm256_mul_pd(rx2, r);
	y = _mm256_add_pd(y, c1p5);
#endif
	r = _mm256_mul_pd(r, y);
	rx = _mm256_mul_pd(r, x);
	rx2 = _mm256_mul_pd(rx, mhalf);
#ifdef __FMA__
	y = _mm256_fmadd_pd(rx2, r, c1p5);
#else
	y = _mm256_mul_pd(rx2, r);
	y = _mm256_add_pd(y, c1p5);
#endif
	r = _mm256_mul_pd(r, y);
	return _mm256_mul_pd(x, r);
#else
	return _mm256_sqrt_pd(x);
#endif
}

inline v4df v1hypotd_core(v4df x, v4df y)
{
	v4df xp = _mm256_max_pd(x, y);
	v4df yp = _mm256_min_pd(x, y);
	v4df zmask = _mm256_cmp_pd(xp, zero, _CMP_EQ_OQ);
	xp = _mm256_blendv_pd(xp, one, zmask);
	v4df arg = _mm256_div_pd(yp, xp);
	xp = _mm256_blendv_pd(xp, zero, zmask);
#ifdef __FMA__
	arg = _mm256_fmadd_pd(arg, arg, one);
#else
	arg = _mm256_mul_pd(arg, arg);
	arg = _mm256_add_pd(arg, one);
#endif
	v4df s = v1sqrtd_core(arg);
	return _mm256_mul_pd(xp, s);
}

#endif // __AVX512F__

VECII_CONST(sqrt_mask1f,0x7fffffff);

#ifdef __AVX2__
VECII_CONST(sqrt_magicf,0x5f375a86);
#else
VECI_CONST(sqrt_magicf,0x5f375a86);
#endif

#ifdef __AVX512F__

inline v16sf v1sqrtf_core(v16sf x)
{
#if 0
	v16si ir = _mm512_castps_si512(x);
	ir = _mm512_srli_epi32(ir, 1);
	ir = _mm512_sub_epi32(sqrt_magicf,ir);
	v16sf r = _mm512_castsi512_ps(ir);
	__mmask16 zmask = _mm512_cmp_ps_mask(x, flt_min, _CMP_LT_OS);
	r = _mm512_mask_load_ps(r, zmask, &zerof);
	// do not precompute x/2.f to avoid underflow to denormalized numbers
	// this may be flushed to zero, but is also very slow....
	v16sf rx = _mm512_mul_ps(r, x);
	v16sf rx2 = _mm512_mul_ps(rx, mhalff);
	v16sf y = _mm512_fmadd_ps(rx2, r, c1p5f);
	r = _mm512_mul_ps(r, y);
	rx = _mm512_mul_ps(r, x);
	rx2 = _mm512_mul_ps(rx, mhalff);
	y = _mm512_fmadd_ps(rx2, r, c1p5f);
	r = _mm512_mul_ps(r, y);
	rx = _mm512_mul_ps(r, x);
	rx2 = _mm512_mul_ps(rx, mhalff);
	y = _mm512_fmadd_ps(rx2, r, c1p5f);
	r = _mm512_mul_ps(r, y);
	return _mm512_mul_ps(x, r);
#else
	return _mm512_sqrt_ps(x);
#endif
}

inline v16sf v1hypotf_core(v16sf x, v16sf y)
{
	v16sf xp = _mm512_max_ps(x, y);
	v16sf yp = _mm512_min_ps(x, y);
	__mmask16 zmask = _mm512_cmp_ps_mask(xp, zerof, _CMP_NEQ_OQ);
	v16sf arg = _mm512_mask_div_ps(zerof, zmask, yp, xp);
	arg = _mm512_fmadd_ps(arg, arg, onef);
	v16sf s = v1sqrtf_core(arg);
	return _mm512_mul_ps(xp, s);
}

#else

inline v8sf v1sqrtf_core(v8sf x)
{
#if 0
#ifdef __AVX2__
	v8si ir = _mm256_castps_si256(x);
	ir = _mm256_srli_epi32(ir, 1);
	ir = _mm256_sub_epi32(sqrt_magicf,ir);
#else
	v4sf xl = _mm256_extractf128_ps(x, 0);
	v4sf xh = _mm256_extractf128_ps(x, 1);
	v4si ixl = _mm_castps_si128(xl);
	v4si ixh = _mm_castps_si128(xh);
	ixl = _mm_srli_epi32(ixl, 1);
	ixh = _mm_srli_epi32(ixh, 1);
	ixl = _mm_sub_epi32(sqrt_magicf,ixl);
	ixh = _mm_sub_epi32(sqrt_magicf,ixh);
	v8si ir = _mm256_setzero_si256();
	ir = _mm256_insertf128_si256(ir, ixl, 0);
	ir = _mm256_insertf128_si256(ir, ixh, 1);
#endif
	v8sf r = _mm256_castsi256_ps(ir);
	v8sf zmask = _mm256_cmp_ps(x, flt_min, _CMP_LT_OQ);
	r = _mm256_blendv_ps(r, zerof, zmask);
	// do not precompute x/2.f to avoid underflow to denormalized numbers
	// this may be flushed to zero, but is also very slow....
	v8sf rx = _mm256_mul_ps(r, x);
	v8sf rx2 = _mm256_mul_ps(rx, mhalff);
#ifdef __FMA__
	v8sf y = _mm256_fmadd_ps(rx2, r, c1p5f);
#else
	v8sf y = _mm256_mul_ps(rx2, r);
	y = _mm256_add_ps(y, c1p5f);
#endif
	r = _mm256_mul_ps(r, y);
	rx = _mm256_mul_ps(r, x);
	rx2 = _mm256_mul_ps(rx, mhalff);
#ifdef __FMA__
	y = _mm256_fmadd_ps(rx2, r, c1p5f);
#else
	y = _mm256_mul_ps(rx2, r);
	y = _mm256_add_ps(y, c1p5f);
#endif
	r = _mm256_mul_ps(r, y);
	rx = _mm256_mul_ps(r, x);
	rx2 = _mm256_mul_ps(rx, mhalff);
#ifdef __FMA__
	y = _mm256_fmadd_ps(rx2, r, c1p5f);
#else
	y = _mm256_mul_ps(rx2, r);
	y = _mm256_add_ps(y, c1p5f);
#endif
	r = _mm256_mul_ps(r, y);
	return _mm256_mul_ps(x, r);
#else
	return _mm256_sqrt_ps(x);
#endif
}

inline v8sf v1hypotf_core(v8sf x, v8sf y)
{
	v8sf xp = _mm256_max_ps(x, y);
	v8sf yp = _mm256_min_ps(x, y);
	v8sf zmask = _mm256_cmp_ps(xp, zerof, _CMP_EQ_OQ);
	xp = _mm256_blendv_ps(xp, onef, zmask);
	v8sf arg = _mm256_div_ps(yp, xp);
	xp = _mm256_blendv_ps(xp, zerof, zmask);
#ifdef __FMA__
	arg = _mm256_fmadd_ps(arg, arg, onef);
#else
	arg = _mm256_mul_ps(arg, arg);
	arg = _mm256_add_ps(arg, onef);
#endif
	v8sf s = v1sqrtf_core(arg);
	return _mm256_mul_ps(xp, s);
}

#endif // __AVX512F__

#endif // __AVX__

#endif
