/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "vectorize.h"
#include "vectorize_math.h"
#include "vectorize_sqrt_core.h"

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

#ifdef __AVX512F__

inline v8df v1sqrtd(v8df x)
{
	__mmask8 invalid1 = _mm512_cmp_pd_mask(x, zero, _CMP_LT_OQ);
	__mmask8 invalid2 = _mm512_cmp_pd_mask(x, dbl_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask8 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1sqrtd", x, invalid) );
	}
	return v1sqrtd_core(x);
}

inline v8df v1hypotd(v8df x, v8df y)
{
	v8di ix = _mm512_castpd_si512(x);
	v8di iy = _mm512_castpd_si512(y);
	ix = _mm512_and_si512(ix, sqrt_mask1);
	iy = _mm512_and_si512(iy, sqrt_mask1);
	x = _mm512_castsi512_pd(ix);
	__mmask8 invalid1 = _mm512_cmp_pd_mask(x, dbl_max, _CMP_NLE_UQ);
	y = _mm512_castsi512_pd(iy);
	__mmask8 invalid2 = _mm512_cmp_pd_mask(y, dbl_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
		DOMAIN_ERROR( DEMsg("v1hypotd", x, invalid1, y, invalid2) );
	return v1hypotd_core(x, y);
}

inline v16sf v1sqrtf(v16sf x)
{
	__mmask16 invalid1 = _mm512_cmp_ps_mask(x, zerof, _CMP_LT_OQ);
	__mmask16 invalid2 = _mm512_cmp_ps_mask(x, flt_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask16 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1sqrtf", x, invalid) );
	}
	return v1sqrtf_core(x);
}

inline v16sf v1hypotf(v16sf x, v16sf y)
{
	v16si ix = _mm512_castps_si512(x);
	v16si iy = _mm512_castps_si512(y);
	ix = _mm512_and_si512(ix, sqrt_mask1f);
	iy = _mm512_and_si512(iy, sqrt_mask1f);
	x = _mm512_castsi512_ps(ix);
	__mmask16 invalid1 = _mm512_cmp_ps_mask(x, flt_max, _CMP_NLE_UQ);
	y = _mm512_castsi512_ps(iy);
	__mmask16 invalid2 = _mm512_cmp_ps_mask(y, flt_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
		DOMAIN_ERROR( DEMsg("v1hypotf", x, invalid1, y, invalid2) );
	return v1hypotf_core(x, y);
}

#else

inline v4df v1sqrtd(v4df x)
{
	v4df invalid1 = _mm256_cmp_pd(x, zero, _CMP_LT_OQ);
	v4df invalid2 = _mm256_cmp_pd(x, dbl_max, _CMP_NLE_UQ);
	v4df invalid = _mm256_or_pd(invalid1, invalid2);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1sqrtd", x, invalid) );
	return v1sqrtd_core(x);
}

inline v4df v1hypotd(v4df x, v4df y)
{
	v4df mask1 = _mm256_castsi256_pd(sqrt_mask1);
	x = _mm256_and_pd(x, mask1);
	v4df invalid1 = _mm256_cmp_pd(x, dbl_max, _CMP_NLE_UQ);
	y = _mm256_and_pd(y, mask1);
	v4df invalid2 = _mm256_cmp_pd(y, dbl_max, _CMP_NLE_UQ);
	v4df invalid = _mm256_or_pd(invalid1, invalid2);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1hypotd", x, invalid1, y, invalid2) );
	return v1hypotd_core(x, y);
}

inline v8sf v1sqrtf(v8sf x)
{
	v8sf invalid1 = _mm256_cmp_ps(x, zerof, _CMP_LT_OQ);
	v8sf invalid2 = _mm256_cmp_ps(x, flt_max, _CMP_NLE_UQ);
	v8sf invalid = _mm256_or_ps(invalid1, invalid2);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1sqrtf", x, invalid) );
	return v1sqrtf_core(x);
}

inline v8sf v1hypotf(v8sf x, v8sf y)
{
	v8sf mask1 = _mm256_castsi256_ps(sqrt_mask1f);
	x = _mm256_and_ps(x, mask1);
	v8sf invalid1 = _mm256_cmp_ps(x, flt_max, _CMP_NLE_UQ);
	y = _mm256_and_ps(y, mask1);
	v8sf invalid2 = _mm256_cmp_ps(y, flt_max, _CMP_NLE_UQ);
	v8sf invalid = _mm256_or_ps(invalid1, invalid2);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1hypotf", x, invalid1, y, invalid2) );
	return v1hypotf_core(x, y);
}

#endif // __AVX512F__

#else

// stub routines, should never be called
inline int v1sqrtd(int) { return 0; }
inline int v1hypotd(int, int) { return 0; }
inline int v1sqrtf(int) { return 0; }
inline int v1hypotf(int, int) { return 0; }

#endif // __AVX__

// wrapper routines to give math functions C++ linkage
// this prevents warnings from the Oracle Studio compiler
inline double wr_sqrtd(double x)
{
	return sqrt(x);
}

inline double wr_hypotd(double x, double y)
{
	return hypot(x, y);
}

inline sys_float wr_sqrtf(sys_float x)
{
	return sqrtf(x);
}

inline sys_float wr_hypotf(sys_float x, sys_float y)
{
	return hypotf(x, y);
}

void vsqrt(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vsqrt()" );

	vecfun( x, y, nlo, nhi, wr_sqrtd, v1sqrtd );
}

void vhypot(const double x1[], const double x2[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vhypot()" );

	vecfun2( x1, x2, y, nlo, nhi, wr_hypotd, v1hypotd );
}

void vsqrt(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vsqrt()" );

	vecfun( x, y, nlo, nhi, wr_sqrtf, v1sqrtf );
}

void vhypot(const sys_float x1[], const sys_float x2[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vhypot()" );

	vecfun2( x1, x2, y, nlo, nhi, wr_hypotf, v1hypotf );
}

void vsqrt(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(sqrt, 1.);
}

void vhypot(double *z, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
	V1FUN2_PD_4(hypot, 1.);
}

void vsqrt(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(sqrt, 1.);
}

void vsqrt(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(sqrt, 1.f);
}

void vhypot(sys_float *z, sys_float x0, sys_float y0, sys_float x1, sys_float y1, sys_float x2, sys_float y2, 
	    sys_float x3, sys_float y3)
{
	V1FUN2_PS_4(hypot, 1.f);
}

void vsqrt(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	   sys_float x6, sys_float x7)
{
	V1FUN_PS_8(sqrt, 1.f);
}

void vhypot(sys_float *z, sys_float x0, sys_float y0, sys_float x1, sys_float y1, sys_float x2, sys_float y2, 
	    sys_float x3, sys_float y3, sys_float x4, sys_float y4, sys_float x5, sys_float y5, sys_float x6,
	    sys_float y6, sys_float x7, sys_float y7)
{
	V1FUN2_PS_8(hypot, 1.f);
}

void vsqrt(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	   sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	   sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(sqrt, 1.f);
}
