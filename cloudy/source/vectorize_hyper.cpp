/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "vectorize.h"
#include "vectorize_math.h"
#include "vectorize_log_core.h"
#include "vectorize_hyper_core.h"

//
// Written by Peter A.M. van Hoof, Royal Observatory of Belgium, Brussels
//
// this file contains vectorized versions of the single and double variants of the asinh()
// function. They are vectorized using AVX instructions, but also make use of AVX2, FMA,
// and AVX512 instructions when available. The basic algorithms for calculating the asinh()
// functions were somewhat simplified from the openlibm library versions available at
// http://openlibm.org/ which is subject to the following copyright:
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice 
// is preserved.
// ====================================================
//
// an optimized version vfast_asinh() is also supplied which can only handle arguments in
// the domain [0,sqrt(FLT_MAX)] or [0,sqrt(DBL_MAX)].
//

#ifdef __AVX__

#ifdef __AVX512F__

inline v8df v1asinhd(v8df x)
{
	v8di ix = _mm512_castpd_si512(x);
	v8di is = _mm512_and_epi64(ix, asinh_mask2);
	ix = _mm512_and_epi64(ix, asinh_mask1);
	x = _mm512_castsi512_pd(ix);
	__mmask8 invalid = _mm512_cmp_pd_mask(x, dbl_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1asinhd", x, invalid) );
	v8df rh = v1logd_core(x);
	rh = _mm512_add_pd(rh, ln_two);
	v8df xl = _mm512_min_pd(x, asinh_2p28);
	v8df rl = v1asinhd_core(xl);
	v8df r = _mm512_setzero_pd();
	__mmask8 mask = _mm512_cmp_pd_mask(xl, x, _CMP_EQ_OQ);
	r = _mm512_mask_add_pd(rh, mask, rl, r);
	v8di ir = _mm512_castpd_si512(r);
	ir = _mm512_or_epi64(ir, is);
	return _mm512_castsi512_pd(ir);
}

inline v8df v1fast_asinhd(v8df x)
{
	__mmask8 invalid1 = _mm512_cmp_pd_mask(x, zero, _CMP_LT_OQ);
	__mmask8 invalid2 = _mm512_cmp_pd_mask(x, sqrt_dbl_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask8 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1fast_asinhd", x, invalid) );
	}
	return v1asinhd_core(x);
}

inline v16sf v1asinhf(v16sf x)
{
	v16si ix = _mm512_castps_si512(x);
	v16si is = _mm512_and_epi32(ix, asinh_mask2f);
	ix = _mm512_and_epi32(ix, asinh_mask1f);
	x = _mm512_castsi512_ps(ix);
	__mmask16 invalid = _mm512_cmp_ps_mask(x, flt_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1asinhf", x, invalid) );
	v16sf rh = v1logf_core(x);
	rh = _mm512_add_ps(rh, ln_twof);
	v16sf xl = _mm512_min_ps(x, asinhf_2p28);
	v16sf rl = v1asinhf_core(xl);
	v16sf r = _mm512_setzero_ps();
	__mmask16 mask = _mm512_cmp_ps_mask(xl, x, _CMP_EQ_OQ);
	r = _mm512_mask_add_ps(rh, mask, rl, r);
	v16si ir = _mm512_castps_si512(r);
	ir = _mm512_or_epi32(ir, is);
	return _mm512_castsi512_ps(ir);
}

inline v16sf v1fast_asinhf(v16sf x)
{
	__mmask16 invalid1 = _mm512_cmp_ps_mask(x, zerof, _CMP_LT_OQ);
	__mmask16 invalid2 = _mm512_cmp_ps_mask(x, sqrt_flt_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask16 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1fast_asinhf", x, invalid) );
	}
	return v1asinhf_core(x);
}

#else

inline v4df v1asinhd(v4df x)
{
	v4df signbit = x;
	v4df mask1 = _mm256_castsi256_pd(asinh_mask1);
	x = _mm256_and_pd(x, mask1);
	v4df invalid = _mm256_cmp_pd(x, dbl_max, _CMP_NLE_UQ);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1asinhd", x, invalid) );
	mask1 = _mm256_castsi256_pd(asinh_mask2);
	signbit = _mm256_and_pd(signbit, mask1);
	// use the approximation for large arguments: asinh(x) ~ log(x) + log(2)
	// this formula can be safely used up to x = DBL_MAX, but is not correct for x < 2^28
	v4df rh = v1logd_core(x);
	rh = _mm256_add_pd(rh, ln_two);
	// use the expression for small arguments: asinh(x) = log1p(x + x^2/(1 + sqrt(1+x^2)))
	// this formula is exact, but will overflow for x > sqrt(DBL_MAX), hence the use of xl
	v4df xl = _mm256_min_pd(x, asinh_2p28);
	v4df rl = v1asinhd_core(xl);
	v4df mask = _mm256_cmp_pd(xl, x, _CMP_EQ_OQ);
	v4df r = _mm256_blendv_pd(rh, rl, mask);
	return _mm256_or_pd(r, signbit);
}

inline v4df v1fast_asinhd(v4df x)
{
	v4df invalid1 = _mm256_cmp_pd(x, zero, _CMP_LT_OQ);
	v4df invalid2 = _mm256_cmp_pd(x, sqrt_dbl_max, _CMP_NLE_UQ);
	v4df invalid = _mm256_or_pd(invalid1, invalid2);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1fast_asinhd", x, invalid) );
	return v1asinhd_core(x);
}

inline v8sf v1asinhf(v8sf x)
{
	v8sf signbit = x;
	v8sf mask1 = _mm256_castsi256_ps(asinh_mask1f);
	x = _mm256_and_ps(x, mask1);
	v8sf invalid = _mm256_cmp_ps(x, flt_max, _CMP_NLE_UQ);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1asinhf", x, invalid) );
	mask1 = _mm256_castsi256_ps(asinh_mask2f);
	signbit = _mm256_and_ps(signbit, mask1);
	v8sf rh = v1logf_core(x);
	rh = _mm256_add_ps(rh, ln_twof);
	v8sf xl = _mm256_min_ps(x, asinhf_2p28);
	v8sf rl = v1asinhf_core(xl);
	v8sf mask = _mm256_cmp_ps(xl, x, _CMP_EQ_OQ);
	v8sf r = _mm256_blendv_ps(rh, rl, mask);
	return _mm256_or_ps(r, signbit);
}

inline v8sf v1fast_asinhf(v8sf x)
{
	v8sf invalid1 = _mm256_cmp_ps(x, zerof, _CMP_LT_OQ);
	v8sf invalid2 = _mm256_cmp_ps(x, sqrt_flt_max, _CMP_NLE_UQ);
	v8sf invalid = _mm256_or_ps(invalid1, invalid2);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1fast_asinhf", x, invalid) );
	return v1asinhf_core(x);
}

#endif // __AVX512F__

#else

// stub routines, should never be called
inline int v1asinhd(int) { return 0; }
inline int v1fast_asinhd(int) { return 0; }
inline int v1asinhf(int) { return 0; }
inline int v1fast_asinhf(int) { return 0; }

#endif // __AVX__

// wrapper routines to give math functions C++ linkage
// this prevents warnings from the Oracle Studio compiler
inline double wr_asinhd(double x)
{
	return asinh(x);
}

inline sys_float wr_asinhf(sys_float x)
{
	return asinhf(x);
}

// Fall-backs for non-AVX hardware -- not actually fast
inline double fast_asinh(double x)
{
	return asinh(x);
}

inline sys_float fast_asinhf(sys_float x)
{
	return asinhf(x);
}

void vasinh(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vasinh()" );

	vecfun( x, y, nlo, nhi, wr_asinhd, v1asinhd );
}

void vfast_asinh(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vfast_asinh()" );

	vecfun( x, y, nlo, nhi, wr_asinhd, v1fast_asinhd );
}

void vasinh(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vasinh()" );

	vecfun( x, y, nlo, nhi, wr_asinhf, v1asinhf );
}

void vfast_asinh(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vfast_asinh()" );

	vecfun( x, y, nlo, nhi, wr_asinhf, v1fast_asinhf );
}

void vasinh(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(asinh, 0.);
}

void vfast_asinh(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(fast_asinh, 0.);
}

void vasinh(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(asinh, 0.);
}

void vfast_asinh(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(fast_asinh, 0.);
}

void vasinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(asinh, 0.f);
}

void vfast_asinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(fast_asinh, 0.f);
}

void vasinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7)
{
	V1FUN_PS_8(asinh, 0.f);
}

void vfast_asinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
		 sys_float x6, sys_float x7)
{
	V1FUN_PS_8(fast_asinh, 0.f);
}

void vasinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(asinh, 0.f);
}

void vfast_asinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
		 sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
		 sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(fast_asinh, 0.f);
}
