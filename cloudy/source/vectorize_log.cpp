/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "vectorize.h"
#include "vectorize_math.h"
#include "vectorize_log_core.h"

//
// Written by Peter A.M. van Hoof, Royal Observatory of Belgium, Brussels
//
// this file contains vectorized versions of the single and double variants of the log()
// function. They are vectorized using AVX instructions, but also make use of AVX2, FMA,
// and AVX512 instructions when available. The basic algorithms for calculating the log()
// functions were slightly simplified from the openlibm library versions available at
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
// This work was also inspired by similar routines written by Giovanni Garberoglio which
// are available at: http://software-lisc.fbk.eu/avx_mathfun/
//

#ifdef __AVX__

#ifdef __AVX512F__

inline v8df v1logd(v8df x)
{
	__mmask8 invalid1 = _mm512_cmp_pd_mask(x, dbl_min, _CMP_LT_OQ);
	__mmask8 invalid2 = _mm512_cmp_pd_mask(x, dbl_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask8 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1logd", x, invalid) );
	}
	return  v1logd_core(x);
}

inline v8df v1log10d(v8df x)
{
	__mmask8 invalid1 = _mm512_cmp_pd_mask(x, dbl_min, _CMP_LT_OQ);
	__mmask8 invalid2 = _mm512_cmp_pd_mask(x, dbl_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask8 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1log10d", x, invalid) );
	}
	x = v1logd_core(x);
	return _mm512_mul_pd(x, log10e);
}

inline v8df v1log1pd(v8df x)
{
	__mmask8 invalid1 = _mm512_cmp_pd_mask(x, mone, _CMP_LE_OQ);
	__mmask8 invalid2 = _mm512_cmp_pd_mask(x, dbl_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask8 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1log1pd", x, invalid) );
	}
	return  v1log1pd_core(x);
}

inline v16sf v1logf(v16sf x)
{
	__mmask16 invalid1 = _mm512_cmp_ps_mask(x, flt_min, _CMP_LT_OQ);
	__mmask16 invalid2 = _mm512_cmp_ps_mask(x, flt_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask16 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1logf", x, invalid) );
	}
	return v1logf_core(x);
}

inline v16sf v1log10f(v16sf x)
{
	__mmask16 invalid1 = _mm512_cmp_ps_mask(x, flt_min, _CMP_LT_OQ);
	__mmask16 invalid2 = _mm512_cmp_ps_mask(x, flt_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask16 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1log10f", x, invalid) );
	}
	x = v1logf_core(x);
	return _mm512_mul_ps(x, log10ef);
}

inline v16sf v1log1pf(v16sf x)
{
	__mmask16 invalid1 = _mm512_cmp_ps_mask(x, monef, _CMP_LE_OQ);
	__mmask16 invalid2 = _mm512_cmp_ps_mask(x, flt_max, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid1, invalid2) )
	{
		__mmask16 invalid = invalid1 | invalid2;
		DOMAIN_ERROR( DEMsg("v1log1pf", x, invalid) );
	}
	return v1log1pf_core(x);
}

#else

inline v4df v1logd(v4df x)
{
	v4df invalid1 = _mm256_cmp_pd(x, dbl_min, _CMP_LT_OQ);
	v4df invalid2 = _mm256_cmp_pd(x, dbl_max, _CMP_NLE_UQ);
	v4df invalid = _mm256_or_pd(invalid1, invalid2);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1logd", x, invalid) );
	return v1logd_core(x);
}

inline v4df v1log10d(v4df x)
{
	v4df invalid1 = _mm256_cmp_pd(x, dbl_min, _CMP_LT_OQ);
	v4df invalid2 = _mm256_cmp_pd(x, dbl_max, _CMP_NLE_UQ);
	v4df invalid = _mm256_or_pd(invalid1, invalid2);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1log10d", x, invalid) );
	x = v1logd_core(x);
	return _mm256_mul_pd(x, log10e);
}

inline v4df v1log1pd(v4df x)
{
	v4df invalid1 = _mm256_cmp_pd(x, mone, _CMP_LE_OQ);
	v4df invalid2 = _mm256_cmp_pd(x, dbl_max, _CMP_NLE_UQ);
	v4df invalid = _mm256_or_pd(invalid1, invalid2);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1log1pd", x, invalid) );
	return v1log1pd_core(x);
}

inline v8sf v1logf(v8sf x)
{
	v8sf invalid1 = _mm256_cmp_ps(x, flt_min, _CMP_LT_OQ);
	v8sf invalid2 = _mm256_cmp_ps(x, flt_max, _CMP_NLE_UQ);
	v8sf invalid = _mm256_or_ps(invalid1, invalid2);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1logf", x, invalid) );
	return v1logf_core(x);
}

inline v8sf v1log10f(v8sf x)
{
	v8sf invalid1 = _mm256_cmp_ps(x, flt_min, _CMP_LT_OQ);
	v8sf invalid2 = _mm256_cmp_ps(x, flt_max, _CMP_NLE_UQ);
	v8sf invalid = _mm256_or_ps(invalid1, invalid2);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1log10f", x, invalid) );
	x = v1logf_core(x);
	return _mm256_mul_ps(x, log10ef);
}

inline v8sf v1log1pf(v8sf x)
{
	v8sf invalid1 = _mm256_cmp_ps(x, monef, _CMP_LE_OQ);
	v8sf invalid2 = _mm256_cmp_ps(x, flt_max, _CMP_NLE_UQ);
	v8sf invalid = _mm256_or_ps(invalid1, invalid2);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1log1pf", x, invalid) );
	return v1log1pf_core(x);
}

#endif // __AVX512F__

#else

// stub routines, should never be called
inline int v1logd(int) { return 0; }
inline int v1log10d(int) { return 0; }
inline int v1log1pd(int) { return 0; }
inline int v1logf(int) { return 0; }
inline int v1log10f(int) { return 0; }
inline int v1log1pf(int) { return 0; }

#endif // __AVX__

// wrapper routines to give math functions C++ linkage
// this prevents warnings from the Oracle Studio compiler
inline double wr_logd(double x)
{
	return log(x);
}

inline sys_float wr_logf(sys_float x)
{
	return logf(x);
}

inline double wr_log10d(double x)
{
	return log10(x);
}

inline sys_float wr_log10f(sys_float x)
{
	return log10f(x);
}

inline double wr_log1pd(double x)
{
	return log1p(x);
}

inline sys_float wr_log1pf(sys_float x)
{
	return log1pf(x);
}


void vlog(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vlog()" );

	vecfun( x, y, nlo, nhi, wr_logd, v1logd );
}

void vlog10(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vlog10()" );

	vecfun( x, y, nlo, nhi, wr_log10d, v1log10d );
}

void vlog1p(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vlog1p()" );

	vecfun( x, y, nlo, nhi, wr_log1pd, v1log1pd );
}

void vlog(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vlog()" );

	vecfun( x, y, nlo, nhi, wr_logf, v1logf );
}

void vlog10(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vlog10()" );

	vecfun( x, y, nlo, nhi, wr_log10f, v1log10f );
}

void vlog1p(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vlog1p()" );

	vecfun( x, y, nlo, nhi, wr_log1pf, v1log1pf );
}

void vlog(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(log, 1.);
}

void vlog10(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(log10, 1.);
}

void vlog1p(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(log1p, 1.);
}

void vlog(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(log, 1.);
}

void vlog10(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(log10, 1.);
}

void vlog1p(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(log1p, 1.);
}

void vlog(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(log, 1.f);
}

void vlog10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(log10, 1.f);
}

void vlog1p(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(log1p, 1.f);
}

void vlog(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7)
{
	V1FUN_PS_8(log, 1.f);
}

void vlog10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7)
{
	V1FUN_PS_8(log10, 1.f);
}

void vlog1p(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7)
{
	V1FUN_PS_8(log1p, 1.f);
}

void vlog(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	  sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(log, 1.f);
}

void vlog10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(log10, 1.f);
}

void vlog1p(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(log1p, 1.f);
}
