/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "vectorize.h"
#include "vectorize_math.h"
#include "vectorize_exp_core.h"

//
// Written by Peter A.M. van Hoof, Royal Observatory of Belgium, Brussels
//
// This file contains vectorized versions of the single and double variants of the exp()
// function. They are vectorized using AVX instructions, but also make use of AVX2, FMA,
// and AVX512 instructions when available. The basic algorithms for calculating the exp()
// functions were slightly modified from versions taken from the Cephes library written by
// Stephen L. Moshier, available at: http://www.netlib.org/cephes/
//
// The algorithms for calculating the expm1() functions are simplified versions of the
// openlibm routines available at http://openlibm.org/ which is subject to the following
// copyright:
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

inline v8df v1expd(v8df x)
{
	__mmask8 invalid = _mm512_cmp_pd_mask(x, exp_max_arg, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expd", x, invalid) );
	x = _mm512_max_pd(exp_min_arg, x);
	v8si m = v1expd_core(x);
	v8df p2n = v1pow2d_core(m);
	return _mm512_mul_pd(x, p2n);
}

inline v8df v1exp10d(v8df x)
{
	__mmask8 invalid = _mm512_cmp_pd_mask(x, exp10_max_arg, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1exp10d", x, invalid) );
	x = _mm512_max_pd(exp10_min_arg, x);
	v8df h = _mm512_mul_pd(x, third);
	v8df y = _mm512_roundscale_pd(h, _MM_FROUND_TRUNC);	
	v8df z = _mm512_mul_pd(y, three);
	x = _mm512_sub_pd(x, z);
	x = _mm512_mul_pd(x, ln_ten);
	x = _mm512_fmadd_pd(z, exp10_c1, x);
	v8si m1 = v1expd_core(x);
	y = _mm512_mul_pd(y, ten);
	v8si m2 = _mm512_cvtpd_epi32(y);
	v8si m = _mm256_add_epi32(m1, m2);
	m = _mm256_max_epi32(m, moff);
	v8df p2n = v1pow2d_core(m);
	return _mm512_mul_pd(x, p2n);
}

inline v8df v1expm1d(v8df x)
{
	__mmask8 invalid = _mm512_cmp_pd_mask(x, expm1_max_arg, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expm1d", x, invalid) );
	x = _mm512_max_pd(expm1_min_arg, x);
	return v1expm1d_core(x);
}

inline v16sf v1expf(v16sf x)
{
	__mmask16 invalid = _mm512_cmp_ps_mask(x, expf_max_arg, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expf", x, invalid) );
	x = _mm512_max_ps(expf_min_arg, x);
	v16si n = v1expf_core(x);
	v16sf p2n = v1pow2f_core(n);
	return _mm512_mul_ps(x, p2n);
}

inline v16sf v1exp10f(v16sf x)
{
	__mmask16 invalid = _mm512_cmp_ps_mask(x, exp10f_max_arg, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1exp10f", x, invalid) );
	x = _mm512_max_ps(exp10f_min_arg, x);
	v16sf h = _mm512_mul_ps(x, thirdf);
	v16sf y = _mm512_roundscale_ps(h, _MM_FROUND_TRUNC);	
	v16sf z = _mm512_mul_ps(y, threef);
	x = _mm512_sub_ps(x, z);
	x = _mm512_mul_ps(x, ln_tenf);
	x = _mm512_fmadd_ps(z, exp10f_c1, x);
	v16si m1 = v1expf_core(x);
	y = _mm512_mul_ps(y, tenf);
	v16si m2 = _mm512_cvtps_epi32(y);
	v16si m = _mm512_add_epi32(m1, m2);
	m = _mm512_max_epi32(m, mofff);
	v16sf p2n = v1pow2f_core(m);
	return _mm512_mul_ps(x, p2n);
}

inline v16sf v1expm1f(v16sf x)
{
	__mmask16 invalid = _mm512_cmp_ps_mask(x, expm1f_max_arg, _CMP_NLE_UQ);
	if( ! _mm512_kortestz(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expm1f", x, invalid) );
	x = _mm512_max_ps(expm1f_min_arg, x);
	return v1expm1f_core(x);
}

#else

inline v4df v1expd(v4df x)
{
	v4df invalid = _mm256_cmp_pd(x, exp_max_arg, _CMP_NLE_UQ);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expd", x, invalid) );
	x = _mm256_max_pd(exp_min_arg, x);
	v4si m = v1expd_core(x);
	v4df p2n = v1pow2d_core(m);
	return _mm256_mul_pd(x, p2n);
}

inline v4df v1exp10d(v4df x)
{
	v4df invalid = _mm256_cmp_pd(x, exp10_max_arg, _CMP_NLE_UQ);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1exp10d", x, invalid) );
	x = _mm256_max_pd(exp10_min_arg, x);
	v4df h = _mm256_mul_pd(x, third);
	v4df y = _mm256_round_pd(h, _MM_FROUND_TRUNC);	
	v4df z = _mm256_mul_pd(y, three);
	x = _mm256_sub_pd(x, z);
	x = _mm256_mul_pd(x, ln_ten);
#ifdef __FMA__
	x = _mm256_fmadd_pd(z, exp10_c1, x);
#else
	h = _mm256_mul_pd(z, exp10_c1);
	x = _mm256_add_pd(h, x);
#endif
	v4si m1 = v1expd_core(x);
	y = _mm256_mul_pd(y, ten);
	v4si m2 = _mm256_cvtpd_epi32(y);
	v4si m = _mm_add_epi32(m1, m2);
	m = _mm_max_epi32(m, moff);
	v4df p2n = v1pow2d_core(m);
	return _mm256_mul_pd(x, p2n);
}

inline v4df v1expm1d(v4df x)
{
	v4df invalid = _mm256_cmp_pd(x, expm1_max_arg, _CMP_NLE_UQ);
	if( ! _mm256_testz_pd(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expm1d", x, invalid) );
	x = _mm256_max_pd(expm1_min_arg, x);
	return v1expm1d_core(x);
}

inline v8sf v1expf(v8sf x)
{
	v8sf invalid = _mm256_cmp_ps(x, expf_max_arg, _CMP_NLE_UQ);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expf", x, invalid) );
	x = _mm256_max_ps(expf_min_arg, x);
	v8si n = v1expf_core(x);
	v8sf p2n = v1pow2f_core(n);
	return _mm256_mul_ps(x, p2n);
}

inline v8sf v1exp10f(v8sf x)
{
	v8sf invalid = _mm256_cmp_ps(x, exp10f_max_arg, _CMP_NLE_UQ);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1exp10f", x, invalid) );
	x = _mm256_max_ps(exp10f_min_arg, x);
	v8sf h = _mm256_mul_ps(x, thirdf);
	v8sf y = _mm256_round_ps(h, _MM_FROUND_TRUNC);	
	v8sf z = _mm256_mul_ps(y, threef);
	x = _mm256_sub_ps(x, z);
	x = _mm256_mul_ps(x, ln_tenf);
#ifdef __FMA__
	x = _mm256_fmadd_ps(z, exp10f_c1, x);
#else
	h = _mm256_mul_ps(z, exp10f_c1);
	x = _mm256_add_ps(h, x);
#endif
	v8si m1 = v1expf_core(x);
	y = _mm256_mul_ps(y, tenf);
	v8si m2 = _mm256_cvtps_epi32(y);
#ifdef __AVX2__
	v8si m = _mm256_add_epi32(m1, m2);
	m = _mm256_max_epi32(m, mofff);
	v8sf p2n = v1pow2f_core(m);
#else
	v4si m1l = _mm256_extractf128_si256(m1, 0);
	v4si m2l = _mm256_extractf128_si256(m2, 0);
	v4si m1u = _mm256_extractf128_si256(m1, 1);
	v4si m2u = _mm256_extractf128_si256(m2, 1);
	m1l = _mm_add_epi32(m1l, m2l);
	m1u = _mm_add_epi32(m1u, m2u);
	m1l = _mm_max_epi32(m1l, mofff);
	m1u = _mm_max_epi32(m1u, mofff);
	m1l = _mm_add_epi32(m1l, offf);
	m1u = _mm_add_epi32(m1u, offf);
	m1l = _mm_slli_epi32(m1l, 23);
	m1u = _mm_slli_epi32(m1u, 23);
	v8si m = _mm256_setzero_si256();
	m = _mm256_insertf128_si256(m, m1l, 0);
	m = _mm256_insertf128_si256(m, m1u, 1);
	v8sf p2n = _mm256_castsi256_ps(m);
#endif
	return _mm256_mul_ps(x, p2n);
}

inline v8sf v1expm1f(v8sf x)
{
	v8sf invalid = _mm256_cmp_ps(x, expm1f_max_arg, _CMP_NLE_UQ);
	if( ! _mm256_testz_ps(invalid, invalid) )
		DOMAIN_ERROR( DEMsg("v1expm1f", x, invalid) );
	x = _mm256_max_ps(expm1f_min_arg, x);
	return v1expm1f_core(x);
}

#endif // __AVX512F__

#else

// stub routines, should never be called
inline int v1expd(int) { return 0; }
inline int v1exp10d(int) { return 0; }
inline int v1expm1d(int) { return 0; }
inline int v1expf(int) { return 0; }
inline int v1exp10f(int) { return 0; }
inline int v1expm1f(int) { return 0; }

#endif // __AVX__

// wrapper routines to give math functions C++ linkage
// this prevents warnings from the Oracle Studio compiler
inline double wr_expd(double x)
{
	return exp(x);
}

inline double wr_expm1d(double x)
{
	return expm1(x);
}

inline sys_float wr_expf(sys_float x)
{
	return expf(x);
}

inline sys_float wr_expm1f(sys_float x)
{
	return expm1f(x);
}

// wrapper routines needed to resolve overloaded functions
// these prevent problems with the Oracle Studio compiler
// it cannot handle overloaded functions as a template argument
inline double wr_exp10d(double x)
{
	return exp10(x);
}

inline sys_float wr_exp10f(sys_float x)
{
	return exp10(x);
}


void vexp(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vexp()" );

	vecfun( x, y, nlo, nhi, wr_expd, v1expd );
}

void vexp10(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vexp10()" );

	vecfun( x, y, nlo, nhi, wr_exp10d, v1exp10d );
}

void vexpm1(const double x[], double y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vexpm1()" );

	vecfun( x, y, nlo, nhi, wr_expm1d, v1expm1d );
}

void vexp(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vexp()" );

	vecfun( x, y, nlo, nhi, wr_expf, v1expf );
}

void vexp10(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vexp10()" );

	vecfun( x, y, nlo, nhi, wr_exp10f, v1exp10f );
}

void vexpm1(const sys_float x[], sys_float y[], long nlo, long nhi)
{
	DEBUG_ENTRY( "vexpm1()" );

	vecfun( x, y, nlo, nhi, wr_expm1f, v1expm1f );
}

void vexp(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(exp, 0.);
}

void vexp10(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(exp10, 0.);
}

void vexpm1(double *y, double x0, double x1, double x2, double x3)
{
	V1FUN_PD_4(expm1, 0.);
}

void vexp(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(exp, 0.);
}

void vexp10(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(exp10, 0.);
}

void vexpm1(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	V1FUN_PD_8(expm1, 0.);
}

void vexp(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(exp, 0.f);
}

void vexp10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(exp10, 0.f);
}

void vexpm1(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3)
{
	V1FUN_PS_4(expm1, 0.f);
}

void vexp(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7)
{
	V1FUN_PS_8(exp, 0.f);
}

void vexp10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7)
{
	V1FUN_PS_8(exp10, 0.f);
}

void vexpm1(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7)
{
	V1FUN_PS_8(expm1, 0.f);
}

void vexp(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	  sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(exp, 0.f);
}

void vexp10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(exp10, 0.f);
}

void vexpm1(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15)
{
	V1FUN_PS_16(expm1, 0.f);
}
