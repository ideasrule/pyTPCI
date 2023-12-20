/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_LOG_CORE_H
#define VECTORIZE_LOG_CORE_H

#include "vectorize_math.h"

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

#if defined(__AVX512F__) || defined(__AVX2__)
VECLL_CONST(log_off,0x3ff);
VECLL_CONST(log_mask2,0x000fffffffffffff);
VECLL_CONST(log_sq2off,0x95f6400000000);
VECLL_CONST(log_mask3,0x0010000000000000);
VECLL_CONST(log_mask4,0x3ff0000000000000);
#else
VECI_CONST(log_off,0x3ff);
VECI_CONST(log_mask2,0x000fffff);
VECI_CONST(log_sq2off,0x95f64);
VECI_CONST(log_mask3,0x00100000);
VECI_CONST(log_mask4,0x3ff00000);
#endif

VECDI_CONST(log_lg1,0xbfe5555555555593); // -6.666666666666735130e-01
VECDI_CONST(log_lg2,0xbfd999999997fa04); // -3.999999999940941908e-01
VECDI_CONST(log_lg3,0xbfd2492494229359); // -2.857142874366239149e-01
VECDI_CONST(log_lg4,0xbfcc71c51d8e78af); // -2.222219843214978396e-01
VECDI_CONST(log_lg5,0xbfc7466496cb03de); // -1.818357216161805012e-01
VECDI_CONST(log_lg6,0xbfc39a09d078c69f); // -1.531383769920937332e-01
VECDI_CONST(log_lg7,0xbfc2f112df3e5244); // -1.479819860511658591e-01

#ifdef __AVX512F__

inline v8df v1logd_core(v8df x)
{
	v8di ix = _mm512_castpd_si512(x);
	v8di k = _mm512_srli_epi64(ix, 52);
	k = _mm512_sub_epi64(k, log_off);
	v8di iy = _mm512_and_epi64(ix, log_mask2);
	v8di iz = _mm512_add_epi64(iy, log_sq2off);
	iz = _mm512_and_epi64(iz, log_mask3);
	v8di ie = _mm512_xor_epi64(iz, log_mask4);
	iy = _mm512_or_epi64(iy, ie);
	iz = _mm512_srli_epi64(iz, 52);
	k = _mm512_add_epi64(k, iz);
	v8df f = _mm512_castsi512_pd(iy);
	f = _mm512_sub_pd(f, one);
	v8df s = _mm512_add_pd(f, two);
	s = _mm512_div_pd(f, s);
#ifdef __AVX512DQ__
	v8df dk = _mm512_cvtepi64_pd(k);
#else
	v8si k2 = _mm512_cvtepi64_epi32(k);
	v8df dk = _mm512_cvtepi32_pd(k2);
#endif
	v8df z = _mm512_mul_pd(s, s);
	v8df R = _mm512_fmadd_pd(log_lg7, z, log_lg6);
	R = _mm512_fmadd_pd(R, z, log_lg5);
	R = _mm512_fmadd_pd(R, z, log_lg4);
	R = _mm512_fmadd_pd(R, z, log_lg3);
	R = _mm512_fmadd_pd(R, z, log_lg2);
	R = _mm512_fmadd_pd(R, z, log_lg1);
	R = _mm512_fmadd_pd(R, z, f);
	R = _mm512_fnmadd_pd(R, s, f);
	return _mm512_fmadd_pd(dk, ln_two, R);
}

inline v8df v1log1pd_core(v8df x)
{
	v8df xarg = x;
	x = _mm512_add_pd(x, one);
	v8df c = _mm512_sub_pd(x, one);
	c = _mm512_sub_pd(xarg, c);
	c = _mm512_div_pd(c, x);
	v8df y = v1logd_core(x);
	return _mm512_add_pd(y, c);
}

#else

inline v4df v1logd_core(v4df x)
{
#ifdef __AVX2__
	v4di ix = _mm256_castpd_si256(x);
	v4di k = _mm256_srli_epi64(ix, 52);
	k = _mm256_sub_epi64(k, log_off);
	v4di iy = _mm256_and_si256(ix, log_mask2);
	v4di iz = _mm256_add_epi64(iy, log_sq2off);
	iz = _mm256_and_si256(iz, log_mask3);
	v4di ie = _mm256_xor_si256(iz, log_mask4);
	iy = _mm256_or_si256(iy, ie);
	iz = _mm256_srli_epi64(iz, 52);
	k = _mm256_add_epi64(k, iz);
	v4di kp = _mm256_permutevar8x32_epi32(k, _mm256_set_epi32(7,5,3,1,6,4,2,0) );
	v4si kk = _mm256_extractf128_si256(kp, 0);
	x = _mm256_castsi256_pd(iy);
	v4df dk = _mm256_cvtepi32_pd(kk);
#else
	v8sf xs =  _mm256_castpd_ps(x);
	xs = _mm256_permute_ps(xs, _MM_SHUFFLE(3, 1, 2, 0));
	v4sf xls = _mm256_extractf128_ps(xs, 0);
	v4si xl = _mm_castps_si128(xls);
	v4sf xhs = _mm256_extractf128_ps(xs, 1);
	v4si xh = _mm_castps_si128(xhs);
	v4si hx = _mm_unpackhi_epi64(xl, xh);
	v4si lx = _mm_unpacklo_epi64(xl, xh);
	v4si k = _mm_srli_epi32(hx, 20);
	k = _mm_sub_epi32(k, log_off);
	hx = _mm_and_si128(hx, log_mask2);
	v4si i = _mm_add_epi32(hx, log_sq2off);
	i = _mm_and_si128(i, log_mask3);
	v4si ii = _mm_xor_si128(i, log_mask4);
	hx = _mm_or_si128(hx, ii);
	v8si xi = _mm256_setzero_si256();
	xl = _mm_unpacklo_epi32(lx,hx);
	xi = _mm256_insertf128_si256(xi, xl, 0);
	xh = _mm_unpackhi_epi32(lx,hx);
	xi = _mm256_insertf128_si256(xi, xh, 1);
	x = _mm256_castsi256_pd(xi);
	i = _mm_srli_epi32(i, 20);
	k = _mm_add_epi32(k, i);
	v4df dk = _mm256_cvtepi32_pd(k);
#endif
	// f and s have opposite sign compared to openlibm!
	v4df f = _mm256_sub_pd(one, x);
	v4df s = _mm256_sub_pd(two, f);
	s = _mm256_div_pd(f, s);
	v4df z = _mm256_mul_pd(s, s);
#ifdef __FMA__
	v4df R = _mm256_fmadd_pd(log_lg7, z, log_lg6);
	R = _mm256_fmadd_pd(R, z, log_lg5);
	R = _mm256_fmadd_pd(R, z, log_lg4);
	R = _mm256_fmadd_pd(R, z, log_lg3);
	R = _mm256_fmadd_pd(R, z, log_lg2);
	R = _mm256_fmadd_pd(R, z, log_lg1);
	R = _mm256_fmsub_pd(R, z, f);
	R = _mm256_fmsub_pd(R, s, f);
	v4df y = _mm256_fmadd_pd(dk, ln_two, R);
#else
	v4df R = _mm256_mul_pd(log_lg7, z);
	R = _mm256_add_pd(R, log_lg6);
	R = _mm256_mul_pd(R, z);
	R = _mm256_add_pd(R, log_lg5);
	R = _mm256_mul_pd(R, z);
	R = _mm256_add_pd(R, log_lg4);
	R = _mm256_mul_pd(R, z);
	R = _mm256_add_pd(R, log_lg3);
	R = _mm256_mul_pd(R, z);
	R = _mm256_add_pd(R, log_lg2);
	R = _mm256_mul_pd(R, z);
	R = _mm256_add_pd(R, log_lg1);
	R = _mm256_mul_pd(R, z);
	R = _mm256_sub_pd(R, f);
	R = _mm256_mul_pd(R, s);
	R = _mm256_sub_pd(R, f);
	v4df R2 = _mm256_mul_pd(dk, ln_two);
	v4df y = _mm256_add_pd(R, R2);
#endif
	return y;
}

inline v4df v1log1pd_core(v4df x)
{
	v4df xarg = x;
	x = _mm256_add_pd(x, one);
	v4df c = _mm256_sub_pd(x, one);
	c = _mm256_sub_pd(xarg, c);
	c = _mm256_div_pd(c, x);
	v4df y = v1logd_core(x);
	return _mm256_add_pd(y, c);
}

#endif // __AVX512F__

#if defined(__AVX512F__) || defined(__AVX2__)
VECII_CONST(log_offf,0x7f);
VECII_CONST(log_mask2f,0x007fffff);
VECII_CONST(log_sq2offf,(0x95f64<<3));
VECII_CONST(log_mask3f,0x00800000);
VECII_CONST(log_mask4f,0x3f800000);
#else
VECI_CONST(log_offf,0x7f);
VECI_CONST(log_mask2f,0x007fffff);
VECI_CONST(log_sq2offf,(0x95f64<<3));
VECI_CONST(log_mask3f,0x00800000);
VECI_CONST(log_mask4f,0x3f800000);
#endif

VECFI_CONST(log_lg1f,0xbf2aaaaa); // -0.66666662693f
VECFI_CONST(log_lg2f,0xbeccce13); // -0.40000972152f
VECFI_CONST(log_lg3f,0xbe91e9ee); // -0.28498786688f
VECFI_CONST(log_lg4f,0xbe789e26); // -0.24279078841f

#ifdef __AVX512F__

inline v16sf v1logf_core(v16sf x)
{
	v16si ix = _mm512_castps_si512(x);
	v16si k = _mm512_srli_epi32(ix, 23);
	k = _mm512_sub_epi32(k, log_offf);
	v16si iy = _mm512_and_epi32(ix, log_mask2f);
	v16si iz = _mm512_add_epi32(iy, log_sq2offf);
	iz = _mm512_and_epi32(iz, log_mask3f);
	v16si ie = _mm512_xor_epi32(iz, log_mask4f);
	iy = _mm512_or_epi32(iy, ie);
	iz = _mm512_srli_epi32(iz, 23);
	k = _mm512_add_epi32(k, iz);
	v16sf f = _mm512_castsi512_ps(iy);
	f = _mm512_sub_ps(f, onef);
	v16sf s = _mm512_add_ps(f, twof);
	s = _mm512_div_ps(f, s);
	v16sf dk = _mm512_cvtepi32_ps(k);
	v16sf z = _mm512_mul_ps(s, s);
	v16sf R = _mm512_fmadd_ps(log_lg4f, z, log_lg3f);
	R = _mm512_fmadd_ps(R, z, log_lg2f);
	R = _mm512_fmadd_ps(R, z, log_lg1f);
	R = _mm512_fmadd_ps(R, z, f);
	R = _mm512_fnmadd_ps(R, s, f);
	return _mm512_fmadd_ps(dk, ln_twof, R);
}

inline v16sf v1log1pf_core(v16sf x)
{
	v16sf xarg = x;
	x = _mm512_add_ps(x, onef);
	v16sf c = _mm512_sub_ps(x, onef);
	c = _mm512_sub_ps(xarg, c);
	c = _mm512_div_ps(c, x);
	v16sf y = v1logf_core(x);
	return _mm512_add_ps(y, c);
}

#else

inline v8sf v1logf_core(v8sf x)
{
	v8si ix = _mm256_castps_si256(x);
#ifdef __AVX2__
	v8si k = _mm256_srli_epi32(ix, 23);
	k = _mm256_sub_epi32(k, log_offf);
	v8si iy = _mm256_and_si256(ix, log_mask2f);
	v8si iz = _mm256_add_epi32(iy, log_sq2offf);
	iz = _mm256_and_si256(iz, log_mask3f);
	v8si ie = _mm256_xor_si256(iz, log_mask4f);
	iy = _mm256_or_si256(iy, ie);
	iz = _mm256_srli_epi32(iz, 23);
	k = _mm256_add_epi32(k, iz);
#else
	v4si ixl = _mm256_extractf128_si256(ix, 0);
	v4si kl = _mm_srli_epi32(ixl, 23);
	kl = _mm_sub_epi32(kl, log_offf);
	v4si iyl = _mm_and_si128(ixl, log_mask2f);
	v4si iz = _mm_add_epi32(iyl, log_sq2offf);
	iz = _mm_and_si128(iz, log_mask3f);
	v4si ie = _mm_xor_si128(iz, log_mask4f);
	iyl = _mm_or_si128(iyl, ie);
	iz = _mm_srli_epi32(iz, 23);
	kl = _mm_add_epi32(kl, iz);
	v4si ixh = _mm256_extractf128_si256(ix, 1);
	v4si kh = _mm_srli_epi32(ixh, 23);
	kh = _mm_sub_epi32(kh, log_offf);
	v4si iyh = _mm_and_si128(ixh, log_mask2f);
	iz = _mm_add_epi32(iyh, log_sq2offf);
	iz = _mm_and_si128(iz, log_mask3f);
	ie = _mm_xor_si128(iz, log_mask4f);
	iyh = _mm_or_si128(iyh, ie);
	iz = _mm_srli_epi32(iz, 23);
	kh = _mm_add_epi32(kh, iz);
	v8si iy = _mm256_setzero_si256();
	iy = _mm256_insertf128_si256(iy, iyl, 0);
	iy = _mm256_insertf128_si256(iy, iyh, 1);
	v8si k = _mm256_setzero_si256();
	k = _mm256_insertf128_si256(k, kl, 0);
	k = _mm256_insertf128_si256(k, kh, 1);
#endif
	x = _mm256_castsi256_ps(iy);
	v8sf dk = _mm256_cvtepi32_ps(k);
	// f and s have opposite sign compared to openlibm!
	v8sf f = _mm256_sub_ps(onef, x);
	v8sf s = _mm256_sub_ps(twof, f);
	s = _mm256_div_ps(f, s);
	v8sf z = _mm256_mul_ps(s, s);
#ifdef __FMA__
	v8sf R = _mm256_fmadd_ps(log_lg4f, z, log_lg3f);
	R = _mm256_fmadd_ps(R, z, log_lg2f);
	R = _mm256_fmadd_ps(R, z, log_lg1f);
	R = _mm256_fmsub_ps(R, z, f);
	R = _mm256_fmsub_ps(R, s, f);
	v8sf y = _mm256_fmadd_ps(dk, ln_twof, R);
#else
	v8sf R = _mm256_mul_ps(log_lg4f, z);
	R = _mm256_add_ps(R, log_lg3f);
	R = _mm256_mul_ps(R, z);
	R = _mm256_add_ps(R, log_lg2f);
	R = _mm256_mul_ps(R, z);
	R = _mm256_add_ps(R, log_lg1f);
	R = _mm256_mul_ps(R, z);
	R = _mm256_sub_ps(R, f);
	R = _mm256_mul_ps(R, s);
	R = _mm256_sub_ps(R, f);
	v8sf R2 = _mm256_mul_ps(dk, ln_twof);
	v8sf y = _mm256_add_ps(R, R2);
#endif
	return y;
}

inline v8sf v1log1pf_core(v8sf x)
{
	v8sf xarg = x;
	x = _mm256_add_ps(x, onef);
	v8sf c = _mm256_sub_ps(x, onef);
	c = _mm256_sub_ps(xarg, c);
	c = _mm256_div_ps(c, x);
	v8sf y = v1logf_core(x);
	return _mm256_add_ps(y, c);
}

#endif // __AVX512F__

#endif // __AVX__

#endif
