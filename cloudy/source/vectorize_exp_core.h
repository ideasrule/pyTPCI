/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_EXP_CORE_H
#define VECTORIZE_EXP_CORE_H

#include "vectorize_math.h"

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

VECDI_CONST(exp_max_arg,0x40862e42fefa39ee); // +709.78271289338385941
VECDI_CONST(exp_min_arg,0xc086232bdd7abcd3); // -708.39641853226419244

VECDI_CONST(exp10_max_arg,0x40734413509f79fe); // +308.25471555991668993
VECDI_CONST(exp10_min_arg,0xc0733a7146f72a42); // -307.65265556858878426

VECDI_CONST(expm1_max_arg,0x40862e42fefa39ef); // +709.78271289338397310
VECDI_CONST(expm1_min_arg,0xc042b708872320e2); // -37.429947750237047899

VECI_CONST(off,0x3ff);
VECI_CONST(moff,0xfffffc01);

VECDI_CONST(exp_c1,0xbfe62e4000000000); // -6.9314575195312500000e-01
VECDI_CONST(exp_c2,0xbeb7f7d1cf79abca); // -1.4286068203094172555e-06

VECDI_CONST(exp_p0,0x3f2089cdd5e44be8); // 1.2617719307481058306e-04
VECDI_CONST(exp_p1,0x3f9f06d10cca2c7e); // 3.0299440770744194562e-02
VECDI_CONST(exp_p2,0x3ff0000000000000); // 1.0000000000000000000e+00
VECDI_CONST(exp_q0,0x3ec92eb6bc365fa0); // 3.0019850513866445972e-06
VECDI_CONST(exp_q1,0x3f64ae39b508b6c0); // 2.5244834034968410830e-03
VECDI_CONST(exp_q2,0x3fcd17099887e074); // 2.2726554820815503266e-01
VECDI_CONST(exp_q3,0x4000000000000000); // 2.000000000000000000000e0

VECDI_CONST(exp10_c1,0xbf8030c37085dc7f); // -7.90550887243868070612e-3

VECDI_CONST(expm1_q1,0xbfa11111111110f4); // -3.33333333333331316428e-02
VECDI_CONST(expm1_q2,0x3f5a01a019fe5585); //  1.58730158725481460165e-03
VECDI_CONST(expm1_q3,0xbf14ce199eaadbb7); // -7.93650757867487942473e-05
VECDI_CONST(expm1_q4,0x3ed0cfca86e65239); //  4.00821782732936239552e-06
VECDI_CONST(expm1_q5,0xbe8afdb76e09c32d); // -2.01099218183624371326e-07

#ifdef __AVX512F__

inline v8df v1pow2d_core(v8si m)
{
	m = _mm256_add_epi32(m, off);
	v8di n = _mm512_cvtepi32_epi64(m);
	n = _mm512_slli_epi64(n, 52);
	return _mm512_castsi512_pd(n);
}

inline v8si v1expd_core(v8df& x)
{
	v8df px = _mm512_mul_pd(x, log2e);
	px = _mm512_floor_pd(px);
	v8si m = _mm512_cvtpd_epi32(px);
	x = _mm512_fmadd_pd(px, exp_c1, x);
	x = _mm512_fmadd_pd(px, exp_c2, x);
	v8df xx = _mm512_mul_pd(x, x);
	px = _mm512_fmadd_pd(exp_p0, xx, exp_p1);
	px = _mm512_fmadd_pd(px, xx, exp_p2);
	v8df qx = _mm512_fmadd_pd(exp_q0, xx, exp_q1);
	qx = _mm512_fmadd_pd(qx, xx, exp_q2);
	qx = _mm512_fmadd_pd(qx, xx, exp_q3);
	px = _mm512_mul_pd(x, px);
	xx = _mm512_sub_pd(qx, px);
	x = _mm512_div_pd(px, xx);
	x = _mm512_fmadd_pd(two, x, one);
	return m;
}

inline v8df v1expm1d_core(v8df x)
{
	v8df px = _mm512_mul_pd(x, log2e);
	px = _mm512_roundscale_pd(px, _MM_FROUND_NINT);
	v8df py = _mm512_min_pd(px, c1023);
	v8si k = _mm512_cvtpd_epi32(py);
	py = _mm512_sub_pd(px, py);
	py = _mm512_add_pd(py, one);
	x = _mm512_fmadd_pd(px, exp_c1, x);
	x = _mm512_fmadd_pd(px, exp_c2, x);
	v8df hfx = _mm512_mul_pd(x, half);
	v8df hxs = _mm512_mul_pd(x, hfx);
	// sign of t is reversed compared to openlibm implementation
	v8df r1 = _mm512_fmadd_pd(expm1_q5, hxs, expm1_q4);
	r1 = _mm512_fmadd_pd(r1, hxs, expm1_q3);
	r1 = _mm512_fmadd_pd(r1, hxs, expm1_q2);
	r1 = _mm512_fmadd_pd(r1, hxs, expm1_q1);
	r1 = _mm512_fmadd_pd(r1, hxs, one);
	v8df t = _mm512_fmadd_pd(r1, hfx, mthree);
	v8df e = _mm512_fmadd_pd(x, t, six);
	v8df h = _mm512_add_pd(r1, t);
	h = _mm512_div_pd(h, e);
	e = _mm512_mul_pd(h, hxs);
	v8df p2n = v1pow2d_core(k);
	v8df p2mn = _mm512_div_pd(one, p2n);
	e = _mm512_mul_pd(e, x);
	e = _mm512_sub_pd(e, hxs);
	e = _mm512_sub_pd(e, x);
	t = _mm512_sub_pd(one, p2mn);
        v8df y = _mm512_sub_pd(t, e);
	y = _mm512_mul_pd(y, py);
	return _mm512_mul_pd(y, p2n);
}

#else

inline v4df v1pow2d_core(v4si m)
{
	m = _mm_add_epi32(m, off);
#ifdef __AVX2__
	v4di n = _mm256_cvtepi32_epi64(m);
	n = _mm256_slli_epi64(n, 52);
#else
	m = _mm_slli_epi32(m, 20);
	v4si z = _mm_set1_epi32(0);
	v4si nl = _mm_unpacklo_epi32(z, m);
	v8si n = _mm256_set1_epi32(0);
	n = _mm256_insertf128_si256(n, nl, 0);
	v4si nu = _mm_unpackhi_epi32(z, m);
	n = _mm256_insertf128_si256(n, nu, 1);
#endif
	return _mm256_castsi256_pd(n);
}

inline v4si v1expd_core(v4df& x)
{
	v4df px = _mm256_mul_pd(x, log2e);
	px = _mm256_floor_pd(px);
	v4si m = _mm256_cvtpd_epi32(px);
#ifdef __FMA__
	x = _mm256_fmadd_pd(px, exp_c1, x);
	x = _mm256_fmadd_pd(px, exp_c2, x);
#else
	v4df y = _mm256_mul_pd(px, exp_c1);
	x = _mm256_add_pd(x, y);
	y = _mm256_mul_pd(px, exp_c2);
	x = _mm256_add_pd(x, y);
#endif
	v4df xx = _mm256_mul_pd(x, x);
#ifdef __FMA__
	px = _mm256_fmadd_pd(exp_p0, xx, exp_p1);
	px = _mm256_fmadd_pd(px, xx, exp_p2);
	v4df qx = _mm256_fmadd_pd(exp_q0, xx, exp_q1);
	qx = _mm256_fmadd_pd(qx, xx, exp_q2);
	qx = _mm256_fmadd_pd(qx, xx, exp_q3);
#else
	px = _mm256_mul_pd(exp_p0, xx);
	px = _mm256_add_pd(px, exp_p1);
	px = _mm256_mul_pd(px, xx);
	px = _mm256_add_pd(px, exp_p2);
	v4df qx = _mm256_mul_pd(exp_q0, xx);
	qx = _mm256_add_pd(qx, exp_q1);
	qx = _mm256_mul_pd(qx, xx);
	qx = _mm256_add_pd(qx, exp_q2);
	qx = _mm256_mul_pd(qx, xx);
	qx = _mm256_add_pd(qx, exp_q3);
#endif
	px = _mm256_mul_pd(x, px);
	xx = _mm256_sub_pd(qx, px);
	x = _mm256_div_pd(px, xx);
#ifdef __FMA__
	x = _mm256_fmadd_pd(two, x, one);
#else
	x = _mm256_mul_pd(x, two);
	x = _mm256_add_pd(x, one);
#endif
	return m;
}

inline v4df v1expm1d_core(v4df x)
{
	v4df px = _mm256_mul_pd(x, log2e);
	px = _mm256_round_pd(px, _MM_FROUND_NINT);
	v4df py = _mm256_min_pd(px, c1023);
	v4si k = _mm256_cvtpd_epi32(py);
	py = _mm256_sub_pd(px, py);
	py = _mm256_add_pd(py, one);
#ifdef __FMA__
	x = _mm256_fmadd_pd(px, exp_c1, x);
	x = _mm256_fmadd_pd(px, exp_c2, x);
#else
	v4df w = _mm256_mul_pd(px, exp_c1);
	x = _mm256_add_pd(x, w);
	w = _mm256_mul_pd(px, exp_c2);
	x = _mm256_add_pd(x, w);
#endif
	v4df hfx = _mm256_mul_pd(x, half);
	v4df hxs = _mm256_mul_pd(x, hfx);
	// sign of t is reversed compared to openlibm implementation
#ifdef __FMA__
	v4df r1 = _mm256_fmadd_pd(expm1_q5, hxs, expm1_q4);
	r1 = _mm256_fmadd_pd(r1, hxs, expm1_q3);
	r1 = _mm256_fmadd_pd(r1, hxs, expm1_q2);
	r1 = _mm256_fmadd_pd(r1, hxs, expm1_q1);
	r1 = _mm256_fmadd_pd(r1, hxs, one);
	v4df t = _mm256_fmadd_pd(r1, hfx, mthree);
	v4df e = _mm256_fmadd_pd(x, t, six);
#else
	v4df r1 = _mm256_mul_pd(expm1_q5, hxs);
	r1 = _mm256_add_pd(r1, expm1_q4);
	r1 = _mm256_mul_pd(r1, hxs);
	r1 = _mm256_add_pd(r1, expm1_q3);
	r1 = _mm256_mul_pd(r1, hxs);
	r1 = _mm256_add_pd(r1, expm1_q2);
	r1 = _mm256_mul_pd(r1, hxs);
	r1 = _mm256_add_pd(r1, expm1_q1);
	r1 = _mm256_mul_pd(r1, hxs);
	r1 = _mm256_add_pd(r1, one);
	v4df t = _mm256_mul_pd(r1, hfx);
	t = _mm256_add_pd(t, mthree);
	v4df e = _mm256_mul_pd(x, t);
	e = _mm256_add_pd(e, six);
#endif
	v4df h = _mm256_add_pd(r1, t);
	h = _mm256_div_pd(h, e);
	e = _mm256_mul_pd(h, hxs);
	v4df p2n = v1pow2d_core(k);
	v4df p2mn = _mm256_div_pd(one, p2n);
	e = _mm256_mul_pd(e, x);
	e = _mm256_sub_pd(e, hxs);
	e = _mm256_sub_pd(e, x);
	t = _mm256_sub_pd(one, p2mn);
        v4df y = _mm256_sub_pd(t, e);
	y = _mm256_mul_pd(y, py);
	return _mm256_mul_pd(y, p2n);
}

#endif // __AVX512F__

VECFI_CONST(expf_max_arg,0x42b17217); // +88.72283173
VECFI_CONST(expf_min_arg,0xc2aeac51); // -87.33655548

VECFI_CONST(exp10f_max_arg,0x421a209a); // +38.53183746
VECFI_CONST(exp10f_min_arg,0xc217b819); // -37.92978287

VECFI_CONST(expm1f_max_arg,0x42b17217); // +88.72283173
VECFI_CONST(expm1f_min_arg,0xc18aa123); // -17.32868004

#if defined(__AVX2__) || defined(__AVX512F__)
VECII_CONST(offf,0x7f);
VECII_CONST(mofff,0xffffff81);
#else
VECI_CONST(offf,0x7f);
VECI_CONST(mofff,0xffffff81);
#endif

VECFI_CONST(expf_c1,0xbf318000); // -0.693359375f
VECFI_CONST(expf_c2,0x395e8083); // 2.12194440e-4f

VECFI_CONST(expf_p0,0x39506967); // 1.9875691500e-4f
VECFI_CONST(expf_p1,0x3ab743ce); // 1.3981999507e-3f
VECFI_CONST(expf_p2,0x3c088908); // 8.3334519073e-3f
VECFI_CONST(expf_p3,0x3d2aa9c1); // 4.1665795894e-2f
VECFI_CONST(expf_p4,0x3e2aaaaa); // 1.6666665459e-1f
VECFI_CONST(expf_p5,0x3f000000); // 5.0000000000e-1f

VECFI_CONST(exp10f_c1,0xbc01861c); // -7.9055088724e-3

VECFI_CONST(expm1f_q1,0xbd088868); // -3.3333212137e-2
VECFI_CONST(expm1f_q2,0x3acf3010); //  1.5807170421e-3

#ifdef __AVX512F__

inline v16sf v1pow2f_core(v16si n)
{
	n = _mm512_add_epi32(n, offf);
	n = _mm512_slli_epi32(n, 23);
	return _mm512_castsi512_ps(n);
}

inline v16si v1expf_core(v16sf& x)
{
	v16sf px = _mm512_mul_ps(x, log2ef);
	px = _mm512_floor_ps(px);
	v16si n = _mm512_cvtps_epi32(px);
	x = _mm512_fmadd_ps(px, expf_c1, x);
	x = _mm512_fmadd_ps(px, expf_c2, x);
	v16sf xx = _mm512_mul_ps(x, x);
	px = _mm512_fmadd_ps(expf_p0, x, expf_p1);
	px = _mm512_fmadd_ps(px, x, expf_p2);
	px = _mm512_fmadd_ps(px, x, expf_p3);
	px = _mm512_fmadd_ps(px, x, expf_p4);
	px = _mm512_fmadd_ps(px, x, expf_p5);
	px = _mm512_mul_ps(px, xx);
	px = _mm512_add_ps(px, x);
	x = _mm512_add_ps(px, onef);
	return n;
}

inline v16sf v1expm1f_core(v16sf x)
{
	v16sf px = _mm512_mul_ps(x, log2ef);
	px = _mm512_roundscale_ps(px, _MM_FROUND_NINT);
	px = _mm512_min_ps(px, c127);
	v16si n = _mm512_cvtps_epi32(px);
	x = _mm512_fmadd_ps(px, expf_c1, x);
	x = _mm512_fmadd_ps(px, expf_c2, x);
	v16sf hfx = _mm512_mul_ps(x, halff);
	v16sf hxs = _mm512_mul_ps(x, hfx);
	// sign of t is reversed compared to openlibm implementation
	v16sf r1 = _mm512_fmadd_ps(expm1f_q2, hxs, expm1f_q1);
	r1 = _mm512_fmadd_ps(r1, hxs, onef);
	v16sf t = _mm512_fmadd_ps(r1, hfx, mthreef);
	v16sf e = _mm512_fmadd_ps(x, t, sixf);
	v16sf h = _mm512_add_ps(r1, t);
	h = _mm512_div_ps(h, e);
	e = _mm512_mul_ps(h, hxs);
	v16sf p2n = v1pow2f_core(n);
	v16sf p2mn = _mm512_div_ps(onef, p2n);
	e = _mm512_mul_ps(e, x);
	e = _mm512_sub_ps(e, hxs);
	e = _mm512_sub_ps(e, x);
	t = _mm512_sub_ps(onef, p2mn);
        v16sf y = _mm512_sub_ps(t, e);
	return _mm512_mul_ps(y, p2n);
}

#else

inline v8sf v1pow2f_core(v8si n)
{
#ifdef __AVX2__
	n = _mm256_add_epi32(n, offf);
	n = _mm256_slli_epi32(n, 23);
#else
	v4si nl = _mm256_extractf128_si256(n, 0);
	nl = _mm_add_epi32(nl, offf);
	v4si nu = _mm256_extractf128_si256(n, 1);
	nu = _mm_add_epi32(nu, offf);
	nl = _mm_slli_epi32(nl, 23);
	nu = _mm_slli_epi32(nu, 23);
	n = _mm256_insertf128_si256(n, nl, 0);
	n = _mm256_insertf128_si256(n, nu, 1);
#endif
	return _mm256_castsi256_ps(n);
}

inline v8si v1expf_core(v8sf& x)
{
	v8sf px = _mm256_mul_ps(x, log2ef);
	px = _mm256_floor_ps(px);
	v8si n = _mm256_cvtps_epi32(px);
#ifdef __FMA__
	x = _mm256_fmadd_ps(px, expf_c1, x);
	x = _mm256_fmadd_ps(px, expf_c2, x);
#else
	v8sf y = _mm256_mul_ps(px, expf_c1);
	x = _mm256_add_ps(x, y);
	y = _mm256_mul_ps(px, expf_c2);
	x = _mm256_add_ps(x, y);
#endif
	v8sf xx = _mm256_mul_ps(x, x);
#ifdef __FMA__
	px = _mm256_fmadd_ps(expf_p0, x, expf_p1);
	px = _mm256_fmadd_ps(px, x, expf_p2);
	px = _mm256_fmadd_ps(px, x, expf_p3);
	px = _mm256_fmadd_ps(px, x, expf_p4);
	px = _mm256_fmadd_ps(px, x, expf_p5);
#else
	px = _mm256_mul_ps(expf_p0, x);
	px = _mm256_add_ps(px, expf_p1);
	px = _mm256_mul_ps(px, x);
	px = _mm256_add_ps(px, expf_p2);
	px = _mm256_mul_ps(px, x);
	px = _mm256_add_ps(px, expf_p3);
	px = _mm256_mul_ps(px, x);
	px = _mm256_add_ps(px, expf_p4);
	px = _mm256_mul_ps(px, x);
	px = _mm256_add_ps(px, expf_p5);
#endif
	px = _mm256_mul_ps(px, xx);
	px = _mm256_add_ps(px, x);
	x = _mm256_add_ps(px, onef);
	return n;
}

inline v8sf v1expm1f_core(v8sf x)
{
	v8sf px = _mm256_mul_ps(x, log2ef);
	px = _mm256_round_ps(px, _MM_FROUND_NINT);
	px = _mm256_min_ps(px, c127);
	v8si n = _mm256_cvtps_epi32(px);
#ifdef __FMA__
	x = _mm256_fmadd_ps(px, expf_c1, x);
	x = _mm256_fmadd_ps(px, expf_c2, x);
#else
	v8sf w = _mm256_mul_ps(px, expf_c1);
	x = _mm256_add_ps(x, w);
	w = _mm256_mul_ps(px, expf_c2);
	x = _mm256_add_ps(x, w);
#endif
	v8sf hfx = _mm256_mul_ps(x, halff);
	v8sf hxs = _mm256_mul_ps(x, hfx);
	// sign of t is reversed compared to openlibm implementation
#ifdef __FMA__
	v8sf r1 = _mm256_fmadd_ps(expm1f_q2, hxs, expm1f_q1);
	r1 = _mm256_fmadd_ps(r1, hxs, onef);
	v8sf t = _mm256_fmadd_ps(r1, hfx, mthreef);
	v8sf e = _mm256_fmadd_ps(x, t, sixf);
#else
	v8sf r1 = _mm256_mul_ps(expm1f_q2, hxs);
	r1 = _mm256_add_ps(r1, expm1f_q1);
	r1 = _mm256_mul_ps(r1, hxs);
	r1 = _mm256_add_ps(r1, onef);
	v8sf t = _mm256_mul_ps(r1, hfx);
	t = _mm256_add_ps(t, mthreef);
	v8sf e = _mm256_mul_ps(x, t);
	e = _mm256_add_ps(e, sixf);
#endif
	v8sf h = _mm256_add_ps(r1, t);
	h = _mm256_div_ps(h, e);
	e = _mm256_mul_ps(h, hxs);
	v8sf p2n = v1pow2f_core(n);
	v8sf p2mn = _mm256_div_ps(onef, p2n);
	e = _mm256_mul_ps(e, x);
	e = _mm256_sub_ps(e, hxs);
	e = _mm256_sub_ps(e, x);
	t = _mm256_sub_ps(onef, p2mn);
        v8sf y = _mm256_sub_ps(t, e);
	return _mm256_mul_ps(y, p2n);
}

#endif // __AVX512F__

#endif // __AVX__

#endif
