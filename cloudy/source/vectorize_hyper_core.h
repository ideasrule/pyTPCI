/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_HYPER_CORE_H
#define VECTORIZE_HYPER_CORE_H

#include "vectorize_math.h"
#include "vectorize_sqrt_core.h"
#include "vectorize_log_core.h"

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

#ifdef __AVX__

VECLL_CONST(asinh_mask1,0x7fffffffffffffff);
VECLL_CONST(asinh_mask2,0x8000000000000000);

VECDI_CONST(asinh_2p28,0x41b0000000000000); // 2^28

#ifdef __AVX512F__

inline v8df v1asinhd_core(v8df x)
{
	v8df x2 = _mm512_mul_pd(x, x);
	v8df arg = _mm512_add_pd(x2, one);
	arg = v1sqrtd_core(arg);
	arg = _mm512_add_pd(arg, one);
	arg = _mm512_div_pd(x2, arg);
	arg = _mm512_add_pd(arg, x);
	return v1log1pd_core(arg);
}

#else

inline v4df v1asinhd_core(v4df x)
{
	v4df x2 = _mm256_mul_pd(x, x);
	v4df arg = _mm256_add_pd(x2, one);
	arg = v1sqrtd_core(arg);
	arg = _mm256_add_pd(arg, one);
	arg = _mm256_div_pd(x2, arg);
	arg = _mm256_add_pd(arg, x);
	return v1log1pd_core(arg);
}

#endif // __AVX512F__

VECII_CONST(asinh_mask1f,0x7fffffff);
VECII_CONST(asinh_mask2f,0x80000000);

VECFI_CONST(asinhf_2p28,0x4d800000); // 2^28

#ifdef __AVX512F__

inline v16sf v1asinhf_core(v16sf x)
{
	v16sf x2 = _mm512_mul_ps(x, x);
	v16sf arg = _mm512_add_ps(x2, onef);
	arg = v1sqrtf_core(arg);
	arg = _mm512_add_ps(arg, onef);
	arg = _mm512_div_ps(x2, arg);
	arg = _mm512_add_ps(arg, x);
	return v1log1pf_core(arg);
}

#else

inline v8sf v1asinhf_core(v8sf x)
{
	v8sf x2 = _mm256_mul_ps(x, x);
	v8sf arg = _mm256_add_ps(x2, onef);
	arg = v1sqrtf_core(arg);
	arg = _mm256_add_ps(arg, onef);
	arg = _mm256_div_ps(x2, arg);
	arg = _mm256_add_ps(arg, x);
	return v1log1pf_core(arg);
}

#endif // __AVX512F__

#endif // __AVX__

#endif
