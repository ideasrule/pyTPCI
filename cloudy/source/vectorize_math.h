/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_MATH_H
#define VECTORIZE_MATH_H

#ifdef __AVX__

// use shorthand for packed types
typedef __m128d v2df;
typedef __m128  v4sf;
typedef __m128i v2di;
typedef __m128i v4si;

typedef __m256d v4df;
typedef __m256  v8sf;
typedef __m256i v4di;
typedef __m256i v8si;

#ifdef __AVX512F__

typedef __m512d v8df;
typedef __m512  v16sf;
typedef __m512i v8di;
typedef __m512i v16si;

// some macros for defining vector constants

#define VECD_CONST(name,x) \
ALIGNED(64) static const double __avx_##name[8] = {x,x,x,x,x,x,x,x}; \
static const v8df& name = *reinterpret_cast<const v8df*>(__avx_##name)

#define VECDI_CONST(name,x) \
ALIGNED(64) static const uint64 __avx_##name[8] = {x,x,x,x,x,x,x,x}; \
static const v8df& name = *reinterpret_cast<const v8df*>(__avx_##name)

#define VECF_CONST(name,x) \
ALIGNED(64) static const sys_float __avx_##name[16] = {x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x}; \
static const v16sf& name = *reinterpret_cast<const v16sf*>(__avx_##name)

#define VECFI_CONST(name,x) \
ALIGNED(64) static const uint32 __avx_##name[16] = {x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x}; \
static const v16sf& name = *reinterpret_cast<const v16sf*>(__avx_##name)

#define VECI_CONST(name,x) \
ALIGNED(32) static const uint32 __avx_##name[8] = {x,x,x,x,x,x,x,x}; \
static const v8si& name = *reinterpret_cast<const v8si*>(__avx_##name)

#define VECII_CONST(name,x) \
ALIGNED(64) static const uint32 __avx_##name[16] = {x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x}; \
static const v16si& name = *reinterpret_cast<const v16si*>(__avx_##name)

#define VECL_CONST(name,x) \
ALIGNED(32) static const uint64 __avx_##name[4] = {x,x,x,x}; \
static const v4di& name = *reinterpret_cast<const v4di*>(__avx_##name)

#define VECLL_CONST(name,x) \
ALIGNED(64) static const uint64 __avx_##name[8] = {x,x,x,x,x,x,x,x}; \
static const v8di& name = *reinterpret_cast<const v8di*>(__avx_##name)

union vpack
{
	ALIGNED(64) v8df pd;
	ALIGNED(64) v8di pl;
	ALIGNED(64) v16sf pf;
	ALIGNED(64) v16si pi;
	ALIGNED(64) double d[8];
	ALIGNED(64) int64 l[8];
	ALIGNED(64) sys_float f[16];
	ALIGNED(64) int32 i[16];
};

inline bool getmaskbit(unsigned int mask, unsigned int n)
{
	return ( ((mask>>n)&1) == 1 );
}

template<class T>
inline string print_pack(const T x[], int npack, unsigned int mask, const string& name)
{
	ostringstream oss;
	for( int i=0; i < npack; ++i )
	{
		if( getmaskbit(mask, i) )
			oss << " => ";
		else
			oss << "    ";
		oss << name << "[" << i << "] = " << x[i] << "\n";
	}
	return oss.str();
}

inline string DEMsg(const string& routine, v8df x, __mmask8 mask)
{
	vpack p;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pd = x;
	oss << print_pack(p.d, 8, mask, "x");
	return oss.str();
}

inline string DEMsg(const string& routine, v16sf x, __mmask16 mask)
{
	vpack p;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pf = x;
	oss << print_pack(p.f, 16, mask, "x");
	return oss.str();
}

inline string DEMsg(const string& routine, v8df x, __mmask8 mask1, v8df y, __mmask8 mask2)
{
	vpack p;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pd = x;
	oss << print_pack(p.d, 8, mask1, "x");
	p.pd = y;
	oss << print_pack(p.d, 8, mask2, "y");
	return oss.str();
}

inline string DEMsg(const string& routine, v16sf x, __mmask16 mask1, v16sf y, __mmask16 mask2)
{
	vpack p;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pf = x;
	oss << print_pack(p.f, 16, mask1, "x");
	p.pf = y;
	oss << print_pack(p.f, 16, mask2, "y");
	return oss.str();
}

#else // __AVX512F__

#define VECD_CONST(name,x) \
ALIGNED(32) static const double __avx_##name[4] = {x,x,x,x}; \
static const v4df& name = *reinterpret_cast<const v4df*>(__avx_##name)

#define VECDI_CONST(name,x) \
ALIGNED(32) static const uint64 __avx_##name[4] = {x,x,x,x}; \
static const v4df& name = *reinterpret_cast<const v4df*>(__avx_##name)

#define VECF_CONST(name,x) \
ALIGNED(32) static const sys_float __avx_##name[8] = {x,x,x,x,x,x,x,x}; \
static const v8sf& name = *reinterpret_cast<const v8sf*>(__avx_##name)

#define VECFI_CONST(name,x) \
ALIGNED(32) static const uint32 __avx_##name[8] = {x,x,x,x,x,x,x,x}; \
static const v8sf& name = *reinterpret_cast<const v8sf*>(__avx_##name)

#define VECI_CONST(name,x) \
ALIGNED(16) static const uint32 __avx_##name[4] = {x,x,x,x}; \
static const v4si& name = *reinterpret_cast<const v4si*>(__avx_##name)

#define VECII_CONST(name,x) \
ALIGNED(32) static const uint32 __avx_##name[8] = {x,x,x,x,x,x,x,x}; \
static const v8si& name = *reinterpret_cast<const v8si*>(__avx_##name)

#define VECL_CONST(name,x) \
ALIGNED(16) static const uint64 __avx_##name[2] = {x,x}; \
static const v2di& name = *reinterpret_cast<const v2di*>(__avx_##name)

#define VECLL_CONST(name,x) \
ALIGNED(32) static const uint64 __avx_##name[4] = {x,x,x,x}; \
static const v4di& name = *reinterpret_cast<const v4di*>(__avx_##name)

union vpack
{
	ALIGNED(32) v4df pd;
	ALIGNED(32) v4di pl;
	ALIGNED(32) v8sf pf;
	ALIGNED(32) v8si pi;
	ALIGNED(32) double d[4];
	ALIGNED(32) int64 l[4];
	ALIGNED(32) sys_float f[8];
	ALIGNED(32) int32 i[8];
};

template<class T, class U>
inline string print_pack(const T x[], const U m[], int npack, const string& name)
{
	ostringstream oss;
	for( int i=0; i < npack; ++i )
	{
		if( m[i] < 0 )
			oss << " => ";
		else
			oss << "    ";
		oss << name << "[" << i << "] = " << x[i] << "\n";
	}
	return oss.str();
}

inline string DEMsg(const string& routine, v4df x, v4df mask)
{
	vpack p, m;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pd = x;
	m.pd = mask;
	oss << print_pack(p.d, m.l, 4, "x");
	return oss.str();
}

inline string DEMsg(const string& routine, v8sf x, v8sf mask)
{
	vpack p, m;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pf = x;
	m.pf = mask;
	oss << print_pack(p.f, m.i, 8, "x");
	return oss.str();
}

inline string DEMsg(const string& routine, v4df x, v4df mask1, v4df y, v4df mask2)
{
	vpack p, m;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pd = x;
	m.pd = mask1;
	oss << print_pack(p.d, m.l, 4, "x");
	p.pd = y;
	m.pd = mask2;
	oss << print_pack(p.d, m.l, 4, "y");
	return oss.str();
}

inline string DEMsg(const string& routine, v8sf x, v8sf mask1, v8sf y, v8sf mask2)
{
	vpack p, m;
	ostringstream oss;
	oss << "domain error in " << routine << ".\n";
	p.pf = x;
	m.pf = mask1;
	oss << print_pack(p.f, m.i, 8, "x");
	p.pf = y;
	m.pf = mask2;
	oss << print_pack(p.f, m.i, 8, "y");
	return oss.str();
}

#endif // __AVX512F__

VECD_CONST(mthree,-3.);
VECD_CONST(mone,-1.);
VECD_CONST(mhalf,-0.5);
VECD_CONST(zero,0.);
VECDI_CONST(third,0x3fd5555555555555); // 0.33333333333333333333
VECD_CONST(half,0.5);
VECD_CONST(one,1.);
VECD_CONST(c1p5,1.5);
VECD_CONST(two,2.);
VECD_CONST(three,3.);
VECD_CONST(six,6.);
VECD_CONST(ten,10.);
VECD_CONST(c1023,1023.);

VECDI_CONST(ln_two,0x3fe62e42fefa39ef); // log(2.)
VECDI_CONST(log2e,0x3ff71547652b82fe);  // log2(e) = 1./log(2.));
VECDI_CONST(ln_ten,0x40026bb1bbb55516); // log(10.)
VECDI_CONST(log10e,0x3fdbcb7b1526e50e); // log10(e) = 1./log(10.)

VECDI_CONST(dbl_min,0x0010000000000000); // DBL_MIN
VECDI_CONST(dbl_max,0x7fefffffffffffff); // DBL_MAX
VECDI_CONST(sqrt_dbl_max,0x5fefffffffffffff); // sqrt(DBL_MAX)

VECF_CONST(mthreef,-3.f);
VECF_CONST(monef,-1.f);
VECF_CONST(mhalff,-0.5f);
VECF_CONST(zerof,0.f);
VECFI_CONST(thirdf,0x3eaaaaab); // 0.3333333333f
VECF_CONST(halff,0.5f);
VECF_CONST(onef,1.f);
VECF_CONST(twof,2.f);
VECF_CONST(c1p5f,1.5f);
VECF_CONST(threef,3.f);
VECF_CONST(sixf,6.f);
VECF_CONST(tenf,10.f);
VECF_CONST(c127,127.f);

VECFI_CONST(ln_twof,0x3f317218); // logf(2.f)
VECFI_CONST(log2ef,0x3fb8aa3b);  // log2f(e) = 1.f/logf(2.f)
VECFI_CONST(ln_tenf,0x40135d8e); // logf(10.f)
VECFI_CONST(log10ef,0x3ede5bd9); // log10f(e) = 1.f/logf(10.f)

VECFI_CONST(flt_min,0x00800000); // FLT_MIN
VECFI_CONST(flt_max,0x7f7fffff); // FLT_MAX
VECFI_CONST(sqrt_flt_max,0x5f7fffff); // sqrt(FLT_MAX)

template<class V>
inline void vecfun_partial(const sys_float x[], sys_float y[], long nlo, long nhi, long npack, V (*vecfun1)(V))
{
	vpack xx, yy;
	for( long i=0; i < npack; ++i )
		xx.f[i] = x[min(nlo+i,nhi-1)];
	yy.pf = vecfun1(xx.pf);
	for( long i=nlo; i < nhi; ++i )
		y[i] = yy.f[i-nlo];
}

template<class V>
inline void vecfun_partial(const double x[], double y[], long nlo, long nhi, long npack, V (*vecfun1)(V))
{
	vpack xx, yy;
	for( long i=0; i < npack; ++i )
		xx.d[i] = x[min(nlo+i,nhi-1)];
	yy.pd = vecfun1(xx.pd);
	for( long i=nlo; i < nhi; ++i )
		y[i] = yy.d[i-nlo];
}

template<class V>
inline void vecfun2_partial(const sys_float x1[], const sys_float x2[], sys_float y[], long nlo, long nhi,
			    long npack, V (*vecfun1)(V, V))
{
	vpack xx1, xx2, yy;
	for( long i=0; i < npack; ++i )
	{
		xx1.f[i] = x1[min(nlo+i,nhi-1)];
		xx2.f[i] = x2[min(nlo+i,nhi-1)];
	}
	yy.pf = vecfun1(xx1.pf, xx2.pf);
	for( long i=nlo; i < nhi; ++i )
		y[i] = yy.f[i-nlo];
}

template<class V>
inline void vecfun2_partial(const double x1[], const double x2[], double y[], long nlo, long nhi, long npack,
			    V (*vecfun1)(V, V))
{
	vpack xx1, xx2, yy;
	for( long i=0; i < npack; ++i )
	{
		xx1.d[i] = x1[min(nlo+i,nhi-1)];
		xx2.d[i] = x2[min(nlo+i,nhi-1)];
	}
	yy.pd = vecfun1(xx1.pd, xx2.pd);
	for( long i=nlo; i < nhi; ++i )
		y[i] = yy.d[i-nlo];
}

// this is the generic wrapper routine for all vectorized math functions
// it makes sure that the vectorized function gets arguments that are properly aligned
// scalfun1 is the normal scalar routine working on a single argument
// vecfun1 is the vectorized routine working on a packed register
template<class T, class V>
inline void vecfun(const T x[], T y[], long nlo, long nhi, T (*)(T), V (*vecfun1)(V))
{
	if( nhi <= nlo )
		return;

	// determine number of elements of type T in a packed register
	long npack = sizeof(V)/sizeof(T);
	if( nhi-nlo < npack )
	{
		vecfun_partial(x, y, nlo, nhi, npack, vecfun1);
		return;
	}
	long i, ilo = nlo;
	long ox = (reinterpret_cast<long>(x)/sizeof(T)+nlo)%npack;
	long oy = (reinterpret_cast<long>(y)/sizeof(T)+nlo)%npack;
	bool lgNonAligned = ( ox != oy );
	T *yl, *ylocal = NULL;
	if( lgNonAligned )
	{
		ylocal = new T[nhi+npack-1];
		long ol = (reinterpret_cast<long>(ylocal)/sizeof(T)+nlo)%npack;
		if( ox >= ol )
			yl = &ylocal[ox-ol];
		else
			yl = &ylocal[ox+npack-ol];
	}
	else
	{
		yl = y;
	}
	// the initial element is not aligned correctly ->
	// use scalar routine until we are correctly aligned
	if( ox > 0 )
	{
		vecfun_partial(x, yl, ilo, min(ilo+npack-ox,nhi), npack, vecfun1);
		ilo = min(ilo+npack-ox,nhi);
	}
	// use vectorized routine as long as there are at least npack evaluations left
	for( i=ilo; i < nhi-npack+1; i += npack )
	{
		const V& xx = *reinterpret_cast<const V*>(&x[i]);
		V& yy = *reinterpret_cast<V*>(&yl[i]);
		yy = vecfun1(xx);
	}
	ilo = i;
	// use partial routine again for the remaining calls (if any)...
	if( ilo < nhi )
		vecfun_partial(x, yl, ilo, nhi, npack, vecfun1);
	if( lgNonAligned )
	{
		for( i=nlo; i < nhi; ++i )
			y[i] = yl[i];
		delete[] ylocal;
	}
}

template<class T, class V>
inline void vecfun2(const T x1[], const T x2[], T y[], long nlo, long nhi, T (*)(T, T), V (*vecfun1)(V, V))
{
	if( nhi <= nlo )
		return;

	// determine number of elements of type T in a packed register
	long npack = sizeof(V)/sizeof(T);
	if( nhi-nlo < npack )
	{
		vecfun2_partial(x1, x2, y, nlo, nhi, npack, vecfun1);
		return;
	}
	long i, ilo = nlo;
	long ox1 = (reinterpret_cast<long>(x1)/sizeof(T)+nlo)%npack;
	long ox2 = (reinterpret_cast<long>(x2)/sizeof(T)+nlo)%npack;
	long oy = (reinterpret_cast<long>(y)/sizeof(T)+nlo)%npack;
	bool lgNonAligned1 = ( ox1 != ox2 );
	bool lgNonAligned2 = ( ox1 != oy );
	const T *x2l;
	T *x2local = NULL;
	if( lgNonAligned1 )
	{
		x2local = new T[nhi+npack-1];
		long ol = (reinterpret_cast<long>(x2local)/sizeof(T)+nlo)%npack;
		T *ptr;
		if( ox1 >= ol )
			ptr = &x2local[ox1-ol];
		else
			ptr = &x2local[ox1+npack-ol];
		memcpy(ptr+nlo, x2+nlo, size_t((nhi-nlo)*sizeof(T)));
		x2l = ptr;
	}
	else
	{
		x2l = x2;
	}
	T *yl, *ylocal = NULL;
	if( lgNonAligned2 )
	{
		ylocal = new T[nhi+npack-1];
		long ol = (reinterpret_cast<long>(ylocal)/sizeof(T)+nlo)%npack;
		if( ox1 >= ol )
			yl = &ylocal[ox1-ol];
		else
			yl = &ylocal[ox1+npack-ol];
	}
	else
	{
		yl = y;
	}
	// the initial element is not aligned correctly ->
	// use scalar routine until we are correctly aligned
	if( ox1 > 0 )
	{
		vecfun2_partial(x1, x2l, yl, ilo, min(ilo+npack-ox1,nhi), npack, vecfun1);
		ilo = min(ilo+npack-ox1,nhi);
	}
	// use vectorized routine as long as there are at least npack evaluations left
	for( i=ilo; i < nhi-npack+1; i += npack )
	{
		const V& xx1 = *reinterpret_cast<const V*>(&x1[i]);
		const V& xx2 = *reinterpret_cast<const V*>(&x2l[i]);
		V& yy = *reinterpret_cast<V*>(&yl[i]);
		yy = vecfun1(xx1, xx2);
	}
	ilo = i;
	// use partial routine again for the remaining calls (if any)...
	if( ilo < nhi )
		vecfun2_partial(x1, x2l, yl, ilo, nhi, npack, vecfun1);
	if( lgNonAligned1 )
		delete[] x2local;
	if( lgNonAligned2 )
	{
		for( i=nlo; i < nhi; ++i )
			y[i] = yl[i];
		delete[] ylocal;
	}
}

#ifdef __AVX512F__
#define V1FUN_PD_4(FUN, V) \
	v8df xx = _mm512_set_pd( V, V, V, V, x3, x2, x1, x0 ); \
	v8df yy = v1##FUN##d(xx); \
	memcpy(y, &yy, 4*sizeof(double))
#else
#define V1FUN_PD_4(FUN, V) \
	v4df xx = _mm256_set_pd( x3, x2, x1, x0 ); \
	v4df yy = v1##FUN##d(xx); \
	memcpy(y, &yy, 4*sizeof(double))
#endif

#ifdef __AVX512F__
#define V1FUN_PD_8(FUN, V) \
	v8df xx = _mm512_set_pd( x7, x6, x5, x4, x3, x2, x1, x0 ); \
	v8df yy = v1##FUN##d(xx); \
	memcpy(y, &yy, 8*sizeof(double))
#else
#define V1FUN_PD_8(FUN, V) \
	v4df yy[2]; \
	v4df xx = _mm256_set_pd( x3, x2, x1, x0 ); \
	yy[0] = v1##FUN##d(xx); \
	xx = _mm256_set_pd( x7, x6, x5, x4 ); \
	yy[1] = v1##FUN##d(xx); \
	memcpy(y, &yy, 8*sizeof(double))
#endif

#ifdef __AVX512F__
#define V1FUN_PS_4(FUN, V) \
	v16sf xx = _mm512_set_ps( V, V, V, V, V, V, V, V, V, V, V, V, x3, x2, x1, x0 ); \
	v16sf yy = v1##FUN##f(xx); \
	memcpy(y, &yy, 4*sizeof(sys_float))
#else
#define V1FUN_PS_4(FUN, V) \
	v8sf xx = _mm256_set_ps( V, V, V, V, x3, x2, x1, x0 ); \
	v8sf yy = v1##FUN##f(xx); \
	memcpy(y, &yy, 4*sizeof(sys_float))
#endif

#ifdef __AVX512F__
#define V1FUN_PS_8(FUN, V) \
	v16sf xx = _mm512_set_ps( V, V, V, V, V, V, V, V, x7, x6, x5, x4, x3, x2, x1, x0 ); \
	v16sf yy = v1##FUN##f(xx); \
	memcpy(y, &yy, 8*sizeof(sys_float))
#else
#define V1FUN_PS_8(FUN, V) \
	v8sf xx = _mm256_set_ps( x7, x6, x5, x4, x3, x2, x1, x0 ); \
	v8sf yy = v1##FUN##f(xx); \
	memcpy(y, &yy, 8*sizeof(sys_float))
#endif

#ifdef __AVX512F__
#define V1FUN_PS_16(FUN, V) \
	v16sf xx = _mm512_set_ps( x15, x14, x13, x12, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0 ); \
	v16sf yy = v1##FUN##f(xx); \
	memcpy(y, &yy, 16*sizeof(sys_float))
#else
#define V1FUN_PS_16(FUN, V) \
	v8sf yy[2]; \
	v8sf xx = _mm256_set_ps( x7, x6, x5, x4, x3, x2, x1, x0 ); \
	yy[0] = v1##FUN##f(xx); \
	xx = _mm256_set_ps( x15, x14, x13, x12, x11, x10, x9, x8 ); \
	yy[1] = v1##FUN##f(xx); \
	memcpy(y, &yy, 16*sizeof(sys_float))
#endif

#ifdef __AVX512F__
#define V1FUN2_PD_4(FUN, V) \
	v8df xx = _mm512_set_pd( V, V, V, V, x3, x2, x1, x0 ); \
	v8df yy = _mm512_set_pd( V, V, V, V, y3, y2, y1, y0 ); \
	v8df zz = v1##FUN##d(xx, yy); \
	memcpy(z, &zz, 4*sizeof(double))
#else
#define V1FUN2_PD_4(FUN, V) \
	v4df xx = _mm256_set_pd( x3, x2, x1, x0 ); \
	v4df yy = _mm256_set_pd( y3, y2, y1, y0 ); \
	v4df zz = v1##FUN##d(xx, yy); \
	memcpy(z, &zz, 4*sizeof(double))
#endif

#ifdef __AVX512F__
#define V1FUN2_PS_4(FUN, V) \
	v16sf xx = _mm512_set_ps( V, V, V, V, V, V, V, V, V, V, V, V, x3, x2, x1, x0 ); \
	v16sf yy = _mm512_set_ps( V, V, V, V, V, V, V, V, V, V, V, V, y3, y2, y1, y0 ); \
	v16sf zz = v1##FUN##f(xx, yy); \
	memcpy(z, &zz, 4*sizeof(sys_float))
#else
#define V1FUN2_PS_4(FUN, V) \
	v8sf xx = _mm256_set_ps( V, V, V, V, x3, x2, x1, x0 ); \
	v8sf yy = _mm256_set_ps( V, V, V, V, y3, y2, y1, y0 ); \
	v8sf zz = v1##FUN##f(xx, yy); \
	memcpy(z, &zz, 4*sizeof(sys_float))
#endif

#ifdef __AVX512F__
#define V1FUN2_PS_8(FUN, V) \
	v16sf xx = _mm512_set_ps( V, V, V, V, V, V, V, V, x7, x6, x5, x4, x3, x2, x1, x0 ); \
	v16sf yy = _mm512_set_ps( V, V, V, V, V, V, V, V, y7, y6, y5, y4, y3, y2, y1, y0 ); \
	v16sf zz = v1##FUN##f(xx, yy); \
	memcpy(z, &zz, 8*sizeof(sys_float))
#else
#define V1FUN2_PS_8(FUN, V) \
	v8sf xx = _mm256_set_ps( x7, x6, x5, x4, x3, x2, x1, x0 ); \
	v8sf yy = _mm256_set_ps( y7, y6, y5, y4, y3, y2, y1, y0 ); \
	v8sf zz = v1##FUN##f(xx, yy); \
	memcpy(z, &zz, 8*sizeof(sys_float))
#endif

#else // __AVX__

// fallback for non-AVX capable hardware
template<class T, class V>
inline void vecfun(const T x[], T y[], long nlo, long nhi, T (*scalfun1)(T), V (*)(V))
{
	for( long i=nlo; i < nhi; ++i )
		y[i] = scalfun1(x[i]);
}

template<class T, class V>
inline void vecfun2(const T x1[], const T x2[], T y[], long nlo, long nhi, T (*scalfun1)(T, T), V (*)(V, V))
{
	for( long i=nlo; i < nhi; ++i )
		y[i] = scalfun1(x1[i], x2[i]);
}

#define V1FUN_4(FUN, V) \
	y[0] = FUN(x0); \
	y[1] = FUN(x1); \
	y[2] = FUN(x2); \
	y[3] = FUN(x3)

#define V1FUN_8(FUN, V) \
	y[0] = FUN(x0); \
	y[1] = FUN(x1); \
	y[2] = FUN(x2); \
	y[3] = FUN(x3); \
	y[4] = FUN(x4); \
	y[5] = FUN(x5); \
	y[6] = FUN(x6); \
	y[7] = FUN(x7)

#define V1FUN_16(FUN, V) \
	y[0] = FUN(x0); \
	y[1] = FUN(x1); \
	y[2] = FUN(x2); \
	y[3] = FUN(x3); \
	y[4] = FUN(x4); \
	y[5] = FUN(x5); \
	y[6] = FUN(x6); \
	y[7] = FUN(x7); \
	y[8] = FUN(x8); \
	y[9] = FUN(x9); \
	y[10] = FUN(x10); \
	y[11] = FUN(x11); \
	y[12] = FUN(x12); \
	y[13] = FUN(x13); \
	y[14] = FUN(x14); \
	y[15] = FUN(x15)

#define V1FUN_PD_4(FUN, V) \
	V1FUN_4(FUN, V)

#define V1FUN_PD_8(FUN, V) \
	V1FUN_8(FUN, V)

#define V1FUN_PS_4(FUN, V) \
	V1FUN_4(FUN##f, V)

#define V1FUN_PS_8(FUN, V) \
	V1FUN_8(FUN##f, V)

#define V1FUN_PS_16(FUN, V) \
	V1FUN_16(FUN##f, V)

#define V1FUN2_4(FUN, V) \
	z[0] = FUN(x0, y0); \
	z[1] = FUN(x1, y1); \
	z[2] = FUN(x2, y2); \
	z[3] = FUN(x3, y3)

#define V1FUN2_8(FUN, V) \
	z[0] = FUN(x0, y0); \
	z[1] = FUN(x1, y1); \
	z[2] = FUN(x2, y2); \
	z[3] = FUN(x3, y3); \
	z[4] = FUN(x4, y4); \
	z[5] = FUN(x5, y5); \
	z[6] = FUN(x6, y6); \
	z[7] = FUN(x7, y7)

#define V1FUN2_PD_4(FUN, V) \
	V1FUN2_4(FUN, V)

#define V1FUN2_PS_4(FUN, V) \
	V1FUN2_4(FUN##f, V)

#define V1FUN2_PS_8(FUN, V) \
	V1FUN2_8(FUN##f, V)

#endif // __AVX__

#endif
