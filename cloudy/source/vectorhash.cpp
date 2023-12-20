/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

//-----------------------------------------------------------------------------
//  This code has been written by:
//
//  Peter A.M. van Hoof
//  Royal Observatory of Belgium
//  p.vanhoof@oma.be
//
//  Vectorhash is a fast, non-cryptographic hash function optimized for large
//  buffers. It has in part been inspired by xoroshiro128+ written by David
//  Blackman and Sebastiano Vigna (for the core of the hashing function) and
//  MurmurHash3 written by Austin Appleby (for the finalization mix and some
//  support routines). This program has been written such that it can be
//  vectorized using SIMD instructions from the SSE2, AVX2, or AVX512f
//  instruction sets. It also has a scalar version for non-Intel platforms.
//  The resulting hash is identical for each of these instruction sets.
//

// this file deliberately does NOT include cddefines.h
// this makes vh128sum.exe a true standalone program and avoids bloating the executable
// with loads of Cloudy stuff; it also avoids lots of spurious recompilations of vh128sum.exe

#include <iostream>
#include <iomanip>
#include <stdint.h>
#include <cstdio>
#include <string>
#include <sstream>
#include <cstdlib>

#ifdef __SSE2__
#include <immintrin.h>
#endif

#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#else
#define _POSIX_MAPPED_FILES 0
#endif

#if _POSIX_MAPPED_FILES > 0
#include <sys/mman.h>
#endif

using namespace std;

//-----------------------------------------------------------------------------
// tunable parameters

// width of the hash (in bits)
// this value MUST be a power of 2 with a minimum of 128
static const size_t vh_hash_width = 128;

// width of largest hardware register that is supported (in bits)
static const size_t vh_hwreg_width = 512;

//-----------------------------------------------------------------------------
// derived quantities

// width of the virtual SIMD register supported in VectorHash (in bits)
// must be at least as wide as the largest hardware register that is used
static const size_t vh_virtreg_width = ( 2*vh_hash_width > vh_hwreg_width ) ? 2*vh_hash_width : vh_hwreg_width;
// width of the hash in uint32_t elements
static const size_t vh_nstate = vh_hash_width/32;
// number of uint32_t's that fit in the virtual register
static const size_t vh_nint = vh_virtreg_width/32;

// the file is read with this blocksize (in bytes)
static const size_t blocksize = 4*vh_nint*sizeof(uint32_t);

//-----------------------------------------------------------------------------
// Platform-specific functions and macros

#if defined(_MSC_VER)

#define ROTL32(x,y)	_rotl(x,y)

#else

inline uint32_t rotl32 ( uint32_t x, int r )
{
	return (x << r) | (x >> (32 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)

#endif

#ifdef _MSC_VER
// posix_memalign not defined on windows
inline int posix_memalign(void **p, size_t a, size_t s)
{
	*p = _aligned_malloc(s, a);
	return ( *p == NULL ) ? errno : 0;
}

inline void posix_memalign_free(void *p)
{
	_aligned_free(p);
}
#else
inline void posix_memalign_free(void *p)
{
	free(p);
}
#endif

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

inline uint32_t getblock32 ( uint32_t p )
{
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
	return (((p & 0xff000000) >> 24) | ((p & 0x00ff0000) >>  8) |
		((p & 0x0000ff00) <<  8) | ((p & 0x000000ff) << 24));
#else
	return p;
#endif
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

inline uint32_t fmix32 ( uint32_t h0, uint32_t h1 = 0xd86b048b )
{
	h1 ^= h0;
	h0 = ROTL32(h0, 11) ^ h1 ^ ROTL32(h1, 13);
	h1 = ROTL32(h1, 19);
	return h0 + h1;
}

//-----------------------------------------------------------------------------

inline void pad_buffer(const char* src, char* buf, size_t len)
{
	for( size_t i=0; i < len; i++ )
		buf[i] = src[i];
	char c = '\0';
	while( len < blocksize )
		buf[len++] = ++c;
}

//-----------------------------------------------------------------------------

inline void VectorHashFold(uint32_t* res, uint32_t* z[], size_t nbits)
{
	uint32_t narr = nbits/32;
	uint32_t nelem = 4*vh_virtreg_width/nbits;
	for( uint32_t i=0; i < narr; i++ )
	{
		res[i] = z[i][0];
		for( size_t j=1; j < nelem; j++ )
			res[i] ^= z[(j+i)%narr][j];
	}
}

#define VEC(X) for( size_t j=0; j < vh_nint; j++ ) { X; }

static void VectorHashFinalize(size_t len, uint32_t* h1, uint32_t* h2, uint32_t* h3, uint32_t* h4, void* out)
{
	uint32_t lenc = fmix32(len & 0x7fffffff);

	VEC( h1[j] ^= lenc );
	VEC( h2[j] ^= lenc );
	VEC( h3[j] ^= lenc );
	VEC( h4[j] ^= lenc );

	for( size_t j=1; j < vh_nint; j++ )
	{
		h1[0] += h2[j];
		h2[0] += h3[j];
		h3[0] += h4[j];
		h4[0] += h1[j];
		h1[j] += h2[0];
		h2[j] += h3[0];
		h3[j] += h4[0];
		h4[j] += h1[0];
	}

	VEC( h1[j] += h2[j] );
	VEC( h1[j] += h3[j] );
	VEC( h1[j] += h4[j] );
	VEC( h2[j] += h1[j] );
	VEC( h3[j] += h1[j] );
	VEC( h4[j] += h1[j] );
	VEC( h1[j] = fmix32(h1[j], h3[j]) );
	VEC( h2[j] = fmix32(h2[j], h4[j]) );
	VEC( h3[j] = fmix32(h3[j], h1[j]) );
	VEC( h4[j] = fmix32(h4[j], h2[j]) );
	VEC( h1[j] += h2[j] );
	VEC( h1[j] += h3[j] );
	VEC( h1[j] += h4[j] );
	VEC( h2[j] += h1[j] );
	VEC( h3[j] += h1[j] );
	VEC( h4[j] += h1[j] );

	for( size_t j=1; j < vh_nint; j++ )
	{
		h1[0] += h4[j];
		h2[0] += h1[j];
		h3[0] += h2[j];
		h4[0] += h3[j];
		h1[j] += h4[0];
		h2[j] += h1[0];
		h3[j] += h2[0];
		h4[j] += h3[0];
	}

	uint32_t* res = (uint32_t*)out;
	uint32_t* z[vh_nstate];
	size_t j = 0;
	size_t nn = vh_nint*4/vh_nstate;
	for( size_t i=0; i < vh_nstate/4; i++ )
		z[j++] = &h1[i*nn];
	for( size_t i=0; i < vh_nstate/4; i++ )
		z[j++] = &h2[i*nn];
	for( size_t i=0; i < vh_nstate/4; i++ )
		z[j++] = &h3[i*nn];
	for( size_t i=0; i < vh_nstate/4; i++ )
		z[j++] = &h4[i*nn];
	VectorHashFold(res, z, vh_hash_width);
}

#undef VEC

static void stateinit(uint32_t st[], uint32_t& seed)
{
	for( size_t i=0; i < vh_nint; i++ )
		st[i] = seed;
	seed = fmix32(seed);
}

#if defined(__AVX512F__)

typedef __m512i v16si;

static const size_t npack = 16;
static const size_t nreg = vh_nint/npack;

#define VEC(X) for( size_t j=0; j < nreg; j++ ) { X; }

inline void VectorHashBody(const v16si* data, v16si h1[], v16si h2[], v16si h3[], v16si h4[])
{
	v16si s[nreg], x1[nreg], x2[nreg];

	VEC( s[j] = _mm512_xor_si512(h1[j], h2[j]) );
	VEC( s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( x1[j] = _mm512_rol_epi32(h1[j], 11) );
	VEC( x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( h1[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( h2[j] = _mm512_rol_epi32(s[j], 19) );

	VEC( s[j] = _mm512_xor_si512(h2[j], h3[j]) );
	VEC( s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( x1[j] = _mm512_rol_epi32(h2[j], 11) );
	VEC( x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( h2[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( h3[j] = _mm512_rol_epi32(s[j], 19) );

	VEC( s[j] = _mm512_xor_si512(h3[j], h4[j]) );
	VEC( s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( x1[j] = _mm512_rol_epi32(h3[j], 11) );
	VEC( x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( h3[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( h4[j] = _mm512_rol_epi32(s[j], 19) );

	VEC( s[j] = _mm512_xor_si512(h4[j], h1[j]) );
	VEC( s[j] = _mm512_xor_si512(s[j], *data++) );
	VEC( x1[j] = _mm512_rol_epi32(h4[j], 11) );
	VEC( x1[j] = _mm512_xor_si512(x1[j], s[j]) );
	VEC( x2[j] = _mm512_slli_epi32(s[j], 14) );
	VEC( h4[j] = _mm512_xor_si512(x1[j], x2[j]) );
	VEC( h1[j] = _mm512_rol_epi32(s[j], 19) );
}

void VectorHash(const void* key, size_t len, uint32_t seed, void* out)
{
	v16si h1[nreg], h2[nreg], h3[nreg], h4[nreg];
	stateinit( (uint32_t*)h1, seed );
	stateinit( (uint32_t*)h2, seed );
	stateinit( (uint32_t*)h3, seed );
	stateinit( (uint32_t*)h4, seed );

	size_t nblocks = len/blocksize;
	const v16si* data = (const v16si*)key;
	for( size_t i=0; i < nblocks; i++ )
	{
		VectorHashBody(data, h1, h2, h3, h4);
		data += 4*nreg;
	}

	// pad the remaining characters and process...
	v16si buf[blocksize/sizeof(v16si)];
	pad_buffer( (const char*)data, (char*)buf, len-nblocks*blocksize );
	VectorHashBody(buf, h1, h2, h3, h4);

	uint32_t* z1 = (uint32_t*)h1;
	uint32_t* z2 = (uint32_t*)h2;
	uint32_t* z3 = (uint32_t*)h3;
	uint32_t* z4 = (uint32_t*)h4;
	VectorHashFinalize(len, z1, z2, z3, z4, out);
}

#elif defined(__AVX2__)

typedef __m256i v8si;

static const size_t npack = 8;
static const size_t nreg = vh_nint/npack;

#define VEC(X) for( size_t j=0; j < nreg; j++ ) { X; }

inline void VectorHashBody(const v8si* data, v8si h1[], v8si h2[], v8si h3[], v8si h4[])
{
	v8si s[nreg], x1[nreg], x2[nreg];

	VEC( s[j] = _mm256_xor_si256(h1[j], h2[j]) );
	VEC( s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( x1[j] = _mm256_slli_epi32(h1[j], 11) );
	VEC( x2[j] = _mm256_srli_epi32(h1[j], 21) );
	VEC( x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( h1[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( h2[j] = _mm256_or_si256(x1[j], x2[j]) );

	VEC( s[j] = _mm256_xor_si256(h2[j], h3[j]) );
	VEC( s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( x1[j] = _mm256_slli_epi32(h2[j], 11) );
	VEC( x2[j] = _mm256_srli_epi32(h2[j], 21) );
	VEC( x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( h2[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( h3[j] = _mm256_or_si256(x1[j], x2[j]) );

	VEC( s[j] = _mm256_xor_si256(h3[j], h4[j]) );
	VEC( s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( x1[j] = _mm256_slli_epi32(h3[j], 11) );
	VEC( x2[j] = _mm256_srli_epi32(h3[j], 21) );
	VEC( x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( h3[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( h4[j] = _mm256_or_si256(x1[j], x2[j]) );

	VEC( s[j] = _mm256_xor_si256(h4[j], h1[j]) );
	VEC( s[j] = _mm256_xor_si256(s[j], *data++) );
	VEC( x1[j] = _mm256_slli_epi32(h4[j], 11) );
	VEC( x2[j] = _mm256_srli_epi32(h4[j], 21) );
	VEC( x1[j] = _mm256_or_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_xor_si256(x1[j], s[j]) );
	VEC( x2[j] = _mm256_slli_epi32(s[j], 14) );
	VEC( h4[j] = _mm256_xor_si256(x1[j], x2[j]) );
	VEC( x1[j] = _mm256_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm256_srli_epi32(s[j], 13) );
	VEC( h1[j] = _mm256_or_si256(x1[j], x2[j]) );
}

void VectorHash(const void* key, size_t len, uint32_t seed, void* out)
{
	v8si h1[nreg], h2[nreg], h3[nreg], h4[nreg];
	stateinit( (uint32_t*)h1, seed );
	stateinit( (uint32_t*)h2, seed );
	stateinit( (uint32_t*)h3, seed );
	stateinit( (uint32_t*)h4, seed );

	size_t nblocks = len/blocksize;
	const v8si* data = (const v8si*)key;		
	for( size_t i=0; i < nblocks; i++ )
	{
		VectorHashBody(data, h1, h2, h3, h4);
		data += 4*nreg;
	}

	// pad the remaining characters and process...
	v8si buf[blocksize/sizeof(v8si)];
	pad_buffer( (const char*)data, (char*)buf, len-nblocks*blocksize );
	VectorHashBody(buf, h1, h2, h3, h4);

	uint32_t* z1 = (uint32_t*)h1;
	uint32_t* z2 = (uint32_t*)h2;
	uint32_t* z3 = (uint32_t*)h3;
	uint32_t* z4 = (uint32_t*)h4;
	VectorHashFinalize(len, z1, z2, z3, z4, out);
}

#elif defined(__SSE2__)

typedef __m128i v4si;

static const size_t npack = 4;
static const size_t nreg = vh_nint/npack;

#define VEC(X) for( size_t j=0; j < nreg; j++ ) { X; }

inline void VectorHashBody(const v4si* data, v4si h1[], v4si h2[], v4si h3[], v4si h4[])
{
	v4si s[nreg], x1[nreg], x2[nreg];

	VEC( s[j] = _mm_xor_si128(h1[j], h2[j]) );
	VEC( s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( x1[j] = _mm_slli_epi32(h1[j], 11) );
	VEC( x2[j] = _mm_srli_epi32(h1[j], 21) );
	VEC( x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( h1[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( h2[j] = _mm_or_si128(x1[j], x2[j]) );

	VEC( s[j] = _mm_xor_si128(h2[j], h3[j]) );
	VEC( s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( x1[j] = _mm_slli_epi32(h2[j], 11) );
	VEC( x2[j] = _mm_srli_epi32(h2[j], 21) );
	VEC( x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( h2[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( h3[j] = _mm_or_si128(x1[j], x2[j]) );

	VEC( s[j] = _mm_xor_si128(h3[j], h4[j]) );
	VEC( s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( x1[j] = _mm_slli_epi32(h3[j], 11) );
	VEC( x2[j] = _mm_srli_epi32(h3[j], 21) );
	VEC( x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( h3[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( h4[j] = _mm_or_si128(x1[j], x2[j]) );

	VEC( s[j] = _mm_xor_si128(h4[j], h1[j]) );
	VEC( s[j] = _mm_xor_si128(s[j], *data++) );
	VEC( x1[j] = _mm_slli_epi32(h4[j], 11) );
	VEC( x2[j] = _mm_srli_epi32(h4[j], 21) );
	VEC( x1[j] = _mm_or_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_xor_si128(x1[j], s[j]) );
	VEC( x2[j] = _mm_slli_epi32(s[j], 14) );
	VEC( h4[j] = _mm_xor_si128(x1[j], x2[j]) );
	VEC( x1[j] = _mm_slli_epi32(s[j], 19) );
	VEC( x2[j] = _mm_srli_epi32(s[j], 13) );
	VEC( h1[j] = _mm_or_si128(x1[j], x2[j]) );
}

void VectorHash(const void* key, size_t len, uint32_t seed, void* out)
{
	v4si h1[nreg], h2[nreg], h3[nreg], h4[nreg];
	stateinit( (uint32_t*)h1, seed );
	stateinit( (uint32_t*)h2, seed );
	stateinit( (uint32_t*)h3, seed );
	stateinit( (uint32_t*)h4, seed );

	size_t nblocks = len/blocksize;
	const v4si* data = (const v4si*)key;
	for( size_t i=0; i < nblocks; i++ )
	{
		VectorHashBody(data, h1, h2, h3, h4);
		data += 4*nreg;
	}

	// pad the remaining characters and process...
	v4si buf[blocksize/sizeof(v4si)];
	pad_buffer( (const char*)data, (char*)buf, len-nblocks*blocksize );
	VectorHashBody(buf, h1, h2, h3, h4);

	uint32_t* z1 = (uint32_t*)h1;
	uint32_t* z2 = (uint32_t*)h2;
	uint32_t* z3 = (uint32_t*)h3;
	uint32_t* z4 = (uint32_t*)h4;
	VectorHashFinalize(len, z1, z2, z3, z4, out);
}

#else

#define VEC(X) for( size_t j=0; j < vh_nint; j++ ) { X; }

inline void VectorHashBody(const uint32_t* data, uint32_t h1[], uint32_t h2[], uint32_t h3[], uint32_t h4[])
{
	uint32_t s[vh_nint];

	VEC( s[j] = h1[j] ^ h2[j] ^ getblock32( *data++ ) );
	VEC( h1[j] = ROTL32(h1[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( h2[j] = ROTL32(s[j], 19) );

	VEC( s[j] = h2[j] ^ h3[j] ^ getblock32( *data++ ) );
	VEC( h2[j] = ROTL32(h2[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( h3[j] = ROTL32(s[j], 19) );

	VEC( s[j] = h3[j] ^ h4[j] ^ getblock32( *data++ ) );
	VEC( h3[j] = ROTL32(h3[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( h4[j] = ROTL32(s[j], 19) );

	VEC( s[j] = h4[j] ^ h1[j] ^ getblock32( *data++ ) );
	VEC( h4[j] = ROTL32(h4[j], 11) ^ s[j] ^ (s[j] << 14) );
	VEC( h1[j] = ROTL32(s[j], 19) );
}

void VectorHash(const void* key, size_t len, uint32_t seed, void* out)
{
	uint32_t h1[vh_nint], h2[vh_nint], h3[vh_nint], h4[vh_nint];
	stateinit( h1, seed );
	stateinit( h2, seed );
	stateinit( h3, seed );
	stateinit( h4, seed );

	size_t nblocks = len/blocksize;
	const uint32_t* data = (const uint32_t*)key;
	for( size_t i=0; i < nblocks; i++ )
	{
		VectorHashBody(data, h1, h2, h3, h4);
		data += 4*vh_nint;
	}

	// pad the remaining characters and process...
	uint32_t buf[blocksize/sizeof(uint32_t)];
	pad_buffer( (const char*)data, (char*)buf, len-nblocks*blocksize );
	VectorHashBody(buf, h1, h2, h3, h4);

	VectorHashFinalize(len, h1, h2, h3, h4, out);
}

#endif

string VHstream(FILE* io)
{
	if( fseek( io, 0, SEEK_END ) != 0 )
		return string();
	long fsize = ftell(io);
	if( fsize < 0 )
		return string();
#if _POSIX_MAPPED_FILES > 0
	int fd = fileno(io);
	char* map = (char*)mmap( NULL, fsize, PROT_READ, MAP_SHARED, fd, 0 );
	if( fsize > 0 && map == MAP_FAILED )
		return string();
	uint32_t state[vh_nstate];
	VectorHash( map, fsize, 0xfd4c799d, state );
	munmap(map, fsize);
#else
	if( fseek( io, 0, SEEK_SET ) != 0 )
		return string();
	void* map = NULL;
	if( fsize > 0 )
	{
		if( posix_memalign( &map, vh_hwreg_width/8, fsize ) != 0 )
			return string();
		if( fread( map, fsize, 1, io ) != 1 )
			return string();
	}
	uint32_t state[vh_nstate];
	VectorHash( map, fsize, 0xfd4c799d, state );
	if( map != NULL )
		posix_memalign_free( map );
#endif

	ostringstream hash;
	for( uint32_t i=0; i < vh_nstate; ++i )
		hash << hex << setfill('0') << setw(8) << state[i];

	return hash.str();
}

string VHstring(const string& s)
{
	uint32_t state[vh_nstate];
	VectorHash( s.data(), s.length(), 0xfd4c799d, state );

	ostringstream hash;
	for( uint32_t i=0; i < vh_nstate; ++i )
		hash << hex << setfill('0') << setw(8) << state[i];

	return hash.str();
}
