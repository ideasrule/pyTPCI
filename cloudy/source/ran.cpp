/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "vectorize.h"
#include "vectorize_math.h"
#include "thirdparty.h"
#include "ran.h"

// aditional PRNG algorithms are here:
// http://xoroshiro.di.unimi.it/

t_ran ran( PRNG_XOSHIRO256STARSTAR );

STATIC void init_seed(uint64 s, int nRANK, uint64 state[], size_t ns, size_t npack, algo_prng algo);
STATIC void convert_double(double* pool, size_t size);
STATIC void convert_float(sys_float* pool, size_t size);
STATIC void convert_zig(double* pool, size_t size);

#if defined(__AVX512F__)

// AVX512f version of xoroshiro128plus() v1.0
void t_ran::p_xoroshiro128plus(uint64* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_xoroshiro128plus()" );

	size_t npack = 8;
	// state and pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(p_state)%64 == 0 );
	ASSERT( reinterpret_cast<unsigned long>(pool)%64 == 0 );
	ASSERT( size%npack == 0 );
	size_t stride = p_ns/SQ_XOROSHIRO128;
	ASSERT( stride >= npack );
	v8di* state0 = (v8di*)(&p_state[0]);
	v8di* state1 = (v8di*)(&p_state[stride]);
	v8di st0 = *state0;
	v8di st1 = *state1;
	v8di* p = (v8di*)pool;
	for( size_t i=0; i < size/npack; ++i )
	{
		const v8di s0 = st0;
		v8di s1 = st1;
		const v8di result = _mm512_add_epi64(s0, s1);

		s1 = _mm512_xor_si512(s0, s1);
		v8di t = _mm512_rol_epi64(s0, 24);
		v8di y = _mm512_slli_epi64(s1, 16);
		y = _mm512_xor_si512(y, s1);
		st0 = _mm512_xor_si512(y, t);
		st1 = _mm512_rol_epi64(s1, 37);

		*p++ = result;
	}
	*state0 = st0;
	*state1 = st1;
}

// AVX512f version of xoshiro256starstar() v1.0
void t_ran::p_xoshiro256starstar(uint64* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_xoshiro256starstar()" );

	size_t npack = 8;
	// state and pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(p_state)%64 == 0 );
	ASSERT( reinterpret_cast<unsigned long>(pool)%64 == 0 );
	ASSERT( size%npack == 0 );
	size_t stride = p_ns/SQ_XOSHIRO256;
	ASSERT( stride >= npack );
	v8di* state0 = (v8di*)(&p_state[0]);
	v8di* state1 = (v8di*)(&p_state[stride]);
	v8di* state2 = (v8di*)(&p_state[2*stride]);
	v8di* state3 = (v8di*)(&p_state[3*stride]);
	v8di st0 = *state0;
	v8di st1 = *state1;
	v8di st2 = *state2;
	v8di st3 = *state3;
	v8di* p = (v8di*)pool;
	for( size_t i=0; i < size/npack; ++i )
	{
		v8di h1 = _mm512_slli_epi64(st1, 2);
		v8di h2 = _mm512_add_epi64(st1, h1); // h2 = st1 * 5
		h2 = _mm512_rol_epi64(h2, 7);
		h1 = _mm512_slli_epi64(h2, 3);
		const v8di result = _mm512_add_epi64(h1, h2); // result = h2 * 9
		const v8di t = _mm512_slli_epi64(st1, 17);

		st2 = _mm512_xor_si512(st2, st0);
		st3 = _mm512_xor_si512(st3, st1);
		st1 = _mm512_xor_si512(st1, st2);
		st0 = _mm512_xor_si512(st0, st3);
		st2 = _mm512_xor_si512(st2, t);
		st3 = _mm512_rol_epi64(st3, 45);
	
		*p++ = result;
	}
	*state0 = st0;
	*state1 = st1;
	*state2 = st2;
	*state3 = st3;
}

#elif defined(__AVX2__)

// AVX2 version of xoroshiro128plus() v1.0
void t_ran::p_xoroshiro128plus(uint64* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_xoroshiro128plus()" );

	size_t npack = 4;
	// state and pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(p_state)%32 == 0 );
	ASSERT( reinterpret_cast<unsigned long>(pool)%32 == 0 );
	ASSERT( size%npack == 0 );
	size_t stride = p_ns/SQ_XOROSHIRO128;
	ASSERT( stride >= npack );
	v4di* state0 = (v4di*)(&p_state[0]);
	v4di* state1 = (v4di*)(&p_state[stride]);
	v4di st0 = *state0;
	v4di st1 = *state1;
	v4di* p = (v4di*)pool;
	for( size_t i=0; i < size/npack; ++i )
	{
		const v4di s0 = st0;
		v4di s1 = st1;
		const v4di result = _mm256_add_epi64(s0, s1);

		s1 = _mm256_xor_si256(s0, s1);
		v4di y = _mm256_slli_epi64(s0, 24);
		v4di z = _mm256_srli_epi64(s0, 40);
		v4di t = _mm256_or_si256(y, z);
		y = _mm256_slli_epi64(s1, 16);
		y = _mm256_xor_si256(y, s1);
		st0 = _mm256_xor_si256(y, t);
		y = _mm256_slli_epi64(s1, 37);
		z = _mm256_srli_epi64(s1, 27);
		st1 = _mm256_or_si256(y, z);

		*p++ = result;
	}
	*state0 = st0;
	*state1 = st1;
}

// AVX2 version of xoshiro256starstar() v1.0
void t_ran::p_xoshiro256starstar(uint64* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_xoshiro256starstar()" );

	size_t npack = 4;
	// state and pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(p_state)%32 == 0 );
	ASSERT( reinterpret_cast<unsigned long>(pool)%32 == 0 );
	ASSERT( size%npack == 0 );
	size_t stride = p_ns/SQ_XOSHIRO256;
	ASSERT( stride >= npack );
	v4di* state0 = (v4di*)(&p_state[0]);
	v4di* state1 = (v4di*)(&p_state[stride]);
	v4di* state2 = (v4di*)(&p_state[2*stride]);
	v4di* state3 = (v4di*)(&p_state[3*stride]);
	v4di st0 = *state0;
	v4di st1 = *state1;
	v4di st2 = *state2;
	v4di st3 = *state3;
	v4di* p = (v4di*)pool;
	for( size_t i=0; i < size/npack; ++i )
	{
		v4di h1 = _mm256_slli_epi64(st1, 2);
		v4di h2 = _mm256_add_epi64(st1, h1); // h2 = st1 * 5
	        h1 = _mm256_slli_epi64(h2, 7);
		h2 = _mm256_srli_epi64(h2, 57);
		h2 = _mm256_or_si256(h1, h2); // h2 = rotl(h2, 7)
		h1 = _mm256_slli_epi64(h2, 3);
		const v4di result = _mm256_add_epi64(h1, h2); // result = h2 * 9
		const v4di t = _mm256_slli_epi64(st1, 17);

		st2 = _mm256_xor_si256(st2, st0);
		st3 = _mm256_xor_si256(st3, st1);
		st1 = _mm256_xor_si256(st1, st2);
		st0 = _mm256_xor_si256(st0, st3);
		st2 = _mm256_xor_si256(st2, t);
		h1 = _mm256_slli_epi64(st3, 45);
		h2 = _mm256_srli_epi64(st3, 19);
		st3 = _mm256_or_si256(h1, h2); // st3 = rotl(st3, 45)
	
		*p++ = result;
	}
	*state0 = st0;
	*state1 = st1;
	*state2 = st2;
	*state3 = st3;
}

// This needs to test for __AVX__ despite the fact that the following routine
// only uses SSE2 features. This is because the SIMD infrastructure (e.g.
// definitions for v2di, etc) is only set up for AVX architecture and higher...
#elif defined(__AVX__)

// SSE2 version of xoroshiro128plus() v1.0
void t_ran::p_xoroshiro128plus(uint64* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_xoroshiro128plus()" );

	size_t npack = 2;
	// state and pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(p_state)%16 == 0 );
	ASSERT( reinterpret_cast<unsigned long>(pool)%16 == 0 );
	ASSERT( size%npack == 0 );
	size_t stride = p_ns/SQ_XOROSHIRO128;
	ASSERT( stride >= npack );
	v2di* state0 = (v2di*)(&p_state[0]);
	v2di* state1 = (v2di*)(&p_state[stride]);
	v2di st0 = *state0;
	v2di st1 = *state1;
	v2di* p = (v2di*)pool;
	for( size_t i=0; i < size/npack; ++i )
	{
		const v2di s0 = st0;
		v2di s1 = st1;
		const v2di result = _mm_add_epi64(s0, s1);

		s1 = _mm_xor_si128(s0, s1);
		v2di y = _mm_slli_epi64(s0, 24);
		v2di z = _mm_srli_epi64(s0, 40);
		v2di t = _mm_or_si128(y, z);
		y = _mm_slli_epi64(s1, 16);
		y = _mm_xor_si128(y, s1);
		st0 = _mm_xor_si128(y, t);
		y = _mm_slli_epi64(s1, 37);
		z = _mm_srli_epi64(s1, 27);
		st1 = _mm_or_si128(y, z);

		*p++ = result;
	}
	*state0 = st0;
	*state1 = st1;
}

// SSE2 version of xoshiro256starstar() v1.0
void t_ran::p_xoshiro256starstar(uint64* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_xoshiro256starstar()" );

	size_t npack = 2;
	// state and pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(p_state)%16 == 0 );
	ASSERT( reinterpret_cast<unsigned long>(pool)%16 == 0 );
	ASSERT( size%npack == 0 );
	size_t stride = p_ns/SQ_XOSHIRO256;
	ASSERT( stride >= npack );
	v2di* state0 = (v2di*)(&p_state[0]);
	v2di* state1 = (v2di*)(&p_state[stride]);
	v2di* state2 = (v2di*)(&p_state[2*stride]);
	v2di* state3 = (v2di*)(&p_state[3*stride]);
	v2di st0 = *state0;
	v2di st1 = *state1;
	v2di st2 = *state2;
	v2di st3 = *state3;
	v2di* p = (v2di*)pool;
	for( size_t i=0; i < size/npack; ++i )
	{
		v2di h1 = _mm_slli_epi64(st1, 2);
		v2di h2 = _mm_add_epi64(st1, h1); // h2 = st1 * 5
	        h1 = _mm_slli_epi64(h2, 7);
		h2 = _mm_srli_epi64(h2, 57);
		h2 = _mm_or_si128(h1, h2); // h2 = rotl(h2, 7)
		h1 = _mm_slli_epi64(h2, 3);
		const v2di result = _mm_add_epi64(h1, h2); // result = h2 * 9
		const v2di t = _mm_slli_epi64(st1, 17);

		st2 = _mm_xor_si128(st2, st0);
		st3 = _mm_xor_si128(st3, st1);
		st1 = _mm_xor_si128(st1, st2);
		st0 = _mm_xor_si128(st0, st3);
		st2 = _mm_xor_si128(st2, t);
		h1 = _mm_slli_epi64(st3, 45);
		h2 = _mm_srli_epi64(st3, 19);
		st3 = _mm_or_si128(h1, h2); // st3 = rotl(st3, 45)
	
		*p++ = result;
	}
	*state0 = st0;
	*state1 = st1;
	*state2 = st2;
	*state3 = st3;
}

#else

void t_ran::p_xoroshiro128plus(uint64* pool, size_t size)
{
	xoroshiro128plus(pool, size, p_state, p_ns);
}

void t_ran::p_xoshiro256starstar(uint64* pool, size_t size)
{
	xoshiro256starstar(pool, size, p_state, p_ns);
}

#endif

// size is size in quadwords (i.e. 64 bits)
void t_ran::p_u64(void* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_u64()" );

	if( p_algo == PRNG_XOROSHIRO128PLUS )
		p_xoroshiro128plus((uint64*)pool, size);
	else if( p_algo == PRNG_XOSHIRO256STARSTAR )
		p_xoshiro256starstar((uint64*)pool, size);
	else
		TotalInsanity();
}

// size is size in quadwords (i.e. 64 bits)
void t_ran::p_dbl(void* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_dbl()" );

	// this routine assumes standard IEEE 64-bit floats
	static_assert( sizeof(double) == 8, "non-standard implementation of type double" );

	p_u64((uint64*)pool, size);
	convert_double((double*)pool, size);
}

// size is size in quadwords (i.e. 64 bits)
void t_ran::p_flt(void* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_flt()" );

	// this routine assumes standard IEEE 32-bit floats
	static_assert( sizeof(sys_float) == 4, "non-standard implementation of type float" );

	p_u64((uint64*)pool, size);
	convert_float((sys_float*)pool, 2*size);
}

// size is size in quadwords (i.e. 64 bits)
void t_ran::p_zig(void* pool, size_t size)
{
	DEBUG_ENTRY( "t_ran::p_zig()" );

	// this routine assumes standard IEEE 64-bit floats
	static_assert( sizeof(double) == 8, "non-standard implementation of type double" );

	p_u64((uint64*)pool, size);
	convert_zig((double*)pool, size);
}

//
// The following code is part of the Ziggurat algorithm for generating normal
// variates centered around 0 with unit variance. This is much faster than the polar
// rejection (a.k.a. Box-Muller) method that can e.g. be found in numerical recipes.
//
// The original algorithm is described by (http://www.jstatsoft.org/v05/i08)
// >>refer	random	generator	Marsaglia G. & Tsang W.W., Journal of Statistical Software 5, 8 (2000)
//
// This algorithm was later slightly improved by
// >>refer	random	generator	Doornik J.A., http://www.doornik.com/research/ziggurat.pdf (2005)
// We implement a somewhat further optimized version in Cloudy
//
// A comparison of the multitude of algorithms for generating normal variates can be found here
//   http://www.cse.cuhk.edu.hk/~phwl/mt/public/archives/papers/grng_acmcs07.pdf
//

static const int ZIGC = 256;

static const double ZigXN[ZIGC+1] = {
	3.9107579595370900e+00, 3.6541528853610088e+00, 3.4492782985609645e+00, 3.3202447338391661e+00,
	3.2245750520470291e+00, 3.1478892895171500e+00, 3.0835261320012330e+00, 3.0278377917686354e+00,
	2.9786032798808448e+00, 2.9343668672078542e+00, 2.8941210536123481e+00, 2.8571387308721325e+00,
	2.8228773968253251e+00, 2.7909211740007858e+00, 2.7609440052788226e+00, 2.7326853590428271e+00,
	2.7059336561218581e+00, 2.6805146432845222e+00, 2.6562830375755024e+00, 2.6331163936303246e+00,
	2.6109105184875485e+00, 2.5895759867069952e+00, 2.5690354526805366e+00, 2.5492215503234608e+00,
	2.5300752321585169e+00, 2.5115444416253423e+00, 2.4935830412696807e+00, 2.4761499396691433e+00,
	2.4592083743333113e+00, 2.4427253181989568e+00, 2.4266709849357260e+00, 2.4110184138996855e+00,
	2.3957431197804806e+00, 2.3808227951706260e+00, 2.3662370567158186e+00, 2.3519672273776600e+00,
	2.3379961487950314e+00, 2.3243080188696230e+00, 2.3108882505998500e+00, 2.2977233489013296e+00,
	2.2848008027229461e+00, 2.2721089902268239e+00, 2.2596370951722178e+00, 2.2473750329458078e+00,
	2.2353133849283280e+00, 2.2234433400909057e+00, 2.2117566428825444e+00, 2.2002455466096480e+00,
	2.1889027716247207e+00, 2.1777214677386416e+00, 2.1666951803526460e+00, 2.1558178198750633e+00,
	2.1450836340462036e+00, 2.1344871828443202e+00, 2.1240233156878157e+00, 2.1136871506849340e+00,
	2.1034740557131468e+00, 2.0933796311370503e+00, 2.0833996939965518e+00, 2.0735302635169788e+00,
	2.0637675478099564e+00, 2.0541079316488648e+00, 2.0445479652157328e+00, 2.0350843537278087e+00,
	2.0257139478620330e+00, 2.0164337349043717e+00, 2.0072408305586849e+00, 1.9981324713565642e+00,
	1.9891060076155713e+00, 1.9801588968985984e+00, 1.9712886979317696e+00, 1.9624930649424619e+00,
	1.9537697423827340e+00, 1.9451165600067539e+00, 1.9365314282737589e+00, 1.9280123340507183e+00,
	1.9195573365912288e+00, 1.9111645637692822e+00, 1.9028322085484464e+00, 1.8945585256687101e+00,
	1.8863418285347764e+00, 1.8781804862909777e+00, 1.8700729210692368e+00, 1.8620176053976323e+00,
	1.8540130597581481e+00, 1.8460578502831198e+00, 1.8381505865807286e+00, 1.8302899196806666e+00,
	1.8224745400917832e+00, 1.8147031759641676e+00, 1.8069745913486934e+00, 1.7992875845475802e+00,
	1.7916409865500100e+00, 1.7840336595472763e+00, 1.7764644955223450e+00, 1.7689324149090779e+00,
	1.7614363653167067e+00, 1.7539753203154551e+00, 1.7465482782794930e+00, 1.7391542612836690e+00,
	1.7317923140507072e+00, 1.7244615029457757e+00, 1.7171609150155407e+00, 1.7098896570690061e+00,
	1.7026468547976139e+00, 1.6954316519322385e+00, 1.6882432094348587e+00, 1.6810807047228233e+00,
	1.6739433309237604e+00, 1.6668302961592867e+00, 1.6597408228557895e+00, 1.6526741470806485e+00,
	1.6456295179023603e+00, 1.6386061967731111e+00, 1.6316034569324220e+00, 1.6246205828305684e+00,
	1.6176568695705342e+00, 1.6107116223673337e+00, 1.6037841560235830e+00, 1.5968737944202613e+00,
	1.5899798700216485e+00, 1.5831017233934714e+00, 1.5762387027333329e+00, 1.5693901634125345e+00,
	1.5625554675284397e+00, 1.5557339834665549e+00, 1.5489250854715355e+00, 1.5421281532263476e+00,
	1.5353425714388431e+00, 1.5285677294350246e+00, 1.5218030207582931e+00, 1.5150478427739924e+00,
	1.5083015962785720e+00, 1.5015636851127065e+00, 1.4948335157777184e+00, 1.4881104970546544e+00,
	1.4813940396253757e+00, 1.4746835556950255e+00, 1.4679784586152309e+00, 1.4612781625074078e+00,
	1.4545820818855233e+00, 1.4478896312776697e+00, 1.4412002248457980e+00, 1.4345132760029464e+00,
	1.4278281970272904e+00, 1.4211443986723231e+00, 1.4144612897724647e+00, 1.4077782768433715e+00,
	1.4010947636762026e+00, 1.3944101509250713e+00, 1.3877238356868846e+00, 1.3810352110727420e+00,
	1.3743436657700305e+00, 1.3676485835943180e+00, 1.3609493430301018e+00, 1.3542453167594306e+00,
	1.3475358711773593e+00, 1.3408203658931521e+00, 1.3340981532160836e+00, 1.3273685776246247e+00,
	1.3206309752177301e+00, 1.3138846731468690e+00, 1.3071289890273539e+00, 1.3003632303274337e+00,
	1.2935866937335176e+00, 1.2867986644897864e+00, 1.2799984157103332e+00, 1.2731852076618437e+00,
	1.2663582870146883e+00, 1.2595168860601442e+00, 1.2526602218912979e+00, 1.2457874955449979e+00,
	1.2388978911020274e+00, 1.2319905747424451e+00, 1.2250646937528080e+00, 1.2181193754817266e+00,
	1.2111537262399112e+00, 1.2041668301405601e+00, 1.1971577478755859e+00, 1.1901255154228016e+00,
	1.1830691426787607e+00, 1.1759876120114898e+00, 1.1688798767268338e+00, 1.1617448594415742e+00,
	1.1545814503558518e+00, 1.1473885054167339e+00, 1.1401648443639958e+00, 1.1329092486483370e+00,
	1.1256204592112944e+00, 1.1182971741150629e+00, 1.1109380460092495e+00, 1.1035416794202682e+00,
	1.0961066278476035e+00, 1.0886313906495142e+00, 1.0811144096988894e+00, 1.0735540657878717e+00,
	1.0659486747575067e+00, 1.0582964833260065e+00, 1.0505956645862071e+00, 1.0428443131393705e+00,
	1.0350404398286053e+00, 1.0271819660307513e+00, 1.0192667174605292e+00, 1.0112924174349784e+00,
	1.0032566795395914e+00, 9.9515699962994308e-01, 9.8699074709384627e-01, 9.7875515528893775e-01,
	9.7044731105886461e-01, 9.6206414321760525e-01, 9.5360240987557265e-01, 9.4505868446257113e-01,
	9.3642934028089686e-01, 9.2771053339623477e-01, 9.1889818364373499e-01, 9.0998795349076900e-01,
	9.0097522445517453e-01, 8.9185507072679238e-01, 8.8262222957891012e-01, 8.7327106808249455e-01,
	8.6379554554682692e-01, 8.5418917100156055e-01, 8.4444495490242366e-01, 8.3455535407951875e-01,
	8.2451220874528863e-01, 8.1430667012806435e-01, 8.0392911698266489e-01, 7.9336905883315278e-01,
	7.8261502329958876e-01, 7.7165442421673935e-01, 7.6047340642208316e-01, 7.4905666200958165e-01,
	7.3738721142583863e-01, 7.2544614090130355e-01, 7.1321228518202273e-01, 7.0066184109758445e-01,
	6.8776789278625772e-01, 6.7449982282743648e-01, 6.6082257423420598e-01, 6.4669571488438893e-01,
	6.3207223637502463e-01, 6.1689698999623555e-01, 6.0110461774394042e-01, 5.8461676609372226e-01,
	5.6733825704047303e-01, 5.4915170231302679e-01, 5.2990972064649511e-01, 5.0942332958593339e-01,
	4.8744396612175434e-01, 4.6363433677176324e-01, 4.3751840218666266e-01, 4.0838913458800075e-01,
	3.7512133285046573e-01, 3.3573751918045946e-01, 2.8617459174726051e-01, 2.1524189591327381e-01,
	0.0000000000000000e+00
};

static const double ZigRN[ZIGC] = {
	9.3438482339457918e-01, 9.4393376707887688e-01, 9.6259114123217282e-01, 9.7118595481318182e-01,
	9.7621833534896407e-01, 9.7955355109525166e-01, 9.8194004595756246e-01, 9.8373938259782689e-01,
	9.8514860539777549e-01, 9.8628466874907117e-01, 9.8722157019173218e-01, 9.8800851576557880e-01,
	9.8867955694410314e-01, 9.8925904142287513e-01, 9.8976486079327719e-01, 9.9021047087164871e-01,
	9.9060619510022785e-01, 9.9096009202198276e-01, 9.9127854840110596e-01, 9.9156669443231094e-01,
	9.9182870051291072e-01, 9.9206799331940876e-01, 9.9228741575504453e-01, 9.9248934712539427e-01,
	9.9267579465715527e-01, 9.9284846405344196e-01, 9.9300881450025391e-01, 9.9315810199357490e-01,
	9.9329741379120706e-01, 9.9342769604768011e-01, 9.9354977616116547e-01, 9.9366438098060894e-01,
	9.9377215174420630e-01, 9.9387365641643144e-01, 9.9396939993917421e-01, 9.9405983279868837e-01,
	9.9414535822376648e-01, 9.9422633826462492e-01, 9.9430309895119195e-01, 9.9437593469006502e-01,
	9.9444511202858621e-01, 9.9451087289022999e-01, 9.9457343736628845e-01, 9.9463300613352923e-01,
	9.9468976255523878e-01, 9.9474387451317581e-01, 9.9479549600995243e-01, 9.9484476857485049e-01,
	9.9489182250074104e-01, 9.9493677793540536e-01, 9.9497974584693893e-01, 9.9502082887992749e-01,
	9.9506012211659289e-01, 9.9509771375503842e-01, 9.9513368571496374e-01, 9.9516811417977458e-01,
	9.9520107008276160e-01, 9.9523261954398701e-01, 9.9526282426362433e-01, 9.9529174187674330e-01,
	9.9531942627388326e-01, 9.9534592789120968e-01, 9.9537129397356761e-01, 9.9539556881334612e-01,
	9.9541879396769928e-01, 9.9544100845638606e-01, 9.9546224894220314e-01, 9.9548254989577112e-01,
	9.9550194374622680e-01, 9.9552046101920322e-01, 9.9553813046331774e-01, 9.9555497916626590e-01,
	9.9557103266149105e-01, 9.9558631502630091e-01, 9.9560084897220869e-01, 9.9561465592819853e-01,
	9.9562775611753773e-01, 9.9564016862870119e-01, 9.9565191148091414e-01, 9.9566300168476796e-01,
	9.9567345529832307e-01, 9.9568328747906887e-01, 9.9569251253207880e-01, 9.9570114395466480e-01,
	9.9570919447780692e-01, 9.9571667610460934e-01, 9.9572360014601191e-01, 9.9572997725396051e-01,
	9.9573581745222939e-01, 9.9574113016506516e-01, 9.9574592424380703e-01, 9.9575020799162972e-01,
	9.9575398918653768e-01, 9.9575727510272860e-01, 9.9576007253043786e-01, 9.9576238779435977e-01,
	9.9576422677074117e-01, 9.9576559490322458e-01, 9.9576649721752453e-01, 9.9576693833499974e-01,
	9.9576692248518850e-01, 9.9576645351736526e-01, 9.9576553491117126e-01, 9.9576416978636662e-01,
	9.9576236091175052e-01, 9.9576011071328807e-01, 9.9575742128147937e-01, 9.9575429437800855e-01,
	9.9575073144169801e-01, 9.9574673359379606e-01, 9.9574230164262523e-01, 9.9573743608760856e-01,
	9.9573213712269715e-01, 9.9572640463921136e-01, 9.9572023822811573e-01, 9.9571363718173422e-01,
	9.9570660049492177e-01, 9.9569912686569617e-01, 9.9569121469533950e-01, 9.9568286208797385e-01,
	9.9567406684961146e-01, 9.9566482648668508e-01, 9.9565513820405349e-01, 9.9564499890248248e-01,
	9.9563440517559698e-01, 9.9562335330629759e-01, 9.9561183926263364e-01, 9.9559985869312606e-01,
	9.9558740692152536e-01, 9.9557447894099416e-01, 9.9556106940769862e-01, 9.9554717263379067e-01,
	9.9553278257976396e-01, 9.9551789284616121e-01, 9.9550249666460933e-01, 9.9548658688815639e-01,
	9.9547015598088462e-01, 9.9545319600676330e-01, 9.9543569861771430e-01, 9.9541765504084812e-01,
	9.9539905606483314e-01, 9.9537989202535504e-01, 9.9536015278961876e-01, 9.9533982773984286e-01,
	9.9531890575569049e-01, 9.9529737519558037e-01, 9.9527522387681022e-01, 9.9525243905442606e-01,
	9.9522900739876285e-01, 9.9520491497157360e-01, 9.9518014720066261e-01, 9.9515468885292668e-01,
	9.9512852400570317e-01, 9.9510163601631507e-01, 9.9507400748969466e-01, 9.9504562024395526e-01,
	9.9501645527377336e-01, 9.9498649271142992e-01, 9.9495571178534647e-01, 9.9492409077593824e-01,
	9.9489160696859402e-01, 9.9485823660357142e-01, 9.9482395482258057e-01, 9.9478873561181080e-01,
	9.9475255174112853e-01, 9.9471537469915738e-01, 9.9467717462391658e-01, 9.9463792022867370e-01,
	9.9459757872262833e-01, 9.9455611572601110e-01, 9.9451349517914489e-01, 9.9446967924496910e-01,
	9.9442462820447775e-01, 9.9437830034447705e-01, 9.9433065183700176e-01, 9.9428163660966262e-01,
	9.9423120620613348e-01, 9.9417930963589485e-01, 9.9412589321226608e-01, 9.9407090037765122e-01,
	9.9401427151481903e-01, 9.9395594374289498e-01, 9.9389585069661779e-01, 9.9383392228723466e-01,
	9.9377008444324044e-01, 9.9370425882895308e-01, 9.9363636253869403e-01, 9.9356630776406951e-01,
	9.9349400143156352e-01, 9.9341934480730776e-01, 9.9334223306551539e-01, 9.9326255481662407e-01,
	9.9318019159069570e-01, 9.9309501727105265e-01, 9.9300689747246695e-01, 9.9291568885747405e-01,
	9.9282123838350755e-01, 9.9272338247255887e-01, 9.9262194609389565e-01, 9.9251674174904159e-01,
	9.9240756834664812e-01, 9.9229420995307360e-01, 9.9217643440235304e-01, 9.9205399174675046e-01,
	9.9192661252615300e-01, 9.9179400583110655e-01, 9.9165585713021331e-01, 9.9151182582775710e-01,
	9.9136154251165642e-01, 9.9120460584495695e-01, 9.9104057904581400e-01, 9.9086898589098571e-01,
	9.9068930616585893e-01, 9.9050097046948593e-01, 9.9030335426539429e-01, 9.9009577104727475e-01,
	9.8987746446202540e-01, 9.8964759919976886e-01, 9.8940525041967620e-01, 9.8914939142954350e-01,
	9.8887887927323703e-01, 9.8859243779956663e-01, 9.8828863768385278e-01, 9.8796587274272862e-01,
	9.8762233171446345e-01, 9.8725596445898278e-01, 9.8686444124682726e-01, 9.8644510343095493e-01,
	9.8599490329659356e-01, 9.8551033021549073e-01, 9.8498731932492467e-01, 9.8442113771148365e-01,
	9.8380624136205808e-01, 9.8313609373663313e-01, 9.8240293339698259e-01, 9.8159747319601653e-01,
	9.8070850631734174e-01, 9.7972238371257558e-01, 9.7862231119118781e-01, 9.7738738919589319e-01,
	9.7599127836112509e-01, 9.7440030911418207e-01, 9.7257074531867649e-01, 9.7044472540755145e-01,
	9.6794407128065607e-01, 9.6496053533935255e-01, 9.6133984665242955e-01, 9.5685442305509605e-01,
	9.5115412025831936e-01, 9.4367126738941931e-01, 9.3342161734665918e-01, 9.1853896462966156e-01,
	8.9501046669157103e-01, 8.5237596455057263e-01, 7.5213489289562085e-01, 0.0000000000000000e+00
};

static const double ZigEX2[ZIGC+1] = {
	4.7746776458665530e-04, 1.2602859304985980e-03, 2.6090727461063629e-03, 4.0379725933718715e-03,
	5.5224032992647540e-03, 7.0508754713921101e-03, 8.6165827694229171e-03, 1.0214971439731100e-02,
	1.1842757857943104e-02, 1.3497450601780807e-02, 1.5177088307982072e-02, 1.6880083152595839e-02,
	1.8605121275783350e-02, 2.0351096230109354e-02, 2.2117062707379922e-02, 2.3902203305873237e-02,
	2.5705804008632656e-02, 2.7527235669693315e-02, 2.9365939758230111e-02, 3.1221417192023690e-02,
	3.3093219458688698e-02, 3.4980941461833073e-02, 3.6884215688691151e-02, 3.8802707404656918e-02,
	4.0736110656078753e-02, 4.2684144916619378e-02, 4.4646552251446536e-02, 4.6623094902089664e-02,
	4.8613553216035145e-02, 5.0617723861121788e-02, 5.2635418276973649e-02, 5.4666461325077916e-02,
	5.6710690106399467e-02, 5.8767952921137984e-02, 6.0838108349751806e-02, 6.2921024437977854e-02,
	6.5016577971470438e-02, 6.7124653828023989e-02, 6.9245144397250269e-02, 7.1377949059141965e-02,
	7.3522973714240991e-02, 7.5680130359194964e-02, 7.7849336702372207e-02, 8.0030515814947509e-02,
	8.2223595813495684e-02, 8.4428509570654661e-02, 8.6645194450867782e-02, 8.8873592068594229e-02,
	9.1113648066700734e-02, 9.3365311913026619e-02, 9.5628536713353335e-02, 9.7903279039215627e-02,
	1.0018949876917202e-01, 1.0248715894230627e-01, 1.0479622562286706e-01, 1.0711666777507288e-01,
	1.0944845714721002e-01, 1.1179156816424558e-01, 1.1414597782825521e-01, 1.1651166562603701e-01,
	1.1888861344334570e-01, 1.2127680548523544e-01, 1.2367622820205140e-01, 1.2608687022065035e-01,
	1.2850872228047364e-01, 1.3094177717412817e-01, 1.3338602969216284e-01, 1.3584147657175735e-01,
	1.3830811644906432e-01, 1.4078594981496831e-01, 1.4327497897404712e-01, 1.4577520800653793e-01,
	1.4828664273312872e-01, 1.5080929068241017e-01, 1.5334316106083767e-01, 1.5588826472506456e-01,
	1.5844461415652022e-01, 1.6101222343811766e-01, 1.6359110823298295e-01, 1.6618128576511007e-01,
	1.6878277480185033e-01, 1.7139559563815562e-01, 1.7401977008249936e-01, 1.7665532144440665e-01,
	1.7930227452353040e-01, 1.8196065560021649e-01, 1.8463049242750454e-01, 1.8731181422451693e-01,
	1.9000465167119307e-01, 1.9270903690432881e-01, 1.9542500351488559e-01, 1.9815258654653811e-01,
	2.0089182249543133e-01, 2.0364274931112150e-01, 2.0640540639867930e-01, 2.0917983462193565e-01,
	2.1196607630785294e-01, 2.1476417525200850e-01, 2.1757417672517837e-01, 2.2039612748101159e-01,
	2.2323007576478959e-01, 2.2607607132326488e-01, 2.2893416541557748e-01, 2.3180441082524852e-01,
	2.3468686187325269e-01, 2.3758157443217368e-01, 2.4048860594144911e-01, 2.4340801542371199e-01,
	2.4633986350223877e-01, 2.4928421241951670e-01, 2.5224112605694377e-01, 2.5521066995567715e-01,
	2.5819291133864802e-01, 2.6118791913376371e-01, 2.6419576399831757e-01, 2.6721651834463184e-01,
	2.7025025636695998e-01, 2.7329705406967580e-01, 2.7635698929678126e-01, 2.7943014176276532e-01,
	2.8251659308484939e-01, 2.8561642681665811e-01, 2.8872972848335393e-01, 2.9185658561828098e-01,
	2.9499708780116257e-01, 2.9815132669790134e-01, 3.0131939610203412e-01, 3.0450139197789627e-01,
	3.0769741250555377e-01, 3.1090755812756371e-01, 3.1413193159763014e-01, 3.1737063803122240e-01,
	3.2062378495823013e-01, 3.2389148237773202e-01, 3.2717384281495859e-01, 3.3047098138053710e-01,
	3.3378301583210851e-01, 3.3711006663841281e-01, 3.4045225704594545e-01, 3.4380971314829134e-01,
	3.4718256395825148e-01, 3.5057094148288120e-01, 3.5397498080156925e-01, 3.5739482014729052e-01,
	3.6083060099117575e-01, 3.6428246813054960e-01, 3.6775056978059623e-01, 3.7123505766982134e-01,
	3.7473608713949141e-01, 3.7825381724723811e-01, 3.8178841087503135e-01, 3.8534003484173396e-01,
	3.8890886002046460e-01, 3.9249506146101076e-01, 3.9609881851754708e-01, 3.9972031498193167e-01,
	4.0335973922286888e-01, 4.0701728433124795e-01, 4.1069314827198322e-01, 4.1438753404270678e-01,
	4.1810064983968459e-01, 4.2183270923135330e-01, 4.2558393133990058e-01, 4.2935454103134152e-01,
	4.3314476911457406e-01, 4.3695485254992927e-01, 4.4078503466776991e-01, 4.4463556539772775e-01,
	4.4850670150921407e-01, 4.5239870686388250e-01, 4.5631185268077357e-01, 4.6024641781492348e-01,
	4.6420268905027884e-01, 4.6818096140782217e-01, 4.7218153846988326e-01, 4.7620473272168379e-01,
	4.8025086591124971e-01, 4.8432026942891160e-01, 4.8841328470771206e-01, 4.9253026364614866e-01,
	4.9667156905479631e-01, 5.0083757512848215e-01, 5.0502866794582879e-01, 5.0924524599813614e-01,
	5.1348772074974303e-01, 5.1775651723220062e-01, 5.2205207467479486e-01, 5.2637484717418670e-01,
	5.3072530440619392e-01, 5.3510393238301956e-01, 5.3951123425954461e-01, 5.4394773119264994e-01,
	5.4841396325792113e-01, 5.5291049042851992e-01, 5.5743789362148632e-01, 5.6199677581727792e-01,
	5.6658776325895177e-01, 5.7121150673807497e-01, 5.7586868297521054e-01, 5.8055999610368347e-01,
	5.8528617926630033e-01, 5.9004799633579197e-01, 5.9484624377099127e-01, 5.9968175262216772e-01,
	6.0455539070054953e-01, 6.0946806492889538e-01, 6.1442072389207680e-01, 6.1941436060903921e-01,
	6.2445001555027424e-01, 6.2952877992812828e-01, 6.3465179929096005e-01, 6.3982027745643899e-01,
	6.4503548082425188e-01, 6.5029874311429459e-01, 6.5561147058322466e-01, 6.6097514778024136e-01,
	6.6639134391238064e-01, 6.7186171990076637e-01, 6.7738803622251309e-01, 6.8297216164879138e-01,
	6.8861608300852706e-01, 6.9432191613003258e-01, 7.0009191814049010e-01, 7.0592850133679741e-01,
	7.1183424888235847e-01, 7.1781193263490140e-01, 7.2386453347288160e-01, 7.2999526456580244e-01,
	7.3620759813126668e-01, 7.4250529634463625e-01, 7.4889244722372672e-01, 7.5537350651175450e-01,
	7.6195334684154647e-01, 7.6863731580333483e-01, 7.7543130498613833e-01, 7.8234183265986190e-01,
	7.8937614357119856e-01, 7.9654233042825462e-01, 8.0384948317638949e-01, 8.1130787431821993e-01,
	8.1892919160941480e-01, 8.2672683395209423e-01, 8.3471629299293038e-01, 8.4291565311844108e-01,
	8.5134625846512368e-01, 8.6003362120300864e-01, 8.6900868804379316e-01, 8.7830965581614684e-01,
	8.8798466076339988e-01, 8.9809592190630405e-01, 9.0872644006056291e-01, 9.1999150504836025e-01,
	9.3206007596899021e-01, 9.4519895345307803e-01, 9.5987909181241593e-01, 9.7710170128273133e-01,
	1.0000000000000000e+00
};

double t_ran::p_ZigTailNormal(bool lgNegative)
{
	DEBUG_ENTRY( "t_ran::p_ZigTailNormal()" );

	static const double R = 3.6541528853610088;
	double x, y;
	do
	{
		x = log(dbl())/R;
		y = log(dbl());
	}
	while( -2.*y < pow2(x) );
	return lgNegative ? x-R : R-x;
}

void t_ran::p_init(uint64 s, int nRANK)
{
	DEBUG_ENTRY( "t_ran::init()" );

	if( p_lgInitialized )
		return;

	p_seed(s, nRANK);
	p_pc.init(8*ND, this, &t_ran::p_u64);
	p_ps.init(4*ND, this, &t_ran::p_u64);
	p_pi.init(2*ND, this, &t_ran::p_u64);
	p_pl.init(ND, this, &t_ran::p_u64);
	p_pd.init(ND, this, &t_ran::p_dbl);
	p_pf.init(2*ND, this, &t_ran::p_flt);
	p_zd.init(ND, this, &t_ran::p_zig);

	p_zigxd = ZigXN;
	p_zigrd = ZigRN;
	p_zige2d = ZigEX2;

	p_lgInitialized = true;
}

void t_ran::p_seed(uint64 s, int nRANK)
{
	DEBUG_ENTRY( "t_ran::p_seed()" );

	// remember the seed for possible reuse in new_rank()
	p_s = s;
	init_seed(s, nRANK, p_state, p_ns, p_npack, p_algo);
	p_pc.reset();
	p_ps.reset();
	p_pi.reset();
	p_pl.reset();
	p_pd.reset();
	p_pf.reset();
	p_zd.reset();
}

uint64 t_ran::p_generate_random_seed()
{
	DEBUG_ENTRY( "t_ran::p_generate_random_seed()" );

#ifdef HAVE_URANDOM
	// generate seed using /dev/urandom -- this is the preferred method since it
	// gives better results if multiple Cloudy runs are started simultaneously in a
	// script (using the local time, they would likely all have the same seed).
	ifstream ran( "/dev/urandom", mode_rb );
	uint64 s;
	ran.read(reinterpret_cast<char*>(&s), sizeof(s));
	return s;
#else
	// generate a random seed using the local time
	time_t rawtime;
	time(&rawtime);
	// now randomize the raw time by taking md5sum of time string
	tm tim;
	localtime_r(&rawtime, &tim);
	const int NCHR = 128;
	char buf[NCHR];
	strftime(buf, NCHR, "%A %F %T %z (%Z)", &tim);
	uint64 s[2];
	MD5string(buf, s);
	return s[0]+s[1];
#endif
}

t_ran::t_ran(algo_prng algo)
{
	DEBUG_ENTRY( "t_ran::t_ran()" );

	p_lgInitialized = false;
	p_algo = algo;
	if( p_algo == PRNG_XOROSHIRO128PLUS )
		p_sq = SQ_XOROSHIRO128;
	else if( p_algo == PRNG_XOSHIRO256STARSTAR )
		p_sq = SQ_XOSHIRO256;
	else
		TotalInsanity();
	p_npack = CD_ALIGN/sizeof(uint64);
	// each quadword in a SIMD register has its own state so that it gets a different PRNG stream
	p_ns = p_sq*p_npack; // p_state[p_ns]
	// make sure the memory is correctly aligned for SIMD instructions
	void* p;
	if( posix_memalign(&p, CD_ALIGN, p_ns*sizeof(uint64)) != 0 )
		throw bad_alloc();
	p_s = 0;
	p_state = (uint64*)p;
}

STATIC void init_seed(uint64 s, int nRANK, uint64 state[], size_t ns, size_t npack, algo_prng algo)
{
	DEBUG_ENTRY( "init_seed()" );

	ASSERT( ns%npack == 0 );
	uint32 sz = ns/npack; // size of the state in quadwords

	// generate the initial state using splitmix64; this works for
	// any size of the state, as long as it is a multiple of 64 bits
	vector<uint64> st(sz);
	for( uint32 i=0; i < sz; ++i )
		st[i] = splitmix64(s);

	for( int r=0; r <= nRANK; ++r )
	{
		for( size_t i=0; i < npack; ++i )
		{
			if( algo == PRNG_XOROSHIRO128PLUS )
			{
				ASSERT( sz == 2 );
				xoroshiro128plus_jump(st[0], st[1]);
			}
			else if( algo == PRNG_XOSHIRO256STARSTAR )
			{
				ASSERT( sz == 4 );
				xoshiro256starstar_jump(st[0], st[1], st[2], st[3]);
			}
			else
			{
				TotalInsanity();
			}
			for( uint32 s=0; s < sz; ++s )
				state[i+s*npack] = st[s];
		}
	}
}

#ifdef __AVX__

VECLL_CONST(dmask1, 0x000fffffffffffff);
VECLL_CONST(dmask2, 0x3ff0000000000000);
VECLL_CONST(dmask3, 0x4000000000000000);
VECDI_CONST(doffset, 0x3fefffffffffffff);
VECD_CONST(doffset2, 3.);
VECDI_CONST(doffset3, 0x3cb0000000000000); // 2^-52

VECII_CONST(fmask1, 0x007fffff);
VECII_CONST(fmask2, 0x3f800000);
VECFI_CONST(foffset, 0x3f7fffff);

#ifdef __AVX512F__

// convert integer variates into double variates in the range (0.,1.)
STATIC void convert_double(double* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_double()" );

	// pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(pool)%64 == 0 );
	v8di* p = (v8di*)pool;
	size_t npack = 8;
	for( size_t i=0; i < size/npack; ++i )
	{
		v8di z = *p;
		z = _mm512_and_si512(z, dmask1);
		z = _mm512_or_si512(z, dmask2);
		v8df x = _mm512_castsi512_pd(z);
		x = _mm512_sub_pd(x, doffset);
		*p++ = _mm512_castpd_si512(x);
	}
}

// convert integer variates into sys_float variates in the range (0.f,1.f)
STATIC void convert_float(sys_float* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_float()" );

	// pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(pool)%64 == 0 );
	v16si* p = (v16si*)pool;
	size_t npack = 16;
	for( size_t i=0; i < size/npack; ++i )
	{
		v16si z = *p;
		z = _mm512_and_si512(z, fmask1);
		z = _mm512_or_si512(z, fmask2);
		v16sf x = _mm512_castsi512_ps(z);
		x = _mm512_sub_ps(x, foffset);
		*p++ = _mm512_castps_si512(x);
	}
}

// convert integer variates into double variates in the range (-1.,1.)
STATIC void convert_zig(double* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_zig()" );

	// pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(pool)%64 == 0 );
	v8di* p = (v8di*)pool;
	size_t npack = 8;
	for( size_t i=0; i < size/npack; ++i )
	{
		v8di z = *p;
		z = _mm512_and_si512(z, dmask1);
		z = _mm512_or_si512(z, dmask3);
		v8df x = _mm512_castsi512_pd(z);
		x = _mm512_sub_pd(x, doffset2);
		x = _mm512_add_pd(x, doffset3);
		*p++ = _mm512_castpd_si512(x);
	}
}

#else

STATIC void convert_double(double* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_double()" );

	// pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(pool)%32 == 0 );
	v4df* p = (v4df*)pool;
	size_t npack = 4;
	for( size_t i=0; i < size/npack; ++i )
	{
		v4df z = *p;
		v4df mask1 = _mm256_castsi256_pd(dmask1);
		z = _mm256_and_pd(z, mask1);
		v4df mask2 = _mm256_castsi256_pd(dmask2);
		z = _mm256_or_pd(z, mask2);
		z = _mm256_sub_pd(z, doffset);
		*p++ = z;
	}
}

STATIC void convert_float(sys_float* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_float()" );

	// pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(pool)%32 == 0 );
	v8sf* p = (v8sf*)pool;
	size_t npack = 8;
	for( size_t i=0; i < size/npack; ++i )
	{
		v8sf z = *p;
		v8sf mask1 = _mm256_castsi256_ps(fmask1);
		z = _mm256_and_ps(z, mask1);
		v8sf mask2 = _mm256_castsi256_ps(fmask2);
		z = _mm256_or_ps(z, mask2);
		z = _mm256_sub_ps(z, foffset);
		*p++ = z;
	}
}

STATIC void convert_zig(double* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_zig()" );

	// pool must be correctly aligned for SIMD access
	ASSERT( reinterpret_cast<unsigned long>(pool)%32 == 0 );
	v4df* p = (v4df*)pool;
	size_t npack = 4;
	for( size_t i=0; i < size/npack; ++i )
	{
		v4df z = *p;
		v4df mask1 = _mm256_castsi256_pd(dmask1);
		z = _mm256_and_pd(z, mask1);
		v4df mask3 = _mm256_castsi256_pd(dmask3);
		z = _mm256_or_pd(z, mask3);
		z = _mm256_sub_pd(z, doffset2);
		z = _mm256_add_pd(z, doffset3);
		*p++ = z;
	}
}

#endif // __AVX512F__

#else

union pund
{
	uint64 i;
	double d;
};

union punf
{
	uint32 i;
	sys_float f;
};

STATIC void convert_double(double* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_double()" );

	pund c = { 0x3fefffffffffffffULL };

	for( size_t i=0; i < size; ++i )
	{
		auto p = reinterpret_cast<pund*>(pool+i);
		// the next two operations create a random FP number in the range [1.,2.)
		p->i &= 0x000fffffffffffffULL;
		p->i |= 0x3ff0000000000000ULL;
		// shift to range (0.,1.)
		p->d -= c.d; // 1 - 2^-53 -- 0x3fefffffffffffff
	}
}

STATIC void convert_float(sys_float* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_float()" );

	punf c = { 0x3f7fffff };

	for( size_t i=0; i < size; ++i )
	{
		auto p = reinterpret_cast<punf*>(pool+i);
		// the next two operations create a random FP number in the range [1.f,2.f)
		p->i &= 0x007fffff;
		p->i |= 0x3f800000;
		// shift to range (0.f,1.f)
		p->f -= c.f; // 1 - 2^-24 -- 0x3f7fffff
	}
}

STATIC void convert_zig(double* pool, size_t size)
{
	DEBUG_ENTRY ( "convert_zig()" );

	pund c = { 0x3cb0000000000000ULL };

	for( size_t i=0; i < size; ++i )
	{
		auto p = reinterpret_cast<pund*>(pool+i);
		// the next two operations create a random FP number in the range [2.,4.)
		p->i &= 0x000fffffffffffffULL;
		p->i |= 0x4000000000000000ULL;
		// shift to range (-1.,1.)
		p->d -= 3.;
		p->d += c.d; // 2^-52 -- 0x3cb0000000000000
	}
}

#endif // __AVX__
