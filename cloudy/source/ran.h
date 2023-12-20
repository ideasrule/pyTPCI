/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef RAN_H
#define RAN_H

class t_ran;

// see: https://isocpp.org/wiki/faq/pointers-to-members
typedef void (t_ran::*t_ran_fun)(void* p, size_t s);

//
// ran_pool holds a pool of PRNG variates, thus allowing efficient initialization using vectorization
//
// init()           sets the size of the pool and the method for filling the pool
// reset()          invalidates the remaining variates in the pool, thus forcing a reinitialization
//                  the next time next() is called; useful when the random seed is changed
// lgInitialized()  true if init() was called
// next()           returns the next variate; this routine is CPU time critical
// p_alloc()        allocates the pool, making sure it is correctly aligned for SIMD access
// p_update_pool()  (re)fill the pool with new random variates
// ~ran_pool()      frees the pool
//
template<class T>
class ran_pool
{
	void* p_pool;     // start of the pool, correctly aligned for SIMD
	T* p_next;        // pointer to next data item that can be returned
	T* p_end;         // points just beyond the end of the pool
	size_t p_size;    // number of elements in the pool
	size_t p_squad;   // size of the pool in quadwords
	t_ran* p_rc;      // pointer to the instantiation of the random number class
	t_ran_fun p_fill; // fills pool, member function of t_ran class

	void p_alloc()
	{
		// make sure the memory is correctly aligned for SIMD instructions
		if( posix_memalign(&p_pool, CD_ALIGN, p_size*sizeof(T)) != 0 )
			throw bad_alloc();
		p_next = p_end = (T*)p_pool + p_size;
	}
	void p_update_pool()
	{
		ASSERT( lgInitialized() );
		if( p_pool == NULL )
			p_alloc();
		(p_rc->*p_fill)(p_pool, p_squad);
		p_next = (T*)p_pool;
	}
public:
	void init(size_t s, t_ran* rc, t_ran_fun f)
	{
		ASSERT( !lgInitialized() );
		ASSERT( s > 0 );
		p_size = s;
		p_squad = (s*sizeof(T))/8;
		// make sure size is multiple of 8 bytes
		ASSERT( p_squad*8 == s*sizeof(T) );
		p_rc = rc;
		p_fill = f;
	}
	void reset()
	{
		p_next = p_end;
	}
	bool lgInitialized() const
	{
		return ( p_size > 0 );
	}
	T next()
	{
		if( UNLIKELY(p_next == p_end) )
			p_update_pool();
		return *p_next++;
	}

	ran_pool()
	{
		p_pool = p_next = p_end = NULL;
		p_size = p_squad = 0;
		p_rc = NULL;
		p_fill = NULL;
	}
	ran_pool(const ran_pool&) = delete;
	ran_pool& operator= (const ran_pool&) = delete;
	~ran_pool()
	{
		posix_memalign_free(p_pool);
	}
};

//
// Two algorithms are supported here: xoroshiro128+ (v1.0) and xoshiro256** (v1.0)
//
// The 2016 version of xoroshiro128+ was the first to be added, and has since been upgraded to v1.0.
// The xoshiro256** algorithm was added later, and is described as an "all-purpose, rock-solid
// generator" with no known deficiencies. This is the default generator in Cloudy, with xoroshiro128+
// being kept as a backup. The latter has the advantage that it is somewhat faster, but has some mild
// deficiencies. To switch, simply change the instantiation of ran at the top of ran.cpp.
//
// For more information see: http://xoshiro.di.unimi.it/
//

enum algo_prng { PRNG_XOROSHIRO128PLUS, PRNG_XOSHIRO256STARSTAR };

//
// t_ran is a global class for efficiently generating pseudo-random numbers. It can generate integer
//       as well as floating point variates (the latter with a uniform or normal distribution)
//
// The variates are stored in pools to allow efficient computation using SIMD instructions. The
// underlying PRNG algorithm to generate the random bits can be freely chosen. Currently supported
// are xoshiro256** and xoroshiro128+. Each quadword in a SIMD register should generate a different,
// non-overlapping stream of random numbers and therefore needs a different state. This is done
// automatically in the seeding process. The algorithm also takes care that each rank in an MPI run
// gets a different, non-overlapping stream of random numbers. The algorithm is not thread-safe in
// openMP runs.
//
//   ND -- size of the pool of variates in quadwords
//   p_algo -- which algorithm should be used as the underlying PRNG
//   p_sq -- size of the state for the underlying PRNG algorithm in quadwords
//   p_npack -- number of quadwords packed into a single SIMD register
//   p_ns -- size of the state for a full SIMD register in quadwords (= p_sq*p_npack)
//   p_s -- the seed that was used to generate the initial state
//   p_state[p_ns] -- the state of the PRNG (needs to be correctly aligned for SIMD access).
//
// i7(), u8(), i15(), u16(), i31(), u32(), i63(), u64() -- return signed or unsigned integer variates
// dbl(), rnm() -- return uniformly distributed floating point variates on (0,1)
// normal() -- return normally distributed floating point variates (double precision only)
//             these numbers are generated using the Ziggurat algorithm
// init() -- does a one time initialization using a randomly generated seed. Subsequent calls will
//           be ignored. There is also a version with a fixed seed.
// new_rank() -- jump ahead in random number sequence for a different rank number. Mainly useful for
//               forked threads.
// get_seed() -- return the seed that was used.
// print_seed() -- return the seed that was used as a string for printing.
//
class t_ran
{
	// pool size in quadwords
	static const size_t ND = 2048;
	// size of the state in quadwords for a single xoroshiro128plus PRNG stream
	static const size_t SQ_XOROSHIRO128 = 2;
	static const size_t SQ_XOSHIRO256 = 4;

	bool p_lgInitialized;
	algo_prng p_algo;
	size_t p_sq;
	size_t p_npack;
	size_t p_ns;
	uint64 p_s;
	uint64* p_state;

	const double* p_zigxd;
	const double* p_zigrd;
	const double* p_zige2d;

	void p_init(uint64 s, int nRANK);
	void p_seed(uint64 s, int nRANK);
	uint64 p_generate_random_seed();
	uint64 p_random_seed();
	void p_xoroshiro128plus(uint64* pool, size_t size);
	void p_xoshiro256starstar(uint64* pool, size_t size);
	double p_ZigTailNormal(bool lgNegative);

	// fill arrays with random variates
	// NB NB -- the pool need to be correctly aligned for SIMD access
	// NB NB -- the size needs to be a multiple of the SIMD vector size
	void p_u64(void* pool, size_t size);
	void p_dbl(void* pool, size_t size);
	void p_flt(void* pool, size_t size);
	void p_zig(void* pool, size_t size);

	ran_pool<uint8> p_pc;
	ran_pool<uint16> p_ps;
	ran_pool<uint32> p_pi;
	ran_pool<uint64> p_pl;
	ran_pool<double> p_pd;
	ran_pool<sys_float> p_pf;
	ran_pool<double> p_zd;
public:
	/* integer random variate on the [0,0x7f] interval */
	int8 i7() { return int8(p_pc.next() >> 1); }
	/* integer random variate on the [0,0xff] interval */
	uint8 u8() { return p_pc.next(); }
	/* integer random variate on the [0,0x7fff] interval */
	int16 i15() { return int16(p_ps.next() >> 1); }
	/* integer random variate on the [0,0xffff] interval */
	uint16 u16() { return p_ps.next(); }
	/* integer random variate on the [0,0x7fffffff] interval */
	int32 i31() { return int32(p_pi.next() >> 1); }
	/* integer random variate on the [0,0xffffffff] interval */
	uint32 u32() { return p_pi.next(); }
	/* integer random variate on the [0,0x7fffffffffffffff] interval */
	int64 i63() { return int64(p_pl.next() >> 1);	}
	/* integer random variate on the [1,0xffffffffffffffff] interval */
	uint64 u64() { return p_pl.next(); }
	/* double precision uniform random variate on the (0,1) interval */
	double dbl() { return p_pd.next(); }
	/* single precision uniform random variate on the (0,1) interval */
	realnum rnm()
	{
#ifdef FLT_IS_DBL
		return p_pd.next();
#else
		return p_pf.next();
#endif
	}
	/* generates a random number with normal distribution and unit standard deviation using Ziggurat */
	double normal()
	{
		while( true )
		{
			double u = p_zd.next();
			uint8 i = u8();
			if( LIKELY(fabs(u) < p_zigrd[i]) )
				return u*p_zigxd[i];
			if( UNLIKELY(i == 0) )
				return p_ZigTailNormal( u < 0. );
			double x = u*p_zigxd[i];
			double e2 = exp( -0.5*pow2(x) );
			if( p_zige2d[i+1] + dbl()*(p_zige2d[i]-p_zige2d[i+1]) < e2 )
				return x;
		}
	}

	void init()
	{
		uint64 s = p_generate_random_seed();
		p_init(s, 0);
	}
	void init(uint64 s, int nRANK)
	{
		p_init(s, nRANK);
	}
	void new_rank(int nRANK)
	{
		ASSERT( p_lgInitialized );
		p_seed(p_s, nRANK);
	}
	uint64 get_seed() const { return p_s; }
	string print_seed() const
	{
		ostringstream oss;
		oss << "PRNG seed: 0x" << setw(16) << setfill('0') << hex << p_s;
		return oss.str();
	}

	explicit t_ran(algo_prng algo);
	t_ran(const t_ran&) = delete;
	t_ran& operator= (const t_ran&) = delete;
	~t_ran()
	{
		posix_memalign_free(p_state);
	}
};

extern t_ran ran;

#endif
