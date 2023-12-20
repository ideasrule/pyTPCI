/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_H
#define VECTORIZE_H

#include "vectorize_reduce.h"
#include "vectorize_exp.h"
#include "vectorize_log.h"
#include "vectorize_sqrt.h"
#include "vectorize_hyper.h"

//
// The t_avx_pool class maintains a pool of scratch arrays that can be used by the vectorization
// routines. The memory is correctly aligned on an AVX boundary.
//
// To use the pool, declare a pointer variable as follows
//
//   #include "vectorize.h"
//   avx_ptr<realnum> a(n), b(nlo, nhi), c(nlo, nhi);
//
// You can index the resulting pointers as you would a normal array
//
//   realnum c = a[0]*b[2];
//
// Indexing the arrays is valid for the following indices:
//
//   a[i] : for 0 <= i < n
//   b[i] : for nlo <= i < nhi
//
// In the case of the array b[], the element b[nlo] will be correctly aligned on an AVX boundary. To
// pass the scratch arrays to the vectorization routine, use the following:
//
//   vexp( b.ptr0(), c.ptr0(), nlo, nhi );
//
// When the avx_ptr variable goes out of scope, the scratch array will be released automatically.
// An avx_ptr supports bounds checking when the BOUNDS_CHECK macro is defined at compile time
// (similar to multi_arr and flex_arr).
//
// The t_avx_pool class maintains a minimum pool of p_min_size arrays. The size of the pool can grow
// when needed. When a scratch array is released with an index number >= p_min_size, the memory will
// be freed immediately to conserve memory.
//
// The default scratch array can hold p_def_size doubles. In the default setup that should be larger
// than rfield.nflux_with_check for efficiency since an array of that size is plausible to be
// requested. If the frequency mesh is larger (e.g. due to the SET CONTINUUM RESOLUTION command) the
// following optimization will be used. The first time an array larger than p_def_size doubles is
// requested, it will have to allocated. When that array is subsequently released, it is swapped
// with an unused smaller array with index < p_min_size, and then the smaller array is freed. This
// gives near optimal performance even for large frequency grids, at the expense of a slight
// overhead when releasing the scratch array.
//

class t_avx_pool
{
	vector<void*> p_ptr;
	vector<size_t> p_size;
	vector<bool> p_used;

	// for efficiency this number should be larger than rfield.nflux_with_check
	// in the default setup -- adjust this constant if that is not the case.
	static const size_t p_def_size = 8500*sizeof(double);
	static const size_t p_min_size = 30;

	void p_alloc(size_t sz, bool lgUsed)
	{
		void *p_ptr_alloc;
		if( posix_memalign(&p_ptr_alloc, CD_ALIGN, sz) != 0 )
			throw bad_alloc();
		p_ptr.push_back(p_ptr_alloc);
		p_size.push_back(sz);
		p_used.push_back(lgUsed);
	}

public:
	t_avx_pool()
	{
		for( size_t i=0; i < p_min_size; ++i )
			p_alloc(p_def_size,false);
	}
	~t_avx_pool()
	{
		for( size_t i=0; i < p_ptr.size(); ++i )
			posix_memalign_free(p_ptr[i]);
	}
	void* avx_alloc(size_t sz)
	{
		void* p_ptr_alloc = NULL;
		for( size_t i=0; i < p_ptr.size(); i++ )
		{
			if( !p_used[i] && sz <= p_size[i] )
			{
				p_ptr_alloc = p_ptr[i];
				p_used[i] = true;
				break;
			}
		}
		if( p_ptr_alloc == NULL )
		{
			p_alloc(sz,true);
			p_ptr_alloc = p_ptr.back();
		}
		return p_ptr_alloc;
	}
	void avx_free(void* p_ptr_alloc)
	{
		size_t i;
		for( i=0; i < p_ptr.size(); i++ )
		{
			if( p_ptr[i] == p_ptr_alloc )
			{
				if( i >= p_min_size )
				{
					// 
					size_t j = p_min_size;
					for( j=0; j < p_min_size; ++j )
					{
						if( !p_used[j] && p_size[j] < p_size[i] )
							break;
					}
					if( j < p_min_size )
					{
						posix_memalign_free(p_ptr[j]);
						p_ptr[j] = p_ptr[i];
						p_size[j] = p_size[i];
					}
					else
					{
						posix_memalign_free(p_ptr[i]);
					}
					p_ptr[i] = NULL;
					p_size[i] = 0;
				}
				p_used[i] = false;
				break;
			}
		}
		// clean up unused entries at the back
		while( p_ptr.back() == NULL )
		{
			p_ptr.pop_back();
			p_size.pop_back();
			p_used.pop_back();
		}
	}
};

extern t_avx_pool avx_pool;

template<class T, bool lgBC=lgBOUNDSCHECKVAL>
class avx_ptr
{
	long p_begin;
	long p_end;
	T* p_ptr_alloc;
	T* p_ptr;

	// make default constructor private since it shouldn't be used
	avx_ptr()
	{
		p_ptr = NULL;
		p_ptr_alloc = NULL;
		p_begin = 0;
		p_end = 0;
	}
	void p_alloc(long begin, long end)
	{
		size_t sz = size_t(max(end-begin,0)*sizeof(T));
		if( sz == 0 )
		{
			p_ptr_alloc = NULL;
			p_ptr = NULL;
			p_begin = 0;
			p_end = 0;
		}		
		else
		{
			p_ptr_alloc = static_cast<T*>(avx_pool.avx_alloc(sz));
			p_ptr = p_ptr_alloc - begin;
			p_begin = begin;
			p_end = end;
		}
	}
public:
	typedef T& reference;
	typedef const T& const_reference;

	explicit avx_ptr(long size)
	{
		p_alloc(0, size);
	}
	avx_ptr(long begin, long end)
	{
		p_alloc(begin, end);
	}
	~avx_ptr()
	{
		if( p_ptr_alloc != NULL )
			avx_pool.avx_free(p_ptr_alloc);
	}
	reference operator[] (long i)
	{
		if( lgBC && ( i < p_begin || i >= p_end ) )
			OUT_OF_RANGE( "avx_ptr::operator[]" );
		return *(p_ptr+i);
	}
	const_reference operator[] (long i) const
	{
		if( lgBC && ( i < p_begin || i >= p_end ) )
			OUT_OF_RANGE( "avx_ptr::operator[]" );
		return *(p_ptr+i);
	}
	T* data()
	{
		return p_ptr_alloc;
	}
	const T* data() const
	{
		return p_ptr_alloc;
	}
	T* ptr0()
	{
		return p_ptr;
	}
	const T* ptr0() const
	{
		return p_ptr;
	}
};

//
// The allocator_avx class can be used to get a vector with internal data aligned on an AVX boundary
//
// Use the class as follows:
//
//   #include "vectorize.h"
//   vector< realnum, allocator_avx<realnum> > arr;
//
// Now that C++11 is in effect, you can also use the following shorthand:
//
//   vector_avx<realnum> arr;
//

template<class T>
class allocator_avx;

template<>
class allocator_avx<void>
{
public:
	typedef size_t      size_type;
	typedef ptrdiff_t   difference_type;
	typedef void*       pointer;
	typedef const void* const_pointer;
	typedef void        value_type;

	template<class U>
	struct rebind { typedef allocator_avx<U> other; };

	typedef true_type propagate_on_container_move_assignment;
};

template<class T>
class allocator_avx : public allocator<T>
{
public:
	typedef size_t     size_type;
	typedef ptrdiff_t  difference_type;
	typedef T*         pointer;
	typedef const T*   const_pointer;
	typedef T&         reference;
	typedef const T&   const_reference;
	typedef T          value_type;
	
	template<class U>
        struct rebind
        {
		typedef allocator_avx<U> other;
	};

	typedef true_type propagate_on_container_move_assignment;

	allocator_avx() throw() {}

	allocator_avx(const allocator_avx& a) throw() : allocator<T>(a) {}

	template<class U>
	allocator_avx(const allocator_avx<U>&) throw() {}

	~allocator_avx() throw() {}

	pointer allocate(size_type n, typename allocator_avx<void>::const_pointer = NULL)
	{
		void* p;
		if( posix_memalign(&p, CD_ALIGN, n*sizeof(T)) != 0 )
			throw bad_alloc();
		return pointer(p);
	}

	void deallocate(pointer p, size_type) throw()
	{
		posix_memalign_free(p);
	}
};

template<class T, class U>
inline bool operator== (const allocator_avx<T>&, const allocator_avx<U>&)
{
	return true;
}

template<class T>
inline bool operator== (const allocator_avx<T>&, const allocator_avx<T>&)
{
	return true;
}

template<class T, class U>
inline bool operator!= (const allocator_avx<T>&, const allocator_avx<U>&)
{
	return false;
}

template<class T>
inline bool operator!= (const allocator_avx<T>&, const allocator_avx<T>&)
{
	return false;
}

template<typename T>
using vector_avx = typename std::vector<T,allocator_avx<T>>;

#endif
