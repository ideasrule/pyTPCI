/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CONTAINER_CLASSES_H_
#define CONTAINER_CLASSES_H_

//! define the supported memory layouts of multi_arr / flex_arr
//!
//! ARPA_TYPE: arrays of pointers to arrays to ...
//! C_TYPE: row-major order, exactly as used in C (last index runs fastest)
//! FLX_TYPE: the one and only layout used by flex_arr
typedef enum { ARPA_TYPE, C_TYPE, FLX_TYPE, ML_TOP } mem_layout;

// magic numbers to identify each memory layout
static const int32 MA_VERS[ML_TOP] = { 120070905, 220070803, 320071126 };

#ifdef USE_C_TYPE
#define MEM_LAYOUT_VAL C_TYPE
#else
#define MEM_LAYOUT_VAL ARPA_TYPE
#endif

//! basic_pntr - base class for generalization of normal pointers that enables bounds checking
//! it comes with the full set of operators that you would expect for a random access pointer
template<class T, bool lgBC>
class basic_pntr
{
	static const int p_nd = lgBC ? 3 : 1;
	T* p_p[p_nd]; // p[0] current pointer, p[1] lower bound, p[2] upper bound

	T* p_index_checked ( const ptrdiff_t n ) const
	{
		T* t = p_p[0]+n;
		if( lgBC ) 
		{
			if( t < p_p[1] || t >= p_p[2] )
				OUT_OF_RANGE( "basic_pntr::p_index_checked()" );
		}
		return t;
	}
	void p_set_vals( T* p0, T* p1, T* p2 )
	{
		if( lgBC )
		{
			p_p[0] = p0; p_p[1] = p1; p_p[2] = p2;
		}
		else
			p_p[0] = p0;
	}
public:
	typedef random_access_iterator_tag iterator_category;
	typedef T            value_type;
	typedef T&           reference;
	typedef const T&     const_reference;
	typedef T*           pointer;
	typedef const T*     const_pointer;
	typedef size_t       size_type;
	typedef ptrdiff_t    difference_type;

	// constructors
	basic_pntr( T* p0, T* p1, T* p2 )
	{
		p_set_vals( p0, p1, p2 );
	}
	basic_pntr( T* p0 )
	{
		p_set_vals( p0, NULL, NULL );
	}
	basic_pntr()
	{
		p_set_vals( NULL, NULL, NULL );
	}
	basic_pntr& operator= ( const basic_pntr& ) = default;
	basic_pntr( const basic_pntr& t ) = default;
	// protected destructor (to prevent direct instantiation or destruction via basic_pntr*, see Alexandrescu p13)
protected:
	~basic_pntr() {} 
public:
	// pre-increment
	basic_pntr& operator++ ()
	{
		++p_p[0];
		return *this;
	}
	// pre-decrement
	basic_pntr& operator-- ()
	{
		--p_p[0];
		return *this;
	}
	// define operators for += and -=, normal arithmetic is defined separately in
	// the derived classes; it cannot be done here since they would require implicit
	// conversion from basic_pntr -> pntr or const_pntr to work; this would also create
	// an implicit and silent conversion from const_pntr -> pntr, which is illegal...
	basic_pntr& operator+= ( const ptrdiff_t n ) { p_p[0] += n; return *this; }
	basic_pntr& operator-= ( const ptrdiff_t n ) { p_p[0] -= n; return *this; }
	// dereference
	T& operator* () const
	{
		return *(p_index_checked(0));
	}
	T* operator-> () const
	{
		return p_index_checked(0);
	}
	T& operator[] ( const ptrdiff_t n ) const
	{
		return *(p_index_checked(n));
	}
	// finally, define the boolean operators...
	bool operator== ( const basic_pntr& t ) const { return p_p[0] == t.p_p[0]; }
	bool operator!= ( const basic_pntr& t ) const { return p_p[0] != t.p_p[0]; }
	bool operator<  ( const basic_pntr& t ) const { return p_p[0] <  t.p_p[0]; }
	bool operator<= ( const basic_pntr& t ) const { return p_p[0] <= t.p_p[0]; }
	bool operator>  ( const basic_pntr& t ) const { return p_p[0] >  t.p_p[0]; }
	bool operator>= ( const basic_pntr& t ) const { return p_p[0] >= t.p_p[0]; }
};

//! pntr - interface class to replace normal pointers
template<class T, bool lgBC>
class pntr : public basic_pntr<T,lgBC>
{
public:
	// constructors are not inherited, so define them again
	pntr( T* p0 ) : basic_pntr<T,lgBC>( p0 ) {}
	pntr( T* p0, T* p1, T* p2 ) : basic_pntr<T,lgBC>( p0, p1, p2 ) {}
	pntr() {}
	// the increment / decrement operators need to be recast...
	// otherwise expressions like p = ++q would be illegal for iterators...
	pntr& operator++ () { return static_cast<pntr&>(basic_pntr<T,lgBC>::operator++()); }
	const pntr operator++ (int) { pntr t = *this; ++(*this); return t; }
	pntr& operator-- () { return static_cast<pntr&>(basic_pntr<T,lgBC>::operator--()); }
	const pntr operator-- (int) { pntr t = *this; --(*this); return t; }
	// define p+n, p-n, p-q
	const pntr operator+ ( const ptrdiff_t n ) const { pntr s = *this; s += n; return s; }
	const pntr operator- ( const ptrdiff_t n ) const { pntr s = *this; s -= n; return s; }
	ptrdiff_t operator- ( const pntr& t ) const { return &(*this[0]) - &t[0]; }
};

// this defines n+p
template<class T, bool lgBC>
inline const pntr<T,lgBC> operator+ ( const ptrdiff_t n, const pntr<T,lgBC>& t )
{
	pntr<T,lgBC> s = t;
	s += n;
	return s;
}

//! const_pntr - same as pntr, except that it replaces const pointers rather than normal pointers
template<class T, bool lgBC>
class const_pntr : public basic_pntr<T,lgBC>
{
public:
	// constructors are not inherited, so define them again
	const_pntr( T* p0 ) : basic_pntr<T,lgBC>( p0 ) {}
	const_pntr( T* p0, T* p1, T* p2 ) : basic_pntr<T,lgBC>( p0, p1, p2 ) {}
	const_pntr() {}
	// make sure we can assign a pntr to a const_pntr by creating an implicit conversion to const_pntr
	const_pntr( const pntr<T,lgBC>& t ) : basic_pntr<T,lgBC>( t ) {}
	// the increment / decrement operators need to be recast...
	// otherwise expressions like *p++ = 1. would be legal for const_iterators...
	const_pntr& operator++ () { return static_cast<const_pntr&>(basic_pntr<T,lgBC>::operator++()); }
	const const_pntr operator++ (int) { const_pntr t = *this; ++(*this); return t; }
	const_pntr& operator-- () { return static_cast<const_pntr&>(basic_pntr<T,lgBC>::operator--()); }
	const const_pntr operator-- (int) { const_pntr t = *this; --(*this); return t; }
	const_pntr& operator+= ( const ptrdiff_t n )
	{
		return static_cast<const_pntr&>(basic_pntr<T,lgBC>::operator+=(n));
	}
	const_pntr& operator-= ( const ptrdiff_t n )
	{
		return static_cast<const_pntr&>(basic_pntr<T,lgBC>::operator-=(n));
	}
	// the dereference operators need to be recast...
	const T& operator* () const { return static_cast<const T&>(basic_pntr<T,lgBC>::operator*()); }
	const T* operator-> () const { return static_cast<const T*>(basic_pntr<T,lgBC>::operator->()); }
	const T& operator[] ( const ptrdiff_t n ) const
	{
		return static_cast<const T&>(basic_pntr<T,lgBC>::operator[](n));
	}
	// define p+n, p-n, p-q
	const const_pntr operator+ ( const ptrdiff_t n ) const { const_pntr s = *this; s += n; return s; }
	const const_pntr operator- ( const ptrdiff_t n ) const { const_pntr s = *this; s -= n; return s; }
	ptrdiff_t operator- ( const const_pntr& t ) const { return &(*this[0]) - &t[0]; }
};

// this defines n+p
template<class T, bool lgBC>
inline const const_pntr<T,lgBC> operator+ ( const ptrdiff_t n, const const_pntr<T,lgBC>& t )
{
	const_pntr<T,lgBC> s = t;
	s += n;
	return s;
}

//! tree_vec - a simple class to store the bounds checking information for multi_arr
struct tree_vec
{
	typedef size_t size_type;

	size_type n;
	tree_vec *d;

private:
	void p_clear0()
	{
		if( d != NULL )
		{
			for( size_type i = 0; i < n; ++i )
				d[i].clear();
			delete[] d;
		}
	}
	void p_clear1()
	{
		n = 0;
		d = NULL;
	}

public:
	tree_vec()
	{
		p_clear1();
	}
	tree_vec(const tree_vec& m)
	{
		p_clear1();
		*this = m;
	}
	~tree_vec()
	{
		p_clear0();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}
	const tree_vec& operator= (const tree_vec& m)
	{
		if( &m != this )
		{
			clear();
			n = m.n;
			if( m.d != NULL )
			{
				d = new tree_vec[n];
				tree_vec *p = d;
				const tree_vec *mp = m.d;
				for( size_type i = 0; i < n; ++i )
					*p++ = *mp++;
			}
		}
		return *this;
	}
	tree_vec& getvec(const size_type i, const size_type index[])
	{
		if( i == 0 )
			return *this;
		else
			return getvec(i-1,index).d[index[i-1]];
	}
	const tree_vec& getvec(const size_type i, const size_type index[]) const
	{
		if( i == 0 )
			return *this;
		else
			return getvec(i-1,index).d[index[i-1]];
	}
};

//! multi_geom - this class maintains all the geometry information for multi_arr
//! keeping it separate makes it easy to clone the information from one multi_arr to another
template<int d,mem_layout ALLOC=MEM_LAYOUT_VAL>
class multi_geom
{
public:
	typedef size_t size_type;

	tree_vec v;

	size_type size;  //! allocated size (number of data elements, pointers are not counted)
	size_type s[d];  //! size of each dimension (only used in C_TYPE layout)
	size_type st[d]; //! stride for each dimension (only used in C_TYPE layout)
	size_type nsl[d];//! sizes of each of the pointer arrays

private:
	void p_clear0()
	{
		v.clear();
	}
	void p_clear1()
	{
		size = 0;
		for( int i=0; i < d; ++i )
		{
			s[i] = 0;
			st[i] = 0;
			nsl[i] = 0;
		}
	}

public:
	multi_geom()
	{
		p_clear1();
	}
	multi_geom(const multi_geom& m)
	{
		p_clear1();
		*this = m;
	}
	~multi_geom()
	{
		p_clear0();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}
	const multi_geom& operator= (const multi_geom& m)
	{
		if( &m != this )
		{
			clear();
			v = m.v;
			size = m.size;
			for( int i=0; i < d; ++i )
			{
				s[i] = m.s[i];
				st[i] = m.st[i];
				nsl[i] = m.nsl[i];
			}
		}
		return *this;
	}
	bool lgInbounds(const size_type n, const size_type index[]) const
	{
		if( n != 0 )
			return ( lgInbounds(n-1,index) && index[n-1] < v.getvec(n-1,index).n );
		else
			return true;
	}
	void reserve(const size_type n, const size_type index[])
	{
		ASSERT( n <= d && index[n-1] > 0 && lgInbounds( n-1, index ) );

		tree_vec& w = v.getvec( n-1, index );
		if( d > n )
		{
			ASSERT( w.d == NULL );
			w.d = new tree_vec[ index[n-1] ];
		}
		w.n = index[n-1];
		s[n-1] = max(s[n-1],index[n-1]);
		nsl[n-1] += index[n-1];
	}
	void reserve_recursive(const size_type n, size_type index[])
	{
		if( n == 0 )
		{
			reserve( n+1, index );
			if( n+1 < d )
				reserve_recursive( n+1, index );
		}
		else
		{
			size_type top = index[n-1];
			for( size_type i=0; i < top; ++i )
			{
				index[n-1] = i;
				reserve( n+1, index );
				if( n+1 < d )
					reserve_recursive( n+1, index );
			}
			index[n-1] = top;
		}
	}
	void finalize(void)
	{
		static_assert( ALLOC == C_TYPE || ALLOC == ARPA_TYPE,
			       "Allocation type must be C_TYPE or ARPA_TYPE" );
		if( ALLOC == ARPA_TYPE )
		{
			size_type n1[d], n2[d];
			for( int dim=0; dim < d; ++dim )
				n1[dim] = n2[dim] = 0L;
			// sanity checks
			p_setupArray( n1, n2, &v, 0 );
			for( int dim=0; dim < d-1; ++dim )
				ASSERT( n1[dim] == nsl[dim] && n2[dim] == nsl[dim+1] );
			size = nsl[d-1];
		}
		else if( ALLOC == C_TYPE )
		{
			st[d-1] = s[d-1]; 
			for( int i = d-2; i >= 0; --i )
				st[i] = st[i+1]*s[i];
			size = st[0];
		}
		else
		{
			TotalInsanity();
		}
	}

private:
	void p_setupArray( size_type n1[], size_type n2[], const tree_vec* w, size_type l )
	{
		for( size_type i=0; i < w->n; ++i )
		{
			n1[l]++;
			// using int(l) is always safe and avoids warnings about the test
			// 'l < d-2' being always false when d == 2 (e.g. with g++ 4.6.0)
			if( int(l) < d-2 )
			{
				p_setupArray( n1, n2, &w->d[i], l+1 );
			}
			n2[l] += w->d[i].n;
		}
	}
};


//
//! The n_pointer and const_n_pointer classes below define indexing into the multi_arr
//!
//! NB NB -- The design of these classes is CRUCIAL for CPU performance! Any redesign
//!          should be thoroughly tested for CPU time regressions on a broad range of
//!          platforms and compilers!
//!
//! The primary goal of the design of these classes should be to make life as easy as
//! possible for the compiler, NOT the programmer! The design should assure that that
//! the methods are FULLY INLINED. Without that the resulting code will be severely
//! crippled!
//!
//! Below are a total of 8 specializations for each of these classes. These could be
//! combined into 2 specialization for each class (one for N > 1 and one for N == 1),
//! but it turns out that the resulting method is too complicated for most compilers
//! to inline. Hence we decided to split the methods up into simpler specializations.
//! As they are written below, nearly every compiler can inline them (at least for the
//! non-bounds_checking versions)...
//!
//! This is the situation on 2008 jan 26:
//!
//! Linux IA32/AMD64   ARPA    C    ARPA/BC  C/BC
//!    g++ 4.3.x        OK     OK     OK      OK   (tested on prerelease)
//!    g++ 4.2.x        OK     OK     OK      OK
//!    g++ 4.1.x        OK     OK     OK      OK
//!    g++ 4.0.x        OK     OK     OK      OK
//!    g++ 3.4.x        ??     OK     ??      ??
//!    g++ 3.3.x        ??     OK     ??      ??
//!
//!    icc 10.1         OK     OK     OK      OK
//!    icc 10.0         OK     OK     OK      OK
//!    icc  9.1         OK     OK     OK      OK
//!    icc  9.0         OK     OK     OK      OK
//!
//!    Sun Studio 12    OK     OK     OK      OK
//!
//!    Pathscale 2.3.1  OK     OK     --      --   (only tested on AMD64)
//!
//!    Portland 6.1-5   --     --     --      --   (only tested on IA32)
//!
//! Linux IA64
//!    g++ 4.2.x        OK     OK     OK      OK
//!    g++ 3.3.x        ??     OK?    ??      ??
//!
//!    icc 10.1         OK     OK     --      --
//!    icc 10.0         OK     OK     --      --
//!    icc  9.1         OK     OK     --      --
//!    icc  9.0         OK     OK     OK      OK
//!    icc  8.1         OK     OK     OK      OK
//!
//! Windows IA32
//!    Visual Studio 8  OK     OK     --      --
//!    icl 10.1         OK     OK     --      --
//!
//! Solaris Ultrasparc
//!    g++ 4.3.x        OK     OK     OK      OK   (tested on prerelease)
//!    g++ 4.2.x        OK     OK     OK      OK
//!    g++ 4.1.x        OK     OK     OK      OK
//!    g++ 4.0.x        OK     OK     OK      OK
//!    g++ 3.4.x        ??     ??     ??      ??
//!    g++ 3.3.x        ??     ??     ??      ??
//!
//!    Sun Studio 12    OK     OK     OK      OK
//!
//! Notes:
//!   - The tests were carried out on the following simple test code:
//!
//!       long test(multi_arr<long,6>& arr)
//!       {
//!         return arr[2][3][4][5][6][7];
//!       }
//!
//!     The tests were done both on ARPA_TYPE and C_TYPE arrays, with and without
//!     bounds-checking enabled (BC indicates bounds-checking enabled). The resulting
//!     (optimized) assembly was inspected by eye.
//!
//!   - "OK" indicates that all methods were inlined and the resulting assembly looked
//!     optimal; "??" indicates that all methods were inlined, but the assembly looked
//!     sub-optimal (i.e. not all overhead was optimized away); "--" indicates that at
//!     least some of the methods were not inlined at all.
//!
//! Indexing a multi_arr works as follows. When the compiler encounters arr[i][j][k],
//! it will first emit a to call multi_arr::operator[] with i as its argument; this
//! operator will construct a temporary n_pointer as follows
//!
//!     arr[i]  -->  n_pointer<T,3>(*p,st[],*v)::operator[](i) -> n_pointer<T,2>
//!
//! here p is the base pointer to the data, st contains the strides of a C_TYPE array
//! (which are not used in ARPA_TYPE arrays) and v points to the tree_vec with the
//! bounds-checking information. The operator[] of the n_pointer will update p and v
//! using the value of i. This will only take care of the first index. At this stage
//! the compiler still needs to emit code for the remainder: n_pointer<T,2>[j][k].
//! The operator[] of the n_pointer<T,2> will take care of the next index:
//!
//!     n_pointer<T,2>::operator[](j)  -->  n_pointer<T,1>
//!
//! So another temporary n_pointer is emitted, and p and v are updated again by the
//! operator[] using j. Now the compiler still needs to emit code for: n_pointer<T,1>[k];
//! n_pointer<T,1> is special since it will absorb the last index, so it should not
//! return yet another n_pointer but a (const) reference to the data item itself:
//!
//!     n_pointer<T,1>::operator[](k)  -->  T& *(p + k)
//!     const_n_pointer<T,1>::operator[](k)  -->  const T& *(p + k)
//


// forward definitions
template<class T, int N, mem_layout ALLOC, bool lgBC> class n_pointer;
template<class T, int N, mem_layout ALLOC, bool lgBC> class const_n_pointer;

template<class T, int N>
class n_pointer<T,N,ARPA_TYPE,false>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
   typedef n_pointer<T,N-1,ARPA_TYPE,false> value_type;
	n_pointer(T* p, const size_t* st=NULL, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		return value_type( *((T**)p_p+i) );
	}
};

template<class T, int N>
class n_pointer<T,N,C_TYPE,false>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef n_pointer<T,N-1,C_TYPE,false> value_type;
	n_pointer(T* p, const size_t* st, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		return value_type( p_p+i*p_st[0], p_st+1 );
	}
};

template<class T>
class n_pointer<T,1,ARPA_TYPE,false>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef T value_type;
	n_pointer(T* p, const size_t* st=NULL, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		return *(p_p + i);
	}
};

template<class T>
class n_pointer<T,1,C_TYPE,false>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef T value_type;
	n_pointer(T* p, const size_t* st, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		return *(p_p + i);
	}
};

template<class T, int N>
class n_pointer<T,N,ARPA_TYPE,true>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef n_pointer<T,N-1,ARPA_TYPE,true> value_type;
	n_pointer(T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "n_pointer::operator[]" );
		return value_type( *((T**)p_p+i), NULL, &p_v->d[i] );
	}
};

template<class T, int N>
class n_pointer<T,N,C_TYPE,true>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef n_pointer<T,N-1,C_TYPE,true> value_type;
	n_pointer(T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "n_pointer::operator[]" );
		return value_type( p_p+i*p_st[0], p_st+1, &p_v->d[i] );
	}
};

template<class T>
class n_pointer<T,1,ARPA_TYPE,true>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef T value_type;
	n_pointer(T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "n_pointer::operator[]" );
		return *(p_p + i);
	}
};

template<class T>
class n_pointer<T,1,C_TYPE,true>
{
	T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef T value_type;
	n_pointer(T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "n_pointer::operator[]" );
		return *(p_p + i);
	}
};

template<class T, int N>
class const_n_pointer<T,N,ARPA_TYPE,false>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const_n_pointer<T,N-1,ARPA_TYPE,false> value_type;
	const_n_pointer(const T* p, const size_t* st=NULL, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		return value_type( *((T**)p_p+i) );
	}
};

template<class T, int N>
class const_n_pointer<T,N,C_TYPE,false>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const_n_pointer<T,N-1,C_TYPE,false> value_type;
	const_n_pointer(const T* p, const size_t* st, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		return value_type( p_p+i*p_st[0], p_st+1 );
	}
};

template<class T>
class const_n_pointer<T,1,ARPA_TYPE,false>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const T value_type;
	const_n_pointer(const T* p, const size_t* st=NULL, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		return *(p_p + i);
	}
};

template<class T>
class const_n_pointer<T,1,C_TYPE,false>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const T value_type;
	const_n_pointer(const T* p, const size_t* st, const tree_vec* v=NULL) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		return *(p_p + i);
	}
};

template<class T, int N>
class const_n_pointer<T,N,ARPA_TYPE,true>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const_n_pointer<T,N-1,ARPA_TYPE,true> value_type;
	const_n_pointer(const T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "const_n_pointer::operator[]" );
		return value_type( *((T**)p_p+i), NULL, &p_v->d[i] );
	}
};

template<class T, int N>
class const_n_pointer<T,N,C_TYPE,true>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const_n_pointer<T,N-1,C_TYPE,true> value_type;
	const_n_pointer(const T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	const value_type operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "const_n_pointer::operator[]" );
		return value_type( p_p+i*p_st[0], p_st+1, &p_v->d[i] );
	}
};

template<class T>
class const_n_pointer<T,1,ARPA_TYPE,true>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const T value_type;
	const_n_pointer(const T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "const_n_pointer::operator[]" );
		return *(p_p + i);
	}
};

template<class T>
class const_n_pointer<T,1,C_TYPE,true>
{
	const T* p_p;
	const size_t* p_st;
	const tree_vec* p_v;
public:
	typedef const T value_type;
	const_n_pointer(const T* p, const size_t* st, const tree_vec* v) : p_p(p), p_st(st), p_v(v) {}
	value_type& operator[] (const size_t i) const
	{
		if( i >= p_v->n )
			OUT_OF_RANGE( "const_n_pointer::operator[]" );
		return *(p_p + i);
	}
};

//
//! multi_arr: generic class for allocating multidimensional arrays.
//! A typical example of its use could be:
//!
//!     multi_arr<double,3> arr; // define a placeholder for the array
//!                       // the first argument is the type of data it holds
//!                       // the second argument is the number of dimensions
//!                       // (between 2 and 6)
//!
//!     arr.alloc(3,4,2); // this will allocate a 3x4x2 block of doubles
//!                       // memory will be allocated as a valarray, so each
//!                       // element will be initialized to T() -- this means
//!                       // that even POD types like double with be zeroed
//! 
//!     multi_arr<double,3> arr(3,4,2); // shorthand for the above
//!
//! The following is an alternative way of allocating the array. It is
//! very similar to the pre-multi_arr way of allocating arrays. In ARPA_TYPE
//! arrays this will help you save memory since only the data elements
//! that are really needed will be allocated. In C_TYPE allocation, the
//! smallest rectangular block will be allocated that can hold all the data.
//! This will use more memory in return for somewhat improved CPU speed.
//! Tests carried out in 2007 showed that the speed advantage of C_TYPE arrays
//! was only 1 to 2 percent. Hence the memory savings were deemed more important
//! and ARPA_TYPE arrays were made the default. However, C_TYPE arrays are
//! guaranteed to be compatible with C code, so these should be used if they
//! are meant to be passed on to C library routines. The example below allocates
//! a triangular matrix.
//!
//!     arr.reserve(3);
//!     for( int i=0, i < 3; ++i )
//!     {
//!       arr.reserve( i, i+1 ); // note that size does not need to be constant!
//!       for( int j=0, j < i+1; ++j )
//!         arr.reserve( i, j, j+1 );
//!     }
//!     arr.alloc();
//!
//! these are plausible ways to use the multi_arr class:
//!        
//!     arr.invalidate(); // this will set float or double arrays to all SNaN
//!                       // it will set any other type array to all 0xff bytes.
//!     arr.zero();       // this will set the array to all zero
//!     arr = -1;         // this will set the array to all -1.
//!
//!     arr[0][0][0] = 1.;
//!     arr[0][0][1] = 2.;
//!     double x = arr[0][0][0]*arr[0][0][1];
//!
//!     multi_arr<double,2,C_TYPE> a(10,10); // allocate C_TYPE array
//!     C_library_routine( a.data(), ... );  // and call C library routine with it
//!
//!     arr.clear();      // this will deallocate the array
//!                       // the destructor will also automatically deallocate
//!
//! the multi_arr class comes with iterators that allow you to speed up memory access
//! even further. using iterators as shown below will generally speed up the code
//! significantly since it avoids calculating the array index over and over inside the
//! body of the loop. especially in tight loops over arrays with high dimension this can
//! become a significant overhead! a const_iterator is also supplied for read-only access,
//! but no reverse_iterators. you can define and initialize an iterator as follows
//!
//!     multi_arr<double,3>::iterator p = arr.begin(n,k);
//!
//! the notation multi_arr<double,3>::iterator is rather cumbersome, so it may be
//! convenient to define something like:
//!
//!     typedef multi_arr<double,3>::iterator md3i;
//!     typedef multi_arr<double,3>::const_iterator md3ci;
//!
//! all the possible combinations for bool, long, realnum and double multi_arr's are
//! predefined below.
//!
//! this is a plausible way to use an iterator:
//!
//!     for( int k=0; i < 4; k++ )
//!     {
//!       for( md3i p = arr.begin(n,k); p != arr.end(n,k); ++p )
//!         *p = 3.;
//!     }
//!
//! however, since many compilers have a hard time figuring out that arr.end() has no
//! side effects, it is better to do the following:
//!
//!     for( int k=0; i < 4; k++ )
//!     {
//!       md3i end = arr.end(n,k);
//!       for( md3i p = arr.begin(n,k); p != end; ++p )
//!         *p = 3.;
//!     }
//!
//! NB NB -- the memory layout may change in future editions, so user code should not
//!          make any assumptions about the layout. the only exception is that the user
//!          may safely assume that for the default memory layout the last index
//!          runs over contiguous memory. this allows for efficient iterator access.
//!          the example above was OK since arr[n][k][0] and arr[n][k][1] are
//!          guaranteed to be adjacent. however the next example is not OK:
//!
//! !! WRONG !!, arr[n][k-1][1] and arr[n][k][0] may NOT be adjacent
//!
//!     md3i p = arr.begin(n,0);
//!     for( int k=0; i < 4; k++ )
//!       for( int i=0; i < 2; i++ )
//!         *p++ = 3.;   // ERROR, this may segfault.
//!                      // bounds checking will catch this (see below for enabling this)
//!
//! you can also use iterators for array-like access via []:
//!
//!     double sum = 0.;
//!     for( int k=0; i < 4; k++ )
//!     {
//!       md3ci p = arr.begin(n,k);
//!       for( int i=0; i < 2; i++ )
//!         sum += p[i];
//!     }
//!
//! last, but not least, the multi_arr class supports array bounds checking, both
//! for direct access through the indexing method, as well as iterator access. To enable
//! bounds checking, simply define the preprocessor macro BOUNDS_CHECK during compilation.
//! the resulting code will be MUCH slower, so this should only be used as a debugging tool.
//

template<class T, int d, mem_layout ALLOC=MEM_LAYOUT_VAL, bool lgBC=lgBOUNDSCHECKVAL>
class multi_arr
{
	// ancillary data describing the memory layout of the multi_arr
	multi_geom<d,ALLOC> p_g;
	T** p_psl[d-1];     // pointer arrays for ARPA structure
	valarray<T> p_dsl;  // this contains the actual data
	T* p_ptr;           // main pointer to allocated structure
	T** p_ptr2;         // used in debugger to get access to internal representation
	T*** p_ptr3;
	T**** p_ptr4;
	T***** p_ptr5;
	T****** p_ptr6;

public:
	typedef random_access_iterator_tag iterator_category;
	typedef T            value_type;
	typedef T&           reference;
	typedef const T&     const_reference;
	typedef T*           pointer;
	typedef const T*     const_pointer;
	typedef size_t       size_type;
	typedef ptrdiff_t    difference_type;
	typedef pntr<T,lgBC>       iterator;
	typedef const_pntr<T,lgBC> const_iterator;

private:
	static const size_type npos = static_cast<size_type>(-1);

	void p_clear0()
	{
		p_g.clear();
		for( int i=0; i < d-1; ++i )
			delete[] p_psl[i];
		p_dsl.resize(0);
	}
	void p_clear1()
	{
		for( int i=0; i < d-1; ++i )
			p_psl[i] = NULL;
		p_ptr = NULL;
		p_ptr2 = NULL;
		p_ptr3 = NULL;
		p_ptr4 = NULL;
		p_ptr5 = NULL;
		p_ptr6 = NULL;
	}

public:
	multi_arr()
	{
		p_clear1();
	}
	multi_arr(const multi_geom<d,ALLOC>& g)
	{
		p_clear1();
		alloc( g );
	}
	multi_arr(size_type d1, size_type d2)
	{
		p_clear1();
		size_type index[] = { d1, d2 };
		alloc( index );
	}
	multi_arr(size_type d1, size_type d2, size_type d3)
	{
		p_clear1();
		size_type index[] = { d1, d2, d3 };
		alloc( index );
	}
	multi_arr(size_type d1, size_type d2, size_type d3, size_type d4)
	{
		p_clear1();
		size_type index[] = { d1, d2, d3, d4 };
		alloc( index );
	}
	multi_arr(size_type d1, size_type d2, size_type d3, size_type d4, size_type d5)
	{
		p_clear1();
		size_type index[] = { d1, d2, d3, d4, d5 };
		alloc( index );
	}
	multi_arr(size_type d1, size_type d2, size_type d3, size_type d4, size_type d5, size_type d6)
	{
		p_clear1();
		size_type index[] = { d1, d2, d3, d4, d5, d6 };
		alloc( index );
	}
	multi_arr(const multi_arr& m)
	{
		p_clear1();
		*this = m;
	}
	~multi_arr()
	{
		p_clear0();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}
	const multi_arr& operator= (const multi_arr& m)
	{
		if( &m != this )
		{
			alloc( m.p_g );
			vals() = m.vals();
		}
		return *this;
	}
	const multi_arr& operator= (const T& val)
	{
		p_dsl = val;
		return *this;
	}
	void zero()
	{
		ASSERT( vals().size() == p_g.size );
		if( p_g.size > 0 )
			memset( data(), 0, p_g.size*sizeof(T) );
	}
	void invalidate()
	{
		ASSERT( vals().size() == p_g.size );
		invalidate_array( data(), p_g.size*sizeof(T) );
	}

	void reserve(size_type i1)
	{
		ASSERT( vals().size() == 0 );
		const size_type index[] = { i1 };
		p_g.reserve( 1, index );
	}
	void reserve(size_type i1, size_type i2)
	{
		ASSERT( vals().size() == 0 );
		const size_type index[] = { i1, i2 };
		p_g.reserve( 2, index );
	}
	void reserve(size_type i1, size_type i2, size_type i3)
	{
		ASSERT( vals().size() == 0 );
		const size_type index[] = { i1, i2, i3 };
		p_g.reserve( 3, index );
	}
	void reserve(size_type i1, size_type i2, size_type i3, size_type i4)
	{
		ASSERT( vals().size() == 0 );
		const size_type index[] = { i1, i2, i3, i4 };
		p_g.reserve( 4, index );
	}
	void reserve(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5)
	{
		ASSERT( vals().size() == 0 );
		const size_type index[] = { i1, i2, i3, i4, i5 };
		p_g.reserve( 5, index );
	}
	void reserve(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5, size_type i6)
	{
		ASSERT( vals().size() == 0 );
		const size_type index[] = { i1, i2, i3, i4, i5, i6 };
		p_g.reserve( 6, index );
	}
	void alloc()
	{
		p_g.finalize();
		static_assert( ALLOC == C_TYPE || ALLOC == ARPA_TYPE,
			       "Allocation type must be C_TYPE or ARPA_TYPE" );

		if( ALLOC == ARPA_TYPE )
		{
			size_type n1[d], n2[d];
			// allocate the pointer arrays ( p_psl[0..d-2] ) and data ( p_dsl )
			for( int dim=0; dim < d; ++dim )
			{
				n1[dim] = n2[dim] = 0L;
				if( dim != d-1 )
				{
					ASSERT( p_psl[dim] == NULL );
					if( p_g.nsl[dim] > 0 )
						p_psl[dim] = new T*[ p_g.nsl[dim] ];
				}
				else
				{
					ASSERT( p_dsl.size() == 0 );
					p_dsl.resize( p_g.nsl[dim] );
				}
			}
			// now initialize all the pointer arrays
			p_setupArray( n1, n2, &p_g.v, 0 );
			p_ptr = (T*)p_psl[0];
		}
		else if( ALLOC == C_TYPE )
		{
			for( int i=0; i < d-1; ++i )
				p_psl[i] = NULL;
			ASSERT( p_dsl.size() == 0 );
			p_dsl.resize( p_g.st[0] );
			p_ptr = &p_dsl[0];
		}
		else
		{
			TotalInsanity();
		}
		p_ptr2 = (T**)p_ptr;
		p_ptr3 = (T***)p_ptr;
		p_ptr4 = (T****)p_ptr;
		p_ptr5 = (T*****)p_ptr;
		p_ptr6 = (T******)p_ptr;
	}
	// clone the geometry from another multi_arr
	void alloc(const multi_geom<d,ALLOC>& g)
	{
		if( &g != &p_g )
		{
			clear();
			p_g = g;
			alloc();
		}
	}
	// set up a rectangular block of data with dimensions d1 x d2 x ....
	void alloc(size_type d1, size_type d2)
	{
		size_type index[] = { d1, d2 };
		alloc( index );
	}
	void alloc(size_type d1, size_type d2, size_type d3)
	{
		size_type index[] = { d1, d2, d3 };
		alloc( index );
	}
	void alloc(size_type d1, size_type d2, size_type d3, size_type d4)
	{
		size_type index[] = { d1, d2, d3, d4 };
		alloc( index );
	}
	void alloc(size_type d1, size_type d2, size_type d3, size_type d4, size_type d5)
	{
		size_type index[] = { d1, d2, d3, d4, d5 };
		alloc( index );
	}
	void alloc(size_type d1, size_type d2, size_type d3, size_type d4, size_type d5, size_type d6)
	{
		size_type index[] = { d1, d2, d3, d4, d5, d6 };
		alloc( index );
	}
	void alloc(size_type index[])
	{
		for( int n=0; n < d; n++ )
			ASSERT( index[n] > 0 );
		clear();
		p_g.reserve_recursive( 0, index );
		alloc();
	}

private:
	// helper routine for alloc(), this fills in the pointer arrays for the ARPA layout
	void p_setupArray( size_type n1[], size_type n2[], const tree_vec* g, size_type l )
	{
		for( size_type i=0; i < g->n; ++i )
		{
			// using int(l) is always safe and avoids warnings about the test
			// 'l < d-2' being always false when d == 2 (e.g. with g++ 4.6.0)
			if( int(l) < d-2 )
			{
				p_psl[l][n1[l]++] = (T*)(p_psl[l+1]+n2[l]);
				p_setupArray( n1, n2, &g->d[i], l+1 );
			}
			else
			{
				p_psl[l][n1[l]++] = &p_dsl[0]+n2[l];
			}
			n2[l] += g->d[i].n;
		}
	}

	// in the p_iterator methods the bound-checking part is split off into a separate
	// routine p_iterator_bc in order to make it easier for compilers to inline the code
	iterator p_iterator(size_type i1, size_type i2) const
	{
		if( lgBC )
			return p_iterator_bc( i1, i2 );
		else
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			return iterator( &(*t)[i1][i2] );
		}
	}
	iterator p_iterator_bc(size_type i1, size_type i2) const
	{
		size_type index[] = { i1 };
		if( p_g.lgInbounds( 1, index ) )
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			size_type n = p_g.v.getvec( 1, index ).n;
			T* s = ( n > 0 ) ? &(*t)[i1][0] : NULL;
			if( i2 == npos )
				return iterator( s+n, s, s+n );
			else
				return iterator( s+i2, s, s+n );
		}
		else
			OUT_OF_RANGE( "multi_arr::p_iterator()" );
	}
	iterator p_iterator(size_type i1, size_type i2, size_type i3) const
	{
		if( lgBC )
			return p_iterator_bc( i1, i2, i3 );
		else
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			return iterator( &(*t)[i1][i2][i3] );
		}
	}
	iterator p_iterator_bc(size_type i1, size_type i2, size_type i3) const
	{
		size_type index[] = { i1, i2 };
		if( p_g.lgInbounds( 2, index ) )
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			size_type n = p_g.v.getvec( 2, index ).n;
			T* s = ( n > 0 ) ? &(*t)[i1][i2][0] : NULL;
			if( i3 == npos )
				return iterator( s+n, s, s+n );
			else
				return iterator( s+i3, s, s+n );
		}
		else
			OUT_OF_RANGE( "multi_arr::p_iterator()" );
	}
	iterator p_iterator(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		if( lgBC )
			return p_iterator_bc(i1, i2, i3, i4);
		else
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			return iterator( &(*t)[i1][i2][i3][i4] );
		}
	}
	iterator p_iterator_bc(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		size_type index[] = { i1, i2, i3 };
		if( p_g.lgInbounds( 3, index ) )
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			size_type n = p_g.v.getvec( 3, index ).n;
			T* s = ( n > 0 ) ? &(*t)[i1][i2][i3][0] : NULL;
			if( i4 == npos )
				return iterator( s+n, s, s+n );
			else
				return iterator( s+i4, s, s+n );
		}
		else
			OUT_OF_RANGE( "multi_arr::p_iterator()" );
	}
	iterator p_iterator(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		if( lgBC )
			return p_iterator_bc(i1, i2, i3, i4, i5);
		else
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			return iterator( &(*t)[i1][i2][i3][i4][i5] );
		}
	}
	iterator p_iterator_bc(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		size_type index[] = { i1, i2, i3, i4 };
		if( p_g.lgInbounds( 4, index ) )
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			size_type n = p_g.v.getvec( 4, index ).n;
			T* s = ( n > 0 ) ? &(*t)[i1][i2][i3][i4][0] : NULL;
			if( i5 == npos )
				return iterator( s+n, s, s+n );
			else
				return iterator( s+i5, s, s+n );
		}
		else
			OUT_OF_RANGE( "multi_arr::p_iterator()" );
	}
	iterator p_iterator(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5, size_type i6) const
	{
		if( lgBC )
			return p_iterator_bc(i1, i2, i3, i4, i5, i6);
		else
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			return iterator( &(*t)[i1][i2][i3][i4][i5][i6] );
		}
	}
	iterator p_iterator_bc(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5, size_type i6) const
	{
		size_type index[] = { i1, i2, i3, i4, i5 };
		if( p_g.lgInbounds( 5, index ) )
		{
			multi_arr<T,d,ALLOC,lgBC>* t = const_cast<multi_arr<T,d,ALLOC,lgBC>*>(this);
			size_type n = p_g.v.getvec( 5, index ).n;
			T* s = ( n > 0 ) ? &(*t)[i1][i2][i3][i4][i5][0] : NULL;
			if( i6 == npos )
				return iterator( s+n, s, s+n );
			else
				return iterator( s+i6, s, s+n );
		}
		else
			OUT_OF_RANGE( "multi_arr::p_iterator()" );
	}

public:
	const n_pointer<T,d,ALLOC,lgBC> n_ptr()
	{
		return n_pointer<T,d,ALLOC,lgBC>( p_ptr, p_g.st+1, &p_g.v );
	}
	const const_n_pointer<T,d,ALLOC,lgBC> n_ptr() const
	{
		return const_n_pointer<T,d,ALLOC,lgBC>( p_ptr, p_g.st+1, &p_g.v );
	}
	typedef n_pointer<T,d-1,ALLOC,lgBC> indexed_type;
	const indexed_type operator[] (size_type i)
	{
		return n_ptr()[i];
	}
	typedef const_n_pointer<T,d-1,ALLOC,lgBC> const_indexed_type;
	const const_indexed_type operator[] (size_type i) const
	{
		return n_ptr()[i];
	}

	reference at(size_type i1, size_type i2)
	{
		size_type index[] = { i1, i2 };
		if( !p_g.lgInbounds( 2, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2];
	}
	const_reference at(size_type i1, size_type i2) const
	{
		size_type index[] = { i1, i2 };
		if( !p_g.lgInbounds( 2, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2];
	}
	reference at(size_type i1, size_type i2, size_type i3)
	{
		size_type index[] = { i1, i2, i3 };
		if( !p_g.lgInbounds( 3, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3];
	}
	const_reference at(size_type i1, size_type i2, size_type i3) const
	{
		size_type index[] = { i1, i2, i3 };
		if( !p_g.lgInbounds( 3, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3];
	}
	reference at(size_type i1, size_type i2, size_type i3, size_type i4)
	{
		size_type index[] = { i1, i2, i3, i4 };
		if( !p_g.lgInbounds( 4, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3][i4];
	}
	const_reference at(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		size_type index[] = { i1, i2, i3, i4 };
		if( !p_g.lgInbounds( 4, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3][i4];
	}
	reference at(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5)
	{
		size_type index[] = { i1, i2, i3, i4, i5 };
		if( !p_g.lgInbounds( 5, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3][i4][i5];
	}
	const_reference at(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		size_type index[] = { i1, i2, i3, i4, i5 };
		if( !p_g.lgInbounds( 5, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3][i4][i5];
	}
	reference at(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5, size_type i6)
	{
		size_type index[] = { i1, i2, i3, i4, i5, i6 };
		if( !p_g.lgInbounds( 6, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3][i4][i5][i6];
	}
	const_reference at(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5, size_type i6) const
	{
		size_type index[] = { i1, i2, i3, i4, i5, i6 };
		if( !p_g.lgInbounds( 6, index ) )
			OUT_OF_RANGE( "multi_arr::at()" );
		return (*this)[i1][i2][i3][i4][i5][i6];
	}

	iterator ptr(size_type i1, size_type i2)
	{
		return p_iterator(i1, i2);
	}
	const_iterator ptr(size_type i1, size_type i2) const
	{
		return p_iterator(i1, i2);
	}
	iterator ptr(size_type i1, size_type i2, size_type i3)
	{
		return p_iterator(i1, i2, i3);
	}
	const_iterator ptr(size_type i1, size_type i2, size_type i3) const
	{
		return p_iterator(i1, i2, i3);
	}
	iterator ptr(size_type i1, size_type i2, size_type i3, size_type i4)
	{
		return p_iterator(i1, i2, i3, i4);
	}
	const_iterator ptr(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		return p_iterator(i1, i2, i3, i4);
	}
	iterator ptr(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5)
	{
		return p_iterator(i1, i2, i3, i4, i5);
	}
	const_iterator ptr(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		return p_iterator(i1, i2, i3, i4, i5);
	}
	iterator ptr(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5, size_type i6)
	{
		return p_iterator(i1, i2, i3, i4, i5, i6);
	}
	const_iterator ptr(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5, size_type i6) const
	{
		return p_iterator(i1, i2, i3, i4, i5, i6);
	}

	iterator begin(size_type i1)
	{
		return p_iterator(i1, 0);
	}
	const_iterator begin(size_type i1) const
	{
		return p_iterator(i1, 0);
	}
	iterator begin(size_type i1, size_type i2)
	{
		return p_iterator(i1, i2, 0);
	}
	const_iterator begin(size_type i1, size_type i2) const
	{
		return p_iterator(i1, i2, 0);
	}
	iterator begin(size_type i1, size_type i2, size_type i3)
	{
		return p_iterator(i1, i2, i3, 0);
	}
	const_iterator begin(size_type i1, size_type i2, size_type i3) const
	{
		return p_iterator(i1, i2, i3, 0);
	}
	iterator begin(size_type i1, size_type i2, size_type i3, size_type i4)
	{
		return p_iterator(i1, i2, i3, i4, 0);
	}
	const_iterator begin(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		return p_iterator(i1, i2, i3, i4, 0);
	}
	iterator begin(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5)
	{
		return p_iterator(i1, i2, i3, i4, i5, 0);
	}
	const_iterator begin(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		return p_iterator(i1, i2, i3, i4, i5, 0);
	}

	iterator end(size_type i1)
	{
		if( lgBC )
			return p_iterator(i1, npos);
		else
			return p_iterator(i1, p_g.v.d[i1].n);
	}
	const_iterator end(size_type i1) const
	{
		if( lgBC )
			return p_iterator(i1, npos);
		else
			return p_iterator(i1, p_g.v.d[i1].n);
	}
	iterator end(size_type i1, size_type i2)
	{
		if( lgBC )
			return p_iterator(i1, i2, npos);
		else
			return p_iterator(i1, i2, p_g.v.d[i1].d[i2].n);
	}
	const_iterator end(size_type i1, size_type i2) const
	{
		if( lgBC )
			return p_iterator(i1, i2, npos);
		else
			return p_iterator(i1, i2, p_g.v.d[i1].d[i2].n);
	}
	iterator end(size_type i1, size_type i2, size_type i3)
	{
		if( lgBC )
			return p_iterator(i1, i2, i3, npos);
		else
			return p_iterator(i1, i2, i3, p_g.v.d[i1].d[i2].d[i3].n);
	}
	const_iterator end(size_type i1, size_type i2, size_type i3) const
	{
		if( lgBC )
			return p_iterator(i1, i2, i3, npos);
		else
			return p_iterator(i1, i2, i3, p_g.v.d[i1].d[i2].d[i3].n);
	}
	iterator end(size_type i1, size_type i2, size_type i3, size_type i4)
	{
		if( lgBC )
			return p_iterator(i1, i2, i3, i4, npos);
		else
			return p_iterator(i1, i2, i3, i4, p_g.v.d[i1].d[i2].d[i3].d[i4].n);
	}
	const_iterator end(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		if( lgBC )
			return p_iterator(i1, i2, i3, i4, npos);
		else
			return p_iterator(i1, i2, i3, i4, p_g.v.d[i1].d[i2].d[i3].d[i4].n);
	}
	iterator end(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5)
	{
		if( lgBC )
			return p_iterator(i1, i2, i3, i4, i5, npos);
		else
			return p_iterator(i1, i2, i3, i4, i5, p_g.v.d[i1].d[i2].d[i3].d[i4].d[i5].n);
	}
	const_iterator end(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		if( lgBC )
			return p_iterator(i1, i2, i3, i4, i5, npos);
		else
			return p_iterator(i1, i2, i3, i4, i5, p_g.v.d[i1].d[i2].d[i3].d[i4].d[i5].n);
	}

	reference front(size_type i1)
	{
		return *begin(i1);
	}
	const_reference front(size_type i1) const
	{
		return *begin(i1);
	}
	reference front(size_type i1, size_type i2)
	{
		return *begin(i1, i2);
	}
	const_reference front(size_type i1, size_type i2) const
	{
		return *begin(i1, i2);
	}
	reference front(size_type i1, size_type i2, size_type i3)
	{
		return *begin(i1, i2, i3);
	}
	const_reference front(size_type i1, size_type i2, size_type i3) const
	{
		return *begin(i1, i2, i3);
	}
	reference front(size_type i1, size_type i2, size_type i3, size_type i4)
	{
		return *begin(i1, i2, i3, i4);
	}
	const_reference front(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		return *begin(i1, i2, i3, i4);
	}
	reference front(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5)
	{
		return *begin(i1, i2, i3, i4, i5);
	}
	const_reference front(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		return *begin(i1, i2, i3, i4, i5);
	}

	reference back(size_type i1)
	{
		return *(end(i1) - 1);
	}
	const_reference back(size_type i1) const
	{
		return *(end(i1) - 1);
	}
	reference back(size_type i1, size_type i2)
	{
		return *(end(i1, i2) - 1);
	}
	const_reference back(size_type i1, size_type i2) const
	{
		return *(end(i1, i2) - 1);
	}
	reference back(size_type i1, size_type i2, size_type i3)
	{
		return *(end(i1, i2, i3) - 1);
	}
	const_reference back(size_type i1, size_type i2, size_type i3) const
	{
		return *(end(i1, i2, i3) - 1);
	}
	reference back(size_type i1, size_type i2, size_type i3, size_type i4)
	{
		return *(end(i1, i2, i3, i4) - 1);
	}
	const_reference back(size_type i1, size_type i2, size_type i3, size_type i4) const
	{
		return *(end(i1, i2, i3, i4) - 1);
	}
	reference back(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5)
	{
		return *(end(i1, i2, i3, i4, i5) - 1);
	}
	const_reference back(size_type i1, size_type i2, size_type i3, size_type i4, size_type i5) const
	{
		return *(end(i1, i2, i3, i4, i5) - 1);
	}

	size_type size() const
	{
		return p_g.size;
	}
	size_type capacity() const
	{
		return p_g.size;
	}
	bool empty() const
	{
		for( int i=0; i < d-1; ++i )
			if( p_psl[i] != NULL )
				return false;
		return ( p_g.size == 0UL && p_dsl.size() == 0 );
	}

	pointer data()
	{
		if( p_g.size > 0 )
			return get_ptr( p_dsl );
		else
			return NULL;
	}
	const_pointer data() const
	{
		if( p_g.size > 0 )
			return get_ptr( p_dsl );
		else
			return NULL;
	}

	const multi_geom<d,ALLOC>& clone() const
	{
		return p_g;
	}

	valarray<T>& vals()
	{
		return p_dsl;
	}
	const valarray<T>& vals() const
	{
		return p_dsl;
	}
};

// predefine commonly used iterators
typedef multi_arr<bool,2>::iterator mb2i;
typedef multi_arr<bool,2>::const_iterator mb2ci;
typedef multi_arr<bool,3>::iterator mb3i;
typedef multi_arr<bool,3>::const_iterator mb3ci;
typedef multi_arr<bool,4>::iterator mb4i;
typedef multi_arr<bool,4>::const_iterator mb4ci;
typedef multi_arr<bool,5>::iterator mb5i;
typedef multi_arr<bool,5>::const_iterator mb5ci;
typedef multi_arr<bool,6>::iterator mb6i;
typedef multi_arr<bool,6>::const_iterator mb6ci;

typedef multi_arr<long,2>::iterator ml2i;
typedef multi_arr<long,2>::const_iterator ml2ci;
typedef multi_arr<long,3>::iterator ml3i;
typedef multi_arr<long,3>::const_iterator ml3ci;
typedef multi_arr<long,4>::iterator ml4i;
typedef multi_arr<long,4>::const_iterator ml4ci;
typedef multi_arr<long,5>::iterator ml5i;
typedef multi_arr<long,5>::const_iterator ml5ci;
typedef multi_arr<long,6>::iterator ml6i;
typedef multi_arr<long,6>::const_iterator ml6ci;

typedef multi_arr<realnum,2>::iterator mr2i;
typedef multi_arr<realnum,2>::const_iterator mr2ci;
typedef multi_arr<realnum,3>::iterator mr3i;
typedef multi_arr<realnum,3>::const_iterator mr3ci;
typedef multi_arr<realnum,4>::iterator mr4i;
typedef multi_arr<realnum,4>::const_iterator mr4ci;
typedef multi_arr<realnum,5>::iterator mr5i;
typedef multi_arr<realnum,5>::const_iterator mr5ci;
typedef multi_arr<realnum,6>::iterator mr6i;
typedef multi_arr<realnum,6>::const_iterator mr6ci;

typedef multi_arr<double,2>::iterator md2i;
typedef multi_arr<double,2>::const_iterator md2ci;
typedef multi_arr<double,3>::iterator md3i;
typedef multi_arr<double,3>::const_iterator md3ci;
typedef multi_arr<double,4>::iterator md4i;
typedef multi_arr<double,4>::const_iterator md4ci;
typedef multi_arr<double,5>::iterator md5i;
typedef multi_arr<double,5>::const_iterator md5ci;
typedef multi_arr<double,6>::iterator md6i;
typedef multi_arr<double,6>::const_iterator md6ci;

// on Mac systems these instantiations need to be extern in order to avoid duplicate symbols
#define INSTANTIATE_MULTI_ARR( TYPE, BC ) \
template class pntr<TYPE,BC>; \
template class const_pntr<TYPE,BC>;

template<class T, bool lgBC=lgBOUNDSCHECKVAL>
class flex_arr
{
	size_t p_size;  // number of elements allocated
	long p_begin;   // first valid array index
	long p_end;     // one beyond last valid array index
	bool p_init;    // set true when alloc() has been called

	T* p_ptr_alloc; // pointer to start of allocated data
	T* p_ptr;       // pointer used for calculating array indices

public:
	typedef random_access_iterator_tag iterator_category;
	typedef T            value_type;
	typedef T&           reference;
	typedef const T&     const_reference;
	typedef T*           pointer;
	typedef const T*     const_pointer;
	typedef long	     size_type;
	typedef ptrdiff_t    difference_type;
	typedef pntr<T,lgBC>       iterator;
	typedef const_pntr<T,lgBC> const_iterator;

private:
	void p_clear0()
	{
		p_free(p_ptr_alloc);
		p_ptr_alloc = NULL;
	}
	void p_clear1()
	{
		p_size = 0;
		p_begin = 0;
		p_end = 0;
		p_init = false;
		p_ptr_alloc = NULL;
		p_ptr = NULL;
	}
	T* p_alloc(size_t size) const
	{
		return new T[size];
	}
	void p_free(T* p) const
	{
		delete[] p;
	}

public:
	flex_arr()
	{
		p_clear1();
	}
	flex_arr(size_type begin, size_type end)
	{
		p_clear1();
		alloc( begin, end );
	}
	flex_arr(const flex_arr& f)
	{
		p_clear1();
		*this = f;
	}
	~flex_arr()
	{
		p_clear0();
	}
	const flex_arr& operator= (const flex_arr& f)
	{
		if( &f != this )
		{
			clear();
			p_size = f.p_size;
			p_begin = f.p_begin;
			p_end = f.p_end;
			p_init = f.p_init;
			if( f.p_ptr_alloc != NULL )
			{
				p_ptr_alloc = p_alloc(p_size);
				pointer p = p_ptr_alloc;
				const_pointer fp = f.p_ptr_alloc;
				for( size_type i=0; i < p_end-p_begin; ++i )
					*p++ = *fp++;
				p_ptr = p_ptr_alloc - p_begin;
			}
		}
		return *this;
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}
	void zero()
	{
		if( p_size > 0 )
			memset( p_ptr_alloc, 0, p_size*sizeof(T) );
	}
	void invalidate()
	{
		invalidate_array( p_ptr_alloc, p_size*sizeof(T) );
	}

	// reserve memory for the array
	void reserve(size_type size)
	{
		// make sure we start with a clean slate...
		clear();
		if( size > 0 )
		{
			ASSERT( p_ptr_alloc == NULL );
			p_ptr_alloc = p_alloc(size);
			p_size = (size_t)size;
		}
	}
	// allocate array with index between begin <= ind < end
	// memory is allocated here, if not already done with reserve() before
	void alloc(size_type begin, size_type end)
	{
		if( (size_t)max(end-begin,0) > p_size )
		{
			clear();

			ASSERT( p_ptr_alloc == NULL );
			p_ptr_alloc = p_alloc(end-begin);
			p_ptr = p_ptr_alloc - begin;
			p_size = (size_t)(end-begin);
		}
		else
		{
			// store was already allocated with reserve()
			p_ptr = p_ptr_alloc - begin;
		}
		p_begin = begin;
		p_end = end;
		p_init = true;
	}
	// adjust upper limit of array, reallocate store if necessary
	void realloc(size_type end)
	{
		ASSERT( p_init );
		if( (size_t)max(end-p_begin,0) > p_size )
		{
			// reallocate the store
			T* nptr_alloc = p_alloc(end-p_begin);
			T* nptr = nptr_alloc - p_begin;
			// copy store over using operator= from T, using memcpy would be a bug!
			// this could trip valgrind / purify since we don't know if this is initialized
			// there is nothing safe we can do here, so the caller should take care of this
			// note that we ignore fields above p_end, we assume nothing of interest is there
			if( p_ptr_alloc != NULL && p_ptr != NULL )
			{
				for( size_type i=p_begin; i < p_end; ++i )
					nptr[i] = p_ptr[i];
				p_free(p_ptr_alloc);
			}
			p_ptr_alloc = nptr_alloc;
			p_ptr = nptr;
			p_size = (size_t)(end-p_begin);
		}
		p_end = end;
	}

private:
	// the p_pointer() method below defines indexing into the flex_arr
	pointer p_pointer(size_type i) const
	{
		return p_ptr+i;
	}

	iterator p_iterator(size_type i) const
	{
		if( lgBC )
			return iterator( p_pointer(i), p_pointer(p_begin), p_pointer(p_end) );
		else
			return iterator( p_pointer(i) );
	}

	bool p_lgInbounds(size_type i) const
	{
		return ( i >= p_begin && i < p_end );
	}

	reference p_index(size_type i) const
	{
		if( lgBC )
		{
			if( ! p_lgInbounds( i ) )
				OUT_OF_RANGE( "flex_arr::p_index()" );
		}
		return *p_pointer(i);
	}

public:
	reference operator[] (size_type i)
	{
		return reference(p_index(i));
	}
	const_reference operator[] (size_type i) const
	{
		return const_reference(p_index(i));
	}

	reference at(size_type i)
	{
		if( ! p_lgInbounds(i) )
			OUT_OF_RANGE( "flex_arr::at()" );
		return (*this)[i];
	}
	const_reference at(size_type i) const
	{
		if( ! p_lgInbounds(i) )
			OUT_OF_RANGE( "flex_arr::at()" );
		return (*this)[i];
	}

	iterator ptr(size_type i)
	{
		return iterator(p_iterator(i));
	}
	const_iterator ptr(size_type i) const
	{
		return const_iterator(p_iterator(i));
	}

	// \todo add: assign, swap?, ... (go over stl_vector.h)

	iterator begin()
	{
		return ptr(p_begin);
	}
	const_iterator begin() const
	{
		return ptr(p_begin);
	}

	iterator end()
	{
		return ptr(p_end);
	}
	const_iterator end() const
	{
		return ptr(p_end);
	}

	reference front()
	{
		return *begin();
	}
	const_reference front() const
	{
		return *begin();
	}

	reference back()
	{
		return *(end()-1);
	}
	const_reference back() const
	{
		return *(end()-1);
	}

	size_type size() const
	{
		return max(p_end-p_begin,0);
	}
	size_type capacity() const
	{
		return p_size;
	}
	bool empty() const
	{
		return ( size() == 0 );
	}

	pointer data()
	{
		return p_ptr_alloc;
	}
	const_pointer data() const
	{
		return p_ptr_alloc;
	}

	pointer ptr0()
	{
		return p_ptr;
	}
	const_pointer ptr0() const
	{
		return p_ptr;
	}
};

// predefine commonly used iterators
typedef flex_arr<bool>::iterator fabi;
typedef flex_arr<bool>::const_iterator fabci;
typedef flex_arr<long>::iterator fali;
typedef flex_arr<long>::const_iterator falci;
typedef flex_arr<realnum>::iterator fari;
typedef flex_arr<realnum>::const_iterator farci;
typedef flex_arr<double>::iterator fadi;
typedef flex_arr<double>::const_iterator fadci;

#endif /* CONTAINER_CLASSES_H_ */
