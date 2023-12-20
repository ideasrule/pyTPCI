/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PROXY_ITERATOR_H_
#define PROXY_ITERATOR_H_
#include <algorithm> // For std::swap

//
// Smart pointer class extends a proxy class on a list to provide
// random_access_iterator functionality
//
// Requirements: base class P has members
//   list_type => typedef for type of base list class
//   m_list => pointer to base list class
//   m_index => offset to current list item
//   ProxyIterator<P,C> is a friend class
//
// Iterator functionality alters the value of m_index
//

template<bool, class T>
struct EnableIf_c
{
	// In general just pass through type
	typedef T type;
};
template<class T>
struct EnableIf_c<false, T>
{
	// unless argument is false, in which case map to private internal
	// type
private:
	class Invalid {};
public:
	typedef Invalid type;
};
template <class Cond, class T>
	struct EnableIf: EnableIf_c<Cond::value, T>
{
};
template <class T, class U>
struct IsSame
{
	static const bool value = false;
};
template <class T>
struct IsSame<T,T>
{
	static const bool value = true;
};
template <class T>
struct Not
{
   static const bool value = ! T::value;
};

template<class P, class C>
class ProxyIterator
{
	P proxy;
	typedef EnableIf<Not<IsSame<P,C> >,
		ProxyIterator<C,C> > const_iterator_type_1;
	typedef typename const_iterator_type_1::type const_iterator_type;
public:
	typedef P value_type;
	typedef int difference_type;
	typedef random_access_iterator_tag iterator_category;
	typedef const value_type &reference;
	typedef const value_type *pointer;
	// Construction
	explicit ProxyIterator(typename P::list_type* list, int index) : 
	proxy(list, index) {}
	ProxyIterator(const ProxyIterator &other) : 
	proxy(other.proxy.m_list, other.proxy.m_index) {}
	explicit ProxyIterator(void) : proxy(NULL, 0) {}
	ProxyIterator &operator=(ProxyIterator other) // Temporary is intentional
	{
		swap(other);
		return *this;
	}
	operator const_iterator_type() const 
	{ 
		return const_iterator_type(proxy.m_list, proxy.m_index); 
	}
	// swap member function (to implement canonical exception-safe assignment)
	void swap(ProxyIterator other)
	{
		std::swap(proxy.m_list,other.proxy.m_list);
		std::swap(proxy.m_index,other.proxy.m_index);
	}
	// Associated test (not required)
	bool associated() const
	{
		return proxy.m_list != NULL;
	}
	// Equality
	bool equals(const ProxyIterator &other) const
	{
		return other.proxy.m_list == proxy.m_list && 
			other.proxy.m_index == proxy.m_index;
	}
	// Dereference
	reference operator*() const
	{
		return proxy;
	}
	pointer operator->() const
	{
		// Return a pointer to the contained proxy, so dereferencing
		// this gives proxy not iterator behaviour.
		return &proxy;
	}
	// Incrementation
	ProxyIterator& operator++()
	{
		++proxy.m_index;
		return *this;
	}
	ProxyIterator operator++(int)
	{
		ProxyIterator old(*this);
		++*this;
		return old;
	}
	ProxyIterator& operator--()
	{
		--proxy.m_index;
		return *this;
	}
	ProxyIterator operator--(int)
	{
		ProxyIterator old(*this);
		--*this;
		return old;
	}
	// Member functions to implement standard out-of-line operators
	// Arithmetic
	const ProxyIterator add(difference_type i) const
	{
		return ProxyIterator(proxy.m_list,proxy.m_index+i);
	}
	// Number of steps between (not required)
 	difference_type diff(const ProxyIterator &other) const
	{
		return proxy.m_index - other.proxy.m_index;
	}
	// Comparison
	int cmp(const ProxyIterator &other) const
	{
		if (proxy.m_index == other.proxy.m_index)
			return 0;
		else if (proxy.m_index > other.proxy.m_index)
			return 1;
		else
			return -1;
	}
	// Compound assignment
	ProxyIterator& operator+=(difference_type i)
	{
		proxy.m_index += i;
		return *this;
	}
	ProxyIterator &operator-=(difference_type i)
	{
		proxy.m_index -= i;
		return *this;
	}
	// Offset dereference
	const value_type operator[](difference_type i) const
	{
		return P(proxy.m_list,proxy.m_index+i);
	}
};
// Identity
template<class P1, class P2, class C>
	inline bool operator==(const ProxyIterator<P1,C> &a, 
								  const ProxyIterator<P2,C> &b)
{
	return ProxyIterator<C,C>(a).equals(ProxyIterator<C,C>(b));
}
template<class P1,class P2, class C>
inline bool operator!=(const ProxyIterator<P1,C> &a, 
							  const ProxyIterator<P2,C> &b)
{
	return !(a == b);
}
// Arithmetic
template<class P, class C>
inline const ProxyIterator<P,C> operator+(
	typename ProxyIterator<P,C>::difference_type i, const ProxyIterator<P,C> &a)
{
	return a.add(i);
}
template<class P, class C>
inline const ProxyIterator<P,C> operator+(
	const ProxyIterator<P,C> &a, typename ProxyIterator<P,C>::difference_type i)
{
	return a.add(i);
}
template<class P, class C>
inline const ProxyIterator<P,C> operator-(
	const ProxyIterator<P,C> &a, typename ProxyIterator<P,C>::difference_type i)
{
	return a.add(-i);
}
// Comparison
template<class P, class C>
inline bool operator>(const ProxyIterator<P,C> &a, 
							 const ProxyIterator<P,C> &b)
{
	return a.cmp(b) > 0;
}
template<class P, class C>
inline bool operator>=(const ProxyIterator<P,C> &a, 
							  const ProxyIterator<P,C> &b)
{
	return a.cmp(b) >= 0;
}
template<class P, class C>
inline bool operator<(const ProxyIterator<P,C> &a, 
							 const ProxyIterator<P,C> &b)
{
	return a.cmp(b) < 0;
}
template<class P, class C>
inline bool operator<=(const ProxyIterator<P,C> &a, 
							  const ProxyIterator<P,C> &b)
{
	return a.cmp(b) <= 0;
}
// Iterator difference (not required)
template<class P, class C>
inline typename ProxyIterator<P,C>::difference_type operator-(
	const ProxyIterator<P,C> &a, const ProxyIterator<P,C> &b)
{
	return a.diff(b);
}

#endif
