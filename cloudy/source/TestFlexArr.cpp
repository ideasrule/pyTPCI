/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "container_classes.h"

namespace {

	struct LongIntFixtureBlank
	{
		flex_arr<long,false> arr;
		LongIntFixtureBlank() {}
		~LongIntFixtureBlank() {}
	};

	template<bool lgBC>
	struct LongIntFixtureGeneric
	{
		flex_arr<long,lgBC> arr;
		LongIntFixtureGeneric()
		{
			arr.alloc(5,10);
			for( int i=5; i < 10; ++i )
				arr[i] = i;
		}
		~LongIntFixtureGeneric() {}
	};

	typedef LongIntFixtureGeneric<false> LongIntFixture;
	typedef LongIntFixtureGeneric<true> LongIntFixtureBC;

	struct LongIntFixtureFill
	{
		flex_arr<long,false> arr;
		LongIntFixtureFill()
		{
			arr.alloc(5,10);
			for( int i=5; i < 10; ++i )
				arr[i] = i;
		}
		~LongIntFixtureFill() {}
		void myfill()
		{
			long i = 0;
			flex_arr<long,false>::iterator p;
			for( p=arr.begin(); p != arr.end(); ++p )
				*p = ++i;
		}
	};

	struct RealNumFixture
	{
		flex_arr<realnum,false> arr;
		RealNumFixture()
		{
			arr.alloc(5,10);
			for( int i=5; i < 10; ++i )
				arr[i] = i;
		}
		~RealNumFixture() {}
	};

	struct DoubleFixture
	{
		flex_arr<double,false> arr;
		DoubleFixture()
		{
			arr.alloc(5,10);
			for( int i=5; i < 10; ++i )
				arr[i] = i;
		}
		~DoubleFixture() {}
	};

	struct StructWithConstructorFixture
	{
		struct a
		{
			long n;
			a() { n = 23; }
			~a() {}
		};
		flex_arr<a,false> arr;
		StructWithConstructorFixture()
		{
			arr.alloc(7,12);
		}
		~StructWithConstructorFixture() {}
	};

	struct TestAllocFixture
	{
		TestAllocFixture() {}
		~TestAllocFixture() {}
		long mytest()
		{
			flex_arr<long,true> a(2,20);

			flex_arr<long,true>::iterator p;

			a.zero();

			long res = 0;
			for (int i=2; i<20; ++i)
			{
				p = a.ptr(i);
				res += ( *p != 0 );
			}

			return res;
		}
	};

	TEST_FIXTURE(LongIntFixture,TestZero)
	{
		arr.zero();
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(0,arr[i]);
		arr.clear();
		arr.zero();
	}

	TEST_FIXTURE(LongIntFixture,TestInvalidate1)
	{
		arr.invalidate();
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(-1,arr[i]);
		arr.clear();
		arr.invalidate();
	}

	TEST_FIXTURE(RealNumFixture,TestInvalidate2)
	{
		arr.invalidate();
		for( int i=5; i < 10; ++i )
			CHECK(isnan(arr[i]));
		arr.clear();
		arr.invalidate();
	}

	TEST_FIXTURE(DoubleFixture,TestInvalidate3)
	{
		arr.invalidate();
		for( int i=5; i < 10; ++i )
			CHECK(isnan(arr[i]));
		arr.clear();
		arr.invalidate();
	}

	TEST_FIXTURE(LongIntFixtureBC,TestAccess1)
	{
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,arr[i]);
		CHECK_THROW(arr[4],out_of_range);
		CHECK_THROW(arr[10],out_of_range);
		const flex_arr<long,true>* carr = &arr;
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,(*carr)[i]);
		CHECK_THROW((*carr)[4],out_of_range);
		CHECK_THROW((*carr)[10],out_of_range);
	}

	TEST_FIXTURE(LongIntFixtureBC,TestAccess2)
	{
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,arr[i]);
		CHECK_THROW(arr[4],out_of_range);
		CHECK_THROW(arr[10],out_of_range);
		const flex_arr<long,true>* carr = &arr;
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,(*carr)[i]);
		CHECK_THROW((*carr)[4],out_of_range);
		CHECK_THROW((*carr)[10],out_of_range);
	}

	TEST_FIXTURE(LongIntFixture,TestAccess3)
	{
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,*arr.ptr(i));
		const flex_arr<long,false>* carr = &arr;
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,*carr->ptr(i));
	}

	TEST_FIXTURE(LongIntFixtureBC,TestAccess4)
	{
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,*arr.ptr(i));
		CHECK_THROW(*arr.ptr(4),out_of_range);
		CHECK_THROW(*arr.ptr(10),out_of_range);
		const flex_arr<long,true>* carr = &arr;
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,*carr->ptr(i));
		CHECK_THROW(*carr->ptr(4),out_of_range);
		CHECK_THROW(*carr->ptr(10),out_of_range);
	}

	TEST_FIXTURE(LongIntFixture,TestAccess5)
	{
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,arr.at(i));
		CHECK_THROW(arr.at(4),out_of_range);
		CHECK_THROW(arr.at(10),out_of_range);
		const flex_arr<long,false>* carr = &arr;
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL(i,carr->at(i));
		CHECK_THROW(carr->at(4),out_of_range);
		CHECK_THROW(carr->at(10),out_of_range);
	}

	TEST_FIXTURE(LongIntFixtureBC,TestBoundaries)
	{
		CHECK( arr.begin() == arr.ptr(5) );
		CHECK( arr.end() == arr.ptr(10) );
		CHECK( &arr.front() == &arr[5] );
		CHECK( &arr.back() == &arr[9] );
		CHECK( arr.data() == &arr[5] );
		CHECK( arr.ptr0() == &arr[5]-5 );
		const flex_arr<long,true>* carr = &arr;
		CHECK( carr->begin() == arr.ptr(5) );
		CHECK( carr->end() == arr.ptr(10) );
		CHECK( &carr->front() == &arr[5] );
		CHECK( &carr->back() == &arr[9] );
		CHECK( carr->data() == &arr[5] );
		CHECK( carr->ptr0() == &arr[5]-5 );
	}

	TEST_FIXTURE(LongIntFixtureFill,TestRealloc)
	{
		myfill();
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL( i-4, arr[i] );
		long *pold = arr.data();
		arr.realloc( 20 );
		long *pnew = arr.data();
		CHECK( pold != pnew );
		for( int i=5; i < 10; ++i )
			CHECK_EQUAL( i-4, arr[i] );
	}

	TEST_FIXTURE(LongIntFixtureBlank,TestAllocationWithoutReservation)
	{
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(0L,arr.capacity());
		CHECK( arr.empty() );
		arr.alloc(5,10);
		CHECK_EQUAL(5L,arr.size());
		CHECK_EQUAL(5L,arr.capacity());
		CHECK( !arr.empty() );
		arr.realloc(30);
		CHECK_EQUAL(25L,arr.size());
		CHECK_EQUAL(25L,arr.capacity());
		CHECK( !arr.empty() );
		arr.realloc(10);
		CHECK_EQUAL(5L,arr.size());
		CHECK_EQUAL(25L,arr.capacity());
		CHECK( !arr.empty() );
		arr.realloc(0);
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(25L,arr.capacity());
		CHECK( arr.empty() );
		arr.clear();
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(0L,arr.capacity());	
		CHECK( arr.empty() );
		arr.alloc(-5,12);
		CHECK_EQUAL(17L,arr.size());
		CHECK_EQUAL(17L,arr.capacity());
		CHECK( !arr.empty() );
		arr.clear();
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(0L,arr.capacity());
		CHECK( arr.empty() );
		arr.alloc(-5,-12);
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(0L,arr.capacity());
		CHECK( arr.empty() );
		arr.realloc(30);
		CHECK_EQUAL(35L,arr.size());
		CHECK_EQUAL(35L,arr.capacity());
		CHECK( !arr.empty() );
	}

	TEST_FIXTURE(LongIntFixtureBlank,TestAllocationWithReservation)
	{
		arr.reserve(100);
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(100L,arr.capacity());
		CHECK( arr.empty() );
		arr.alloc(5,10);
		CHECK_EQUAL(5L,arr.size());
		CHECK_EQUAL(100L,arr.capacity());
		CHECK( !arr.empty() );
		arr.realloc(30);
		CHECK_EQUAL(25L,arr.size());
		CHECK_EQUAL(100L,arr.capacity());
		CHECK( !arr.empty() );
		arr.realloc(200);
		CHECK_EQUAL(195L,arr.size());
		CHECK_EQUAL(195L,arr.capacity());
		CHECK( !arr.empty() );
		arr.reserve(5);
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(5L,arr.capacity());
		CHECK( arr.empty() );
		arr.alloc(-5,-10);
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(5L,arr.capacity());
		CHECK( arr.empty() );
		arr.reserve(-5);
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(0L,arr.capacity());
		CHECK( arr.empty() );
		arr.reserve(10);
		CHECK_EQUAL(0L,arr.size());
		CHECK_EQUAL(10L,arr.capacity());
		CHECK( arr.empty() );
		arr.alloc(-10,10);
		CHECK_EQUAL(20L,arr.size());
		CHECK_EQUAL(20L,arr.capacity());
		CHECK( !arr.empty() );
	}

	TEST_FIXTURE(StructWithConstructorFixture,TestAllocationWithConstructor)
	{
		for( int i=7; i < 12; ++i )
			CHECK_EQUAL(23,arr[i].n);
	}

	// check whether the variant form for allocating works correctly
	// this also tests p_iterator in bounds-checking mode
	TEST_FIXTURE(TestAllocFixture,TestVariantAlloc)
	{
		CHECK_EQUAL(0,mytest());
	}

	TEST_FIXTURE(LongIntFixtureBC,TestCopyOperator)
	{
		flex_arr<long,true> arr2(12,234);
		CHECK( arr.size() != arr2.size() );
		arr2.zero();
		CHECK_EQUAL(0,arr2[114]);
		arr2 = arr;
		CHECK( arr.size() == arr2.size() );
		// check that copies are distinct
		CHECK( &arr[5] != &arr2[5] );
		for (int i=5; i<10; ++i)
			CHECK_EQUAL(i,arr2[i]);

		CHECK_THROW(arr2[4],out_of_range);
		CHECK_THROW(arr2[10],out_of_range);

		// is it safe to copy to oneself?
		auto copy = &arr2;
		arr2 = *copy;
		// have the contents been preserved?
		CHECK_EQUAL(9,arr2[9]);

		// now copy using constructor
		flex_arr<long,true> arr3 = arr;
		CHECK( arr.size() == arr3.size() );
		// check that copies are distinct
		CHECK( &arr[5] != &arr3[5] );
		for (int i=5; i<10; ++i)
			CHECK_EQUAL(i,arr3[i]);

		arr.clear();
		arr2 = arr;
		CHECK_EQUAL(0L,arr2.size());
		arr2.reserve( 100 );
	}

}
