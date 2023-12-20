/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "container_classes.h"
#include <functional>

namespace {

	// use this for out-of-bounds array access tests to avoid comiler warnings
	long m1 = -1L;

	template<mem_layout ALLOC,bool lgBC>
	struct LongInt2DFixtureGeneric
	{
		multi_arr<long,2,ALLOC,lgBC> arr;
		LongInt2DFixtureGeneric()
		{
			arr.alloc(10,9);
			for (int i=0; i<10; ++i)
				for (int j=0; j<9; ++j)
					arr[i][j] = i*10+j;
		}
		~LongInt2DFixtureGeneric() {}
	};

	typedef LongInt2DFixtureGeneric<ARPA_TYPE,false> LongInt2DFixture;
	typedef LongInt2DFixtureGeneric<ARPA_TYPE,true> LongInt2DFixtureBC;
	typedef LongInt2DFixtureGeneric<C_TYPE,false> LongInt2DFixtureCType;
	typedef LongInt2DFixtureGeneric<C_TYPE,true> LongInt2DFixtureCTypeBC;

	template<mem_layout ALLOC,bool lgBC>
	struct LongInt3DFixtureGeneric
	{
		multi_arr<long,3,ALLOC,lgBC> arr;
		LongInt3DFixtureGeneric()
		{
			arr.alloc(10,9,8);
			for (int i=0; i<10; ++i)
				for (int j=0; j<9; ++j)
					for (int k=0; k<8; ++k)
						arr[i][j][k] = i*100+j*10+k;
		}
		~LongInt3DFixtureGeneric() {}
	};

	typedef LongInt3DFixtureGeneric<ARPA_TYPE,false> LongInt3DFixture;
	typedef LongInt3DFixtureGeneric<ARPA_TYPE,true> LongInt3DFixtureBC;
	typedef LongInt3DFixtureGeneric<C_TYPE,false> LongInt3DFixtureCType;
	typedef LongInt3DFixtureGeneric<C_TYPE,true> LongInt3DFixtureCTypeBC;

	struct RealNum3DFixture
	{
		multi_arr<realnum,3,ARPA_TYPE,false> arr;
		RealNum3DFixture()
		{
			arr.alloc(10,9,8);
		}
		~RealNum3DFixture() {}
	};

	struct Double3DFixture
	{
		multi_arr<double,3,ARPA_TYPE,false> arr;
		Double3DFixture()
		{
			arr.alloc(10,9,8);
		}
		~Double3DFixture() {}
	};

	template<mem_layout ALLOC>
	struct StructWithConstructor3DFixtureGeneric
	{
		struct a
		{
			long n;
			a() { n = 23; }
			~a() {}
		};
		multi_arr<a,3,ALLOC,false> arr;
		StructWithConstructor3DFixtureGeneric()
		{
			arr.alloc(2,3,4);
		}
		~StructWithConstructor3DFixtureGeneric() {}
	};

	typedef StructWithConstructor3DFixtureGeneric<ARPA_TYPE> StructWithConstructor3DFixture;
	typedef StructWithConstructor3DFixtureGeneric<C_TYPE> StructWithConstructor3DFixtureCType;

	template<mem_layout ALLOC,bool lgBC>
	struct LongInt4DFixtureGeneric
	{
		multi_arr<long,4,ALLOC,lgBC> arr;
		LongInt4DFixtureGeneric()
		{
			arr.alloc(10,9,8,7);
			for (int i=0; i<10; ++i)
				for (int j=0; j<9; ++j)
					for (int k=0; k<8; ++k)
						for (int l=0; l<7; ++l)
							arr[i][j][k][l] = i*1000+j*100+k*10+l;
		}
		~LongInt4DFixtureGeneric() {}
	};

	typedef LongInt4DFixtureGeneric<ARPA_TYPE,false> LongInt4DFixture;
	typedef LongInt4DFixtureGeneric<ARPA_TYPE,true> LongInt4DFixtureBC;
	typedef LongInt4DFixtureGeneric<C_TYPE,false> LongInt4DFixtureCType;
	typedef LongInt4DFixtureGeneric<C_TYPE,true> LongInt4DFixtureCTypeBC;

	template<mem_layout ALLOC,bool lgBC>
	struct LongInt5DFixtureGeneric
	{
		multi_arr<long,5,ALLOC,lgBC> arr;
		LongInt5DFixtureGeneric()
		{
			arr.alloc(10,9,8,7,6);
			for (int i=0; i<10; ++i)
				for (int j=0; j<9; ++j)
					for (int k=0; k<8; ++k)
						for (int l=0; l<7; ++l)
							for (int m=0; m<6; ++m)
								arr[i][j][k][l][m] = i*10000+j*1000+k*100+l*10+m;
		}
		~LongInt5DFixtureGeneric() {}
	};

	typedef LongInt5DFixtureGeneric<ARPA_TYPE,false> LongInt5DFixture;
	typedef LongInt5DFixtureGeneric<ARPA_TYPE,true> LongInt5DFixtureBC;
	typedef LongInt5DFixtureGeneric<C_TYPE,false> LongInt5DFixtureCType;
	typedef LongInt5DFixtureGeneric<C_TYPE,true> LongInt5DFixtureCTypeBC;

	template<mem_layout ALLOC,bool lgBC>
	struct LongInt6DFixtureGeneric
	{
		multi_arr<long,6,ALLOC,lgBC> arr;
		LongInt6DFixtureGeneric()
		{
			arr.alloc(10,9,8,7,6,5);
			for (int i=0; i<10; ++i)
				for (int j=0; j<9; ++j)
					for (int k=0; k<8; ++k)
						for (int l=0; l<7; ++l)
							for (int m=0; m<6; ++m)
								for (int n=0; n<5; ++n)
									arr[i][j][k][l][m][n] =
										i*100000+j*10000+k*1000+l*100+m*10+n;
		}
		~LongInt6DFixtureGeneric() {}
	};

	typedef LongInt6DFixtureGeneric<ARPA_TYPE,false> LongInt6DFixture;
	typedef LongInt6DFixtureGeneric<ARPA_TYPE,true> LongInt6DFixtureBC;
	typedef LongInt6DFixtureGeneric<C_TYPE,false> LongInt6DFixtureCType;
	typedef LongInt6DFixtureGeneric<C_TYPE,true> LongInt6DFixtureCTypeBC;

	struct LongInt6DFixtureExplicitReserve
	{
		multi_arr<long,6,ARPA_TYPE,true> arr;
		LongInt6DFixtureExplicitReserve()
		{
			arr.reserve(10);
			for (int i=0; i<10; ++i)
			{
				arr.reserve(i,9);
				for (int j=0; j<9; ++j)
				{
					arr.reserve(i,j,8);
					for (int k=0; k<8; ++k)
					{
						arr.reserve(i,j,k,7);
						for (int l=0; l<7; ++l)
						{
							arr.reserve(i,j,k,l,6);
							for (int m=0; m<6; ++m)
							{
								arr.reserve(i,j,k,l,m,5);
							}
						}
					}
				}
			}
			arr.alloc();
		}
		~LongInt6DFixtureExplicitReserve() {}
	};

	struct TestAllocFixture
	{
		TestAllocFixture() {}
		~TestAllocFixture() {}
		long mytest()
		{
			multi_arr<long,2,ARPA_TYPE,true> a2(2,3);
			multi_arr<long,3,ARPA_TYPE,true> a3(2,3,4);
			multi_arr<long,4,ARPA_TYPE,true> a4(2,3,4,5);
			multi_arr<long,5,ARPA_TYPE,true> a5(2,3,4,5,6);
			multi_arr<long,6,ARPA_TYPE,true> a6(2,3,4,5,6,7);

			multi_arr<long,2,ARPA_TYPE,true>::iterator p2;
			multi_arr<long,3,ARPA_TYPE,true>::iterator p3;
			multi_arr<long,4,ARPA_TYPE,true>::iterator p4;
			multi_arr<long,5,ARPA_TYPE,true>::iterator p5;
			multi_arr<long,6,ARPA_TYPE,true>::iterator p6;

			a2.zero();
			a3.zero();
			a4.zero();
			a5.zero();
			a6.zero();

			long res = 0;
			for (int i=0; i<2; ++i)
				for (int j=0; j<3; ++j)
				{
					p2 = a2.ptr(i,j);
					res += ( *p2 != 0 );
					for (int k=0; k<4; ++k)
					{
						p3 = a3.ptr(i,j,k);
						res += ( *p3 != 0 );
						for (int l=0; l<5; ++l)
						{
							p4 = a4.ptr(i,j,k,l);
							res += ( *p4 != 0 );
							for (int m=0; m<6; ++m)
							{
								p5 = a5.ptr(i,j,k,l,m);
								res += ( *p5 != 0 );
								for (int n=0; n<7; ++n)
								{
									p6 = a6.ptr(i,j,k,l,m,n);
									res += ( *p6 != 0 );
								}
							}
						}
					}
				}

			return res;
		}
	};

	struct LongInt3DCLayoutFixture
	{
		multi_arr<long,3,C_TYPE,false> arr;
		LongInt3DCLayoutFixture()
		{
			arr.alloc(10,10,10);
			for (int i=0; i<10; ++i)
				for (int j=0; j<10; ++j)
					for (int k=0; k<10; ++k)
						arr[i][j][k] = 100*i+10*j+k;
		}
		~LongInt3DCLayoutFixture() {}
		long mytest(const long a[][10][10])
		{
			long res = 0;
			for (int i=0; i<10; ++i)
				for (int j=0; j<10; ++j)
					for (int k=0; k<10; ++k)
						res += ( a[i][j][k] != 100*i+10*j+k );

			return res;
		}
	};

	template<mem_layout ALLOC>
	struct LongInt3DCloneFixtureGeneric
	{
		multi_arr<long,3,ALLOC,true> arr;
		LongInt3DCloneFixtureGeneric()
		{
			arr.reserve(10);
			for (int i=0; i<10; ++i)
			{
				arr.reserve(i,i+1);
				for (int j=0; j<i+1; ++j)
					arr.reserve(i,j,j+1);
			}
			arr.alloc();
		}
		~LongInt3DCloneFixtureGeneric() {}
	};

	typedef LongInt3DCloneFixtureGeneric<ARPA_TYPE> LongInt3DCloneFixture;
	typedef LongInt3DCloneFixtureGeneric<C_TYPE> LongInt3DCloneFixtureCType;

	template<mem_layout ALLOC>
	struct LongInt2DEmptyDimGeneric
	{
		multi_arr<long,2,ALLOC,true> arr;
		LongInt2DEmptyDimGeneric()
		{
			arr.reserve(2);
			for (int i=1; i<2; ++i)
			{
				arr.reserve(i,2);
			}
			arr.alloc();
		}
		~LongInt2DEmptyDimGeneric() {}
	};

	typedef LongInt2DEmptyDimGeneric<ARPA_TYPE> LongInt2DEmptyDim;
	typedef LongInt2DEmptyDimGeneric<C_TYPE> LongInt2DEmptyDimCType;

	template<mem_layout ALLOC>
	struct LongInt3DEmptyDimGeneric
	{
		multi_arr<long,3,ALLOC,true> arr;
		LongInt3DEmptyDimGeneric()
		{
			arr.reserve(2);
			for (int i=0; i<2; ++i)
			{
				arr.reserve(i,2);
				for (int j=1; j<2; ++j)
				{
					arr.reserve(i,j,2);
				}
			}
			arr.alloc();
		}
		~LongInt3DEmptyDimGeneric() {}
	};

	typedef LongInt3DEmptyDimGeneric<ARPA_TYPE> LongInt3DEmptyDim;
	typedef LongInt3DEmptyDimGeneric<C_TYPE> LongInt3DEmptyDimCType;

	template<mem_layout ALLOC>
	struct LongInt4DEmptyDimGeneric
	{
		multi_arr<long,4,ALLOC,true> arr;
		LongInt4DEmptyDimGeneric()
		{
			arr.reserve(2);
			for (int i=0; i<2; ++i)
			{
				arr.reserve(i,2);
				for (int j=0; j<2; ++j)
				{
					arr.reserve(i,j,2);
					for (int k=1; k<2; ++k)
					{
						arr.reserve(i,j,k,2);
					}
				}
			}
			arr.alloc();
		}
		~LongInt4DEmptyDimGeneric() {}
	};

	typedef LongInt4DEmptyDimGeneric<ARPA_TYPE> LongInt4DEmptyDim;
	typedef LongInt4DEmptyDimGeneric<C_TYPE> LongInt4DEmptyDimCType;

	template<mem_layout ALLOC>
	struct LongInt5DEmptyDimGeneric
	{
		multi_arr<long,5,ALLOC,true> arr;
		LongInt5DEmptyDimGeneric()
		{
			arr.reserve(2);
			for (int i=0; i<2; ++i)
			{
				arr.reserve(i,2);
				for (int j=0; j<2; ++j)
				{
					arr.reserve(i,j,2);
					for (int k=0; k<2; ++k)
					{
						arr.reserve(i,j,k,2);
						for (int l=1; l<2; ++l)
						{
							arr.reserve(i,j,k,l,2);
						}
					}
				}
			}
			arr.alloc();
		}
		~LongInt5DEmptyDimGeneric() {}
	};

	typedef LongInt5DEmptyDimGeneric<ARPA_TYPE> LongInt5DEmptyDim;
	typedef LongInt5DEmptyDimGeneric<C_TYPE> LongInt5DEmptyDimCType;

	template<mem_layout ALLOC>
	struct LongInt6DEmptyDimGeneric
	{
		multi_arr<long,6,ALLOC,true> arr;
		LongInt6DEmptyDimGeneric()
		{
			arr.reserve(2);
			for (int i=0; i<2; ++i)
			{
				arr.reserve(i,2);
				for (int j=0; j<2; ++j)
				{
					arr.reserve(i,j,2);
					for (int k=0; k<2; ++k)
					{
						arr.reserve(i,j,k,2);
						for (int l=0; l<2; ++l)
						{
							arr.reserve(i,j,k,l,2);
							for (int m=1; m<2; ++m)
							{
								arr.reserve(i,j,k,l,m,2);
							}
						}
					}
				}
			}
			arr.alloc();
		}
		~LongInt6DEmptyDimGeneric() {}
	};

	typedef LongInt6DEmptyDimGeneric<ARPA_TYPE> LongInt6DEmptyDim;
	typedef LongInt6DEmptyDimGeneric<C_TYPE> LongInt6DEmptyDimCType;

	// first test all the forms of iterator arithmetic on the pntr class
	TEST_FIXTURE(LongInt2DFixture,TestIteratorPostIncrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p;
		p = arr.ptr(5,6);
		p++;
		CHECK_EQUAL(57,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorPreIncrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p;
		p = arr.ptr(5,6);
		++p;
		CHECK_EQUAL(57,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorPostDecrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p;
		p = arr.ptr(5,6);
		p--;
		CHECK_EQUAL(55,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorPreDecrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p;
		p = arr.ptr(5,6);
		--p;
		CHECK_EQUAL(55,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorAddition1)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p;
		p = arr.ptr(5,6);
		p += 2;
		CHECK_EQUAL(58,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorSubtraction1)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p;
		p = arr.ptr(5,6);
		p -= 2;
		CHECK_EQUAL(54,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorAddition2)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = p + 2;
		CHECK_EQUAL(58,*q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorSubtraction2)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = p - 2;
		CHECK_EQUAL(54,*q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorAddition3)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = 2 + p;
		CHECK_EQUAL(58,*q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorSubtraction3)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = p - 2;
		CHECK_EQUAL(2,p-q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorBrackets)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p;
		p = arr.ptr(5,6);
		CHECK_EQUAL(57,p[1]);
	}

	// now test all 6 forms of logical comparison on the pntr class
	TEST_FIXTURE(LongInt2DFixture,TestIteratorComparison1)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = ++p;
		CHECK(p == q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorComparison2)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = p++;
		CHECK(p != q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorComparison3)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = p--;
		CHECK(p < q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorComparison4)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = --p;
		CHECK(p <= q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorComparison5)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = p++;
		CHECK(p > q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestIteratorComparison6)
	{
		multi_arr<long,2,ARPA_TYPE,false>::iterator p, q;
		p = arr.ptr(5,6);
		q = --p;
		CHECK(p >= q);
	}

	// test out-of-bounds access via a pntr
	TEST_FIXTURE(LongInt2DFixtureBC,TestIteratorBoundsCheck)
	{
		multi_arr<long,2,ARPA_TYPE,true>::iterator p;
		p = arr.ptr(5,8);
		CHECK_THROW(UNUSED long i = *++p,out_of_range);
	}

	// now do all the tests from above again on the const_pntr class
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorPostIncrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p;
		p = arr.ptr(5,6);
		p++;
		CHECK_EQUAL(57,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorPreIncrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p;
		p = arr.ptr(5,6);
		++p;
		CHECK_EQUAL(57,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorPostDecrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p;
		p = arr.ptr(5,6);
		p--;
		CHECK_EQUAL(55,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorPreDecrement)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p;
		p = arr.ptr(5,6);
		--p;
		CHECK_EQUAL(55,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorAddition1)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p;
		p = arr.ptr(5,6);
		p += 2;
		CHECK_EQUAL(58,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorSubtraction1)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p;
		p = arr.ptr(5,6);
		p -= 2;
		CHECK_EQUAL(54,*p);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorAddition2)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = p + 2;
		CHECK_EQUAL(58,*q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorSubtraction2)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = p - 2;
		CHECK_EQUAL(54,*q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorAddition3)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = 2 + p;
		CHECK_EQUAL(58,*q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorSubtraction3)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = p - 2;
		CHECK_EQUAL(2,p-q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorBrackets)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p;
		p = arr.ptr(5,6);
		CHECK_EQUAL(57,p[1]);
	}

	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorComparison1)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = ++p;
		CHECK(p == q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorComparison2)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = p++;
		CHECK(p != q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorComparison3)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = p--;
		CHECK(p < q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorComparison4)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = --p;
		CHECK(p <= q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorComparison5)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = p++;
		CHECK(p > q);
	}
	TEST_FIXTURE(LongInt2DFixture,TestConstIteratorComparison6)
	{
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p, q;
		p = arr.ptr(5,6);
		q = --p;
		CHECK(p >= q);
	}

	TEST_FIXTURE(LongInt2DFixtureBC,TestConstIteratorBoundsCheck)
	{
		multi_arr<long,2,ARPA_TYPE,true>::const_iterator p;
		p = arr.ptr(5,8);
		CHECK_THROW(UNUSED long i = *++p,out_of_range);
	}

	TEST_FIXTURE(LongInt3DFixture,TestDataEmpty)
	{
		arr.clear();
		// calling alloc() should be safe...
		arr.alloc();
		CHECK_EQUAL( 0UL, arr.size() );
		CHECK( NULL == arr.data() );
	}

	// now we test the zero(), invalidate(), and value assignment methods
	TEST_FIXTURE(LongInt3DFixture,TestFill)
	{
		arr.zero();
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK_EQUAL(0,arr[i][j][k]);
		arr.invalidate();
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK_EQUAL(-1,arr[i][j][k]);
		arr = -7;
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK_EQUAL(-7,arr[i][j][k]);
		arr.clear();
		// these operations should be safe and do nothing
		arr.zero();
		arr.invalidate();
		arr = -7;
	}
	// same for C_TYPE arrays, also checks whether any legal index
	// accidentally throws an out_of_range exception
	TEST_FIXTURE(LongInt3DFixtureCTypeBC,TestFillCType)
	{
		arr.zero();
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK_EQUAL(0.,arr[i][j][k]);
		arr.invalidate();
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK_EQUAL(-1,arr[i][j][k]);
		arr = -9;
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK_EQUAL(-9,arr[i][j][k]);
	}
	// now test whether realnum and double are properly invalidated
	TEST_FIXTURE(RealNum3DFixture,TestInvalidRealNum)
	{
		arr.invalidate();
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK(isnan(arr[i][j][k]));
		arr.clear();
		// these operations should be safe and do nothing
		arr.invalidate();
	}
	TEST_FIXTURE(Double3DFixture,TestInvalidDouble)
	{
		arr.invalidate();
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					CHECK(isnan(arr[i][j][k]));		
		arr.clear();
		// these operations should be safe and do nothing
		arr.invalidate();
	}

	// test access through both indexing and ptr methods for all dimensions
	// some of this will already have been tested implicitly above
	TEST_FIXTURE(LongInt2DFixture,Test2DIndexedValue)
	{
		CHECK_EQUAL(98,arr[9][8]);
		CHECK_EQUAL(37,*arr.ptr(3,7));
		const multi_arr<long,2,ARPA_TYPE,false>* carr = &arr;
		multi_arr<long,2,ARPA_TYPE,false>::const_iterator p = carr->ptr(4,6);
		CHECK_EQUAL(46,*p);
		const long& q = (*carr)[7][5];
		CHECK_EQUAL(75,q);
	}
	TEST_FIXTURE(LongInt3DFixture,Test3DIndexedValue)
	{
		CHECK_EQUAL(987,arr[9][8][7]);
		CHECK_EQUAL(137,*arr.ptr(1,3,7));
		const multi_arr<long,3,ARPA_TYPE,false>* carr = &arr;
		multi_arr<long,3,ARPA_TYPE,false>::const_iterator p = carr->ptr(4,6,2);
		CHECK_EQUAL(462,*p);
		const long& q = (*carr)[7][5][6];
		CHECK_EQUAL(756,q);
	}
	TEST_FIXTURE(LongInt4DFixture,Test4DIndexedValue)
	{
		CHECK_EQUAL(9876,arr[9][8][7][6]);
		CHECK_EQUAL(1342,*arr.ptr(1,3,4,2));
		const multi_arr<long,4,ARPA_TYPE,false>* carr = &arr;
		multi_arr<long,4,ARPA_TYPE,false>::const_iterator p = carr->ptr(4,6,2,0);
		CHECK_EQUAL(4620,*p);
		const long& q = (*carr)[7][5][6][3];
		CHECK_EQUAL(7563,q);
	}
	TEST_FIXTURE(LongInt5DFixture,Test5DIndexedValue)
	{
		CHECK_EQUAL(98765,arr[9][8][7][6][5]);
		CHECK_EQUAL(10342,*arr.ptr(1,0,3,4,2));
		const multi_arr<long,5,ARPA_TYPE,false>* carr = &arr;
		multi_arr<long,5,ARPA_TYPE,false>::const_iterator p = carr->ptr(4,6,2,0,1);
		CHECK_EQUAL(46201,*p);
		const long& q = (*carr)[7][5][6][3][2];
		CHECK_EQUAL(75632,q);
	}
	TEST_FIXTURE(LongInt6DFixture,Test6DIndexedValue)
	{
		CHECK_EQUAL(987654,arr[9][8][7][6][5][4]);
		CHECK_EQUAL(106423,*arr.ptr(1,0,6,4,2,3));
		const multi_arr<long,6,ARPA_TYPE,false>* carr = &arr;
		multi_arr<long,6,ARPA_TYPE,false>::const_iterator p = carr->ptr(4,6,2,0,1,3);
		CHECK_EQUAL(462013,*p);
		const long& q = (*carr)[7][5][6][3][2][1];
		CHECK_EQUAL(756321,q);
	}

	// check whether out-of-bounds access is detected
	// create two exceptions for each dimension,
	// just below and above the valid range
	TEST_FIXTURE(LongInt2DFixture,Test2DAtOutOfBounds)
	{
		CHECK_THROW(arr.at(m1,1),out_of_range);
		CHECK_THROW(arr.at(10,1),out_of_range);
		CHECK_THROW(arr.at(7,m1),out_of_range);
		CHECK_THROW(arr.at(7,9),out_of_range);
		const multi_arr<long,2,ARPA_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1),out_of_range);
		CHECK_THROW(carr->at(10,1),out_of_range);
		CHECK_THROW(carr->at(7,m1),out_of_range);
		CHECK_THROW(carr->at(7,9),out_of_range);
	}
	TEST_FIXTURE(LongInt3DFixture,Test3DAtOutOfBounds)
	{
		CHECK_THROW(arr.at(m1,1,2),out_of_range);
		CHECK_THROW(arr.at(10,1,2),out_of_range);
		CHECK_THROW(arr.at(7,m1,3),out_of_range);
		CHECK_THROW(arr.at(7,9,3),out_of_range);
		CHECK_THROW(arr.at(7,3,m1),out_of_range);
		CHECK_THROW(arr.at(7,3,8),out_of_range);
		const multi_arr<long,3,ARPA_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2),out_of_range);
		CHECK_THROW(carr->at(10,1,2),out_of_range);
		CHECK_THROW(carr->at(7,m1,3),out_of_range);
		CHECK_THROW(carr->at(7,9,3),out_of_range);
		CHECK_THROW(carr->at(7,3,m1),out_of_range);
		CHECK_THROW(carr->at(7,3,8),out_of_range);
	}
	TEST_FIXTURE(LongInt4DFixture,Test4DAtOutOfBounds)
	{
		CHECK_THROW(arr.at(m1,1,2,3),out_of_range);
		CHECK_THROW(arr.at(10,1,2,3),out_of_range);
		CHECK_THROW(arr.at(7,m1,3,2),out_of_range);
		CHECK_THROW(arr.at(7,9,3,2),out_of_range);
		CHECK_THROW(arr.at(7,3,m1,1),out_of_range);
		CHECK_THROW(arr.at(7,3,8,1),out_of_range);
		CHECK_THROW(arr.at(7,3,1,m1),out_of_range);
		CHECK_THROW(arr.at(7,3,1,7),out_of_range);
		const multi_arr<long,4,ARPA_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2,3),out_of_range);
		CHECK_THROW(carr->at(10,1,2,3),out_of_range);
		CHECK_THROW(carr->at(7,m1,3,2),out_of_range);
		CHECK_THROW(carr->at(7,9,3,2),out_of_range);
		CHECK_THROW(carr->at(7,3,m1,1),out_of_range);
		CHECK_THROW(carr->at(7,3,8,1),out_of_range);
		CHECK_THROW(carr->at(7,3,1,m1),out_of_range);
		CHECK_THROW(carr->at(7,3,1,7),out_of_range);
	}
	TEST_FIXTURE(LongInt5DFixture,Test5DAtOutOfBounds)
	{
		CHECK_THROW(arr.at(m1,1,2,3,4),out_of_range);
		CHECK_THROW(arr.at(10,1,2,3,4),out_of_range);
		CHECK_THROW(arr.at(7,m1,3,2,4),out_of_range);
		CHECK_THROW(arr.at(7,9,3,2,4),out_of_range);
		CHECK_THROW(arr.at(7,3,m1,1,0),out_of_range);
		CHECK_THROW(arr.at(7,3,8,1,0),out_of_range);
		CHECK_THROW(arr.at(7,3,1,m1,2),out_of_range);
		CHECK_THROW(arr.at(7,3,1,7,2),out_of_range);
		CHECK_THROW(arr.at(7,3,1,2,m1),out_of_range);
		CHECK_THROW(arr.at(7,3,1,2,6),out_of_range);
		const multi_arr<long,5,ARPA_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2,3,4),out_of_range);
		CHECK_THROW(carr->at(10,1,2,3,4),out_of_range);
		CHECK_THROW(carr->at(7,m1,3,2,4),out_of_range);
		CHECK_THROW(carr->at(7,9,3,2,4),out_of_range);
		CHECK_THROW(carr->at(7,3,m1,1,0),out_of_range);
		CHECK_THROW(carr->at(7,3,8,1,0),out_of_range);
		CHECK_THROW(carr->at(7,3,1,m1,2),out_of_range);
		CHECK_THROW(carr->at(7,3,1,7,2),out_of_range);
		CHECK_THROW(carr->at(7,3,1,2,m1),out_of_range);
		CHECK_THROW(carr->at(7,3,1,2,6),out_of_range);
	}
	TEST_FIXTURE(LongInt6DFixture,Test6DAtOutOfBounds)
	{
		CHECK_THROW(arr.at(m1,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.at(10,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.at(7,m1,3,2,4,1),out_of_range);
		CHECK_THROW(arr.at(7,9,3,2,4,1),out_of_range);
		CHECK_THROW(arr.at(7,3,m1,1,0,2),out_of_range);
		CHECK_THROW(arr.at(7,3,8,1,0,2),out_of_range);
		CHECK_THROW(arr.at(7,3,1,m1,2,0),out_of_range);
		CHECK_THROW(arr.at(7,3,1,7,2,0),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,m1,3),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,6,3),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,3,m1),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,3,5),out_of_range);
		const multi_arr<long,6,ARPA_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2,3,4,0),out_of_range);
		CHECK_THROW(carr->at(10,1,2,3,4,0),out_of_range);
		CHECK_THROW(carr->at(7,m1,3,2,4,1),out_of_range);
		CHECK_THROW(carr->at(7,9,3,2,4,1),out_of_range);
		CHECK_THROW(carr->at(7,3,m1,1,0,2),out_of_range);
		CHECK_THROW(carr->at(7,3,8,1,0,2),out_of_range);
		CHECK_THROW(carr->at(7,3,1,m1,2,0),out_of_range);
		CHECK_THROW(carr->at(7,3,1,7,2,0),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,m1,3),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,6,3),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,3,m1),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,3,5),out_of_range);
	}

	TEST_FIXTURE(LongInt2DFixtureBC,Test2DOutOfBounds)
	{
		CHECK_THROW(arr[m1][1],out_of_range);
		CHECK_THROW(arr[10][1],out_of_range);
		CHECK_THROW(arr[7][m1],out_of_range);
		CHECK_THROW(arr[7][9],out_of_range);
		CHECK_THROW(arr.ptr(m1,1),out_of_range);
		CHECK_THROW(arr.ptr(10,1),out_of_range);
		const multi_arr<long,2,ARPA_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1],out_of_range);
		CHECK_THROW((*carr)[10][1],out_of_range);
		CHECK_THROW((*carr)[7][m1],out_of_range);
		CHECK_THROW((*carr)[7][9],out_of_range);
	}
	TEST_FIXTURE(LongInt3DFixtureBC,Test3DOutOfBounds)
	{
		CHECK_THROW(arr[m1][1][2],out_of_range);
		CHECK_THROW(arr[10][1][2],out_of_range);
		CHECK_THROW(arr[7][m1][3],out_of_range);
		CHECK_THROW(arr[7][9][3],out_of_range);
		CHECK_THROW(arr[7][3][m1],out_of_range);
		CHECK_THROW(arr[7][3][8],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3),out_of_range);
		const multi_arr<long,3,ARPA_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2],out_of_range);
		CHECK_THROW((*carr)[10][1][2],out_of_range);
		CHECK_THROW((*carr)[7][m1][3],out_of_range);
		CHECK_THROW((*carr)[7][9][3],out_of_range);
		CHECK_THROW((*carr)[7][3][m1],out_of_range);
		CHECK_THROW((*carr)[7][3][8],out_of_range);
	}
	TEST_FIXTURE(LongInt4DFixtureBC,Test4DOutOfBounds)
	{
		CHECK_THROW(arr[m1][1][2][3],out_of_range);
		CHECK_THROW(arr[10][1][2][3],out_of_range);
		CHECK_THROW(arr[7][m1][3][2],out_of_range);
		CHECK_THROW(arr[7][9][3][2],out_of_range);
		CHECK_THROW(arr[7][3][m1][1],out_of_range);
		CHECK_THROW(arr[7][3][8][1],out_of_range);
		CHECK_THROW(arr[7][3][1][m1],out_of_range);
		CHECK_THROW(arr[7][3][1][7],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2,3),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2,3),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3,2),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,m1,1),out_of_range);
		CHECK_THROW(arr.ptr(7,3,8,1),out_of_range);
		const multi_arr<long,4,ARPA_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2][3],out_of_range);
		CHECK_THROW((*carr)[10][1][2][3],out_of_range);
		CHECK_THROW((*carr)[7][m1][3][2],out_of_range);
		CHECK_THROW((*carr)[7][9][3][2],out_of_range);
		CHECK_THROW((*carr)[7][3][m1][1],out_of_range);
		CHECK_THROW((*carr)[7][3][8][1],out_of_range);
		CHECK_THROW((*carr)[7][3][1][m1],out_of_range);
		CHECK_THROW((*carr)[7][3][1][7],out_of_range);
	}
	TEST_FIXTURE(LongInt5DFixtureBC,Test5DOutOfBounds)
	{
		CHECK_THROW(arr[m1][1][2][3][4],out_of_range);
		CHECK_THROW(arr[10][1][2][3][4],out_of_range);
		CHECK_THROW(arr[7][m1][3][2][4],out_of_range);
		CHECK_THROW(arr[7][9][3][2][4],out_of_range);
		CHECK_THROW(arr[7][3][m1][1][0],out_of_range);
		CHECK_THROW(arr[7][3][8][1][0],out_of_range);
		CHECK_THROW(arr[7][3][1][m1][2],out_of_range);
		CHECK_THROW(arr[7][3][1][7][2],out_of_range);
		CHECK_THROW(arr[7][3][1][2][m1],out_of_range);
		CHECK_THROW(arr[7][3][1][2][6],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2,3,4),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2,3,4),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3,2,4),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3,2,4),out_of_range);
		CHECK_THROW(arr.ptr(7,3,m1,1,0),out_of_range);
		CHECK_THROW(arr.ptr(7,3,8,1,0),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,m1,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,7,2),out_of_range);
		const multi_arr<long,5,ARPA_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2][3][4],out_of_range);
		CHECK_THROW((*carr)[10][1][2][3][4],out_of_range);
		CHECK_THROW((*carr)[7][m1][3][2][4],out_of_range);
		CHECK_THROW((*carr)[7][9][3][2][4],out_of_range);
		CHECK_THROW((*carr)[7][3][m1][1][0],out_of_range);
		CHECK_THROW((*carr)[7][3][8][1][0],out_of_range);
		CHECK_THROW((*carr)[7][3][1][m1][2],out_of_range);
		CHECK_THROW((*carr)[7][3][1][7][2],out_of_range);
		CHECK_THROW((*carr)[7][3][1][2][m1],out_of_range);
		CHECK_THROW((*carr)[7][3][1][2][6],out_of_range);
	}
	TEST_FIXTURE(LongInt6DFixtureBC,Test6DOutOfBounds)
	{
		CHECK_THROW(arr[m1][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr[10][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr[7][m1][3][2][4][1],out_of_range);
		CHECK_THROW(arr[7][9][3][2][4][1],out_of_range);
		CHECK_THROW(arr[7][3][m1][1][0][2],out_of_range);
		CHECK_THROW(arr[7][3][8][1][0][2],out_of_range);
		CHECK_THROW(arr[7][3][1][m1][2][0],out_of_range);
		CHECK_THROW(arr[7][3][1][7][2][0],out_of_range);
		CHECK_THROW(arr[7][4][1][2][m1][3],out_of_range);
		CHECK_THROW(arr[7][4][1][2][6][3],out_of_range);
		CHECK_THROW(arr[7][4][1][2][3][m1],out_of_range);
		CHECK_THROW(arr[7][4][1][2][3][5],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3,2,4,1),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3,2,4,1),out_of_range);
		CHECK_THROW(arr.ptr(7,3,m1,1,0,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,8,1,0,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,m1,2,0),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,7,2,0),out_of_range);
		CHECK_THROW(arr.ptr(7,4,1,2,m1,3),out_of_range);
		CHECK_THROW(arr.ptr(7,4,1,2,6,3),out_of_range);
		const multi_arr<long,6,ARPA_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2][3][4][0],out_of_range);
		CHECK_THROW((*carr)[10][1][2][3][4][0],out_of_range);
		CHECK_THROW((*carr)[7][m1][3][2][4][1],out_of_range);
		CHECK_THROW((*carr)[7][9][3][2][4][1],out_of_range);
		CHECK_THROW((*carr)[7][3][m1][1][0][2],out_of_range);
		CHECK_THROW((*carr)[7][3][8][1][0][2],out_of_range);
		CHECK_THROW((*carr)[7][3][1][m1][2][0],out_of_range);
		CHECK_THROW((*carr)[7][3][1][7][2][0],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][m1][3],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][6][3],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][3][m1],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][3][5],out_of_range);
	}

	// now repeat access and out-of-bounds tests for C_TYPE arrays
	TEST_FIXTURE(LongInt2DFixtureCType,Test2DIndexedValueCType)
	{
		CHECK_EQUAL(98,arr[9][8]);
		CHECK_EQUAL(37,*arr.ptr(3,7));
		const multi_arr<long,2,C_TYPE,false>* carr = &arr;
		multi_arr<long,2,C_TYPE,false>::const_iterator p = carr->ptr(4,6);
		CHECK_EQUAL(46,*p);
		const long& q = (*carr)[7][5];
		CHECK_EQUAL(75,q);
	}
	TEST_FIXTURE(LongInt3DFixtureCType,Test3DIndexedValueCType)
	{
		CHECK_EQUAL(987,arr[9][8][7]);
		CHECK_EQUAL(137,*arr.ptr(1,3,7));
		const multi_arr<long,3,C_TYPE,false>* carr = &arr;
		multi_arr<long,3,C_TYPE,false>::const_iterator p = carr->ptr(4,6,2);
		CHECK_EQUAL(462,*p);
		const long& q = (*carr)[7][5][6];
		CHECK_EQUAL(756,q);
	}
	TEST_FIXTURE(LongInt4DFixtureCType,Test4DIndexedValueCType)
	{
		CHECK_EQUAL(9876,arr[9][8][7][6]);
		CHECK_EQUAL(1342,*arr.ptr(1,3,4,2));
		const multi_arr<long,4,C_TYPE,false>* carr = &arr;
		multi_arr<long,4,C_TYPE,false>::const_iterator p = carr->ptr(4,6,2,0);
		CHECK_EQUAL(4620,*p);
		const long& q = (*carr)[7][5][6][3];
		CHECK_EQUAL(7563,q);
	}
	TEST_FIXTURE(LongInt5DFixtureCType,Test5DIndexedValueCType)
	{
		CHECK_EQUAL(98765,arr[9][8][7][6][5]);
		CHECK_EQUAL(10342,*arr.ptr(1,0,3,4,2));
		const multi_arr<long,5,C_TYPE,false>* carr = &arr;
		multi_arr<long,5,C_TYPE,false>::const_iterator p = carr->ptr(4,6,2,0,1);
		CHECK_EQUAL(46201,*p);
		const long& q = (*carr)[7][5][6][3][2];
		CHECK_EQUAL(75632,q);
	}
	TEST_FIXTURE(LongInt6DFixtureCType,Test6DIndexedValueCType)
	{
		CHECK_EQUAL(987654,arr[9][8][7][6][5][4]);
		CHECK_EQUAL(106423,*arr.ptr(1,0,6,4,2,3));
		const multi_arr<long,6,C_TYPE,false>* carr = &arr;
		multi_arr<long,6,C_TYPE,false>::const_iterator p = carr->ptr(4,6,2,0,1,3);
		CHECK_EQUAL(462013,*p);
		const long& q = (*carr)[7][5][6][3][2][1];
		CHECK_EQUAL(756321,q);
	}

	TEST_FIXTURE(LongInt2DFixtureCType,Test2DAtOutOfBoundsCType)
	{
		CHECK_THROW(arr.at(m1,1),out_of_range);
		CHECK_THROW(arr.at(10,1),out_of_range);
		CHECK_THROW(arr.at(7,m1),out_of_range);
		CHECK_THROW(arr.at(7,9),out_of_range);
		const multi_arr<long,2,C_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1),out_of_range);
		CHECK_THROW(carr->at(10,1),out_of_range);
		CHECK_THROW(carr->at(7,m1),out_of_range);
		CHECK_THROW(carr->at(7,9),out_of_range);
	}
	TEST_FIXTURE(LongInt3DFixtureCType,Test3DAtOutOfBoundsCType)
	{
		CHECK_THROW(arr.at(m1,1,2),out_of_range);
		CHECK_THROW(arr.at(10,1,2),out_of_range);
		CHECK_THROW(arr.at(7,m1,3),out_of_range);
		CHECK_THROW(arr.at(7,9,3),out_of_range);
		CHECK_THROW(arr.at(7,3,m1),out_of_range);
		CHECK_THROW(arr.at(7,3,8),out_of_range);
		const multi_arr<long,3,C_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2),out_of_range);
		CHECK_THROW(carr->at(10,1,2),out_of_range);
		CHECK_THROW(carr->at(7,m1,3),out_of_range);
		CHECK_THROW(carr->at(7,9,3),out_of_range);
		CHECK_THROW(carr->at(7,3,m1),out_of_range);
		CHECK_THROW(carr->at(7,3,8),out_of_range);
	}
	TEST_FIXTURE(LongInt4DFixtureCType,Test4DAtOutOfBoundsCType)
	{
		CHECK_THROW(arr.at(m1,1,2,3),out_of_range);
		CHECK_THROW(arr.at(10,1,2,3),out_of_range);
		CHECK_THROW(arr.at(7,m1,3,2),out_of_range);
		CHECK_THROW(arr.at(7,9,3,2),out_of_range);
		CHECK_THROW(arr.at(7,3,m1,1),out_of_range);
		CHECK_THROW(arr.at(7,3,8,1),out_of_range);
		CHECK_THROW(arr.at(7,3,1,m1),out_of_range);
		CHECK_THROW(arr.at(7,3,1,7),out_of_range);
		const multi_arr<long,4,C_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2,3),out_of_range);
		CHECK_THROW(carr->at(10,1,2,3),out_of_range);
		CHECK_THROW(carr->at(7,m1,3,2),out_of_range);
		CHECK_THROW(carr->at(7,9,3,2),out_of_range);
		CHECK_THROW(carr->at(7,3,m1,1),out_of_range);
		CHECK_THROW(carr->at(7,3,8,1),out_of_range);
		CHECK_THROW(carr->at(7,3,1,m1),out_of_range);
		CHECK_THROW(carr->at(7,3,1,7),out_of_range);
	}
	TEST_FIXTURE(LongInt5DFixtureCType,Test5DAtOutOfBoundsCType)
	{
		CHECK_THROW(arr.at(m1,1,2,3,4),out_of_range);
		CHECK_THROW(arr.at(10,1,2,3,4),out_of_range);
		CHECK_THROW(arr.at(7,m1,3,2,4),out_of_range);
		CHECK_THROW(arr.at(7,9,3,2,4),out_of_range);
		CHECK_THROW(arr.at(7,3,m1,1,0),out_of_range);
		CHECK_THROW(arr.at(7,3,8,1,0),out_of_range);
		CHECK_THROW(arr.at(7,3,1,m1,2),out_of_range);
		CHECK_THROW(arr.at(7,3,1,7,2),out_of_range);
		CHECK_THROW(arr.at(7,3,1,2,m1),out_of_range);
		CHECK_THROW(arr.at(7,3,1,2,6),out_of_range);
		const multi_arr<long,5,C_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2,3,4),out_of_range);
		CHECK_THROW(carr->at(10,1,2,3,4),out_of_range);
		CHECK_THROW(carr->at(7,m1,3,2,4),out_of_range);
		CHECK_THROW(carr->at(7,9,3,2,4),out_of_range);
		CHECK_THROW(carr->at(7,3,m1,1,0),out_of_range);
		CHECK_THROW(carr->at(7,3,8,1,0),out_of_range);
		CHECK_THROW(carr->at(7,3,1,m1,2),out_of_range);
		CHECK_THROW(carr->at(7,3,1,7,2),out_of_range);
		CHECK_THROW(carr->at(7,3,1,2,m1),out_of_range);
		CHECK_THROW(carr->at(7,3,1,2,6),out_of_range);
	}
	TEST_FIXTURE(LongInt6DFixtureCType,Test6DAtOutOfBoundsCType)
	{
		CHECK_THROW(arr.at(m1,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.at(10,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.at(7,m1,3,2,4,1),out_of_range);
		CHECK_THROW(arr.at(7,9,3,2,4,1),out_of_range);
		CHECK_THROW(arr.at(7,3,m1,1,0,2),out_of_range);
		CHECK_THROW(arr.at(7,3,8,1,0,2),out_of_range);
		CHECK_THROW(arr.at(7,3,1,m1,2,0),out_of_range);
		CHECK_THROW(arr.at(7,3,1,7,2,0),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,m1,3),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,6,3),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,3,m1),out_of_range);
		CHECK_THROW(arr.at(7,4,1,2,3,5),out_of_range);
		const multi_arr<long,6,C_TYPE,false>* carr = &arr;
		CHECK_THROW(carr->at(m1,1,2,3,4,0),out_of_range);
		CHECK_THROW(carr->at(10,1,2,3,4,0),out_of_range);
		CHECK_THROW(carr->at(7,m1,3,2,4,1),out_of_range);
		CHECK_THROW(carr->at(7,9,3,2,4,1),out_of_range);
		CHECK_THROW(carr->at(7,3,m1,1,0,2),out_of_range);
		CHECK_THROW(carr->at(7,3,8,1,0,2),out_of_range);
		CHECK_THROW(carr->at(7,3,1,m1,2,0),out_of_range);
		CHECK_THROW(carr->at(7,3,1,7,2,0),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,m1,3),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,6,3),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,3,m1),out_of_range);
		CHECK_THROW(carr->at(7,4,1,2,3,5),out_of_range);
	}

	TEST_FIXTURE(LongInt2DFixtureCTypeBC,Test2DOutOfBoundsCType)
	{
		CHECK_THROW(arr[m1][1],out_of_range);
		CHECK_THROW(arr[10][1],out_of_range);
		CHECK_THROW(arr[7][m1],out_of_range);
		CHECK_THROW(arr[7][9],out_of_range);
		CHECK_THROW(arr.ptr(m1,1),out_of_range);
		CHECK_THROW(arr.ptr(10,1),out_of_range);
		const multi_arr<long,2,C_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1],out_of_range);
		CHECK_THROW((*carr)[10][1],out_of_range);
		CHECK_THROW((*carr)[7][m1],out_of_range);
		CHECK_THROW((*carr)[7][9],out_of_range);
	}
	TEST_FIXTURE(LongInt3DFixtureCTypeBC,Test3DOutOfBoundsCType)
	{
		CHECK_THROW(arr[m1][1][2],out_of_range);
		CHECK_THROW(arr[10][1][2],out_of_range);
		CHECK_THROW(arr[7][m1][3],out_of_range);
		CHECK_THROW(arr[7][9][3],out_of_range);
		CHECK_THROW(arr[7][3][m1],out_of_range);
		CHECK_THROW(arr[7][3][8],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3),out_of_range);
		const multi_arr<long,3,C_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2],out_of_range);
		CHECK_THROW((*carr)[10][1][2],out_of_range);
		CHECK_THROW((*carr)[7][m1][3],out_of_range);
		CHECK_THROW((*carr)[7][9][3],out_of_range);
		CHECK_THROW((*carr)[7][3][m1],out_of_range);
		CHECK_THROW((*carr)[7][3][8],out_of_range);
	}
	TEST_FIXTURE(LongInt4DFixtureCTypeBC,Test4DOutOfBoundsCType)
	{
		CHECK_THROW(arr[m1][1][2][3],out_of_range);
		CHECK_THROW(arr[10][1][2][3],out_of_range);
		CHECK_THROW(arr[7][m1][3][2],out_of_range);
		CHECK_THROW(arr[7][9][3][2],out_of_range);
		CHECK_THROW(arr[7][3][m1][1],out_of_range);
		CHECK_THROW(arr[7][3][8][1],out_of_range);
		CHECK_THROW(arr[7][3][1][m1],out_of_range);
		CHECK_THROW(arr[7][3][1][7],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2,3),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2,3),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3,2),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,m1,1),out_of_range);
		CHECK_THROW(arr.ptr(7,3,8,1),out_of_range);
		const multi_arr<long,4,C_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2][3],out_of_range);
		CHECK_THROW((*carr)[10][1][2][3],out_of_range);
		CHECK_THROW((*carr)[7][m1][3][2],out_of_range);
		CHECK_THROW((*carr)[7][9][3][2],out_of_range);
		CHECK_THROW((*carr)[7][3][m1][1],out_of_range);
		CHECK_THROW((*carr)[7][3][8][1],out_of_range);
		CHECK_THROW((*carr)[7][3][1][m1],out_of_range);
		CHECK_THROW((*carr)[7][3][1][7],out_of_range);
	}
	TEST_FIXTURE(LongInt5DFixtureCTypeBC,Test5DOutOfBoundsCType)
	{
		CHECK_THROW(arr[m1][1][2][3][4],out_of_range);
		CHECK_THROW(arr[10][1][2][3][4],out_of_range);
		CHECK_THROW(arr[7][m1][3][2][4],out_of_range);
		CHECK_THROW(arr[7][9][3][2][4],out_of_range);
		CHECK_THROW(arr[7][3][m1][1][0],out_of_range);
		CHECK_THROW(arr[7][3][8][1][0],out_of_range);
		CHECK_THROW(arr[7][3][1][m1][2],out_of_range);
		CHECK_THROW(arr[7][3][1][7][2],out_of_range);
		CHECK_THROW(arr[7][3][1][2][m1],out_of_range);
		CHECK_THROW(arr[7][3][1][2][6],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2,3,4),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2,3,4),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3,2,4),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3,2,4),out_of_range);
		CHECK_THROW(arr.ptr(7,3,m1,1,0),out_of_range);
		CHECK_THROW(arr.ptr(7,3,8,1,0),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,m1,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,7,2),out_of_range);
		const multi_arr<long,5,C_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2][3][4],out_of_range);
		CHECK_THROW((*carr)[10][1][2][3][4],out_of_range);
		CHECK_THROW((*carr)[7][m1][3][2][4],out_of_range);
		CHECK_THROW((*carr)[7][9][3][2][4],out_of_range);
		CHECK_THROW((*carr)[7][3][m1][1][0],out_of_range);
		CHECK_THROW((*carr)[7][3][8][1][0],out_of_range);
		CHECK_THROW((*carr)[7][3][1][m1][2],out_of_range);
		CHECK_THROW((*carr)[7][3][1][7][2],out_of_range);
		CHECK_THROW((*carr)[7][3][1][2][m1],out_of_range);
		CHECK_THROW((*carr)[7][3][1][2][6],out_of_range);
	}
	TEST_FIXTURE(LongInt6DFixtureCTypeBC,Test6DOutOfBoundsCType)
	{
		CHECK_THROW(arr[m1][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr[10][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr[7][m1][3][2][4][1],out_of_range);
		CHECK_THROW(arr[7][9][3][2][4][1],out_of_range);
		CHECK_THROW(arr[7][3][m1][1][0][2],out_of_range);
		CHECK_THROW(arr[7][3][8][1][0][2],out_of_range);
		CHECK_THROW(arr[7][3][1][m1][2][0],out_of_range);
		CHECK_THROW(arr[7][3][1][7][2][0],out_of_range);
		CHECK_THROW(arr[7][4][1][2][m1][3],out_of_range);
		CHECK_THROW(arr[7][4][1][2][6][3],out_of_range);
		CHECK_THROW(arr[7][4][1][2][3][m1],out_of_range);
		CHECK_THROW(arr[7][4][1][2][3][5],out_of_range);
		CHECK_THROW(arr.ptr(m1,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.ptr(10,1,2,3,4,0),out_of_range);
		CHECK_THROW(arr.ptr(7,m1,3,2,4,1),out_of_range);
		CHECK_THROW(arr.ptr(7,9,3,2,4,1),out_of_range);
		CHECK_THROW(arr.ptr(7,3,m1,1,0,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,8,1,0,2),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,m1,2,0),out_of_range);
		CHECK_THROW(arr.ptr(7,3,1,7,2,0),out_of_range);
		CHECK_THROW(arr.ptr(7,4,1,2,m1,3),out_of_range);
		CHECK_THROW(arr.ptr(7,4,1,2,6,3),out_of_range);
		const multi_arr<long,6,C_TYPE,true>* carr = &arr;
		CHECK_THROW((*carr)[m1][1][2][3][4][0],out_of_range);
		CHECK_THROW((*carr)[10][1][2][3][4][0],out_of_range);
		CHECK_THROW((*carr)[7][m1][3][2][4][1],out_of_range);
		CHECK_THROW((*carr)[7][9][3][2][4][1],out_of_range);
		CHECK_THROW((*carr)[7][3][m1][1][0][2],out_of_range);
		CHECK_THROW((*carr)[7][3][8][1][0][2],out_of_range);
		CHECK_THROW((*carr)[7][3][1][m1][2][0],out_of_range);
		CHECK_THROW((*carr)[7][3][1][7][2][0],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][m1][3],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][6][3],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][3][m1],out_of_range);
		CHECK_THROW((*carr)[7][4][1][2][3][5],out_of_range);
	}

	// check whether the constructor is executed for multi_arr members
	// also checks whether indirection works for pntr and const_pntr
	TEST_FIXTURE(StructWithConstructor3DFixture,TestConstructorExecuted)
	{
		multi_arr<a,3,ARPA_TYPE,false>::iterator p1 = arr.ptr(1,0,3);
		CHECK_EQUAL(23,p1->n);
		multi_arr<a,3,ARPA_TYPE,false>::const_iterator p2 = arr.ptr(1,2,0);
		CHECK_EQUAL(23,p2->n);
	}
	TEST_FIXTURE(StructWithConstructor3DFixtureCType,TestConstructorExecutedCType)
	{
		multi_arr<a,3,C_TYPE,false>::iterator p1 = arr.ptr(1,0,3);
		CHECK_EQUAL(23,p1->n);
		multi_arr<a,3,C_TYPE,false>::const_iterator p2 = arr.ptr(1,2,0);
		CHECK_EQUAL(23,p2->n);
	}

	// check whether the size(), capacity(), and empty() methods work correctly
	TEST_FIXTURE(LongInt3DFixtureBC,TestSize)
	{
		CHECK_EQUAL((size_t)(10*9*8),arr.size());
		arr.clear();
		CHECK(arr.empty());
		// check whether we can allocate a new array
		arr.alloc(10,11,12);
		CHECK_EQUAL((size_t)(10*11*12),arr.capacity());
		arr.zero();
		// check whether all the strides are properly set up
		// this is done through bounds checking
		// this also tests bounds checking const_iterator
		for (int i=0; i<10; ++i)
			for (int j=0; j<11; ++j)
			{
				multi_arr<long,3,ARPA_TYPE,true>::const_iterator p = arr.begin(i,j);
				for (int k=0; k<12; ++k)
					CHECK_EQUAL(0,p[k]);
			}
	}
		
	// check whether the size(), capacity(), and empty() methods work correctly
	// and iterator is compatible with STL algorithms
	TEST_FIXTURE(LongInt3DFixtureBC,TestSizeAlgorithm)
	{
		arr.clear();
		arr.alloc(10,11,12);
		arr.zero();

		for (int i=0; i<10; ++i)
			for (int j=0; j<11; ++j)
			{
				CHECK_EQUAL(0,count_if(arr.begin(i,j),arr.end(i,j),
					bind1st(not_equal_to<long>(),0)));
			}
	}

	TEST_FIXTURE(LongInt3DFixtureCTypeBC,TestSizeCType)
	{
		CHECK_EQUAL((size_t)(10*9*8),arr.size());
		arr.clear();
		CHECK(arr.empty());
		// check whether we can allocate a new array
		arr.alloc(10,11,12);
		CHECK_EQUAL((size_t)(10*11*12),arr.capacity());
		arr.zero();
		// check whether all the strides are properly set up
		// this is done through bounds checking
		// this also tests bounds checking const_iterator
		for (int i=0; i<10; ++i)
			for (int j=0; j<11; ++j)
			{
				multi_arr<long,3,C_TYPE,true>::const_iterator p = arr.begin(i,j);
				for (int k=0; k<12; ++k)
					CHECK_EQUAL(0,p[k]);
			}
	}

	// check whether explicit space reservation works
	TEST_FIXTURE(LongInt6DFixtureExplicitReserve,TestExplicitReserve)
	{
		arr.zero();
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					for (int l=0; l<7; ++l)
						for (int m=0; m<6; ++m)
							for (int n=0; n<5; ++n)
								CHECK_EQUAL(0,arr[i][j][k][l][m][n]);
		// it is not safe to call reserve without clearing first
		CHECK_THROW(arr.reserve(2),bad_assert);
	}

	// check whether the variant form for allocating works correctly
	// this also tests p_iterator in bounds-checking mode
	TEST_FIXTURE(TestAllocFixture,TestVariantAlloc)
	{
		CHECK_EQUAL(0,mytest());
	}

	// check whether the data() method yields valid pointer
	// for both ARPA_TYPE and C_TYPE arrays
	TEST_FIXTURE(LongInt3DFixture,TestData)
	{
		CHECK_EQUAL(&arr[0][0][0],arr.data());
		arr[0][0][0] = 1234;
		CHECK_EQUAL(1234,*arr.data());
		const multi_arr<long,3,ARPA_TYPE,false>* carr = &arr;
		CHECK_EQUAL(&arr[0][0][0],carr->data());
	}
	TEST_FIXTURE(LongInt3DFixtureCType,TestDataCType)
	{
		CHECK_EQUAL(&arr[0][0][0],arr.data());
		arr[0][0][0] = 1234;
		CHECK_EQUAL(1234,*arr.data());
		const multi_arr<long,3,C_TYPE,false>* carr = &arr;
		CHECK_EQUAL(&arr[0][0][0],carr->data());
	}

	TEST_FIXTURE(LongInt6DFixture,TestCopyOperator)
	{
		multi_arr<long,6,ARPA_TYPE,false> arr2(1,2,3,4,5,6);
		CHECK( arr.size() != arr2.size() );
		arr2.zero();
		CHECK_EQUAL(0,arr2[0][1][2][3][4][5]);
		arr2 = arr;
		CHECK( arr.size() == arr2.size() );
		// check that the copies are distinct
		CHECK( &arr[0][1][2][3][4][5] != &arr2[0][1][2][3][4][5] );
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					for (int l=0; l<7; ++l)
						for (int m=0; m<6; ++m)
							for (int n=0; n<5; ++n)
							{
								long z = i*100000+j*10000+k*1000+l*100+m*10+n;
								CHECK_EQUAL(z,arr2[i][j][k][l][m][n]);
							}

		// is it safe to copy to oneself?
		auto copy = &arr2;
		arr2 = *copy;
		// have the contents been preserved?
		CHECK_EQUAL(712344,arr2[7][1][2][3][4][4]);

		// now copy using the constructor
		multi_arr<long,6,ARPA_TYPE,false> arr3 = arr;
		CHECK( arr.size() == arr3.size() );
		// check that the copies are distinct
		CHECK( &arr[7][1][2][3][4][4] != &arr3[7][1][2][3][4][4] );
		CHECK_EQUAL(arr[7][1][2][3][4][4],arr3[7][1][2][3][4][4]);

		arr.clear();
		// copying an empty arr should clear arr2
		arr2 = arr;
		CHECK(arr2.empty());
		// also check the copy constructor
		multi_arr<long,6,ARPA_TYPE,false> arr4( arr );
		CHECK(arr4.empty());
	}

	TEST_FIXTURE(LongInt6DFixtureBC,TestCopyOperatorBC)
	{
		multi_arr<long,6,ARPA_TYPE,true> arr2(1,2,3,4,5,6);
		CHECK( arr.size() != arr2.size() );
		arr2.zero();
		CHECK_EQUAL(0,arr2[0][1][2][3][4][5]);
		arr2 = arr;
		CHECK( arr.size() == arr2.size() );
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					for (int l=0; l<7; ++l)
						for (int m=0; m<6; ++m)
							for (int n=0; n<5; ++n)
							{
								long z = i*100000+j*10000+k*1000+l*100+m*10+n;
								CHECK_EQUAL(z,arr2[i][j][k][l][m][n]);
							}

		CHECK_THROW(arr2[m1][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr2[10][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr2[7][m1][3][2][4][1],out_of_range);
		CHECK_THROW(arr2[7][9][3][2][4][1],out_of_range);
		CHECK_THROW(arr2[7][3][m1][1][0][2],out_of_range);
		CHECK_THROW(arr2[7][3][8][1][0][2],out_of_range);
		CHECK_THROW(arr2[7][3][1][m1][2][0],out_of_range);
		CHECK_THROW(arr2[7][3][1][7][2][0],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][m1][3],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][6][3],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][3][m1],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][3][5],out_of_range);

		// is it safe to copy to oneself?
		auto copy = &arr2;
		arr2 = *copy;
		// have the contents been preserved?
		CHECK_EQUAL(712344,arr2[7][1][2][3][4][4]);

		arr.clear();
		arr2 = arr;
		CHECK(arr2.empty());
		// also check the copy constructor
		multi_arr<long,6,ARPA_TYPE,true> arr4( arr );
		CHECK(arr4.empty());
	}

	TEST_FIXTURE(LongInt6DFixtureCType,TestCopyOperatorCType)
	{
		multi_arr<long,6,C_TYPE,false> arr2(1,2,3,4,5,6);
		CHECK( arr.size() != arr2.size() );
		arr2.zero();
		CHECK_EQUAL(0,arr2[0][1][2][3][4][5]);
		arr2 = arr;
		CHECK( arr.size() == arr2.size() );
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					for (int l=0; l<7; ++l)
						for (int m=0; m<6; ++m)
							for (int n=0; n<5; ++n)
							{
								long z = i*100000+j*10000+k*1000+l*100+m*10+n;
								CHECK_EQUAL(z,arr2[i][j][k][l][m][n]);
							}

		// is it safe to copy to oneself?
		auto copy = &arr2;
		arr2 = *copy;
		// have the contents been preserved?
		CHECK_EQUAL(712344,arr2[7][1][2][3][4][4]);

		arr.clear();
		arr2 = arr;
		CHECK(arr2.empty());
		// also check the copy constructor
		multi_arr<long,6,C_TYPE,false> arr4( arr );
		CHECK(arr4.empty());
	}

	TEST_FIXTURE(LongInt6DFixtureCTypeBC,TestCopyOperatorCTypeBC)
	{
		multi_arr<long,6,C_TYPE,true> arr2(1,2,3,4,5,6);
		CHECK( arr.size() != arr2.size() );
		arr2.zero();
		CHECK_EQUAL(0,arr2[0][1][2][3][4][5]);
		arr2 = arr;
		CHECK( arr.size() == arr2.size() );
		for (int i=0; i<10; ++i)
			for (int j=0; j<9; ++j)
				for (int k=0; k<8; ++k)
					for (int l=0; l<7; ++l)
						for (int m=0; m<6; ++m)
							for (int n=0; n<5; ++n)
							{
								long z = i*100000+j*10000+k*1000+l*100+m*10+n;
								CHECK_EQUAL(z,arr2[i][j][k][l][m][n]);
							}

		CHECK_THROW(arr2[m1][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr2[10][1][2][3][4][0],out_of_range);
		CHECK_THROW(arr2[7][m1][3][2][4][1],out_of_range);
		CHECK_THROW(arr2[7][9][3][2][4][1],out_of_range);
		CHECK_THROW(arr2[7][3][m1][1][0][2],out_of_range);
		CHECK_THROW(arr2[7][3][8][1][0][2],out_of_range);
		CHECK_THROW(arr2[7][3][1][m1][2][0],out_of_range);
		CHECK_THROW(arr2[7][3][1][7][2][0],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][m1][3],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][6][3],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][3][m1],out_of_range);
		CHECK_THROW(arr2[7][4][1][2][3][5],out_of_range);

		// is it safe to copy to oneself?
		auto copy = &arr2;
		arr2 = *copy;
		// have the contents been preserved?
		CHECK_EQUAL(712344,arr2[7][1][2][3][4][4]);

		arr.clear();
		arr2 = arr;
		CHECK(arr2.empty());
		// also check the copy constructor
		multi_arr<long,6,C_TYPE,true> arr4( arr );
		CHECK(arr4.empty());
	}

	// check that C_TYPE actually has standard C layout
	TEST_FIXTURE(LongInt3DCLayoutFixture,TestCLayout)
	{
		CHECK_EQUAL(0L,mytest((long(*)[10][10])arr.data()));
	}

	// test that the geometry of cloned arrays is OK
	TEST_FIXTURE(LongInt3DCloneFixture,TestCloning)
	{
		// the types of arr and dolly need not match!
		multi_arr<bool,3,ARPA_TYPE,true> dolly( arr.clone() );
		CHECK_THROW(dolly[10][0][0],out_of_range);
		for (int i=0; i<10; ++i)
		{
			CHECK_THROW(dolly[i][i+1][0],out_of_range);
			for (int j=0; j<i+1; ++j)
				CHECK_THROW(dolly[i][j][j+1],out_of_range);
		}
		// check that cloning and destroying an uninitialized multi_arr is safe
		multi_arr<long,4,ARPA_TYPE,false> arr2;
		multi_arr<bool,4,ARPA_TYPE,false> dolly2( arr2.clone() );
		CHECK(dolly2.empty());
		// check that cloning oneself is safe (not very useful though...)
		multi_arr<long,5,ARPA_TYPE,false> arr3(3,3,3,3,3);
		arr3.alloc( arr3.clone() );
		// check that the array was not cleared
		CHECK_EQUAL(243UL,arr3.size());
	}

	// same as above, but for C_TYPE arrays
	TEST_FIXTURE(LongInt3DCloneFixtureCType,TestCloningCType)
	{
		// the types of arr and dolly need not match!
		multi_arr<bool,3,C_TYPE,true> dolly( arr.clone() );
		CHECK_THROW(dolly[10][0][0],out_of_range);
		for (int i=0; i<10; ++i)
		{
			CHECK_THROW(dolly[i][i+1][0],out_of_range);
			for (int j=0; j<i+1; ++j)
				CHECK_THROW(dolly[i][j][j+1],out_of_range);
		}
		// check that cloning and destroying an uninitialized multi_arr is safe
		multi_arr<long,4,C_TYPE,false> arr2;
		multi_arr<bool,4,C_TYPE,false> dolly2( arr2.clone() );
		CHECK(dolly2.empty());
		// check that cloning oneself is safe (not very useful though...)
		multi_arr<long,5,ARPA_TYPE,false> arr3(3,3,3,3,3);
		arr3.alloc( arr3.clone() );
		// check that the array was not cleared
		CHECK_EQUAL(243UL,arr3.size());
	}

	TEST_FIXTURE(LongInt2DFixture,Test2DBeginEnd)
	{
		CHECK( arr.begin(2) == arr.ptr(2,0) );
		CHECK( arr.end(2) == arr.ptr(2,9) );
		CHECK_EQUAL(20L,arr.front(2));
		CHECK_EQUAL(28L,arr.back(2));
		const multi_arr<long,2,ARPA_TYPE,false>* carr = &arr;
		CHECK( carr->begin(2) == carr->ptr(2,0) );
		CHECK( carr->end(2) == carr->ptr(2,9) );
		CHECK_EQUAL(20L,carr->front(2));
		CHECK_EQUAL(28L,carr->back(2));
	}
	TEST_FIXTURE(LongInt3DFixture,Test3DBeginEnd)
	{
		CHECK( arr.begin(2,4) == arr.ptr(2,4,0) );
		CHECK( arr.end(2,4) == arr.ptr(2,4,8) );
		CHECK_EQUAL(240L,arr.front(2,4));
		CHECK_EQUAL(247L,arr.back(2,4));
		const multi_arr<long,3,ARPA_TYPE,false>* carr = &arr;
		CHECK( carr->begin(2,4) == carr->ptr(2,4,0) );
		CHECK( carr->end(2,4) == carr->ptr(2,4,8) );
		CHECK_EQUAL(240L,carr->front(2,4));
		CHECK_EQUAL(247L,carr->back(2,4));
	}
	TEST_FIXTURE(LongInt4DFixture,Test4DBeginEnd)
	{
		CHECK( arr.begin(2,4,7) == arr.ptr(2,4,7,0) );
		CHECK( arr.end(2,4,7) == arr.ptr(2,4,7,7) );
		CHECK_EQUAL(2470L,arr.front(2,4,7));
		CHECK_EQUAL(2476L,arr.back(2,4,7));
		const multi_arr<long,4,ARPA_TYPE,false>* carr = &arr;
		CHECK( carr->begin(2,4,7) == carr->ptr(2,4,7,0) );
		CHECK( carr->end(2,4,7) == carr->ptr(2,4,7,7) );
		CHECK_EQUAL(2470L,carr->front(2,4,7));
		CHECK_EQUAL(2476L,carr->back(2,4,7));
	}
	TEST_FIXTURE(LongInt5DFixture,Test5DBeginEnd)
	{
		CHECK( arr.begin(2,4,7,5) == arr.ptr(2,4,7,5,0) );
		CHECK( arr.end(2,4,7,5) == arr.ptr(2,4,7,5,6) );
		CHECK_EQUAL(24750L,arr.front(2,4,7,5));
		CHECK_EQUAL(24755L,arr.back(2,4,7,5));
		const multi_arr<long,5,ARPA_TYPE,false>* carr = &arr;
		CHECK( carr->begin(2,4,7,5) == carr->ptr(2,4,7,5,0) );
		CHECK( carr->end(2,4,7,5) == carr->ptr(2,4,7,5,6) );
		CHECK_EQUAL(24750L,carr->front(2,4,7,5));
		CHECK_EQUAL(24755L,carr->back(2,4,7,5));
	}
	TEST_FIXTURE(LongInt6DFixture,Test6DBeginEnd)
	{
		CHECK( arr.begin(2,4,7,5,1) == arr.ptr(2,4,7,5,1,0) );
		CHECK( arr.end(2,4,7,5,1) == arr.ptr(2,4,7,5,1,5) );
		CHECK_EQUAL(247510L,arr.front(2,4,7,5,1));
		CHECK_EQUAL(247514L,arr.back(2,4,7,5,1));
		const multi_arr<long,6,ARPA_TYPE,false>* carr = &arr;
		CHECK( carr->begin(2,4,7,5,1) == carr->ptr(2,4,7,5,1,0) );
		CHECK( carr->end(2,4,7,5,1) == carr->ptr(2,4,7,5,1,5) );
		CHECK_EQUAL(247510L,carr->front(2,4,7,5,1));
		CHECK_EQUAL(247514L,carr->back(2,4,7,5,1));
	}

	TEST_FIXTURE(LongInt2DFixtureCTypeBC,Test2DBeginEndCTypeBC)
	{
		CHECK( arr.begin(2) == arr.ptr(2,0) );
		CHECK( arr.end(2) == arr.ptr(2,9) );
		CHECK_EQUAL(20L,arr.front(2));
		CHECK_EQUAL(28L,arr.back(2));
		const multi_arr<long,2,C_TYPE,true>* carr = &arr;
		CHECK( carr->begin(2) == carr->ptr(2,0) );
		CHECK( carr->end(2) == carr->ptr(2,9) );
		CHECK_EQUAL(20L,carr->front(2));
		CHECK_EQUAL(28L,carr->back(2));
	}
	TEST_FIXTURE(LongInt3DFixtureCTypeBC,Test3DBeginEndCTypeBC)
	{
		CHECK( arr.begin(2,4) == arr.ptr(2,4,0) );
		CHECK( arr.end(2,4) == arr.ptr(2,4,8) );
		CHECK_EQUAL(240L,arr.front(2,4));
		CHECK_EQUAL(247L,arr.back(2,4));
		const multi_arr<long,3,C_TYPE,true>* carr = &arr;
		CHECK( carr->begin(2,4) == carr->ptr(2,4,0) );
		CHECK( carr->end(2,4) == carr->ptr(2,4,8) );
		CHECK_EQUAL(240L,carr->front(2,4));
		CHECK_EQUAL(247L,carr->back(2,4));
	}
	TEST_FIXTURE(LongInt4DFixtureCTypeBC,Test4DBeginEndCTypeBC)
	{
		CHECK( arr.begin(2,4,7) == arr.ptr(2,4,7,0) );
		CHECK( arr.end(2,4,7) == arr.ptr(2,4,7,7) );
		CHECK_EQUAL(2470L,arr.front(2,4,7));
		CHECK_EQUAL(2476L,arr.back(2,4,7));
		const multi_arr<long,4,C_TYPE,true>* carr = &arr;
		CHECK( carr->begin(2,4,7) == carr->ptr(2,4,7,0) );
		CHECK( carr->end(2,4,7) == carr->ptr(2,4,7,7) );
		CHECK_EQUAL(2470L,carr->front(2,4,7));
		CHECK_EQUAL(2476L,carr->back(2,4,7));
	}
	TEST_FIXTURE(LongInt5DFixtureCTypeBC,Test5DBeginEndCTypeBC)
	{
		CHECK( arr.begin(2,4,7,5) == arr.ptr(2,4,7,5,0) );
		CHECK( arr.end(2,4,7,5) == arr.ptr(2,4,7,5,6) );
		CHECK_EQUAL(24750L,arr.front(2,4,7,5));
		CHECK_EQUAL(24755L,arr.back(2,4,7,5));
		const multi_arr<long,5,C_TYPE,true>* carr = &arr;
		CHECK( carr->begin(2,4,7,5) == carr->ptr(2,4,7,5,0) );
		CHECK( carr->end(2,4,7,5) == carr->ptr(2,4,7,5,6) );
		CHECK_EQUAL(24750L,carr->front(2,4,7,5));
		CHECK_EQUAL(24755L,carr->back(2,4,7,5));
	}
	TEST_FIXTURE(LongInt6DFixtureCTypeBC,Test6DBeginEndCTypeBC)
	{
		CHECK( arr.begin(2,4,7,5,1) == arr.ptr(2,4,7,5,1,0) );
		CHECK( arr.end(2,4,7,5,1) == arr.ptr(2,4,7,5,1,5) );
		CHECK_EQUAL(247510L,arr.front(2,4,7,5,1));
		CHECK_EQUAL(247514L,arr.back(2,4,7,5,1));
		const multi_arr<long,6,C_TYPE,true>* carr = &arr;
		CHECK( carr->begin(2,4,7,5,1) == carr->ptr(2,4,7,5,1,0) );
		CHECK( carr->end(2,4,7,5,1) == carr->ptr(2,4,7,5,1,5) );
		CHECK_EQUAL(247510L,carr->front(2,4,7,5,1));
		CHECK_EQUAL(247514L,carr->back(2,4,7,5,1));
	}

	// can an indexed array element be used in variable length argument lists?
	TEST_FIXTURE(LongInt3DFixture,Test3DVarLengthArgument)
	{
		char buf[100];
		sprintf( buf, "%ld", arr[2][3][4] );
		long res;
		sscanf( buf, "%ld", &res );
		CHECK_EQUAL(234L,res);
	}

	TEST_FIXTURE(LongInt3DFixtureCType,Test3DVarLengthArgumentCType)
	{
		char buf[100];
		sprintf( buf, "%ld", arr[2][3][4] );
		long res;
		sscanf( buf, "%ld", &res );
		CHECK_EQUAL(234L,res);
	}

	TEST_FIXTURE(LongInt2DEmptyDim,Test2DEmptyDimIterator)
	{
		// this should not crash
		multi_arr<long,2,ARPA_TYPE,true>::const_iterator p = arr.begin(0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt3DEmptyDim,Test3DEmptyDimIterator)
	{
		// this should not crash
		multi_arr<long,3,ARPA_TYPE,true>::const_iterator p = arr.begin(0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt4DEmptyDim,Test4DEmptyDimIterator)
	{
		// this should not crash
		multi_arr<long,4,ARPA_TYPE,true>::const_iterator p = arr.begin(0,0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt5DEmptyDim,Test5DEmptyDimIterator)
	{
		// this should not crash
		multi_arr<long,5,ARPA_TYPE,true>::const_iterator p = arr.begin(0,0,0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt6DEmptyDim,Test6DEmptyDimIterator)
	{
		// this should not crash
		multi_arr<long,6,ARPA_TYPE,true>::const_iterator p = arr.begin(0,0,0,0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt2DEmptyDimCType,Test2DEmptyDimIteratorCType)
	{
		// this should not crash
		multi_arr<long,2,C_TYPE,true>::const_iterator p = arr.begin(0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt3DEmptyDimCType,Test3DEmptyDimIteratorCType)
	{
		// this should not crash
		multi_arr<long,3,C_TYPE,true>::const_iterator p = arr.begin(0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt4DEmptyDimCType,Test4DEmptyDimIteratorCType)
	{
		// this should not crash
		multi_arr<long,4,C_TYPE,true>::const_iterator p = arr.begin(0,0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt5DEmptyDimCType,Test5DEmptyDimIteratorCType)
	{
		// this should not crash
		multi_arr<long,5,C_TYPE,true>::const_iterator p = arr.begin(0,0,0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}

	TEST_FIXTURE(LongInt6DEmptyDimCType,Test6DEmptyDimIteratorCType)
	{
		// this should not crash
		multi_arr<long,6,C_TYPE,true>::const_iterator p = arr.begin(0,0,0,0,0);
		// bogus test so that variable gets used
		CHECK( p == p );
	}
}
