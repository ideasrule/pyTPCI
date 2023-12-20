#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "ran.h"
#include "vectorize.h"
#include "vectorhash.h"

namespace {

	TEST(TestAvxPool)
	{
		// there is not a lot we can check here since the interesting bits are hidden
		// behind the interface, so the real testing will be done with valgrind-test
		t_avx_pool pool;
		for( int i=0; i < 2048; ++i )
		{
			size_t sz = (ran.u32()%200000) + 1;
			char *p = static_cast<char*>(pool.avx_alloc(sz));
			size_t ip = reinterpret_cast<size_t>(p);
			// test for correct alignment
			CHECK( p != NULL && (ip/CD_ALIGN)*CD_ALIGN == ip );
			for( size_t j=0; j < sz; ++j )
				p[j] = 'c';
			pool.avx_free(p);
		}
		// check that many simultaneous requests work correctly
		const int ntest = 200;
		char* pp[ntest];
		for( int i=0; i < ntest; ++i )
		{
			pp[i] = static_cast<char*>(pool.avx_alloc(4000));
			size_t ip = reinterpret_cast<size_t>(pp[i]);
			CHECK( pp[i] != NULL && (ip/CD_ALIGN)*CD_ALIGN == ip );
			for( size_t j=0; j < 4000; ++j )
				pp[i][j] = 'c';
		}
		for( int i=0; i < ntest; ++i )
			pool.avx_free(pp[i]);
	}

	TEST(TestAvxPtr)
	{
		avx_ptr<double,true> p1(100);
		CHECK( p1.data() != NULL && p1.data() == p1.ptr0() );
		for( int i=0; i < 100; ++i )
			p1[i] = 1.;
		CHECK_THROW(( p1[-1] = 1. ),out_of_range);
		CHECK_THROW(( p1[100] = 1. ),out_of_range);
		avx_ptr<double> p2(0);
		CHECK( p2.data() == NULL );
		avx_ptr<double> p3(-10);
		CHECK( p3.data() == NULL );
		avx_ptr<double,true> p4(-10,100);
		CHECK( p4.data() != NULL && p4.data() == p4.ptr0()-10 );
		for( int i=-10; i < 100; ++i )
			p4[i] = 1.;
		CHECK_THROW(( p4[-11] = 1. ),out_of_range);
		CHECK_THROW(( p1[100] = 1. ),out_of_range);
		avx_ptr<double> p5(10,10);
		CHECK( p5.data() == NULL );
		avx_ptr<double> p6(10,-10);
		CHECK( p6.data() == NULL );
	}

	TEST(TestAllocatorAvx)
	{
		// there is not a lot we can check here since the interesting bits are hidden
		// behind the interface, so the real testing will be done with valgrind-test
		for( int i=0; i < 2048; ++i )
		{
			size_t sz = (ran.u32()%200000) + 1;
			vector_avx<int> a(sz);
			size_t ip = reinterpret_cast<size_t>(get_ptr(a));
			CHECK( ip != 0 && (ip/CD_ALIGN)*CD_ALIGN == ip );
			for( size_t j=0; j < sz; ++j )
				a[j] = 3;
		}
	}

	TEST(TestReduceAd)
	{
		double a[32];
		for( int i=0; i < 32; ++i )
			a[i] = double(i);
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val = reduce_a( a, ilo, ihi );
				double expect = double(ihi*(ihi-1)/2 - ilo*(ilo-1)/2);
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceAf)
	{
		sys_float a[32];
		for( int i=0; i < 32; ++i )
			a[i] = sys_float(i);
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				sys_float val = reduce_a( a, ilo, ihi );
				sys_float expect = sys_float(ihi*(ihi-1)/2 - ilo*(ilo-1)/2);
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceABd)
	{
		double a[32], b[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = double(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val = reduce_ab( a, b, ilo, ihi );
				double expect = double(ihi*(ihi-1)*(2*ihi-1)/6 - ilo*(ilo-1)*(2*ilo-1)/6);
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceABdf)
	{
		double a[32];
		sys_float b[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = sys_float(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val = reduce_ab( a, b, ilo, ihi );
				double expect = double(ihi*(ihi-1)*(2*ihi-1)/6 - ilo*(ilo-1)*(2*ilo-1)/6);
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceABf)
	{
		sys_float a[32], b[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = sys_float(i);
			b[i] = sys_float(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				sys_float val = reduce_ab( a, b, ilo, ihi );
				sys_float expect = sys_float(ihi*(ihi-1)*(2*ihi-1)/6 - ilo*(ilo-1)*(2*ilo-1)/6);
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceABCd)
	{
		double a[32], b[32], c[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = double(i);
			c[i] = double(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val = reduce_abc( a, b, c, ilo, ihi );
				double expect = double(pow2(ihi*(ihi-1)/2) - pow2(ilo*(ilo-1)/2));
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceABCddf)
	{
		double a[32], b[32];
		sys_float c[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = double(i);
			c[i] = sys_float(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val = reduce_abc( a, b, c, ilo, ihi );
				double expect = double(pow2(ihi*(ihi-1)/2) - pow2(ilo*(ilo-1)/2));
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceABCdff)
	{
		double a[32];
		sys_float b[32], c[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = sys_float(i);
			c[i] = sys_float(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val = reduce_abc( a, b, c, ilo, ihi );
				double expect = double(pow2(ihi*(ihi-1)/2) - pow2(ilo*(ilo-1)/2));
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceABCf)
	{
		sys_float a[32], b[32], c[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = sys_float(i);
			b[i] = sys_float(i);
			c[i] = sys_float(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				sys_float val = reduce_abc( a, b, c, ilo, ihi );
				sys_float expect = sys_float(pow2(ihi*(ihi-1)/2) - pow2(ilo*(ilo-1)/2));
				CHECK( fp_equal( val, expect ) );
			}
		}
	}

	TEST(TestReduceAB_Ad)
	{
		double a[32], b[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = double(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val2, val = reduce_ab_a( a, b, ilo, ihi, &val2 );
				double expect = double(ihi*(ihi-1)*(2*ihi-1)/6 - ilo*(ilo-1)*(2*ilo-1)/6);
				double expect2 = double(ihi*(ihi-1)/2 - ilo*(ilo-1)/2);
				CHECK( fp_equal( val, expect ) );
				CHECK( fp_equal( val2, expect2 ) );
			}
		}
	}

	TEST(TestReduceAB_Adf)
	{
		double a[32];
		sys_float b[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = sys_float(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val2, val = reduce_ab_a( a, b, ilo, ihi, &val2 );
				double expect = double(ihi*(ihi-1)*(2*ihi-1)/6 - ilo*(ilo-1)*(2*ilo-1)/6);
				double expect2 = double(ihi*(ihi-1)/2 - ilo*(ilo-1)/2);
				CHECK( fp_equal( val, expect ) );
				CHECK( fp_equal( val2, expect2 ) );
				val = reduce_ab_a( b, a, ilo, ihi, &val2 );
				CHECK( fp_equal( val, expect ) );
				CHECK( fp_equal( val2, expect2 ) );
			}
		}
	}

	TEST(TestReduceAB_Af)
	{
		double a[32], b[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = double(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val2, val = reduce_ab_a( a, b, ilo, ihi, &val2 );
				double expect = double(ihi*(ihi-1)*(2*ihi-1)/6 - ilo*(ilo-1)*(2*ilo-1)/6);
				double expect2 = double(ihi*(ihi-1)/2 - ilo*(ilo-1)/2);
				CHECK( fp_equal( val, expect ) );
				CHECK( fp_equal( val2, expect2 ) );
			}
		}
	}

	TEST(TestReduceABC_ABd)
	{
		double a[32], b[32], c[32];
		for( int i=0; i < 32; ++i )
		{
			a[i] = double(i);
			b[i] = double(i);
			c[i] = double(i);
		}
		for( int ilo=0; ilo < 32; ilo++ )
		{
			for( int ihi=ilo+1; ihi <= 32; ihi++ )
			{
				double val2, val = reduce_abc_ab( a, b, c, ilo, ihi, &val2 );
				double expect = double(pow2(ihi*(ihi-1)/2) - pow2(ilo*(ilo-1)/2));
				double expect2 = double(ihi*(ihi-1)*(2*ihi-1)/6 - ilo*(ilo-1)*(2*ilo-1)/6);
				CHECK( fp_equal( val, expect ) );
				CHECK( fp_equal( val2, expect2 ) );
			}
		}
	}

	TEST(TestExp10d)
	{
		// exp10(double) is defined in cddefines.h
		double maxarg = log10(DBL_MAX);
		for( int i=0; i < 2048; ++i )
		{
			double arg = ran.dbl()*633. - 324.;
			if( arg < maxarg )
			{
				double val1 = exp10(arg);
				double val2 = pow(10.,arg);
				// inaccuracies in pow(10.,x) will contribute to the error as well...
				// so this test may differ a bit from one glibc version to another.
				CHECK( fp_equal( val1, val2, 10 ) );
			}
		}
	}

	TEST(TestExp10f)
	{
		// exp10(sys_float) is defined in cddefines.h
		sys_float maxarg = log10f(FLT_MAX);
		for( int i=0; i < 2048; ++i )
		{
			sys_float arg = ran.rnm()*84.f - 45.f;
			if( arg < maxarg )
			{
				sys_float val1 = exp10(arg);
				sys_float val2 = sys_float(pow(10.,double(arg)));
				CHECK( fp_equal( val1, val2, 6 ) );
			}
		}
	}

	TEST(TestVexpd)
	{
		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];
#ifdef __AVX__
		double minarg = log(DBL_MIN);
		// check correct behavior near the underflow margin
		// our implementation doesn't support gradual underflow
		double a1 = nextafter(minarg,-DBL_MAX);
		vexp(val,a1,a1,a1,a1);
		CHECK( val[0] == 0. );
		double a2 = nextafter(minarg,DBL_MAX);
		vexp(val,a2,a2,a2,a2);
		CHECK( val[0] > 0. );
		// check correct behavior near the overflow margin
		double maxarg = log(DBL_MAX);
		a1 = nextafter(maxarg,-DBL_MAX);
		vexp(val,a1,a1,a1,a1);
		CHECK( val[0] < DBL_MAX );
		a2 = nextafter(maxarg,DBL_MAX);
		CHECK_THROW( vexp(val,a2,a2,a2,a2), domain_error );
		a2 = -numeric_limits<double>().infinity();
		vexp(val,a2,a2,a2,a2);
		CHECK( val[0] == 0. );
		a2 = numeric_limits<double>().infinity();
		CHECK_THROW( vexp(val,a2,a2,a2,a2), domain_error );
		a2 = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vexp(val,a2,a2,a2,a2), domain_error );
		// now check results over the entire range
		for( int i=0; i < sz; )
		{
			double x = (ran.dbl()-0.5)*1420.;
			if( x < maxarg )
				arg[i++] = x;
		}
		vexp( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
		{
			if( arg[i] < minarg )
				CHECK( val[i] == 0. );
			else
				CHECK( fp_equal( val[i], exp(arg[i]) ) );
		}
#endif
		for( int i=0; i < 32; ++i )
			arg[i] = double(i);
		// finally check that non-aligned arrays are treated correctly
		for( int i=0; i < 32; ++i )
		{
			for( int j=i+1; j <= 32; ++j )
			{
				invalidate_array( val, 32*sizeof(double) );
				vexp( arg, val, i, j );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], exp(arg[k]) ) );
				}
				invalidate_array( val, 32*sizeof(double) );
				vexp( &arg[i], val, 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k >= j-i )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], exp(arg[k+i]) ) );
				}
			}
		}
	}

	TEST(TestVexp10d)
	{
#ifdef __AVX__
		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];

		double minarg = log10(DBL_MIN);
		double a1 = nextafter(minarg,-DBL_MAX);
		vexp10(val,a1,a1,a1,a1);
		CHECK( val[0] == 0. );
		double a2 = nextafter(minarg,DBL_MAX);
		vexp10(val,a2,a2,a2,a2);
		CHECK( val[0] > 0. );
		double maxarg = log10(DBL_MAX);
		a1 = nextafter(maxarg,-DBL_MAX);
		vexp10(val,a1,a1,a1,a1);
		CHECK( val[0] < DBL_MAX );
		a2 = nextafter(maxarg,DBL_MAX);
		CHECK_THROW( vexp10(val,a2,a2,a2,a2), domain_error );
		a2 = -numeric_limits<double>().infinity();
		vexp10(val,a2,a2,a2,a2);
		CHECK( val[0] == 0. );
		a2 = numeric_limits<double>().infinity();
		CHECK_THROW( vexp10(val,a2,a2,a2,a2), domain_error );
		a2 = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vexp10(val,a2,a2,a2,a2), domain_error );
		for( int i=0; i < sz; )
		{
			double x = (ran.dbl()-0.5)*617.;
			if( x < maxarg )
				arg[i++] = x;
		}
		vexp10( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
		{
			if( arg[i] < minarg )
				CHECK( val[i] == 0. );
			else
				CHECK( fp_equal( val[i], pow(10.,arg[i]), 10 ) );
		}
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVexpm1d)
	{
#ifdef __AVX__
		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];

		double maxarg = log(DBL_MAX);
		double a1 = nextafter(maxarg,-DBL_MAX);
		vexpm1(val,a1,a1,a1,a1);
		CHECK( val[0] < DBL_MAX );
		double a2 = nextafter(maxarg,DBL_MAX);
		CHECK_THROW( vexpm1(val,a2,a2,a2,a2), domain_error );
		a2 = -numeric_limits<double>().infinity();
		vexpm1(val,a2,a2,a2,a2);
		CHECK( val[0] == -1. );
		a2 = numeric_limits<double>().infinity();
		CHECK_THROW( vexpm1(val,a2,a2,a2,a2), domain_error );
		a2 = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vexpm1(val,a2,a2,a2,a2), domain_error );
		double y[8];
		vexpm1( y, 1.e-10, 8.e-8, 4.e-5, 6.e-3, 2.e-1, 1., 1., 1. );
		// constants below were derived from Abramowitz & Stegun, Handbook of Mathematical Functions
		CHECK( fp_equal( y[0], 1.000000000050000e-10 ) );
		CHECK( fp_equal( y[1], 8.00000032000000853e-8 ) );
		CHECK( fp_equal( y[2], 4.00008000106667733342e-5 ) );
		CHECK( fp_equal( y[3], 6.0180360540648648555845e-3 ) );
		CHECK( fp_equal( y[4], 2.214027581601698339210720e-1 ) );
		CHECK( fp_equal( y[5], 1.7182818284590452353602875 ) );		
		for( int i=0; i < sz; )
		{
			double x = ran.dbl()*749.-39.;
			if( x < maxarg )
				arg[i++] = x;
		}
		vexpm1( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], expm1(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVexpf)
	{
		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
#ifdef __AVX__
		sys_float minarg = logf(FLT_MIN);
		sys_float a1 = nextafterf(minarg,-FLT_MAX);
		vexp(val,a1,a1,a1,a1);
		CHECK( val[0] == 0.f );
		sys_float a2 = nextafterf(minarg,FLT_MAX);
		vexp(val,a2,a2,a2,a2);
		CHECK( val[0] > 0.f );
		sys_float maxarg = logf(FLT_MAX);
		a1 = nextafterf(maxarg,-FLT_MAX);
		vexp(val,a1,a1,a1,a1);
		CHECK( val[0] < FLT_MAX );
		a2 = nextafterf(maxarg,FLT_MAX);
		CHECK_THROW( vexp(val,a2,a2,a2,a2), domain_error );
		a2 = -numeric_limits<sys_float>().infinity();
		vexp(val,a2,a2,a2,a2);
		CHECK( val[0] == 0.f );
		a2 = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vexp(val,a2,a2,a2,a2), domain_error );
		a2 = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vexp(val,a2,a2,a2,a2), domain_error );
		for( int i=0; i < sz; )
		{
			sys_float x = (ran.rnm()-0.5f)*178.f;
			if( x < maxarg )
				arg[i++] = x;
		}
		vexp( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
		{
			if( arg[i] < minarg )
				CHECK( val[i] == 0.f );
			else
				CHECK( fp_equal( val[i], expf(arg[i]) ) );
		}
#endif
		for( int i=0; i < 32; ++i )
			arg[i] = sys_float(i);
		for( int i=0; i < 32; ++i )
		{
			for( int j=i+1; j <= 32; ++j )
			{
				invalidate_array( val, 32*sizeof(sys_float) );
				vexp( arg, val, i, j );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], expf(arg[k]) ) );
				}
				invalidate_array( val, 32*sizeof(sys_float) );
				vexp( &arg[i], val, 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k >= j-i )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], expf(arg[k+i]) ) );
				}
			}
		}
	}

	TEST(TestVexp10f)
	{
#ifdef __AVX__
		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];

		sys_float minarg = log10f(FLT_MIN);
		sys_float a1 = nextafterf(minarg,-FLT_MAX);
		vexp10(val,a1,a1,a1,a1);
		CHECK( val[0] == 0.f );
		sys_float a2 = nextafterf(minarg,FLT_MAX);
		vexp10(val,a2,a2,a2,a2);
		CHECK( val[0] > 0.f );
		sys_float maxarg = log10f(FLT_MAX);
		a1 = nextafterf(maxarg,-FLT_MAX);
		vexp10(val,a1,a1,a1,a1);
		CHECK( val[0] < FLT_MAX );
		a2 = nextafterf(maxarg,FLT_MAX);
		CHECK_THROW( vexp10(val,a2,a2,a2,a2), domain_error );
		a2 = -numeric_limits<sys_float>().infinity();
		vexp10(val,a2,a2,a2,a2);
		CHECK( val[0] == 0.f );
		a2 = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vexp10(val,a2,a2,a2,a2), domain_error );
		a2 = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vexp10(val,a2,a2,a2,a2), domain_error );
		for( int i=0; i < sz; )
		{
			sys_float x = (ran.rnm()-0.5f)*78.f;
			if( x < maxarg )
				arg[i++] = x;
		}
		vexp10( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
		{
			if( arg[i] < minarg )
				CHECK( val[i] == 0.f );
			else
				CHECK( fp_equal( val[i], powf(10.f,arg[i]), 8 ) );
		}
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVexpm1f)
	{
#ifdef __AVX__
		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];

		sys_float maxarg = logf(FLT_MAX);
		sys_float a1 = nextafterf(maxarg,-FLT_MAX);
		vexpm1(val,a1,a1,a1,a1);
		CHECK( val[0] < FLT_MAX );
		sys_float a2 = nextafterf(maxarg,FLT_MAX);
		CHECK_THROW( vexpm1(val,a2,a2,a2,a2), domain_error );
		a2 = -numeric_limits<sys_float>().infinity();
		vexpm1(val,a2,a2,a2,a2);
		CHECK( val[0] == -1.f );
		a2 = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vexpm1(val,a2,a2,a2,a2), domain_error );
		a2 = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vexpm1(val,a2,a2,a2,a2), domain_error );
		sys_float y[8];
		vexpm1( y, 1.e-10f, 8.e-8f, 4.e-5f, 6.e-3f, 2.e-1f, 1.f, 1.f, 1.f );
		// constants below were derived from Abramowitz & Stegun, Handbook of Mathematical Functions
		CHECK( fp_equal( y[0], 1.000000000050000e-10f ) );
		CHECK( fp_equal( y[1], 8.00000032000000853e-8f ) );
		CHECK( fp_equal( y[2], 4.00008000106667733342e-5f ) );
		CHECK( fp_equal( y[3], 6.0180360540648648555845e-3f ) );
		CHECK( fp_equal( y[4], 2.214027581601698339210720e-1f ) );
		CHECK( fp_equal( y[5], 1.7182818284590452353602875f ) );		
		for( int i=0; i < sz; )
		{
			sys_float x = ran.rnm()*108.f-19.f;
			if( x < maxarg )
				arg[i++] = x;
		}
		vexpm1( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], expm1f(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVexpdx4)
	{
		double y[4];
		vexp(y,1.1,1.2,1.3,1.4);
		CHECK( fp_equal( y[0], exp(1.1) ) );
		CHECK( fp_equal( y[1], exp(1.2) ) );
		CHECK( fp_equal( y[2], exp(1.3) ) );
		CHECK( fp_equal( y[3], exp(1.4) ) );
	}

	TEST(TestVexp10dx4)
	{
		double y[4];
		vexp10(y,1.1,1.2,1.3,1.4);
		CHECK( fp_equal( y[0], pow(10.,1.1) ) );
		CHECK( fp_equal( y[1], pow(10.,1.2) ) );
		CHECK( fp_equal( y[2], pow(10.,1.3) ) );
		CHECK( fp_equal( y[3], pow(10.,1.4) ) );
	}

	TEST(TestVexpm1dx4)
	{
		double y[4];
		vexpm1(y,1.e-10,1.2,1.3,1.4);
		CHECK( fp_equal( y[0], expm1(1.e-10) ) );
		CHECK( fp_equal( y[1], expm1(1.2) ) );
		CHECK( fp_equal( y[2], expm1(1.3) ) );
		CHECK( fp_equal( y[3], expm1(1.4) ) );
	}

	TEST(TestVexpdx8)
	{
		double y[8];
		vexp(y,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8);
		CHECK( fp_equal( y[0], exp(1.1) ) );
		CHECK( fp_equal( y[1], exp(1.2) ) );
		CHECK( fp_equal( y[2], exp(1.3) ) );
		CHECK( fp_equal( y[3], exp(1.4) ) );
		CHECK( fp_equal( y[4], exp(1.5) ) );
		CHECK( fp_equal( y[5], exp(1.6) ) );
		CHECK( fp_equal( y[6], exp(1.7) ) );
		CHECK( fp_equal( y[7], exp(1.8) ) );
	}

	TEST(TestVexp10dx8)
	{
		double y[8];
		vexp10(y,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8);
		CHECK( fp_equal( y[0], pow(10.,1.1) ) );
		CHECK( fp_equal( y[1], pow(10.,1.2) ) );
		CHECK( fp_equal( y[2], pow(10.,1.3) ) );
		CHECK( fp_equal( y[3], pow(10.,1.4) ) );
		CHECK( fp_equal( y[4], pow(10.,1.5) ) );
		CHECK( fp_equal( y[5], pow(10.,1.6) ) );
		CHECK( fp_equal( y[6], pow(10.,1.7) ) );
		CHECK( fp_equal( y[7], pow(10.,1.8) ) );
	}

	TEST(TestVexpm1dx8)
	{
		double y[8];
		vexpm1(y,1.e-10,1.2,1.3,1.4,1.5,1.6,1.7,1.8);
		CHECK( fp_equal( y[0], expm1(1.e-10) ) );
		CHECK( fp_equal( y[1], expm1(1.2) ) );
		CHECK( fp_equal( y[2], expm1(1.3) ) );
		CHECK( fp_equal( y[3], expm1(1.4) ) );
		CHECK( fp_equal( y[4], expm1(1.5) ) );
		CHECK( fp_equal( y[5], expm1(1.6) ) );
		CHECK( fp_equal( y[6], expm1(1.7) ) );
		CHECK( fp_equal( y[7], expm1(1.8) ) );
	}

	TEST(TestVexpfx4)
	{
		sys_float y[4];
		vexp(y,1.1f,1.2f,1.3f,1.4f);
		CHECK( fp_equal( y[0], expf(1.1f) ) );
		CHECK( fp_equal( y[1], expf(1.2f) ) );
		CHECK( fp_equal( y[2], expf(1.3f) ) );
		CHECK( fp_equal( y[3], expf(1.4f) ) );
	}

	TEST(TestVexp10fx4)
	{
		sys_float y[4];
		vexp10(y,1.1f,1.2f,1.3f,1.4f);
		CHECK( fp_equal( y[0], powf(10.f,1.1f) ) );
		CHECK( fp_equal( y[1], powf(10.f,1.2f) ) );
		CHECK( fp_equal( y[2], powf(10.f,1.3f) ) );
		CHECK( fp_equal( y[3], powf(10.f,1.4f) ) );
	}

	TEST(TestVexpm1fx4)
	{
		sys_float y[4];
		vexpm1(y,1.e-10f,1.2f,1.3f,1.4f);
		CHECK( fp_equal( y[0], expm1f(1.e-10f) ) );
		CHECK( fp_equal( y[1], expm1f(1.2f) ) );
		CHECK( fp_equal( y[2], expm1f(1.3f) ) );
		CHECK( fp_equal( y[3], expm1f(1.4f) ) );
	}

	TEST(TestVexpfx8)
	{
		sys_float y[8];
		vexp(y,1.1f,1.2f,1.3f,1.4f,1.5f,1.6f,1.7f,1.8f);
		CHECK( fp_equal( y[0], expf(1.1f) ) );
		CHECK( fp_equal( y[1], expf(1.2f) ) );
		CHECK( fp_equal( y[2], expf(1.3f) ) );
		CHECK( fp_equal( y[3], expf(1.4f) ) );
		CHECK( fp_equal( y[4], expf(1.5f) ) );
		CHECK( fp_equal( y[5], expf(1.6f) ) );
		CHECK( fp_equal( y[6], expf(1.7f) ) );
		CHECK( fp_equal( y[7], expf(1.8f) ) );
	}

	TEST(TestVexp10fx8)
	{
		sys_float y[8];
		vexp10(y,1.1f,1.2f,1.3f,1.4f,1.5f,1.6f,1.7f,1.8f);
		CHECK( fp_equal( y[0], powf(10.f,1.1f) ) );
		CHECK( fp_equal( y[1], powf(10.f,1.2f) ) );
		CHECK( fp_equal( y[2], powf(10.f,1.3f) ) );
		CHECK( fp_equal( y[3], powf(10.f,1.4f) ) );
		CHECK( fp_equal( y[4], powf(10.f,1.5f) ) );
		CHECK( fp_equal( y[5], powf(10.f,1.6f) ) );
		CHECK( fp_equal( y[6], powf(10.f,1.7f) ) );
		CHECK( fp_equal( y[7], powf(10.f,1.8f) ) );
	}

	TEST(TestVexpm1fx8)
	{
		sys_float y[8];
		vexpm1(y,1.e-10f,1.2f,1.3f,1.4f,1.5f,1.6f,1.7f,1.8f);
		CHECK( fp_equal( y[0], expm1f(1.e-10f) ) );
		CHECK( fp_equal( y[1], expm1f(1.2f) ) );
		CHECK( fp_equal( y[2], expm1f(1.3f) ) );
		CHECK( fp_equal( y[3], expm1f(1.4f) ) );
		CHECK( fp_equal( y[4], expm1f(1.5f) ) );
		CHECK( fp_equal( y[5], expm1f(1.6f) ) );
		CHECK( fp_equal( y[6], expm1f(1.7f) ) );
		CHECK( fp_equal( y[7], expm1f(1.8f) ) );
	}

	TEST(TestVexpfx16)
	{
		sys_float y[16];
		vexp(y,1.1f,1.2f,1.3f,1.4f,1.5f,1.6f,1.7f,1.8f,2.11f,2.22f,2.33f,2.44f,2.55f,2.66f,2.77f,2.88f);
		CHECK( fp_equal( y[0], expf(1.1f) ) );
		CHECK( fp_equal( y[1], expf(1.2f) ) );
		CHECK( fp_equal( y[2], expf(1.3f) ) );
		CHECK( fp_equal( y[3], expf(1.4f) ) );
		CHECK( fp_equal( y[4], expf(1.5f) ) );
		CHECK( fp_equal( y[5], expf(1.6f) ) );
		CHECK( fp_equal( y[6], expf(1.7f) ) );
		CHECK( fp_equal( y[7], expf(1.8f) ) );
		CHECK( fp_equal( y[8], expf(2.11f) ) );
		CHECK( fp_equal( y[9], expf(2.22f) ) );
		CHECK( fp_equal( y[10], expf(2.33f) ) );
		CHECK( fp_equal( y[11], expf(2.44f) ) );
		CHECK( fp_equal( y[12], expf(2.55f) ) );
		CHECK( fp_equal( y[13], expf(2.66f) ) );
		CHECK( fp_equal( y[14], expf(2.77f) ) );
		CHECK( fp_equal( y[15], expf(2.88f) ) );
	}

	TEST(TestVexp10fx16)
	{
		sys_float y[16];
		vexp10(y,1.1f,1.2f,1.3f,1.4f,1.5f,1.6f,1.7f,1.8f,2.11f,2.22f,2.33f,2.44f,2.55f,2.66f,2.77f,2.88f);
		CHECK( fp_equal( y[0], powf(10.f,1.1f) ) );
		CHECK( fp_equal( y[1], powf(10.f,1.2f) ) );
		CHECK( fp_equal( y[2], powf(10.f,1.3f) ) );
		CHECK( fp_equal( y[3], powf(10.f,1.4f) ) );
		CHECK( fp_equal( y[4], powf(10.f,1.5f) ) );
		CHECK( fp_equal( y[5], powf(10.f,1.6f) ) );
		CHECK( fp_equal( y[6], powf(10.f,1.7f) ) );
		CHECK( fp_equal( y[7], powf(10.f,1.8f) ) );
		CHECK( fp_equal( y[8], powf(10.f,2.11f) ) );
		CHECK( fp_equal( y[9], powf(10.f,2.22f) ) );
		CHECK( fp_equal( y[10], powf(10.f,2.33f) ) );
		CHECK( fp_equal( y[11], powf(10.f,2.44f) ) );
		CHECK( fp_equal( y[12], powf(10.f,2.55f) ) );
		CHECK( fp_equal( y[13], powf(10.f,2.66f) ) );
		CHECK( fp_equal( y[14], powf(10.f,2.77f) ) );
		CHECK( fp_equal( y[15], powf(10.f,2.88f) ) );
	}

	TEST(TestVexpm1fx16)
	{
		sys_float y[16];
		vexpm1(y,1.e-10f,1.2f,1.3f,1.4f,1.5f,1.6f,1.7f,1.8f,2.11f,2.22f,2.33f,2.44f,2.55f,2.66f,2.77f,2.88f);
		CHECK( fp_equal( y[0], expm1f(1.e-10f) ) );
		CHECK( fp_equal( y[1], expm1f(1.2f) ) );
		CHECK( fp_equal( y[2], expm1f(1.3f) ) );
		CHECK( fp_equal( y[3], expm1f(1.4f) ) );
		CHECK( fp_equal( y[4], expm1f(1.5f) ) );
		CHECK( fp_equal( y[5], expm1f(1.6f) ) );
		CHECK( fp_equal( y[6], expm1f(1.7f) ) );
		CHECK( fp_equal( y[7], expm1f(1.8f) ) );
		CHECK( fp_equal( y[8], expm1f(2.11f) ) );
		CHECK( fp_equal( y[9], expm1f(2.22f) ) );
		CHECK( fp_equal( y[10], expm1f(2.33f) ) );
		CHECK( fp_equal( y[11], expm1f(2.44f) ) );
		CHECK( fp_equal( y[12], expm1f(2.55f) ) );
		CHECK( fp_equal( y[13], expm1f(2.66f) ) );
		CHECK( fp_equal( y[14], expm1f(2.77f) ) );
		CHECK( fp_equal( y[15], expm1f(2.88f) ) );
	}

	TEST(TestVlogd)
	{
#ifdef __AVX__
		double x, y[4];
		x = DBL_MIN/10.;
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = 0.;
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = -1.;
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().infinity();
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = -numeric_limits<double>().infinity();
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*1418.17-708.39);
		vlog( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], log(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVlog10d)
	{
#ifdef __AVX__
		double x, y[4];
		x = DBL_MIN/10.;
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = 0.;
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = -1.;
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().infinity();
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = -numeric_limits<double>().infinity();
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*1418.17-708.39);
		vlog10( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], log10(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVlog1pd)
	{
#ifdef __AVX__
		double x, y[4];
		x = -1.;
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = -2.;
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().infinity();
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = -numeric_limits<double>().infinity();
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];
		// first test arguments close to zero
		for( int i=0; i < sz; ++i )
		{
			arg[i] = exp(ran.dbl()*258.-260.);
			if( (i%2) == 1 )
				arg[i] *= -1.;
		}
		vlog1p( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], log1p(arg[i]) ) );
		// next test full range
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*1418.17-708.39) - 0.9999999999;
		vlog1p( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], log1p(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVlogf)
	{
#ifdef __AVX__
		sys_float x, y[4];
		x = FLT_MIN/10.f;
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = 0.f;
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = -1.f;
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = -numeric_limits<sys_float>().infinity();
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vlog(y,x,x,x,x), domain_error );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*176.05-87.33);
		vlog( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], logf(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVlog10f)
	{
#ifdef __AVX__
		sys_float x, y[4];
		x = FLT_MIN/10.f;
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = 0.f;
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = -1.f;
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = -numeric_limits<sys_float>().infinity();
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vlog10(y,x,x,x,x), domain_error );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*176.05-87.33);
		vlog10( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], log10f(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVlog1pf)
	{
#ifdef __AVX__
		sys_float x, y[4];
		x = -1.f;
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = -2.f;
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = -numeric_limits<sys_float>().infinity();
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vlog1p(y,x,x,x,x), domain_error );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
		// first test arguments close to zero
		for( int i=0; i < sz; ++i )
		{
			arg[i] = exp(ran.dbl()*78.-80.);
			if( (i%2) == 1 )
				arg[i] *= -1.f;
		}
		vlog1p( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], log1pf(arg[i]) ) );
		// next test full range
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*176.05-87.33) - 0.99999;
		vlog1p( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], log1pf(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVlogdx4)
	{
		double y[4];
		vlog(y,11.,22.,33.,44.);
		CHECK( fp_equal( y[0], log(11.) ) );
		CHECK( fp_equal( y[1], log(22.) ) );
		CHECK( fp_equal( y[2], log(33.) ) );
		CHECK( fp_equal( y[3], log(44.) ) );
	}

	TEST(TestVlog10dx4)
	{
		double y[4];
		vlog10(y,11.,22.,33.,44.);
		CHECK( fp_equal( y[0], log10(11.) ) );
		CHECK( fp_equal( y[1], log10(22.) ) );
		CHECK( fp_equal( y[2], log10(33.) ) );
		CHECK( fp_equal( y[3], log10(44.) ) );
	}

	TEST(TestVlog1pdx4)
	{
		double y[4];
		vlog1p(y,11.,22.,33.,44.);
		CHECK( fp_equal( y[0], log1p(11.) ) );
		CHECK( fp_equal( y[1], log1p(22.) ) );
		CHECK( fp_equal( y[2], log1p(33.) ) );
		CHECK( fp_equal( y[3], log1p(44.) ) );
	}

	TEST(TestVlogdx8)
	{
		double y[8];
		vlog(y,11.,22.,33.,44.,55.,66.,77.,88.);
		CHECK( fp_equal( y[0], log(11.) ) );
		CHECK( fp_equal( y[1], log(22.) ) );
		CHECK( fp_equal( y[2], log(33.) ) );
		CHECK( fp_equal( y[3], log(44.) ) );
		CHECK( fp_equal( y[4], log(55.) ) );
		CHECK( fp_equal( y[5], log(66.) ) );
		CHECK( fp_equal( y[6], log(77.) ) );
		CHECK( fp_equal( y[7], log(88.) ) );
	}

	TEST(TestVlog10dx8)
	{
		double y[8];
		vlog10(y,11.,22.,33.,44.,55.,66.,77.,88.);
		CHECK( fp_equal( y[0], log10(11.) ) );
		CHECK( fp_equal( y[1], log10(22.) ) );
		CHECK( fp_equal( y[2], log10(33.) ) );
		CHECK( fp_equal( y[3], log10(44.) ) );
		CHECK( fp_equal( y[4], log10(55.) ) );
		CHECK( fp_equal( y[5], log10(66.) ) );
		CHECK( fp_equal( y[6], log10(77.) ) );
		CHECK( fp_equal( y[7], log10(88.) ) );
	}

	TEST(TestVlog1pdx8)
	{
		double y[8];
		vlog1p(y,11.,22.,33.,44.,55.,66.,77.,88.);
		CHECK( fp_equal( y[0], log1p(11.) ) );
		CHECK( fp_equal( y[1], log1p(22.) ) );
		CHECK( fp_equal( y[2], log1p(33.) ) );
		CHECK( fp_equal( y[3], log1p(44.) ) );
		CHECK( fp_equal( y[4], log1p(55.) ) );
		CHECK( fp_equal( y[5], log1p(66.) ) );
		CHECK( fp_equal( y[6], log1p(77.) ) );
		CHECK( fp_equal( y[7], log1p(88.) ) );
	}

	TEST(TestVlogfx4)
	{
		sys_float y[4];
		vlog(y,11.f,22.f,33.f,44.f);
		CHECK( fp_equal( y[0], logf(11.f) ) );
		CHECK( fp_equal( y[1], logf(22.f) ) );
		CHECK( fp_equal( y[2], logf(33.f) ) );
		CHECK( fp_equal( y[3], logf(44.f) ) );
	}

	TEST(TestVlog10fx4)
	{
		sys_float y[4];
		vlog10(y,11.f,22.f,33.f,44.f);
		CHECK( fp_equal( y[0], log10f(11.f) ) );
		CHECK( fp_equal( y[1], log10f(22.f) ) );
		CHECK( fp_equal( y[2], log10f(33.f) ) );
		CHECK( fp_equal( y[3], log10f(44.f) ) );
	}

	TEST(TestVlog1pfx4)
	{
		sys_float y[4];
		vlog1p(y,11.f,22.f,33.f,44.f);
		CHECK( fp_equal( y[0], log1pf(11.f) ) );
		CHECK( fp_equal( y[1], log1pf(22.f) ) );
		CHECK( fp_equal( y[2], log1pf(33.f) ) );
		CHECK( fp_equal( y[3], log1pf(44.f) ) );
	}

	TEST(TestVlogfx8)
	{
		sys_float y[8];
		vlog(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f);
		CHECK( fp_equal( y[0], logf(11.f) ) );
		CHECK( fp_equal( y[1], logf(22.f) ) );
		CHECK( fp_equal( y[2], logf(33.f) ) );
		CHECK( fp_equal( y[3], logf(44.f) ) );
		CHECK( fp_equal( y[4], logf(55.f) ) );
		CHECK( fp_equal( y[5], logf(66.f) ) );
		CHECK( fp_equal( y[6], logf(77.f) ) );
		CHECK( fp_equal( y[7], logf(88.f) ) );
	}

	TEST(TestVlog10fx8)
	{
		sys_float y[8];
		vlog10(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f);
		CHECK( fp_equal( y[0], log10f(11.f) ) );
		CHECK( fp_equal( y[1], log10f(22.f) ) );
		CHECK( fp_equal( y[2], log10f(33.f) ) );
		CHECK( fp_equal( y[3], log10f(44.f) ) );
		CHECK( fp_equal( y[4], log10f(55.f) ) );
		CHECK( fp_equal( y[5], log10f(66.f) ) );
		CHECK( fp_equal( y[6], log10f(77.f) ) );
		CHECK( fp_equal( y[7], log10f(88.f) ) );
	}

	TEST(TestVlog1pfx8)
	{
		sys_float y[8];
		vlog1p(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f);
		CHECK( fp_equal( y[0], log1pf(11.f) ) );
		CHECK( fp_equal( y[1], log1pf(22.f) ) );
		CHECK( fp_equal( y[2], log1pf(33.f) ) );
		CHECK( fp_equal( y[3], log1pf(44.f) ) );
		CHECK( fp_equal( y[4], log1pf(55.f) ) );
		CHECK( fp_equal( y[5], log1pf(66.f) ) );
		CHECK( fp_equal( y[6], log1pf(77.f) ) );
		CHECK( fp_equal( y[7], log1pf(88.f) ) );
	}

	TEST(TestVlogfx16)
	{
		sys_float y[16];
		vlog(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f,111.f,122.f,133.f,144.f,155.f,166.f,177.f,188.f);
		CHECK( fp_equal( y[0], logf(11.f) ) );
		CHECK( fp_equal( y[1], logf(22.f) ) );
		CHECK( fp_equal( y[2], logf(33.f) ) );
		CHECK( fp_equal( y[3], logf(44.f) ) );
		CHECK( fp_equal( y[4], logf(55.f) ) );
		CHECK( fp_equal( y[5], logf(66.f) ) );
		CHECK( fp_equal( y[6], logf(77.f) ) );
		CHECK( fp_equal( y[7], logf(88.f) ) );
		CHECK( fp_equal( y[8], logf(111.f) ) );
		CHECK( fp_equal( y[9], logf(122.f) ) );
		CHECK( fp_equal( y[10], logf(133.f) ) );
		CHECK( fp_equal( y[11], logf(144.f) ) );
		CHECK( fp_equal( y[12], logf(155.f) ) );
		CHECK( fp_equal( y[13], logf(166.f) ) );
		CHECK( fp_equal( y[14], logf(177.f) ) );
		CHECK( fp_equal( y[15], logf(188.f) ) );
	}

	TEST(TestVlog10fx16)
	{
		sys_float y[16];
		vlog10(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f,111.f,122.f,133.f,144.f,155.f,166.f,177.f,188.f);
		CHECK( fp_equal( y[0], log10f(11.f) ) );
		CHECK( fp_equal( y[1], log10f(22.f) ) );
		CHECK( fp_equal( y[2], log10f(33.f) ) );
		CHECK( fp_equal( y[3], log10f(44.f) ) );
		CHECK( fp_equal( y[4], log10f(55.f) ) );
		CHECK( fp_equal( y[5], log10f(66.f) ) );
		CHECK( fp_equal( y[6], log10f(77.f) ) );
		CHECK( fp_equal( y[7], log10f(88.f) ) );
		CHECK( fp_equal( y[8], log10f(111.f) ) );
		CHECK( fp_equal( y[9], log10f(122.f) ) );
		CHECK( fp_equal( y[10], log10f(133.f) ) );
		CHECK( fp_equal( y[11], log10f(144.f) ) );
		CHECK( fp_equal( y[12], log10f(155.f) ) );
		CHECK( fp_equal( y[13], log10f(166.f) ) );
		CHECK( fp_equal( y[14], log10f(177.f) ) );
		CHECK( fp_equal( y[15], log10f(188.f) ) );
	}

	TEST(TestVlog1pfx16)
	{
		sys_float y[16];
		vlog1p(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f,111.f,122.f,133.f,144.f,155.f,166.f,177.f,188.f);
		CHECK( fp_equal( y[0], log1pf(11.f) ) );
		CHECK( fp_equal( y[1], log1pf(22.f) ) );
		CHECK( fp_equal( y[2], log1pf(33.f) ) );
		CHECK( fp_equal( y[3], log1pf(44.f) ) );
		CHECK( fp_equal( y[4], log1pf(55.f) ) );
		CHECK( fp_equal( y[5], log1pf(66.f) ) );
		CHECK( fp_equal( y[6], log1pf(77.f) ) );
		CHECK( fp_equal( y[7], log1pf(88.f) ) );
		CHECK( fp_equal( y[8], log1pf(111.f) ) );
		CHECK( fp_equal( y[9], log1pf(122.f) ) );
		CHECK( fp_equal( y[10], log1pf(133.f) ) );
		CHECK( fp_equal( y[11], log1pf(144.f) ) );
		CHECK( fp_equal( y[12], log1pf(155.f) ) );
		CHECK( fp_equal( y[13], log1pf(166.f) ) );
		CHECK( fp_equal( y[14], log1pf(177.f) ) );
		CHECK( fp_equal( y[15], log1pf(188.f) ) );
	}

	TEST(TestVsqrtd)
	{
#ifdef __AVX__
		double x, y[4];
		x = -1.;
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = -numeric_limits<double>().infinity();
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().infinity();
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = 0.;
		vsqrt(y,x,x,x,x);
		CHECK( y[0] == 0. );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*1418.17-708.39);
		vsqrt( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], sqrt(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVhypotd)
	{
		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg1[sz];
		ALIGNED(CD_ALIGN) double arg2[sz];
		ALIGNED(CD_ALIGN) double val[sz];
#ifdef __AVX__
		double y[4], x = -numeric_limits<double>().infinity();
		CHECK_THROW( vhypot(y,x,1.,x,1.,x,1.,x,1.), domain_error );
		CHECK_THROW( vhypot(y,1.,x,1.,x,1.,x,1.,x), domain_error );
		x = numeric_limits<double>().infinity();
		CHECK_THROW( vhypot(y,x,1.,x,1.,x,1.,x,1.), domain_error );
		CHECK_THROW( vhypot(y,1.,x,1.,x,1.,x,1.,x), domain_error );
		x = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vhypot(y,x,1.,x,1.,x,1.,x,1.), domain_error );
		CHECK_THROW( vhypot(y,1.,x,1.,x,1.,x,1.,x), domain_error );
		vhypot(y,0.,0.,1.,0.,0.,1.,DBL_MAX/2.,DBL_MAX/2.);
		CHECK( fp_equal( y[0], 0. ) );
		CHECK( fp_equal( y[1], 1. ) );
		CHECK( fp_equal( y[2], 1. ) );
		CHECK( fp_equal( y[3], DBL_MAX/sqrt(2.) ) );
		// now check results over the entire range
		for( int i=0; i < sz; ++i )
		{
			arg1[i] = exp(ran.dbl()*1417.82-708.39);
			if( (i%4) < 2 )
				arg1[i] = -arg1[i];
			arg2[i] = exp(ran.dbl()*1417.82-708.39);
			if( (i%2) == 1 )
				arg2[i] = -arg2[i];
		}
		vhypot( arg1, arg2, val, 0, sz );
		for( int i=0; i < sz; ++i )
		{
			CHECK( fp_equal( val[i], hypot(arg1[i], arg2[i]) ) );
		}
#endif
		for( int i=0; i < 32; ++i )
		{
			arg1[i] = double(i);
			arg2[i] = double(i);
		}
		// finally check that non-aligned arrays are treated correctly
		for( int i=0; i < 32; ++i )
		{
			for( int j=i+1; j <= 32; ++j )
			{
				invalidate_array( val, 32*sizeof(double) );
				vhypot( arg1, arg2, val, i, j );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], double(k)*sqrt(2.) ) );
				}
				invalidate_array( val, 32*sizeof(double) );
				vhypot( &arg1[i], &arg2[i], val, 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k >= j-i )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], double(k+i)*sqrt(2.) ) );
				}
				invalidate_array( val, 32*sizeof(double) );
				vhypot( &arg1[i], arg2, &val[i], 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], sqrt(double((k-i)*(k-i)+k*k)) ) );
				}
				invalidate_array( val, 32*sizeof(double) );
				vhypot( arg1, &arg2[i], &val[i], 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnan(val[k]) );
					else
						CHECK( fp_equal( val[k], sqrt(double((k-i)*(k-i)+k*k)) ) );
				}
				invalidate_array( val, 32*sizeof(double) );
				vhypot( &arg1[i], &arg2[i/2], val, 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k >= j-i )
						CHECK( isnan(val[k]) );
					else
					{
						double ref2 = double((k+i)*(k+i) + (k+i/2)*(k+i/2));
						CHECK( fp_equal( val[k], sqrt(ref2) ) );
					}
				}
			}
		}
	}

	TEST(TestVsqrtf)
	{
#ifdef __AVX__
		sys_float x, y[4];
		x = -1.f;
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = -numeric_limits<sys_float>().infinity();
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vsqrt(y,x,x,x,x), domain_error );
		x = 0.f;
		vsqrt(y,x,x,x,x);
		CHECK( y[0] == 0.f );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*176.05-87.33);
		vsqrt( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], sqrtf(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVhypotf)
	{
		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg1[sz];
		ALIGNED(CD_ALIGN) sys_float arg2[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
#ifdef __AVX__
		sys_float y[4], x = -numeric_limits<sys_float>().infinity();
		CHECK_THROW( vhypot(y,x,1.f,x,1.f,x,1.f,x,1.f), domain_error );
		CHECK_THROW( vhypot(y,1.f,x,1.f,x,1.f,x,1.f,x), domain_error );
		x = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vhypot(y,x,1.f,x,1.f,x,1.f,x,1.f), domain_error );
		CHECK_THROW( vhypot(y,1.f,x,1.f,x,1.f,x,1.f,x), domain_error );
		x = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vhypot(y,x,1.f,x,1.f,x,1.f,x,1.f), domain_error );
		CHECK_THROW( vhypot(y,1.f,x,1.f,x,1.f,x,1.f,x), domain_error );
		vhypot(y,0.f,0.f,1.f,0.f,0.f,1.f,FLT_MAX/2.f,FLT_MAX/2.f);
		CHECK( fp_equal( y[0], 0.f ) );
		CHECK( fp_equal( y[1], 1.f ) );
		CHECK( fp_equal( y[2], 1.f ) );
		CHECK( fp_equal( y[3], FLT_MAX/sqrtf(2.f) ) );
		// now check results over the entire range
		for( int i=0; i < sz; ++i )
		{
			arg1[i] = exp(ran.dbl()*175.70-87.33);
			if( (i%4) < 2 )
				arg1[i] = -arg1[i];
			arg2[i] = exp(ran.dbl()*175.70-87.33);
			if( (i%2) == 1 )
				arg2[i] = -arg2[i];
		}
		vhypot( arg1, arg2, val, 0, sz );
		for( int i=0; i < sz; ++i )
		{
			CHECK( fp_equal( val[i], hypotf(arg1[i], arg2[i]) ) );
		}
#endif
		for( int i=0; i < 32; ++i )
		{
			arg1[i] = sys_float(i);
			arg2[i] = sys_float(i);
		}
		// finally check that non-aligned arrays are treated correctly
		for( int i=0; i < 32; ++i )
		{
			for( int j=i+1; j <= 32; ++j )
			{
				invalidate_array( val, 32*sizeof(sys_float) );
				vhypot( arg1, arg2, val, i, j );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnanf(val[k]) );
					else
						CHECK( fp_equal( val[k], sys_float(k)*sqrtf(2.f) ) );
				}
				invalidate_array( val, 32*sizeof(sys_float) );
				vhypot( &arg1[i], &arg2[i], val, 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k >= j-i )
						CHECK( isnanf(val[k]) );
					else
						CHECK( fp_equal( val[k], sys_float(k+i)*sqrtf(2.f) ) );
				}
				invalidate_array( val, 32*sizeof(sys_float) );
				vhypot( &arg1[i], arg2, &val[i], 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnanf(val[k]) );
					else
						CHECK( fp_equal( val[k], sqrtf(sys_float((k-i)*(k-i)+k*k)) ) );
				}
				invalidate_array( val, 32*sizeof(sys_float) );
				vhypot( arg1, &arg2[i], &val[i], 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k < i || k >= j )
						CHECK( isnanf(val[k]) );
					else
						CHECK( fp_equal( val[k], sqrtf(sys_float((k-i)*(k-i)+k*k)) ) );
				}
				invalidate_array( val, 32*sizeof(sys_float) );
				vhypot( &arg1[i], &arg2[i/2], val, 0, j-i );
				for( int k=0; k < 32; ++k )
				{
					if( k >= j-i )
						CHECK( isnanf(val[k]) );
					else
					{
						sys_float ref2 = sys_float((k+i)*(k+i) + (k+i/2)*(k+i/2));
						CHECK( fp_equal( val[k], sqrtf(ref2) ) );
					}
				}
			}
		}
	}

	TEST(TestVsqrtdx4)
	{
		double y[4];
		vsqrt(y,11.,22.,33.,44.);
		CHECK( fp_equal( y[0], sqrt(11.) ) );
		CHECK( fp_equal( y[1], sqrt(22.) ) );
		CHECK( fp_equal( y[2], sqrt(33.) ) );
		CHECK( fp_equal( y[3], sqrt(44.) ) );
	}

	TEST(TestVhypotdx4)
	{
		double y[4];
		vhypot(y,11.,22.,33.,44.,55.,66.,77.,88.);
		CHECK( fp_equal( y[0], hypot(11.,22.) ) );
		CHECK( fp_equal( y[1], hypot(33.,44.) ) );
		CHECK( fp_equal( y[2], hypot(55.,66.) ) );
		CHECK( fp_equal( y[3], hypot(77.,88.) ) );
	}

	TEST(TestVsqrtdx8)
	{
		double y[8];
		vsqrt(y,11.,22.,33.,44.,55.,66.,77.,88.);
		CHECK( fp_equal( y[0], sqrt(11.) ) );
		CHECK( fp_equal( y[1], sqrt(22.) ) );
		CHECK( fp_equal( y[2], sqrt(33.) ) );
		CHECK( fp_equal( y[3], sqrt(44.) ) );
		CHECK( fp_equal( y[4], sqrt(55.) ) );
		CHECK( fp_equal( y[5], sqrt(66.) ) );
		CHECK( fp_equal( y[6], sqrt(77.) ) );
		CHECK( fp_equal( y[7], sqrt(88.) ) );
	}

	TEST(TestVsqrtfx4)
	{
		sys_float y[4];
		vsqrt(y,11.f,22.f,33.f,44.f);
		CHECK( fp_equal( y[0], sqrtf(11.f) ) );
		CHECK( fp_equal( y[1], sqrtf(22.f) ) );
		CHECK( fp_equal( y[2], sqrtf(33.f) ) );
		CHECK( fp_equal( y[3], sqrtf(44.f) ) );
	}

	TEST(TestVhypotfx4)
	{
		sys_float y[4];
		vhypot(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f);
		CHECK( fp_equal( y[0], hypotf(11.f,22.f) ) );
		CHECK( fp_equal( y[1], hypotf(33.f,44.f) ) );
		CHECK( fp_equal( y[2], hypotf(55.f,66.f) ) );
		CHECK( fp_equal( y[3], hypotf(77.f,88.f) ) );
	}

	TEST(TestVsqrtfx8)
	{
		sys_float y[8];
		vsqrt(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f);
		CHECK( fp_equal( y[0], sqrtf(11.f) ) );
		CHECK( fp_equal( y[1], sqrtf(22.f) ) );
		CHECK( fp_equal( y[2], sqrtf(33.f) ) );
		CHECK( fp_equal( y[3], sqrtf(44.f) ) );
		CHECK( fp_equal( y[4], sqrtf(55.f) ) );
		CHECK( fp_equal( y[5], sqrtf(66.f) ) );
		CHECK( fp_equal( y[6], sqrtf(77.f) ) );
		CHECK( fp_equal( y[7], sqrtf(88.f) ) );
	}

	TEST(TestVhypotfx8)
	{
		sys_float y[8];
		vhypot(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f,111.f,122.f,133.f,144.f,155.f,166.f,177.f,188.f);
		CHECK( fp_equal( y[0], hypotf(11.f,22.f) ) );
		CHECK( fp_equal( y[1], hypotf(33.f,44.f) ) );
		CHECK( fp_equal( y[2], hypotf(55.f,66.f) ) );
		CHECK( fp_equal( y[3], hypotf(77.f,88.f) ) );
		CHECK( fp_equal( y[4], hypotf(111.f,122.f) ) );
		CHECK( fp_equal( y[5], hypotf(133.f,144.f) ) );
		CHECK( fp_equal( y[6], hypotf(155.f,166.f) ) );
		CHECK( fp_equal( y[7], hypotf(177.f,188.f) ) );
	}

	TEST(TestVsqrtfx16)
	{
		sys_float y[16];
		vsqrt(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f,111.f,122.f,133.f,144.f,155.f,166.f,177.f,188.f);
		CHECK( fp_equal( y[0], sqrtf(11.f) ) );
		CHECK( fp_equal( y[1], sqrtf(22.f) ) );
		CHECK( fp_equal( y[2], sqrtf(33.f) ) );
		CHECK( fp_equal( y[3], sqrtf(44.f) ) );
		CHECK( fp_equal( y[4], sqrtf(55.f) ) );
		CHECK( fp_equal( y[5], sqrtf(66.f) ) );
		CHECK( fp_equal( y[6], sqrtf(77.f) ) );
		CHECK( fp_equal( y[7], sqrtf(88.f) ) );
		CHECK( fp_equal( y[8], sqrtf(111.f) ) );
		CHECK( fp_equal( y[9], sqrtf(122.f) ) );
		CHECK( fp_equal( y[10], sqrtf(133.f) ) );
		CHECK( fp_equal( y[11], sqrtf(144.f) ) );
		CHECK( fp_equal( y[12], sqrtf(155.f) ) );
		CHECK( fp_equal( y[13], sqrtf(166.f) ) );
		CHECK( fp_equal( y[14], sqrtf(177.f) ) );
		CHECK( fp_equal( y[15], sqrtf(188.f) ) );
	}

	TEST(TestVasinhd)
	{
#ifdef __AVX__
		double x, y[4];
		x = -numeric_limits<double>().infinity();
		CHECK_THROW( vasinh(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().infinity();
		CHECK_THROW( vasinh(y,x,x,x,x), domain_error );
		x = -numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vasinh(y,x,x,x,x), domain_error );
		x = 0.;
		vasinh(y,x,x,x,x);
		CHECK( y[0] == 0. );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];
		for( int i=0; i < sz; ++i )
		{
			arg[i] = exp(ran.dbl()*1418.17-708.39);
			if( (i%2) == 1 )
				arg[i] = -arg[i];
		}
		vasinh( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], asinh(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVasinhdFast)
	{
#ifdef __AVX__
		double x, y[4];
		x = -numeric_limits<double>().infinity();
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().infinity();
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = numeric_limits<double>().quiet_NaN();
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = -1.;
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		double z = sqrt(DBL_MAX);
		x = nextafter(z,DBL_MAX);
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = nextafter(z,-DBL_MAX);
		vfast_asinh(y,x,x,x,x);
		CHECK( y[0] > 0. && y[0] < DBL_MAX );
		x = 0.;
		vfast_asinh(y,x,x,x,x);
		CHECK( y[0] == 0. );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) double arg[sz];
		ALIGNED(CD_ALIGN) double val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*1063.28-708.39);
		vfast_asinh( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], asinh(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVasinhf)
	{
#ifdef __AVX__
		sys_float x, y[4];
		x = -numeric_limits<sys_float>().infinity();
		CHECK_THROW( vasinh(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vasinh(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vasinh(y,x,x,x,x), domain_error );
		x = 0.f;
		vasinh(y,x,x,x,x);
		CHECK( y[0] == 0.f );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
		for( int i=0; i < sz; ++i )
		{
			arg[i] = exp(ran.dbl()*176.05-87.33);
			if( (i%2) == 1 )
				arg[i] = -arg[i];
		}
		vasinh( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], asinhf(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVasinhfFast)
	{
#ifdef __AVX__
		sys_float x, y[4];
		x = -numeric_limits<sys_float>().infinity();
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().infinity();
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = numeric_limits<sys_float>().quiet_NaN();
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = -1.f;
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		sys_float z = sqrt(FLT_MAX);
		x = nextafterf(z,FLT_MAX);
		CHECK_THROW( vfast_asinh(y,x,x,x,x), domain_error );
		x = nextafterf(z,-FLT_MAX);
		vfast_asinh(y,x,x,x,x);
		CHECK( y[0] > 0.f && y[0] < FLT_MAX );
		x = 0.f;
		vfast_asinh(y,x,x,x,x);
		CHECK( y[0] == 0.f );

		const int sz = 2048;
		ALIGNED(CD_ALIGN) sys_float arg[sz];
		ALIGNED(CD_ALIGN) sys_float val[sz];
		for( int i=0; i < sz; ++i )
			arg[i] = exp(ran.dbl()*131.69-87.33);
		vfast_asinh( arg, val, 0, sz );
		for( int i=0; i < sz; ++i )
			CHECK( fp_equal( val[i], asinhf(arg[i]) ) );
#else
		CHECK( fp_equal( 1., 1. ) );
#endif
	}

	TEST(TestVasinhdx4)
	{
		double y[4];
		vasinh(y,11.,22.,33.,44.);
		CHECK( fp_equal( y[0], asinh(11.) ) );
		CHECK( fp_equal( y[1], asinh(22.) ) );
		CHECK( fp_equal( y[2], asinh(33.) ) );
		CHECK( fp_equal( y[3], asinh(44.) ) );
	}

	TEST(TestVasinhdFastx4)
	{
		double y[4];
		vfast_asinh(y,11.,22.,33.,44.);
		CHECK( fp_equal( y[0], asinh(11.) ) );
		CHECK( fp_equal( y[1], asinh(22.) ) );
		CHECK( fp_equal( y[2], asinh(33.) ) );
		CHECK( fp_equal( y[3], asinh(44.) ) );
	}

	TEST(TestVasinhdx8)
	{
		double y[8];
		vasinh(y,11.,22.,33.,44.,55.,66.,77.,88.);
		CHECK( fp_equal( y[0], asinh(11.) ) );
		CHECK( fp_equal( y[1], asinh(22.) ) );
		CHECK( fp_equal( y[2], asinh(33.) ) );
		CHECK( fp_equal( y[3], asinh(44.) ) );
		CHECK( fp_equal( y[4], asinh(55.) ) );
		CHECK( fp_equal( y[5], asinh(66.) ) );
		CHECK( fp_equal( y[6], asinh(77.) ) );
		CHECK( fp_equal( y[7], asinh(88.) ) );
	}

	TEST(TestVasinhdFastx8)
	{
		double y[8];
		vfast_asinh(y,11.,22.,33.,44.,55.,66.,77.,88.);
		CHECK( fp_equal( y[0], asinh(11.) ) );
		CHECK( fp_equal( y[1], asinh(22.) ) );
		CHECK( fp_equal( y[2], asinh(33.) ) );
		CHECK( fp_equal( y[3], asinh(44.) ) );
		CHECK( fp_equal( y[4], asinh(55.) ) );
		CHECK( fp_equal( y[5], asinh(66.) ) );
		CHECK( fp_equal( y[6], asinh(77.) ) );
		CHECK( fp_equal( y[7], asinh(88.) ) );
	}

	TEST(TestVasinhfx4)
	{
		sys_float y[4];
		vasinh(y,11.f,22.f,33.f,44.f);
		CHECK( fp_equal( y[0], asinhf(11.f) ) );
		CHECK( fp_equal( y[1], asinhf(22.f) ) );
		CHECK( fp_equal( y[2], asinhf(33.f) ) );
		CHECK( fp_equal( y[3], asinhf(44.f) ) );
	}

	TEST(TestVasinhfFastx4)
	{
		sys_float y[4];
		vfast_asinh(y,11.f,22.f,33.f,44.f);
		CHECK( fp_equal( y[0], asinhf(11.f) ) );
		CHECK( fp_equal( y[1], asinhf(22.f) ) );
		CHECK( fp_equal( y[2], asinhf(33.f) ) );
		CHECK( fp_equal( y[3], asinhf(44.f) ) );
	}

	TEST(TestVasinhfx8)
	{
		sys_float y[8];
		vasinh(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f);
		CHECK( fp_equal( y[0], asinhf(11.f) ) );
		CHECK( fp_equal( y[1], asinhf(22.f) ) );
		CHECK( fp_equal( y[2], asinhf(33.f) ) );
		CHECK( fp_equal( y[3], asinhf(44.f) ) );
		CHECK( fp_equal( y[4], asinhf(55.f) ) );
		CHECK( fp_equal( y[5], asinhf(66.f) ) );
		CHECK( fp_equal( y[6], asinhf(77.f) ) );
		CHECK( fp_equal( y[7], asinhf(88.f) ) );
	}

	TEST(TestVasinhfFastx8)
	{
		sys_float y[8];
		vfast_asinh(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f);
		CHECK( fp_equal( y[0], asinhf(11.f) ) );
		CHECK( fp_equal( y[1], asinhf(22.f) ) );
		CHECK( fp_equal( y[2], asinhf(33.f) ) );
		CHECK( fp_equal( y[3], asinhf(44.f) ) );
		CHECK( fp_equal( y[4], asinhf(55.f) ) );
		CHECK( fp_equal( y[5], asinhf(66.f) ) );
		CHECK( fp_equal( y[6], asinhf(77.f) ) );
		CHECK( fp_equal( y[7], asinhf(88.f) ) );
	}

	TEST(TestVasinhfx16)
	{
		sys_float y[16];
		vasinh(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f,111.f,122.f,133.f,144.f,155.f,166.f,177.f,188.f);
		CHECK( fp_equal( y[0], asinhf(11.f) ) );
		CHECK( fp_equal( y[1], asinhf(22.f) ) );
		CHECK( fp_equal( y[2], asinhf(33.f) ) );
		CHECK( fp_equal( y[3], asinhf(44.f) ) );
		CHECK( fp_equal( y[4], asinhf(55.f) ) );
		CHECK( fp_equal( y[5], asinhf(66.f) ) );
		CHECK( fp_equal( y[6], asinhf(77.f) ) );
		CHECK( fp_equal( y[7], asinhf(88.f) ) );
		CHECK( fp_equal( y[8], asinhf(111.f) ) );
		CHECK( fp_equal( y[9], asinhf(122.f) ) );
		CHECK( fp_equal( y[10], asinhf(133.f) ) );
		CHECK( fp_equal( y[11], asinhf(144.f) ) );
		CHECK( fp_equal( y[12], asinhf(155.f) ) );
		CHECK( fp_equal( y[13], asinhf(166.f) ) );
		CHECK( fp_equal( y[14], asinhf(177.f) ) );
		CHECK( fp_equal( y[15], asinhf(188.f) ) );
	}

	TEST(TestVasinhfFastx16)
	{
		sys_float y[16];
		vfast_asinh(y,11.f,22.f,33.f,44.f,55.f,66.f,77.f,88.f,111.f,122.f,133.f,144.f,155.f,166.f,177.f,188.f);
		CHECK( fp_equal( y[0], asinhf(11.f) ) );
		CHECK( fp_equal( y[1], asinhf(22.f) ) );
		CHECK( fp_equal( y[2], asinhf(33.f) ) );
		CHECK( fp_equal( y[3], asinhf(44.f) ) );
		CHECK( fp_equal( y[4], asinhf(55.f) ) );
		CHECK( fp_equal( y[5], asinhf(66.f) ) );
		CHECK( fp_equal( y[6], asinhf(77.f) ) );
		CHECK( fp_equal( y[7], asinhf(88.f) ) );
		CHECK( fp_equal( y[8], asinhf(111.f) ) );
		CHECK( fp_equal( y[9], asinhf(122.f) ) );
		CHECK( fp_equal( y[10], asinhf(133.f) ) );
		CHECK( fp_equal( y[11], asinhf(144.f) ) );
		CHECK( fp_equal( y[12], asinhf(155.f) ) );
		CHECK( fp_equal( y[13], asinhf(166.f) ) );
		CHECK( fp_equal( y[14], asinhf(177.f) ) );
		CHECK( fp_equal( y[15], asinhf(188.f) ) );
	}

	TEST(TestVHstring)
	{
		string test;
		// vh128sum of an empty file...
		CHECK( VHstring( test ) == "30a574ab9824cb4358b310345eb60ab0" );
		CHECK( VHstring( test ).length() == 32 );
		// check that leading zeros are printed correctly
		test = "aaaaaaaaaaaaaaaaaaaaaa\n";
		CHECK( VHstring( test ) == "0ba6d0f084b09afac9413c8c59f4b2ae" );
		CHECK( VHstring( test ).length() == 32 );
	}
}
