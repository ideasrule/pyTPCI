#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "ran.h"
#include "service.h"

namespace {

	TEST(TestFFmtRead)
	{
		bool lgEol;
		char buf[128];
		for( int i=0; i < 4096; ++i )
		{
			double x = exp10(ran.dbl()*600. - 300.);
			if( ran.u8()&1 )
				x = -x;
			sprintf( buf, " %.16e", x );
			long j = 1;
			double y = FFmtRead( buf, &j, 128, &lgEol );
			CHECK( fp_equal( x, y, 2 ) );
			CHECK( !lgEol );
		}
		// test edge cases
		long j = 1;
		double x = FFmtRead( "HYDROGEN\t1", &j, 10, &lgEol );
		CHECK( !lgEol && x == 1. );
		j = 10;
		x = FFmtRead( "HYDROGEN\t1", &j, 10, &lgEol );
		CHECK( !lgEol && x == 1. );
		// these are mainly intended to check agains buffer overruns
		j = 1;
		x = FFmtRead( "+", &j, 1, &lgEol );
		CHECK( lgEol && x == 0. );
		j = 1;
		x = FFmtRead( "-.", &j, 2, &lgEol );
		CHECK( lgEol && x == 0. );
		j = 1;
		x = FFmtRead( ".", &j, 1, &lgEol );
		CHECK( lgEol && x == 0. );
	}

	TEST(TestPowi)
	{
		for( int i=0; i < 2048; ++i )
		{
			double arg1 = exp10(ran.dbl()*4. - 2.);
			if( ran.u8()&1 )
				arg1 = -arg1;
			long arg2 = (ran.i15()%200) - 100;
			CHECK( fp_equal( powi(arg1,arg2), pow(arg1,(double)arg2), max(abs(arg2),3) ) );
		}
	}

	TEST(TestPowpq)
	{
		for( int i=0; i < 1024; ++i )
		{
			double arg1 = ran.dbl()*1.e50;
			int arg2 = (ran.i7()%9) + 1;
			int arg3 = (ran.i7()%8) + 2;
			// powpq() may be more accurate than pow(), e.g. for powpq(x,1,3)...
			CHECK( fp_equal( powpq(arg1,arg2,arg3), pow(arg1,(double)arg2/(double)arg3), 128 ) );
		}
		for( int i=0; i < 1024; ++i )
		{
			double arg1 = ran.dbl()*1.e50;
			int arg2 = -((ran.i7()%9) + 1);
			int arg3 = (ran.i7()%8) + 2;
			CHECK( fp_equal( powpq(arg1,arg2,arg3), pow(arg1,(double)arg2/(double)arg3), 128 ) );
		}
	}
}
