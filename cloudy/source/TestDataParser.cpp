/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "lines.h"
#include "parser.h"

namespace {
	TEST(TestFPReadDouble)
	{
		DataParser d;
		d.setline("1000 \t 2e3 +2 -3. -3.125  2.2e-12  .3e+03");
		double x;
		d.getToken(x);
		CHECK( x == 1000. );
		d.getToken(x);
		CHECK( x == 2000. );
		d.getToken(x);
		CHECK( x == 2. );
		d.getToken(x);
		CHECK( x == -3. );
		d.getToken(x);
		CHECK( x == -3.125 );
		d.getToken(x);
		CHECK( fp_equal(x, 2.2e-12) );
		CHECK( d.getTokenOptional(x) );
		CHECK( x == 300. );
		CHECK( !d.getTokenOptional(x) );
		CHECK( x == 0. );

		// test special values
		d.setline(" -0.0e10 1.7976931348623157e+308 2.225073858507202e-308 ");
		d.getToken(x);
		CHECK( signbit(x) && x == 0. );
		d.getToken(x);
		CHECK( fp_equal(x, DBL_MAX) );
		d.getToken(x);
		CHECK( fp_equal(x, DBL_MIN) );
		d.setline("-1.7976931348623157e+308 -2.225073858507202e-308");
		d.getToken(x);
		CHECK( fp_equal(x, -DBL_MAX) );
		d.getToken(x);
		CHECK( fp_equal(x, -DBL_MIN) );
		d.setline("1.e-500 1.e500");
		d.getToken(x);
		CHECK( fp_equal(x, 0.) );
		d.getToken(x);
		CHECK( std::isinf(x) );

		// test failure modes, for this redirect output to temporary file
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline(".");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("+");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline(".e12");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1.2e+");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1.2e+12e-12");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-+1.2e+12");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-1.2e+-12");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-1..2e+12");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-1.2+3e-1+2");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("a");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1.234u");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("inf");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("nan");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		// and clean up the mess, fclose() will automatically remove the file
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestFPReadFloat)
	{
		DataParser d;
		d.setline("1000 \t 2e3 +2 -3. -3.5  2.2e-12  .3e+03");
		sys_float x;
		d.getToken(x);
		CHECK( x == 1000.f );
		d.getToken(x);
		CHECK( x == 2000.f );
		d.getToken(x);
		CHECK( x == 2.f );
		d.getToken(x);
		CHECK( x == -3.f );
		d.getToken(x);
		CHECK( x == -3.5f );
		d.getToken(x);
		CHECK( fp_equal(x, 2.2e-12f) );
		CHECK( d.getTokenOptional(x) );
		CHECK( x == 300.f );
		CHECK( !d.getTokenOptional(x) );
		CHECK( x == 0.f );

		// test special values
		d.setline(" -0.0e10 3.4028234663852886e+38 1.1754943508222875e-38");
		d.getToken(x);
		CHECK( signbit(x) && x == 0.f );
		d.getToken(x);
		CHECK( fp_equal(x, FLT_MAX) );
		d.getToken(x);
		CHECK( fp_equal(x, FLT_MIN) );
		d.setline("-3.4028234663852886e+38 -1.1754943508222875e-38");
		d.getToken(x);
		CHECK( fp_equal(x, -FLT_MAX) );
		d.getToken(x);
		CHECK( fp_equal(x, -FLT_MIN) );
	}

	TEST(TestIntReadInt64)
	{
		DataParser d;
		d.setline("1000 \t +2 -3");
		int64 x;
		d.getToken(x);
		CHECK( x == 1000LL );
		d.getToken(x);
		CHECK( x == 2LL );
		CHECK( d.getTokenOptional(x) );
		CHECK( x == -3LL );
		CHECK( !d.getTokenOptional(x) );
		CHECK( x == 0LL );

		// test special values
		d.setline("-9223372036854775808 9223372036854775807");
		d.getToken(x);
		CHECK( x == INT64_MIN );
		d.getToken(x);
		CHECK( x == INT64_MAX );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline("+");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("a");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1+2");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1234u");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestIntReadUnsignedInt64)
	{
		DataParser d;
		d.setline("1000 \t +2");
		uint64 x;
		d.getToken(x);
		CHECK( x == 1000ULL );
		CHECK( d.getTokenOptional(x) );
		CHECK( x == 2ULL );
		CHECK( !d.getTokenOptional(x) );
		CHECK( x == 0ULL );

		// test special values
		d.setline("18446744073709551615");
		d.getToken(x);
		CHECK( x == UINT64_MAX );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline("+");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-2");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("a");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1+2");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1234u");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestIntReadInt32)
	{
		DataParser d;
		d.setline("1000 \t +2 -3");
		int32 x;
		d.getToken(x);
		CHECK( x == 1000 );
		d.getToken(x);
		CHECK( x == 2 );
		CHECK( d.getTokenOptional(x) );
		CHECK( x == -3 );
		CHECK( !d.getTokenOptional(x) );
		CHECK( x == 0 );

		// test special values
		d.setline("-2147483648 2147483647");
		d.getToken(x);
		CHECK( x == INT32_MIN );
		d.getToken(x);
		CHECK( x == INT32_MAX );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline("+");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("a");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1+2");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1234u");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestIntReadUnsignedInt32)
	{
		DataParser d;
		d.setline("1000 \t +2");
		uint32 x;
		d.getToken(x);
		CHECK( x == 1000U );
		CHECK( d.getTokenOptional(x) );
		CHECK( x == 2U );
		CHECK( !d.getTokenOptional(x) );
		CHECK( x == 0U );

		// test special values
		d.setline("4294967295");
		d.getToken(x);
		CHECK( x == UINT32_MAX );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline("+");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("-2");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("a");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1+2");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		d.setline("1234u");
		CHECK_THROW( d.getToken(x), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestReadChar)
	{
		DataParser d;
		d.setline("  a \t b  7 ");
		char c;
		d.getToken(c);
		CHECK( c == 'a' );
		d.getToken(c);
		CHECK( c == 'b' );
		d.getToken(c);
		CHECK( c == '7' );
		CHECK( !d.getTokenOptional(c) );
	}


	TEST(TestReadWord)
	{
		DataParser d;
		d.setline("some fancy   text");
		string s;
		d.getToken(s);
		CHECK( s == "some" );
		d.getToken(s);
		CHECK( s == "fancy" );
		CHECK( d.getTokenOptional(s) );
		CHECK( s == "text" );
		CHECK( !d.getTokenOptional(s) );
		CHECK( s == "" );
	}

	TEST(TestReadQuoted)
	{
		DataParser d;
		d.setline(" \"some fancy\" text");
		string s;
		d.getQuote(s);
		CHECK( s == "some fancy" );
		CHECK( !d.getQuoteOptional(s) );
		CHECK( s == "" );
		d.setline("\"some fancy\"");
		d.getQuote(s);
		CHECK( s == "some fancy" );
		d.setline("1.23");
		double x;
		d.getToken(x);
		CHECK( fp_equal(x, 1.23) );
		CHECK( !d.getQuoteOptional(s) );
		CHECK( s == "" );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline("some fancy text");
		CHECK_THROW( d.getQuote(s), cloudy_exit );
		d.setline(" \"some fancy text" );
		CHECK_THROW( d.getQuote(s), cloudy_exit );
		d.setline(" \"some fancy\"\"text\"" );
		CHECK_THROW( d.getQuote(s), cloudy_exit );
		d.setline(" \"some fancy\"text" );
		CHECK_THROW( d.getQuote(s), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestReadLineID)
	{
		DataParser d;
		LineID line;
		d.setline("H  1  1216");
		d.getLineID(line);
		CHECK( line.chLabel == "H  1" && line.wave == 1216_r && line.indLo < 0 && line.indHi < 0 );
		d.setline("H  1 1216A");
		d.getLineID(line);
		CHECK( line.chLabel == "H  1" && line.wave == 1216_r && line.ELo < 0_r );
		d.setline("\"Fe 2b\"  12.00m");
		d.getLineID(line);
		CHECK( line.chLabel == "Fe 2b" && line.wave == 12.00e4_r );
		d.setline("\"Fe 2b  \"  12.00c # comment");
		d.getLineID(line);
		CHECK( line.chLabel == "Fe 2b" && line.wave == 12.00e8_r );
		d.setline("CO  12.00C # comment");
		d.getLineID(line);
		CHECK( line.chLabel == "CO" && line.wave == 12.00e8_r );
		d.setline("Al 2 1670. index=1,5");
		d.getLineID(line);
		CHECK( line.chLabel == "Al 2" && line.wave== 1670_r && line.indLo == 1 && line.indHi == 5 && line.ELo < 0_r );
		d.setline("Al 2 1670. Elow=1");
		d.getLineID(line);
		CHECK( line.chLabel == "Al 2" && line.wave == 1670_r && line.indLo < 0 && line.indHi < 0 && line.ELo == 1_r );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		// not really a failure, but a warning
		d.setline("H 1   1216a");
		d.getLineID(line);
		CHECK( line.chLabel == "H  1" && line.wave == 1216._r );
		d.setline("\"H 1\"   1216M");
		d.getLineID(line);
		CHECK( line.chLabel == "H  1" && line.wave == 1216.e4_r );
		// the rest are all real errors
		d.setline("monitor line \"H  1\" 1216");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("H 1");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("H  1");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("TOTL12.00m");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2b 12.00m");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("\"Fe 2b  12.00m");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  comment  12.00m");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00m  index");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00m  index=3");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00m  index=-1,2");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00m  index=3,2");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00m  elow");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00m  elow=-1");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00m  keyword");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		d.setline("Fe 2  12.00 m");
		CHECK_THROW( d.getLineID(line), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestReadArray)
	{
		DataParser d;
		d.setline("1 2 3 4 5 6");
		long x[10];
		d.getToken(x, 6);
		for( long i=0; i < 6; ++i )
			CHECK( x[i] == i+1 );
		d.setline("1 2 3 4 5 6");
		CHECK( d.getTokenOptional(x, 10) == 6 );
		for( long i=0; i < 6; ++i )
			CHECK( x[i] == i+1 );
		d.setline("   \t ");
		CHECK( d.getTokenOptional(x, 10) == 0 );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline("1 2 3 4 5 6");
		CHECK_THROW( d.getToken(x, 7), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestReadMagic)
	{
		DataParser d;
		d.setline("20170812   ");
		CHECK( !d.lgEOL() );
		d.checkMagic(20170812);
		d.checkEOL();
		CHECK( d.lgEOL() );
		d.setline("   2017 812 ");
		d.checkMagic(2017, 812);
		d.checkEOL();
		CHECK( d.lgEOL() );
		d.setline("2017 8 12");
		d.checkMagic(2017, 8, 12);
		d.checkEOL();
		CHECK( d.lgEOL() );
		d.setline("2017\t 8 \t12\t53");
		d.checkMagic(2017, 8, 12, 53);
		d.checkEOL();
		CHECK( d.lgEOL() );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.setline("20170812");
		CHECK_THROW( d.checkMagic(20170811), cloudy_exit );
		d.setline("2017 812");
		CHECK_THROW( d.checkMagic(2016, 812), cloudy_exit );
		d.setline("2017 812");
		CHECK_THROW( d.checkMagic(2017, 811), cloudy_exit );
		d.setline("2017 8 12");
		CHECK_THROW( d.checkMagic(2016, 8, 12), cloudy_exit );
		d.setline("2017 8 12");
		CHECK_THROW( d.checkMagic(2017, 7, 12), cloudy_exit );
		d.setline("2017 8 12");
		CHECK_THROW( d.checkMagic(2017, 8, 11), cloudy_exit );
		d.setline("2017 8\t 12 53");
		CHECK_THROW( d.checkMagic(2016, 8, 12, 53), cloudy_exit );
		d.setline("2017 \t8 12 53");
		CHECK_THROW( d.checkMagic(2017, 7, 12, 53), cloudy_exit );
		d.setline("2017 8 12 \t 53");
		CHECK_THROW( d.checkMagic(2017, 8, 11, 53), cloudy_exit );
		d.setline("2017 8\t12 \t 53");
		CHECK_THROW( d.checkMagic(2017, 8, 12, 52), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestSkip)
	{
		DataParser d;
		d.setline("123 456 789");
		long x;
		d.getToken(x);
		CHECK( x == 123 );
		d.skipTo(7);
		d.getToken(x);
		CHECK( x == 789 );

		// test failure mode
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		CHECK_THROW( d.skipTo(4), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;

		d.setline("123 Jlo 1 Jhi 2 Jlo 4 Jhi 9");
		d.getToken(x);
		CHECK( x == 123 );
		d.skipAfter("Jhi");
		d.getToken(x);
		CHECK( x == 2 );
		d.skipAfter("Jhi ");
		d.getToken(x);
		CHECK( x == 9 );

		// test failure mode
		bak = ioQQQ;
		tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		CHECK_THROW( d.skipAfter("Jhi"), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestWarning)
	{
		// test warning -- there is nothing we can check in the macros
		// but valgrind may churn something up when this is executed
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		DataParser d;
		d.setline("123 456 789");
		long x;
		d.getToken(x);
		if( x > 100 )
			d.warning( "value may be too large" );
		d.getToken(x);
		CHECK( x == 456 );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestFileAccess)
	{
		char name[] = "temp.XXXXXX";
		int fd = mkstemp(name);
		FILE *io = fdopen(fd, "w");
		fprintf(io, "# comment line\n");
		fprintf(io, "2.567 # 1.234 comment\n");
		fprintf(io, "*********\n");
		fprintf(io, " *********\n");
		fprintf(io, "\n");
		fprintf(io, "  # comment\n");
		fprintf(io, " \t \n");
		fclose(io);
		DataParser d(name, ES_NONE);
		CHECK( d.isOpen() );
		CHECK( d.getline() );
		double x;
		d.getToken(x);
		CHECK( fp_equal(x, 2.567) );
		CHECK( !d.getTokenOptional(x) );
		d.checkEOL();
		CHECK( d.getline() );
		CHECK( !d.lgEODMarker() );
		CHECK( d.getline() );
		CHECK( !d.lgEODMarker() );
		CHECK( !d.getline() );
		CHECK( d.lgEOF() );
		d.rewind();
		CHECK( d.getline() );
		d.getToken(x);
		CHECK( fp_equal(x, 2.567) );
		d.open(name, ES_STARS_ONLY);
		CHECK( d.getline() );
		d.getToken(x);
		CHECK( fp_equal(x, 2.567) );
		CHECK( d.getline() );
		CHECK( d.lgEODMarker() );
		CHECK( d.getline() );
		CHECK( !d.lgEODMarker() );
		CHECK( !d.getline() );
		d.open(name, ES_STARS_AND_BLANKS);
		CHECK( d.getline() );
		d.getToken(x);
		CHECK( fp_equal(x, 2.567) );
		CHECK( d.getline() );
		CHECK( d.lgEODMarker() );
		CHECK( d.getline() );
		CHECK( !d.lgEODMarker() );
		CHECK( d.getline() );
		CHECK( d.lgEODMarker() );
		CHECK( d.getline() );
		CHECK( d.lgEODMarker() );
		CHECK( d.getline() );
		CHECK( d.lgEODMarker() );
		CHECK( !d.getline() );

		// test failure modes
		d.open("nonexistent.aG5Qdp", ES_NONE, AS_TRY);
		CHECK( !d.isOpen() );

		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		d.open(name, ES_NONE);
		CHECK( d.getline() );
		CHECK_THROW( d.checkEOD(), cloudy_exit );
		d.open(name, ES_STARS_ONLY);
		CHECK( d.getline() );
		// check that d.checkEOD() does not throw an exception
		// since it does not return any value, we need this funny syntax...
		CHECK( (d.checkEOD(), true) );
		d.open(name, ES_INVALID);
		CHECK_THROW( d.getline(), cloudy_exit );
		d.close();
		remove(name);
		d.setline( "****" );
		CHECK_THROW( d.lgEODMarker(), cloudy_exit );
		CHECK_THROW( d.open("nonexistent.aG5Qdp", ES_NONE), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}

	TEST(TestFileAccessDOS)
	{
		char name[] = "temp.XXXXXX";
		int fd = mkstemp(name);
		FILE *io = fdopen(fd, "w");
		// deliberately put \r in the wrong place so that the final number is stripped
		// this way we can easily see the effect of stripping ^M from files with DOS EOL
		fprintf(io, "2.567 \r 1.234\n");
		// deliberately omit final EOL marker, this happens from time to time in data files...
		fprintf(io, "1.4778");
		fclose(io);
		DataParser d(name, ES_NONE);
		CHECK( d.getline() );
		double x;
		d.getToken(x);
		CHECK( fp_equal(x, 2.567) );
		CHECK( !d.getTokenOptional(x) );
		d.checkEOL();
		CHECK( d.getline() );
		d.getToken(x);
		CHECK( fp_equal(x, 1.4778) );
		CHECK( !d.getline() );
		CHECK( d.lgEOF() );
		d.close();
		remove(name);
	}
}
