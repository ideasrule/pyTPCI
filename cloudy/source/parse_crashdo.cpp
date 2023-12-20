/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCrashDo any of several tests to check that the code can crash */
#include "cddefines.h"
#include "parser.h"
#include "grid.h"
#include "container_classes.h"
#include "vectorize.h"

#ifdef _MSC_VER
	/* disable warning about undefined vars being used - one of the tests shall do exactly that */
#	pragma warning( disable : 4700 )
	/* disable warning about division by zero */
#	pragma warning( disable : 4756 )
#endif

#ifdef __INTEL_COMPILER
#	pragma warning( disable : 592 )
#endif

#ifdef __clang__
#	pragma clang diagnostic ignored "-Wuninitialized"
#endif

#ifdef __GNUC_EXCL__
#	pragma GCC diagnostic ignored "-Wuninitialized"
#	if ( __GNUC__ > 4 ) || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 7 )
#	pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#	endif
#endif

#if __SUNPRO_CC >= 20800
#	pragma error_messages (off,SEC_UNINITIALIZED_MEM_READ)
#endif

/* this is size of array used in array bounds exceeded crash test */
const int ARR_SIZE = 10;

/* static variable used in undefined and bounds tests */
static double ar2[ARR_SIZE];

// force optimization off; any level of optimization will kill the 
// functionality of this routine
#if defined(_MSC_VER) || defined(__ICC)
#pragma optimize("", off)
#elif defined(__PGI)
#pragma global opt=0
#endif

/*ParseCrashDo any of several tests to check that the code can crash */
void ParseCrashDo(Parser &p)
{
	double ar1, br1;
	bool lgCrash = false;

	DEBUG_ENTRY( "ParseCrashDo()" );

	/* div by 0 to get crash as check on FP environment */
	if( p.nMatch("ZERO") )
	{
		fprintf(ioQQQ," I will now div by 0 to get crash.  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then"
			" there are problems.\n");
		fflush(ioQQQ);
		ar1 = 1. / ZeroNum;
		fprintf(ioQQQ," I am still alive - something is wrong, result is %e\n",
			ar1);
		lgCrash = true;
	}

	/* use an undefined number */
	else if( p.nMatch("UNDE") )
	{
		double A_variable_which_SHOULD_be_used_uninitialized;
		fprintf(ioQQQ," Now I will now use an undefined variable off the stack.  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then"
				" there are problems.\n");
		fflush(ioQQQ);
		/*lint -e530 a not initialized */
		A_variable_which_SHOULD_be_used_uninitialized *= 1e-10f;
		/*lint +e530 a not initialized */

		fprintf(ioQQQ," I am still alive - something is wrong, the result of the multiplication of"
				" undefined by 1e-10 is %e\n", A_variable_which_SHOULD_be_used_uninitialized );
		fflush(ioQQQ);
		lgCrash = true;
	}

	/* make overflow to get crash as check on FP environment */
	else if( p.nMatch("OVER") && p.nMatch("LONG") )
	{ 
		long lng;
		fprintf(ioQQQ," I will now make long overflow to get crash.  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then"
			" there are problems.\n");
		fflush(ioQQQ);
		lng = (long)((double)LONG_MAX*sqrt(1e6));
		fprintf(ioQQQ," I am still alive - something is wrong, the result was %li\n",
			lng);
		lgCrash = true;
	}

	/* make overflow to get crash as check on FP environment */
	else if( p.nMatch("OVER") )
	{ 
		ar1 = 1e-20; 
		fprintf(ioQQQ," I will now make floating point overflow to get crash.  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then"
			" there are problems.\n");
		fflush(ioQQQ);
		br1 = DBL_MAX / ar1;
		fprintf(ioQQQ," I am still alive - something is wrong, the result was %e\n",
			br1);
		lgCrash = true;
	}

	/* dereference NULL pointer to get a segmentation fault */
	else if( p.nMatch("SEGF") )
	{
		double* p = (double*)ZeroPtr;
		fprintf(ioQQQ," I will now dereference a NULL pointer.  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then there are problems.\n");
		fflush(ioQQQ);
		br1 = *p;
		fprintf(ioQQQ," I am still alive - something is wrong, the result was %e\n",
			br1);
		lgCrash = true;
	}

	/* assert false test to get crash as check on environment */
	else if( p.nMatch("ASSE") )
	{ 
		fprintf(ioQQQ," I will now assert that a false statement is true to get a crash.\n\n");
		fprintf(ioQQQ," The correct behavior is for the statement \"PROBLEM DISASTER An assert"
			" has been thrown, this is bad\" to be printed, followed by lots more scary"
			" looking messages.\n\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - the assert macro is not working ....\" then"
			" there are problems.\n\n");
		fflush(ioQQQ);
		ASSERT( DBL_MAX <  ZeroNum );
		fprintf(ioQQQ," I am still alive - the assert macro is not working in this executable.\n");
		lgCrash = true;
	}

	/* assert ratios of zeros (NaN) to get crash as check on environment */
	else if( p.nMatch(" NAN") )
	{ 
		ar1 = 0.;
		fprintf(ioQQQ," I will now make invalid operation (div 0 by 0) to get crash.  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then"
			" there are problems.\n");
		fflush(ioQQQ);
		br1 = ar1 / ZeroNum;
		fprintf(ioQQQ," I am still alive - something is wrong, the result was %e\n",
			br1);
		lgCrash = true;
	}

	/* assert that the set_NaN routine works properly for floats */
	else if( p.nMatch("SETN") && p.nMatch("FLOA") )
	{
		sys_float f;
		fprintf(ioQQQ," I will now initialize a float to a signaling NaN. This should never crash!\n");
		set_NaN(f);
		fprintf(ioQQQ," Initialization finished. I will now perform an operation on this variable."
			"  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then"
			" there are problems.\n");
		fflush(ioQQQ);
		f *= 2.f;
		fprintf(ioQQQ," I am still alive - something is wrong, the result was %e\n",
			f);
		lgCrash = true;
	}

	/* assert that the set_NaN routine works properly for doubles */
	else if( p.nMatch("SETN") )
	{
		double d;
		fprintf(ioQQQ," I will now initialize a double to a signaling NaN. This should never crash!\n");
		set_NaN(d);
		fprintf(ioQQQ," Initialization finished. I will now perform an operation on this variable."
			"  Hold on.\n");
		fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong ....\" then"
			" there are problems.\n");
		fflush(ioQQQ);
		d *= 2.;
		fprintf(ioQQQ," I am still alive - something is wrong, the result was %e\n",
			d);
		lgCrash = true;
	}

	/* test what happens with an array index out of bounds
	 * two options, low for [<0] and high for [>limit] */
	else if( p.nMatch("BOUN") )
	{
		double x;

		/* read offset */
		x = p.FFmtRead();
		if( p.lgEOL() && p.nMatch(" LOW" ) )
			x = -2.;
		if( p.lgEOL() && p.nMatch("HIGH" ) )
			x = 2.;

		/* if x >= 0 (which includes default case where x is not entered)
		 * i will be x beyond the end of the array, or x before the start */
		long int i = ( x >= 0. ) ? (long)(x+0.5) + ARR_SIZE : (long)(x-0.5);

		/* must turn off PCLint detection of logical errors in this block */
		if( p.nMatch("STAT") )
		{
			fprintf(ioQQQ," I will now access static array element ar2[%ld].  Hold on.\n", i );
			fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
				" there are problems.\n");
			fflush(ioQQQ);
			ar2[i] = 1e-10;

			fprintf(ioQQQ," I am still alive - something is wrong\n" );
			fflush(ioQQQ);
		}
		else if( p.nMatch("STAC") || p.nMatch("AUTO") )
		{
			double a[ARR_SIZE];
			fprintf(ioQQQ," I will now access automatic array element a[%ld].  Hold on.\n", i );
			fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
				" there are problems.\n");
			fflush(ioQQQ);
			a[i] = 1e-10;

			fprintf(ioQQQ," I am still alive - something is wrong, return value was %.2e\n", a[i] );
			fflush(ioQQQ);
		}
		else if( p.nMatch("ARRA") )
		{
			array<int,ARR_SIZE> ibound;		
			fprintf(ioQQQ," I will now access array element ibound[%ld].  Hold on.\n", i );
			fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
				" there are problems.\n");
			fflush(ioQQQ);
			ibound[i] = 1;
			fprintf(ioQQQ," I am still alive - something is wrong, return value is %i\n" , ibound[i] );
			fflush(ioQQQ);
		}
		else if( p.nMatch("VECT") )
		{
			vector<int> ibound(ARR_SIZE);		
			fprintf(ioQQQ," I will now access vector array element ibound[%ld].  Hold on.\n", i );
			fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
				" there are problems.\n");
			fflush(ioQQQ);
			ibound[i] = 1;
			fprintf(ioQQQ," I am still alive - something is wrong, return value is %i\n" , ibound[i] );
			fflush(ioQQQ);
		}
		else if( p.nMatch("MULT") )
		{
			/* this tests the multi_arr class testing which occurs if the 
			 * macro BOUNDS_CHECK is set at compile time */
			multi_arr<double,2> b;
			b.reserve(3);
			for( int j=0; j < 3; j++ )
				b.reserve(j,ARR_SIZE+j);
			b.alloc();
			if( p.nMatch("ITER") )
			{
				fprintf(ioQQQ," I will now access multi_arr array element *b.ptr(0,%ld)."
					"  Hold on.\n", i );
				fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
					" there are problems.\n\n");
				fflush(ioQQQ);
				md2i p = b.ptr(0,i);
				*p = 2.;
				fprintf(ioQQQ," I am still alive - something is wrong, return value is %g\n", *p );
				fflush(ioQQQ);
			}
			else
			{
				fprintf(ioQQQ," I will now access multi_arr array element b[0][%ld].  Hold on.\n", i );
				fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
					" there are problems.\n\n");
				fflush(ioQQQ);
				b[0][i] = 2.;
				fprintf(ioQQQ," I am still alive - something is wrong, return value is %g\n", b[0][i] );
				fflush(ioQQQ);
			}
			b.clear();
		}
		else if( p.nMatch("AVXP") )
		{
			avx_ptr<int> ibound(ARR_SIZE);
			fprintf(ioQQQ," I will now access avx_ptr element ibound[%ld].  Hold on.\n", i );
			fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
				" there are problems.\n");
			fflush(ioQQQ);
			ibound[i] = 1;
			fprintf(ioQQQ," I am still alive - something is wrong, return value is %i\n" , ibound[i] );
			fflush(ioQQQ);			
		}
		else
		{
			fprintf(ioQQQ," The CRASH BOUNDS command has five different tests.  One must be specified\n" );
			fprintf(ioQQQ," The STATIC option tests a static declared array, and the STACK or AUTO option"
				" tests an automatic array - these test pgcc and g++.\n");
			fprintf(ioQQQ," The ARRAY or VECTOR option tests the STL container - this tests g++.\n");
			fprintf(ioQQQ," The MULTI option tests if bounds checking is enabled in the multi_arr class"
				" (i.e., if the preprocessor macro BOUNDS_CHECK has been set).\n" );
			fprintf(ioQQQ," The AVXPTR option tests if bounds checking is enabled in the avx_ptr class.\n");
			fprintf(ioQQQ," All have a number as an optional argument, the array index to be accessed.\n");
			fflush(ioQQQ);
		}
		lgCrash = true;
	}

	/* test the isnan function */
	else if( p.nMatch("ISNA") )
	{
		if( p.nMatch("FLOA") )
		{
			sys_float ff;
			fprintf(ioQQQ," I will now set a float to SNaN. This should never crash!\n" );
			set_NaN( ff );
			fprintf(ioQQQ," I will now test this variable with the isnan function\n" );
			fprintf(ioQQQ," The correct behavior is for the statement \"PROBLEM DISASTER An assert"
				" has been thrown, this is bad\" to be printed, followed by lots more scary"
				" looking messages.\n\n");
			fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
				" there are problems.\n");
			ASSERT( !isnan( ff ) );
			fprintf(ioQQQ," I am still alive - something is wrong, value is %e\n", ff );
		}
		else
		{
			double dd;
			fprintf(ioQQQ," I will now set a double to SNaN. This should never crash!\n" );
			set_NaN( dd );
			fprintf(ioQQQ," I will now test this variable with the isnan function\n" );
			fprintf(ioQQQ," The correct behavior is for the statement \"PROBLEM DISASTER An assert"
				" has been thrown, this is bad\" to be printed, followed by lots more scary"
				" looking messages.\n\n");
			fprintf(ioQQQ," If the next line says \"I am still alive - something is wrong\" then"
				" there are problems.\n");
			ASSERT( !isnan( dd ) );
			fprintf(ioQQQ," I am still alive - something is wrong, value is %e\n", dd );
		}
		lgCrash = true;
	}

	/* test if a C++ exception is caught */
	else if( p.nMatch("EXCE") )
	{
		fprintf(ioQQQ," I will now throw a C++ exception of type out_of_range()\n" );
		fprintf(ioQQQ," The correct behavior is for the statement \"DISASTER - An out_of_range"
			" exception was caught, what() = Cloudy Test. Bailing out...\" to be printed.\n\n");
		fprintf(ioQQQ," If you get any other message, the exception was not caught correctly.\n\n");
		throw out_of_range( "Cloudy Test" );
		fprintf(ioQQQ," If you see this statement, the exception did not terminate the program.\n" );
		lgCrash = true;
	}

	/* test if domain errors from vectorized math routines are caught correctly */
	else if( p.nMatch("DOMA") )
	{
#if __AVX__
		double x1, x2, y[4];
		x1 = numeric_limits<double>().infinity();
		x2 = numeric_limits<double>().quiet_NaN();
		fprintf(ioQQQ," I will now invoke vhypot() with invalid arguments.\n" );
		fprintf(ioQQQ," The correct behavior is for the statement \"DISASTER - A vectorized math"
			" routine threw a domain_error. Bailing out...\" to be printed.\n\n" );
		fprintf(ioQQQ," If you get any other message, the exception was not caught correctly.\n\n");
		vhypot(y,0.,0.,-1.,x1,3.,3.,x2,2.);
		fprintf(ioQQQ," If you see this statement, the exception did not terminate the program.\n" );
		lgCrash = true;
#else
		fprintf(ioQQQ," AVX vectorization is not enabled, skipping crash test...\n" );
#endif
	}

	/* test if TotalInsanity is caught */
	else if( p.nMatch("INSA") )
	{
		fprintf(ioQQQ,
				  " I will now call TotalInsanity(), which is used to report when an internal\n"
				  " inconsistency has been found & the code must exit.\n"
				  " The correct behavior is to print the following statement:\n\n"
				  "   \"Something that cannot happen, has happened.\n"
				  "    This is TotalInsanity, I live in service.cpp.\"\n\n"
				  " and exit reporting \"PROBLEM DISASTER\".\n"
				  " If the code continues, the TotalInsanity() function is broken.\n\n");
		fprintf(ioQQQ," Calling TotalInsanity()...\n\n");
		TotalInsanity();
		fprintf(ioQQQ," If you see this statement, then TotalInsanity() did not terminate the program.\n" );
		lgCrash = true;
	}
	/* randomly crash some models in a grid, this will throw all sorts of exceptions */
	else if( p.nMatch("GRID") )
	{
		fprintf( ioQQQ, " When the grid starts, certain grid points will randomly crash\n"
			" due to random errors. However, the grid should still run to\n"
			" completion normally and the save files should contain placeholder\n"
			" output for all crashed grid points.\n" );
		grid.lgCrash = true;
	}
	/* simulate aborting the code */
	else if( p.nMatch("ABOR") )
	{
		fprintf( ioQQQ, " This command will cause Cloudy to abort. The correct behavior is to print:\n"
			 " \"ABORT DISASTER PROBLEM - Cloudy aborted, reason: failed due to CRASH ABORT\"\n"
			 " If you get any other message, then the abort did not work correctly.\n"
			 " I will now abort the code....\n\n" );
		throw cloudy_abort( "failed due to CRASH ABORT" );
		fprintf( ioQQQ, " If you see this statement, the abort did not terminate the program.\n" );
		lgCrash = true;
	}

	else
	{
		fprintf(ioQQQ,
			"Crash option not found - valid options are ZERO, UNDEfined, OVERflow, SEGFault, ASSErt,"
			" _NAN, SETNan, BOUNds, ISNAn, EXCEption, DOMAin, INSAnity, GRID, and ABORt.\nSorry.\n");
		lgCrash = true;
	}

	if( lgCrash )
	{
		cdEXIT(EXIT_FAILURE);
	}
}
