/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that calls cloudy with 2 different models, called twice 
to confirm that code is properly initialized */
#include "cddefines.h"
#include "cddrive.h"
/*
This model runs two sets of models that have identical boundary conditions.  
The output and results should also be identical.  This tests that all 
variables have been properly reset.  If the results of the second 
calculation differ from the first then the code has some memory of 
previous calculations, a bad feature if it is to be used to compute 
a grid of calculations.

There are two sets of output.  The files comp4a.out and comp4b.out have the main
output from the code.  The files file1.txt and file2.txt have the punched results
files.  These files should be identical if all is well.
*/

int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
		long int i;
		char chCard[200];

		for( i=0; i<2; ++i )
		{
			/* initialize the code for this run */
			cdInit();

			if( i==0 )
				cdOutput( "comp4a.out" );
			else
				cdOutput( "comp4b.out" );

			/*this is a very simple constant temp model */
			cdRead("test " );
			/* must have some grains to malloc their space in this grid */
			cdRead("grains ism abundance -10 " );
			cdRead("no times " );
			/* write results to either file1.txt or file2.txt */
			cdRead("print lines column  " );
			sprintf(chCard , "punch results column \"file%li.txt\" hide ",i+1);
			cdRead( chCard );
			if( cdDrive() )
				exit_status = ES_FAILURE;
			/* end of the first model */

			/* start of the second model, fully molecular */
			cdInit();
			cdRead("blackbody 5000 " );
			cdRead("luminosity total solar linear 2 " );
			cdRead("brems 6 " );
			cdRead("luminosity total solar log -2.7 " );
			cdRead("hden 10 " );
			cdRead("grains ism " );
			cdRead("metals deplete ");
			cdRead("stop temperature 10K linear " );
			cdRead("radius 15.8  " );
			cdRead("stop zone 1 " );
			cdRead("no times " );
			cdRead("print lines column  " );
			/* write results to either file1.txt or file2.txt */
			sprintf(chCard , "punch results column \"file%li.txt\" hide ",i+1);
			if( cdDrive() )
				exit_status = ES_FAILURE;
			/* end of the second model */

			/* start of the third model */
			cdInit();
			/* inputs */
			cdRead( "ioniz -1 "  );
			cdRead( "sphere "  );
			cdRead( "grains ism function sublimation"  );
			cdRead( "metals deplete ");
			cdRead( "table agn "  );
			cdRead( "hden 11 " );
			cdRead( "stop column density 19 " );
			cdRead( "stop zone 2 " );
			cdRead( "no times " );
			cdRead( "print lines column  " );
			/* write results to either file1.txt or file2.txt */
			sprintf(chCard , "punch results column \"file%li.txt\" hide ",i+1);
			/* actually call the code */
			if( cdDrive() )
				exit_status = ES_FAILURE;
		}

		// redirect to stdout
		cdOutput();

		printf("calculations complete, now compare comp4a.out and comp4b.out, and file1.txt and file2.txt\n");
		printf("vsdiff comp4a.out comp4b.out\n");
		printf("vsdiff file1.txt file2.txt\n");

		cdEXIT(exit_status);
	}
	catch( bad_alloc & )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Most likely your computer "
			 "ran out of memory.\n Try monitoring the memory use of your run. Bailing out...\n" );
		exit_status = ES_BAD_ALLOC;
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_OUT_OF_RANGE;
	}
	catch( domain_error& e )
	{
		fprintf( ioQQQ, " DISASTER - A vectorized math routine threw a domain_error. Bailing out...\n" );
		fprintf( ioQQQ, " What() = %s", e.what() );
		exit_status = ES_DOMAIN_ERROR;
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line() , e.comment() );
		exit_status = ES_BAD_ASSERT;
	}
	catch( bad_signal& e )
	{
		if( e.sig() == SIGILL )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
			exit_status = ES_ILLEGAL_INSTRUCTION;
		}
		else if( e.sig() == SIGFPE )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
			exit_status = ES_FP_EXCEPTION;
		}
		else if( e.sig() == SIGSEGV )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
			exit_status = ES_SEGFAULT;
		}
#		ifdef SIGBUS
		else if( e.sig() == SIGBUS )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
			exit_status = ES_BUS_ERROR;
		}
#		endif
		else
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );
			exit_status = ES_UNKNOWN_SIGNAL;
		}
	}
	catch( cloudy_abort& e )
	{
		fprintf( ioQQQ, " ABORT DISASTER PROBLEM - Cloudy aborted, reason: %s\n", e.comment() );
		exit_status = ES_CLOUDY_ABORT;
	}
	catch( cloudy_exit& e )
	{
		if( ioQQQ != NULL )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		exit_status = e.exit_status();
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}

	cdPrepareExit(exit_status);

	return exit_status;
}

