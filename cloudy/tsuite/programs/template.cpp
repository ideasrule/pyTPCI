/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* Template main program for calling Cloudy as a subroutine
 * Routine returns 0 if model is ok, and 1 if problems occurred. */
#include "cddefines.h"
#include "cddrive.h"

int main( int argc, char *argv[] )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" ); // do not remove this!

	try {
		// ============================================================================
		// START ENTERING YOUR CODE AFTER THIS LINE
		// inside this try block you can enter your code to call Cloudy as a subroutine
		// replace the code between here and the next ===== line with your own program

		// the code always needs to be initialized first
		cdInit();

		// to write output to a file, you MUST open it with open_data()
		// no need to check for a NULL pointer on return as open_data() will have done that already
		FILE* io = open_data("sample.out", "w");

		// you can call Cloudy in a loop if you want...
		for( int i=0; i < 1; ++i )
		{
			// the code always needs to be initialized first
			cdInit();
			// replace this with a series of calls to cdRead to define the input script
			// this particular command line exercises the smoke test
			cdRead( "test" );
			// this calls Cloudy to execute the input script you defined above
			if( cdDrive() )
				exit_status = ES_FAILURE;
		}

		cdEXIT(exit_status); // always exit with cdEXIT, this assures files are properly closed.
		// ============================================================================
		// THIS IS THE END OF THE PROGRAM YOU WROTE
		// this ends the try block - your code calling Cloudy ends above this line
	}
	// here we catch all the possible exceptions that the code can throw, do not remove this!
	catch( bad_alloc )
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
