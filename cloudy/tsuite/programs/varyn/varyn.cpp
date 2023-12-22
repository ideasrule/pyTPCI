/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* test case to vary density at constant temperature, to show how line ratios change */
#include "cddefines.h"
#include "cddrive.h"
/*int main( int argc, char *argv[] )*/
int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
		double hden , hdeninc;
		char chLine[100];
		
		/* initialize the code for this run */
		cdInit();

		/* calculation's results */
		FILE *ioRES = open_data("varyn.txt","w");
		fprintf(ioRES,"density\t[OII] 3726\t[O II] 3729\t[O II] 2471\t[O II] 7323\t[O II] 7332\n");

		/* the increment in the temperature */
		hdeninc = 0.5;
		/* the hydrogen density that will be used */
		hden = -1.;

		/* this is limit on 32 bit double */
		while( hden < 8. )
		{
			double relint , absint ;
			/* initialize the code for this run */
			cdInit();
			cdTalk(false);
			/*cdNoExec( );*/
			printf("hden %g\n",hden);
			fprintf(ioRES,"%g",hden);

			/* inputs */
			cdRead( "title vary density at constant temperature"  );

			sprintf(chLine,"hden %f ",hden);
			cdRead( chLine  );

			cdRead( "table agn "  );
			cdRead( "stop zone 1 "  );
			cdRead( "ionization parameter -3 "  );
			cdRead( "constant temperature 4 "  );
			cdRead( "normalize to \"O  2\" 3726.03 "  );

			/* actually call the code */
			if( cdDrive() )
				exit_status = ES_FAILURE;

			/* now get the lines */
			cdLine( "O  2" , 3726.03 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "O  2" , 3728.81 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "O  2" , 2470.34 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "O  2" , 7319.99 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			cdLine( "O  2" , 7330.73 , &relint , & absint  );
			fprintf(ioRES, "\t%e", relint );
			fprintf(ioRES, "\n");

			hden += hdeninc;
		}

		cdEXIT(exit_status);
	}
	catch( bad_alloc& )
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

