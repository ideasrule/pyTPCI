/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* test case to show behavior of very high metals environmet going to lte
 * cloud irradiated by STE blackbody, should equilibriate at continuum temperature
 * output hizlte.txt gives input continuum temperature, deduced electron temperature, error */
#include "cddefines.h"
#include "cddrive.h"
/*int main( int argc, char *argv[] )*/
int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
		double telog , teinc;
		double hden , te, error,BigError;

		FILE *ioRES ;
		char chLine[100];

		/* initialize the code for this run */
		cdInit();

		/* calculation's results */
		ioRES = open_data("hizlte.txt","w");
		fprintf(ioRES,"T(BB)\tT(equil)\terror\n");

		/* the first temperature */
		telog = 4.;
		/* the increment in the temperature */
		teinc = 0.5;
		/* the hydrogen density that will be used */
		hden = 10.;

		BigError = 0.;
		/* this is limit on 32 bit double */
		while( telog < 7.99 )
		{
			/* initialize the code for this run */
			cdInit();
			cdTalk(false);
			/*cdNoExec( );*/
			printf("te %g ",telog);

			/* inputs */
			cdRead( "title high metal clouds into lte"  );

			cdRead( "iterate"  );
			cdRead( "set dr 0"  );
			cdRead( "abundances starburst 10 "  );

			sprintf(chLine,"blackbody lte %f ",telog);
			cdRead( chLine  );

			sprintf(chLine,"hden %f ",hden);
			cdRead( chLine  );
			cdRead( "stop zone 1 "  );

			/* actually call the code */
			if( cdDrive() )
				exit_status = ES_FAILURE;

			te = cdTemp_last();
			error = (pow(10.,telog)-te)/te;
			BigError = MAX2( BigError , fabs(error) );
			/*remember largest error */
			fprintf(ioRES,"%e\t%e\t%e\n", telog , log10(te) ,error);
			printf(" %g %g\n",log10(te),error );
			telog += teinc;
		}
		printf("biggest error was %g\n", BigError );
		fprintf(ioRES,"biggest error was %g\n", BigError );
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
#ifdef CATCH_SIGNAL
	catch( bad_signal& e )
	{
		if( ioQQQ != NULL )
		{
			if( e.sig() == SIGINT || e.sig() == SIGQUIT )
			{
				fprintf( ioQQQ, " User interrupt request. Bailing out...\n" );
				exit_status = ES_USER_INTERRUPT;
			}
			else if( e.sig() == SIGTERM )
			{
				fprintf( ioQQQ, " Termination request. Bailing out...\n" );
				exit_status = ES_TERMINATION_REQUEST;
			}
			else if( e.sig() == SIGILL )
			{
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
				exit_status = ES_ILLEGAL_INSTRUCTION;
			}
			else if( e.sig() == SIGFPE )
			{
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
				exit_status = ES_FP_EXCEPTION;
			}
			else if( e.sig() == SIGSEGV )
			{
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
				exit_status = ES_SEGFAULT;
			}
#			ifdef SIGBUS
			else if( e.sig() == SIGBUS )
			{
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
				exit_status = ES_BUS_ERROR;
			}
#			endif
			else
			{
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );
				exit_status = ES_UNKNOWN_SIGNAL;
			}

		}
	}
#endif
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
