/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that calls cloudy when used as a stand-alone program */
#include "cddefines.h"
#include "cddrive.h"

/*int main( int argc, char *argv[] )*/
int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
		double hdenLimit , hdenInit , hden, TeInc ,
			r5007,r4363,r1665,ro3_88,ro3_52, r1661 ,
			absolute , temp, TeLimit, TeInit , hdenInc;

		FILE *ioDATA ;
		char chLine[100];

		/* initialize the code for this run */
		cdInit();

		cdOutput( "vary_nete.out" );

		/* calculation's results will go to this file*/
		ioDATA = open_data("vary_nete.txt","w");

		/* the range of density, and the increment in density, for this grid */
		hdenInit = 0.;
		hdenLimit = 7.;
		hdenInc = 1.;

		/* the range of log temperatures for the grid */
		TeInit = 3.5;
		TeLimit = 4.5;
		/* multiplicative inc - log */
		TeInc = .2;

		/* set the density and temperature to the initial values */
		hden = hdenInit;
		temp = TeInit;

		/* print the header for the data file */
		fprintf(ioDATA,
			"density\ttemp\t[O III] 5007\t[O III] 4363\t[O III] 1665\t[O III] 88\t[O III] 52\t[O III] 1661\n");

		while( hden < 1.01 * hdenLimit )
		{
			while( temp < 1.01*TeLimit )
			{
				/* initialize the code for this run */
				cdInit();

				/* option to not execute the code, uncomment when debugging setup */
				/*cdNoexec( );*/

				/* gas temperature for this calculation */
				sprintf(chLine,"constant temperature %f", temp);
				cdRead( chLine );

				/* gas density for this calculation */
				sprintf(chLine,"hden %f", hden);
				cdRead( chLine );

				/* only want the first zone */
				cdRead( "stop zone 1" );

				/* speed up calculation by not included some elements and lines */
				cdRead( "init file \"fast.ini\"" );

				cdRead( "normalize to \"o  3\" 5006.84" );

				/* an incident continuum must be specified to get the code
				 * to run at all - not very important since we will set
				 * the gas temperature */
				cdRead( "blackbody 40000" );
				cdRead( "ionization parameter -3" );

				/* actually call the code */
				if( cdDrive() )
					exit_status = ES_FAILURE;
				/* flush the output so we see it on the screen */
				fflush(ioQQQ);

				fprintf(ioDATA,"%.3e\t%.5f",	hden, temp );
				fprintf(stderr,"%.3e\t%.5f",	hden, temp );

				/************************ O III ******************/
				if( cdLine( "O  3", 5006.84, &r5007 , &absolute ) < 0 )
				{
					printf("did not find O  3 5007\n");
					exit_status = ES_FAILURE;
				}
				/* this should be unity */
				fprintf(ioDATA,"\t%.3e",	r5007 );
				fprintf(stderr,"\t%.3e",	r5007 );

				if( cdLine( "O  3", 4363.21, &r4363 , &absolute ) < 0 )
				{
					printf("did not find totl 4363\n");
					exit_status = ES_FAILURE;
				}
				fprintf(ioDATA,"\t%.3e",	r4363 );
				fprintf(stderr,"\t%.3e",	r4363 );

				if( cdLine( "blnd", 1666., &r1665 , &absolute ) < 0 )
				{
					printf("did not find totl 1665\n");
					exit_status = ES_FAILURE;
				}
				fprintf(ioDATA,"\t%.3e",	r1665 );
				fprintf(stderr,"\t%.3e",	r1665 );

				if( cdLine( "O  3", 883323., &ro3_88 , &absolute ) < 0 )
				{
					printf("did not find O  3 88\n");
					exit_status = ES_FAILURE;
				}
				fprintf(ioDATA,"\t%.3e",	ro3_88 );
				fprintf(stderr,"\t%.3e",	ro3_88 );

				if( cdLine( "O  3", 518004., &ro3_52 , &absolute ) < 0 )
				{
					printf("did not find O  3 52\n");
					exit_status = ES_FAILURE;
				}
				fprintf(ioDATA,"\t%.3e",	ro3_52 );
				fprintf(stderr,"\t%.3e",	ro3_52 );

				if( cdLine( "O  3", 1660.81, &r1661 , &absolute ) < 0 )
				{
					printf("did not find O  3 1661\n");
					exit_status = ES_FAILURE;
				}
				fprintf(ioDATA,"\t%.3e",	r1661 );
				fprintf(stderr,"\t%.3e",	r1661 );


				/************************* end lines with lf and flush it ***************/
				fprintf(ioDATA,"\n" );
				fflush(ioDATA );
				fprintf(stderr,"\n" );

				temp += TeInc;
			}
			temp = TeInit;
			hden += hdenInc;
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

