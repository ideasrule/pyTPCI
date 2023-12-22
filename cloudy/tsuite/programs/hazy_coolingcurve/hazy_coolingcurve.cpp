/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* runs pure collisional models at range of temperatures, prints cooling */
#include "cddefines.h"
#include "cddrive.h"

/*int main( int argc, char *argv[] )*/
int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
		double telog , cooling, hden , tehigh , dTelog;

		FILE *ioRES;
		char chLine[100];

		/* initialize the code for this run */
		cdInit();

		/* file to output short form of calculation's results */
		ioRES = open_data("hazy_coolingcurve.txt","w");

		/* the log of the hydrogen density cm-3 */
		hden = 0.;

		/* the log of the initial and final temperatures 
		 * this calc has no ionization at very
		 * low temperatures, except for cosmic ray collisions and the cosmic background */
		telog = 1.;
		tehigh = 9.;

		/* increment in log of T, normally 0.1 */
		dTelog = 0.1;

		/* print simple header */
		fprintf(ioRES, "log T\tlog cool/n2\n" ); 
		fprintf(stderr,"log T\tlog cool/n2\n" ); 

		/* we do not want to generate any output */
		cdOutput( "hazy_coolingcurve.out" );

		/* loop over all densities */
		while( telog <= (tehigh+0.0001) )
		{
			/* initialize the code for this run */
			cdInit();

			cdTalk( true );

			/* if this is uncommented the calculation will not be done,
			 * but all parameters will be generated, as a quick way to see
			 * that grid is set up properly */
			/*cdNoExec( );*/

			/* cosmic background, microwave and hard parts */
			cdRead( "background 0 "  );

			/* cosmic ray background, this is included because drives
			 * the chemistry in molecular gas */
			cdRead( "cosmic ray background "  );

			/* this is a pure collisional model to turn off photoionization 
			   cdRead( "no photoionization "  );*/

			/* set a continuum shape, even though not used 
			   cdRead( "table agn "  );
			   cdRead( "phi(h) -4 "  );*/

			/* do only one zone */
			cdRead( "stop zone 1 "  );

			/* set the hydrogen density */
			sprintf(chLine,"hden %f ",hden);
			cdRead( chLine  );

			/* sets the gas kinetic temperature */
			sprintf(chLine,"constant temper %f ",telog);
			cdRead( chLine  );

			/* identify sources of heating and cooling */
			cdRead( "punch heating \"hazy_coolingcurve.het\" last no hash no clobber "  );
			cdRead( "punch cooling \"hazy_coolingcurve.col\" last no hash no clobber "  );

			/* actually call the code */
			if( cdDrive() )
				exit_status = ES_FAILURE;

			/* get cooling for last zone */
			cooling = cdCooling_last();

			/* want to print cooling over density squared */
			cooling = cooling / pow(10.,hden*hden);

			fprintf(ioRES, "%.2f\t%.3f", telog , log10( cooling ) ); 
			fprintf(stderr,"%.2f\t%.3f", telog , log10( cooling ) ); 

			if( exit_status == ES_FAILURE )
			{
				fprintf(ioRES ,"\t problems!!\n");
				fprintf(stderr,"\t problems!!\n");
			}
			else
			{
				fprintf(ioRES ,"\n");
				fprintf(stderr,"\n");
			}

			telog += dTelog;
		}

		fclose(ioRES);

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

