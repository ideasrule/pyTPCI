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
		bool lgFirst=true;
		double flux , hden , temp;
		double TotalPressure;		
		/* gas pressure */
		double GasPressure;				
		/* radiation pressure */
		double RadiationPressure;
		/* the ionization parameter */
		double u;


		FILE *ioRES ;
		char chLine[100];

		/* initialize the code for this run */
		cdInit();

		/* calculation's results */
		ioRES = open_data("hazy_kmt.txt","w");

		flux = 18.;
		for( hden=1.; hden<12.01; hden += 0.25 )
		{
			/* initialize the code for this run */
			cdInit();
			cdTalk(false);
			/*cdNoExec( );*/

			/* inputs */
			cdRead( "background 0 .0000000001"  );
			cdRead( "table agn "  );
			cdRead( "stop zone 1 "  );
			cdRead( "# set radiative recombination badnell "  );
			cdRead( "# set dielectronic recombination badnell "  );
			cdRead( "# set dielectronic recombination kludge off "  );

			sprintf(chLine,"phi(h) %f ",flux);
			cdRead( chLine  );

			sprintf(chLine,"hden %f ",hden);
			cdRead( chLine  );

			/* log of the ionization parameter */
			u = flux - 10.4771 - hden;

			/* actually call the code */
			if( cdDrive() )
				exit_status = ES_FAILURE;
			temp = cdTemp_last();
			cdPressure_last(
				/* total pressure, all forms*/
				&TotalPressure,			
				/* gas pressure */
				&GasPressure,				
				/* radiation pressure */
				&RadiationPressure);	

			if( lgFirst )
			{
				/* print a header */
				lgFirst = false;
				fprintf( ioRES , "hden\tlogU\tlogU/T\tlog T\tPgas\n");
				fprintf( stderr, "hden\tlogU\tlogU/T\tlog T\tPgas\n");
			}
			fprintf(ioRES,"%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
				hden, u, u-log10(temp) , log10(temp) , GasPressure ); 

			fprintf(stderr,"%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
				hden, u, u-log10(temp) , log10(temp) , GasPressure ); 
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

