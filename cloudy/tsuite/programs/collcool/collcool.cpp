/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program calling cloudy to produce a table giving ionization vs temperature */
#include "cddefines.h"
#include "cddrive.h"

int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
#		define NMAX 100
		/* NELEMI and NELEMF define the range of elements for which the cooling will be calculated */
#		define NELEMI 0
#		define NELEMF LIMELM
		double xCoolSave[NMAX][LIMELM], tesave[NMAX];
		double telog , teinc, te_last, te_first;
		double hden, element_abund ;
		long int nte , i;
		long int nelem;
		/* this is the list of element names used to make printout */
		char chElementName[LIMELM][12] =
			{ "Hydrogen\0" ,"Helium\0" ,"Lithium\0" ,"Beryllium\0" ,"Boron\0" ,
			  "Carbon\0" ,"Nitrogen\0" ,"Oxygen\0" ,"Fluorine\0" ,"Neon\0" ,
			  "Sodium\0" ,"Magnesium\0" ,"Aluminium\0" ,"Silicon\0" ,"Phosphorus\0" ,
			  "Sulphur\0" ,"Chlorine\0" ,"Argon\0" ,"Potassium\0" ,"Calcium\0" ,
			  "Scandium\0" ,"Titanium\0" ,"Vanadium\0" ,"Chromium\0" ,"Manganese\0" ,
			  "Iron\0" ,"Cobalt\0" ,"Nickel\0" ,"Copper\0" ,"Zinc\0"};

		FILE *ioRES ;
		char chLine[100];

		/* initialize the code for this run */
		cdInit();

		/* calculation's results are saved here */
		ioRES = open_data("collcool.txt","w");

		for( i = 0; i < NMAX; i++ )
		{
			for( nelem=0; nelem<LIMELM; ++nelem)
			{
				xCoolSave[i][nelem] = -1.;
			}
		}

		/* the first temperature */
		te_first = 4.;
		te_last = 9.;
		/* the increment in the temperature */
		teinc = 0.2;
		/* the log of the hydrogen density that will be used */
		hden = 0.;
		/* the log of the ratio of the abundance of the element in question to hydrogen */
		element_abund = 5.;

		int npoints = (te_last-te_first)/teinc + 1;
		

		for( nelem=NELEMI; nelem<NELEMF; ++nelem)
		{
			nte = 0;
			telog = te_first;

			while( telog <= te_last+0.01 && nte<NMAX)
			{
				/* initialize the code for this run */
				cdInit();
				cdTalk(false);
				/*cdNoExec( );*/
				printf("%s  Te = %g\n",chElementName[nelem],telog);
				
				sprintf(chLine,"coronal %3.1f ",telog);
				cdRead( chLine  );

				/* just do the first zone - only want ionization distribution */
				cdRead( "stop zone 1 "  );
				cdRead( "set dr 0 "  );

				/* the hydrogen density entered as a log */
				sprintf(chLine,"hden %i ",int(hden));
				cdRead( chLine  );

				cdRead( "set dr 0"  );

				// cannot turn elements off and on during one run, so set them
				// all to small abundances, relative to H
				cdRead( "element helium abundance -9"  );
				cdRead( "metals -9"  );

				sprintf(chLine,"element %.11s abundance %f ",chElementName[nelem],element_abund);
				cdRead( chLine  );

				cdRead( "set eden 0"  );

				/* actually call the code */
				if( cdDrive() )
					exit_status = ES_FAILURE;

				xCoolSave[nte][nelem] = cdCooling_last();
				if( nelem != 0 )
				{
					xCoolSave[nte][nelem] /= pow(10.,element_abund);
				}
				//printf("%6.2e\t\%6.2e\n",xCoolSave[nte][nelem],cdCooling_last());

				tesave[nte] = telog;
				telog += teinc;
				++nte;
			}
		}
			
		/* this generates large printout */
		fprintf(ioRES,"Te");
		for( nelem=NELEMI; nelem<NELEMF; ++nelem)
		{
			fprintf(ioRES,"\t%s",chElementName[nelem]);
		}
		fprintf(ioRES,"\n");
		
		for( i = 0; i < npoints; i++ )
		{
			fprintf(ioRES,"%6.2e",pow(10.,tesave[i]));
			for( nelem=NELEMI; nelem<NELEMF; ++nelem)
			{
				fprintf(ioRES,"\t%6.2e",xCoolSave[i][nelem]);
			}
			fprintf(ioRES,"\n");
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
