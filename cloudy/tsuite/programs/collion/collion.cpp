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
		/* faintest ionization fraction to print */
#		define FAINT 1e-9
		/* following is number of ion stages per line */
#		define NELEM 15
#		define NMAX 100
		double xIonSave[NMAX][LIMELM][LIMELM+1] , tesave[NMAX];
		double telog , teinc, te_last;
		double hden ;
		long int nte , i;
		long int nelem , ion;
		/* this is the list of element names used to query code results */
		char chElementNameShort[LIMELM][5] = { "HYDR" , "HELI" ,
						       "LITH" , "BERY" , "BORO" , "CARB" , "NITR" , "OXYG" , "FLUO" ,
						       "NEON" , "SODI" , "MAGN" , "ALUM" , "SILI" , "PHOS" , "SULP" ,
						       "CHLO" , "ARGO" , "POTA" , "CALC" , "SCAN" , "TITA" , "VANA" ,
						       "CHRO" , "MANG" , "IRON" , "COBA" , "NICK" , "COPP" , "ZINC" };
		/* this is the list of element names used to make printout */
		char chElementName[LIMELM][11] =
			{ "Hydrogen  " ,"Helium    " ,"Lithium   " ,"Beryllium " ,"Boron     " ,
			  "Carbon    " ,"Nitrogen  " ,"Oxygen    " ,"Fluorine  " ,"Neon      " ,
			  "Sodium    " ,"Magnesium " ,"Aluminium " ,"Silicon   " ,"Phosphorus" ,
			  "Sulphur   " ,"Chlorine  " ,"Argon     " ,"Potassium " ,"Calcium   " ,
			  "Scandium  " ,"Titanium  " ,"Vanadium  " ,"Chromium  " ,"Manganese " ,
			  "Iron      " ,"Cobalt    " ,"Nickel    " ,"Copper    " ,"Zinc      "};

		FILE *ioRES ;
		char chLine[100];

		/* initialize the code for this run */
		cdInit();

		/* calculation's results are saved here */
		ioRES = open_data("collion.txt","w");
		fprintf(ioRES,"  log fractional ionization for species with abundances > %.2e\n",
			FAINT );

		/* the first temperature */
		telog = 3.;
		te_last = 9.;
		/* the increment in the temperature */
		teinc = 0.1;
		/* the log of the hydrogen density that will be used */
		hden = 0.;

		nte = 0;
		while( telog <= te_last+0.01 && nte<NMAX)
		{
			/* initialize the code for this run */
			cdInit();
			cdTalk(false);
			/*cdNoExec( );*/
			printf("te %g\n",telog);

			/* input continuum is very faint cosmit background - this
			 * should be negligible */
			cdRead( "background 0 .0000000001"  );
 
			/* just do the first zone - only want ionization distribution */
			cdRead( "stop zone 1 "  );

			/* the hydrogen density entered as a log */
			sprintf(chLine,"hden %f ",hden);
			cdRead( chLine  );

			/* this says to compute very small stages of ionization - we normally trim up
			 * the ionizaton so that only important stages are done */
			cdRead( "set trim -20 "  );

			/* the log of the gas kinetic temperature */
			sprintf(chLine,"constant temper %f ",telog);
			cdRead( chLine  );

			/* actually call the code */
			if( cdDrive() )
				exit_status = ES_FAILURE;

			/* now save ionization distribution for later printout */
			for( nelem=0; nelem<LIMELM; ++nelem)
			{
				for( ion=1; ion<=nelem+2;++ion)
				{
					if( cdIonFrac(chElementNameShort[nelem],
						      ion, &xIonSave[nte][nelem][ion],
						      "radius",false) )
					{
						fprintf(ioRES,"\t problems!!\n");
						fprintf(stderr,"\t problems!!\n");
					}
				}
			}

			tesave[nte] = telog;
			telog += teinc;
			++nte;
		}

		/* this generates large printout */
		for( nelem=0; nelem<LIMELM; ++nelem)
		{
			fprintf(ioRES,"\n   Element %li %s\n",
				nelem+1,chElementName[nelem]);
			fprintf(ioRES,"   Te");
			for(i=1; i<MIN2(NELEM+1,nelem+3);++i)
			{
				fprintf(ioRES,"%6li",i);
			}
			fprintf(ioRES,"\n");
			for(i=0;i<nte; ++i)
			{
				fprintf(ioRES,"  %5.2f",tesave[i]);
				for(ion=1;ion<MIN2(NELEM+1,nelem+3);++ion)
				{
					if( xIonSave[i][nelem][ion]>FAINT )
					{
						fprintf(ioRES,"%6.2f",log10(xIonSave[i][nelem][ion]) );
					}
					else
					{
						fprintf(ioRES,"%6s","--  ");
					}
				}
				fprintf(ioRES,"\n");
			}
			/* nelem is on the C scale */
			if( nelem>=NELEM-1 )
			{
				fprintf(ioRES,"\n   Element %li %s\n",
					nelem+1,chElementName[nelem]);
				fprintf(ioRES,"   Te");
				for(i=NELEM+1; i<nelem+3;++i)
				{
					fprintf(ioRES,"%6li",i);
				}
				fprintf(ioRES,"\n");
				for(i=0;i<nte;++i)
				{
					fprintf(ioRES,"  %5.2f",tesave[i]);
					for(ion=NELEM+1;ion<nelem+3;++ion)
					{
						if( xIonSave[i][nelem][ion]>FAINT )
						{
							fprintf(ioRES,"%6.2f",log10(xIonSave[i][nelem][ion]) );
						}
						else
						{
							fprintf(ioRES,"%6s","--  ");
						}
					}
					fprintf(ioRES,"\n");
				}
			}
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

