/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseMetal parse parameters on metal command */
#include "cddefines.h"
#include "input.h"
#include "optimize.h"
#include "elementnames.h"
#include "grainvar.h"
#include "called.h"
#include "abund.h"
#include "parser.h"


/* GetMetalsDeplete data file used by METALS DEPLETE command */
STATIC void GetMetalsDeplete( Parser &p, const bool lgPrintMetalsDeplete )
{
	/* save depletions read in from external file */
	static int lgFirst=true;
	if( lgFirst )
	{
		lgFirst = false;
		string chFile;	/*file name for table read */
		if( p.nMatch( "\"" ) )
		{
			/*
			 * if a quote occurs on the line then get the ini file name
			 * this will also set the name in chCard and OrgCard to spaces
			 * so later keywords do not key off it
			 */
			if( p.GetQuote( chFile ) )
				p.StringError();
		}
		else
		{
			/* no quote appeared, so this is the default name, cloudy.ini */
			chFile = "ISM_CloudyClassic.dpl";
		}
		string chPath = "abundances" + cpu.i().chDirSeparator() + chFile;
		FILE *ioDATA = open_data( chPath, "r" );	// will abort if not found

		if( lgPrintMetalsDeplete )
			fprintf(ioQQQ," First call, GetMetalsDeplete opened file %s \n", chPath.c_str() );

		// init with no depletion set, equal to 1
		for(int nelem=0; nelem<LIMELM; ++nelem)
			abund.DepletionScaleFactor[nelem] = 1.;

		string chLine;
		while( read_whole_line( chLine, ioDATA ) )
		{
			if( lgPrintMetalsDeplete )
				fprintf(ioQQQ, "line: %s", chLine.c_str() );

			/* field of stars or empty line end data */
			if( chLine.size() == 0 || chLine[0]=='\n' || chLine[0]=='\r' || chLine[0]=='*' )
				break;

			/* skip comment */
			if( chLine[0]=='#' )
				continue;

			size_t pp;
			/* erase EOL character */
			if( (pp = chLine.find_first_of("\n\r")) != string::npos )
				chLine.erase(pp);

			string chCAPS = chLine;
			caps(chCAPS);

			bool lgFound = false;
			for(int nelem=0; nelem<LIMELM; nelem++)
			{
				if( chCAPS.find(elementnames.chElementNameShort[nelem]) != string::npos )
				{
					lgFound = true;
					long int i = 1;
					bool lgEOL;

					abund.DepletionScaleFactor[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					if( lgPrintMetalsDeplete )
						fprintf(ioQQQ, " Parse:\t%s\t%.3f\n",
								elementnames.chElementNameShort[nelem] , abund.DepletionScaleFactor[nelem] );
					if( abund.DepletionScaleFactor[nelem]<0. || abund.DepletionScaleFactor[nelem]>1. )
					{
						fprintf(ioQQQ," The grain depletion for element %s was %.3e, it must be between 0 and 1\n",
								elementnames.chElementNameShort[nelem] , abund.DepletionScaleFactor[nelem] );
						cdEXIT(EXIT_FAILURE);
					}

					//we shouldn't need to continue once an element name is found on the line...
					break;
				}
			}
			if( !lgFound )
			{
				fprintf(ioQQQ, "PROBLEM in METALS DEPLETE: did not identify element name on this line: %s\n",
					chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
		}

		fclose( ioDATA );
	}
}

static double AX[LIMELM] , BX[LIMELM] , ZX[LIMELM];
static int lgSetJenkins09[LIMELM];

/* EvalJenkins convert Fstar into Jenkins 2009 depletion factors */
STATIC void EvalJenkins( const double Fstar , const double DxLimit )
{
	for( int nelem=0; nelem<LIMELM; ++nelem )
	{
		if( lgSetJenkins09[nelem] )
		{
			abund.DepletionScaleFactor[nelem] = pow(10., BX[nelem] + AX[nelem]*(Fstar-ZX[nelem]) );
			abund.DepletionScaleFactor[nelem] = MAX2( 0., abund.DepletionScaleFactor[nelem] );
			abund.DepletionScaleFactor[nelem] = MIN2( DxLimit, abund.DepletionScaleFactor[nelem] );
		}
		else
			abund.DepletionScaleFactor[nelem] = 1.;
	}
}

/* GetJenkins09 parse the Jenkins 2009 metals data file */
STATIC void GetJenkins09( Parser &p, const bool lgPrtJenkins09, const double Fstar, const double DxLimit )
{
	static int lgFirst=true;
	if( lgFirst )
	{
		enum { DEBUG_ENTRY = false };

		lgFirst = false;
		string chFile;	/*file name for table read */
		if( p.nMatch( "\"" ) )
		{
			/*
			 * if a quote occurs on the line then get the ini file name
			 * this will also set the name in chCard and OrgCard to spaces
			 * so later keywords do not key off it
			 */
			if( p.GetQuote( chFile ) )
				p.StringError();
		}
		else
		{
			/* no quote appeared, so this is the default name, cloudy.ini */
			chFile = "ISM_Jenkins09_Gunasekera21.dpl";
		}
		string chPath = "abundances" + cpu.i().chDirSeparator() + chFile;
		FILE *ioDATA = open_data( chPath, "r" );	// will abort if not found
		if( DEBUG_ENTRY && lgPrtJenkins09 )
			fprintf(ioQQQ," First call, opened file %s \n\n", chPath.c_str() );

		// init with no depletion set
		for(int nelem=0; nelem<LIMELM; ++nelem)
			lgSetJenkins09[nelem] = false;

		string chLine;
		while( read_whole_line( chLine, ioDATA ) )
		{
			if( DEBUG_ENTRY && lgPrtJenkins09 )
				fprintf(ioQQQ, "line: %s", chLine.c_str() );

			/* field of stars or empty line end data */
			if( chLine.size() == 0 || chLine[0]=='\n' || chLine[0]=='\r' || chLine[0]=='*' )
				break;

			/* skip comment */
			if( chLine[0]=='#' )
				continue;

			size_t pp;
			/* erase EOL character */
			if( (pp = chLine.find_first_of("\n\r")) != string::npos )
				chLine.erase(pp);

			string chCAPS = chLine;
			caps(chCAPS);

			bool lgFound = false;
			if( DEBUG_ENTRY && lgPrtJenkins09 )
				fprintf(ioQQQ, " GetJenkins09 reports contents of data file:\n" );
			for(int nelem=0; nelem<LIMELM; nelem++)
			{
				if( chCAPS.find(elementnames.chElementNameShort[nelem]) != string::npos )
				{
					lgFound = true;
					long int i = 1;
					bool lgEOL;

					lgSetJenkins09[nelem] = true;
					AX[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					BX[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					ZX[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					if( DEBUG_ENTRY && lgPrtJenkins09 )
						fprintf(ioQQQ, " GetJenkins09:\t%s\t%.3f\t%.3f\t%.3f\n",
								elementnames.chElementNameShort[nelem] , AX[nelem] , BX[nelem] , ZX[nelem] );

					//we shouldn't need to continue once an element name is found on the line...
					break;
				}
			}
			if( !lgFound )
			{
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES: did not identify element name on this line: %s\n",
					chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
		}

		fclose( ioDATA );
	}

	/* option to print table of depletions */
	if( lgPrtJenkins09 )
	{
		fprintf(ioQQQ,"\n GetJenkins09: report of range of depletion scale factors follows:\n Fstar");
		for( int nelem=0; nelem<LIMELM; ++nelem )
			fprintf(ioQQQ,"\t%s" , elementnames.chElementNameShort[nelem] );
		fprintf(ioQQQ,"\tSum depleted\n");

		/* print range of Fstar, do not change original value */
		for( double FstarLoc=0; FstarLoc<1.01; FstarLoc+=0.1 )
		{
			fprintf(ioQQQ,"%.3f",FstarLoc);
			EvalJenkins( FstarLoc , DxLimit );
			for( int nelem=0; nelem<LIMELM; ++nelem )
				fprintf(ioQQQ,"\t%.3e", abund.DepletionScaleFactor[nelem]);
			fprintf(ioQQQ,"\t%.3e", abund.SumDepletedAtoms());
			fprintf(ioQQQ,"\n");
		}
	}

	/* return correct depletions for specified Fstar */
	EvalJenkins( Fstar , DxLimit );
}

void ParseMetal(Parser &p)
{
	DEBUG_ENTRY( "ParseMetal()" );

	/* default upper limit to depletion scale factor */	
	double DxLimit = 1e38;

	/* parse the metals command */
	abund.lgAbnReference = false;	
	if( p.nMatch("DEPL") )
	{
		/* keyword depletion is present
		 * deplete by set of scale factors */
		abund.lgDepln = true;

		/* option to use Jenkins 2009ApJ...700.1299J fits to ISM depletion
		 * syntax is METALS DEPLETE JENKINS 2009 FSTAR=0.5 */
		if( p.nMatch("JENKINS") && p.nMatch("2009") )
		{
			/* read Fstar, Jenkins' basic parameter, and confirm it is in range */
			realnum Fstar = 0;
			Fstar = (realnum)p.FFmtRead();
			if( fabs( Fstar-2009) > 0.1 )
			{
				fprintf(ioQQQ," The first number of METALS DEPLETE JENKINS 2009 "
						"must be 2009, it was %.3f\n Sorry\n\n ", Fstar);
				cdEXIT( EXIT_FAILURE );
			}
			/* Fstar depletion strength */
			Fstar = (realnum)p.FFmtRead();

			// check for limit specifying upper bound on depletion scale factor
			if( p.nMatch( "LIMIT" ) )
				DxLimit = p.FFmtRead();

			/* sort out whether log */
			bool lgLogOn = false;
			if( p.nMatch(" LOG") )
				lgLogOn = true;
			else if( p.nMatch("LINE") )
				lgLogOn = false;

			/* save log of Fstar since needed for vary option, convert Fstar to linear */
			double FsLog , DxLimitLog;
			if( Fstar <0. || lgLogOn )
			{
				FsLog = Fstar;
				Fstar = exp10(Fstar);
				DxLimitLog = DxLimit;
				DxLimit = exp10(DxLimit);
			}
			else
			{
				FsLog = log10(MAX2(1e-38,Fstar));
				DxLimitLog = log10(DxLimit);
			}

			/* check range of final Fstar and DxLimit */
			if(Fstar<0 || Fstar>1. )
			{
				fprintf(ioQQQ, "Fstar on METALS DEPLETE JENKINS command must be "
						"between 0 and 1 and was %.3f\n Sorry\n\n", Fstar );
				cdEXIT( EXIT_FAILURE );
			}

			if( DxLimit < 0 )
			{
				fprintf(ioQQQ,
					" The upper limit to the depletion scale"
					" factor must be greater than zero, it was %.2e\n",
					DxLimit);
				cdEXIT( EXIT_FAILURE );
			}

			/* print option, to report details and print table of values */
			bool lgPrtJenkins09 = false;
			if( p.nMatch("PRINT")  )
				lgPrtJenkins09 = true;

			if( lgPrtJenkins09 )
				fprintf(ioQQQ,"\n Jenkins 2009, print set, found Fstar = %.3e limit = %.3e\n" , Fstar , DxLimit );

			GetJenkins09( p, lgPrtJenkins09, Fstar, DxLimit );

			/* vary option */
			if( optimize.lgVarOn )
			{
				strcpy( optimize.chVarFmt[optimize.nparm], "METALS DEPLETE JENKINS 2009=%f LOG LIMIT=%f" );
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				optimize.vparm[0][optimize.nparm] = (realnum)FsLog;
				optimize.vparm[2][optimize.nparm] = (realnum)DxLimitLog;
				optimize.nvarxt[optimize.nparm] = 2;
				/* initial increment is 0.1 and allowed range i 0 to 1 */
				optimize.vincr[optimize.nparm] = 0.1;
				optimize.varang[optimize.nparm][0] = 0.f;
				optimize.varang[optimize.nparm][1] = 1.f;

				++optimize.nparm;
			}
		}
		else
		{
			realnum scale = realnum(p.FFmtRead());
			if( !p.lgEOL() )
			{
				fprintf( ioQQQ, "\nPROBLEM -- No number allowed with bare 'metals deplete',"
			        		" but found '%g'\n", scale );
				cdEXIT(EXIT_FAILURE);
			}

			/* metals deplete command - not Jenkins */
			bool lgPrintMetalsDeplete = false;
			if( p.nMatch("PRINT")  )
				lgPrintMetalsDeplete = true;

			/* no keyword - use stored abundances */
			GetMetalsDeplete( p, lgPrintMetalsDeplete );

			if( lgPrintMetalsDeplete )
			{
				for( long int i=0; i < LIMELM; i++ )
				{
					fprintf(ioQQQ,"DEBUGGG depnew %s\t%.3e\n", elementnames.chElementName[i], abund.DepletionScaleFactor[i] );
				}
				fprintf(ioQQQ,"Sum depleted \t%.3e\n", abund.SumDepletedAtoms());
			}
		}
	}
	else
	{
		/* metal depletion factor, if negative then it is the log */
		abund.ScaleMetals = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("the depletion factor");

		/* sort out whether log */
		bool lgLogOn = false;
		if( p.nMatch(" LOG") )
		{
			lgLogOn = true;
		}
		else if( p.nMatch("LINE") )
		{
			lgLogOn = false;
		}

		double dmlog;
		if( abund.ScaleMetals <= 0. || lgLogOn )
		{
			dmlog = abund.ScaleMetals;
			abund.ScaleMetals = exp10(abund.ScaleMetals);
		}
		else
		{
			dmlog = log10(abund.ScaleMetals);
		}

		/* option to vary grain abundance as well */
		bool lgGrains;
		if( p.nMatch("GRAI") )
		{
			lgGrains = true;
			gv.GrainMetal = abund.ScaleMetals;
		}
		else
		{
			lgGrains = false;
			gv.GrainMetal = 1.;
		}

		/* vary option */
		if( optimize.lgVarOn )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "METALS= %f LOG" );
			if( lgGrains )
				strcat( optimize.chVarFmt[optimize.nparm], " GRAINS" );
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			optimize.vparm[0][optimize.nparm] = (realnum)dmlog;
			optimize.vincr[optimize.nparm] = 0.5;
			optimize.nvarxt[optimize.nparm] = 1;
			++optimize.nparm;
		}
	}

	return;
}
