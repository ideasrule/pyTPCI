/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*grid_do called by cdDrive */
#include "cddefines.h"
#include "conv.h"
#include "input.h"
#include "called.h"
#include "version.h"
#include "init.h"
#include "prt.h"
#include "trace.h"
#include "grainvar.h"
#include "parse.h"
#include "grid.h"
#include "atmdat.h"
#include "flux.h"

/* grid_do called by cdDrive, calls gridXspec or optimize_do */
void grid_do()
{
	char chNote[8];
	long int i, 
	  ii, 
	  j;
	realnum ptem[LIMPAR]; 

	DEBUG_ENTRY( "grid_do()" );

	/* main driver for optimization runs
	 * Drives cloudy to grid variables;*/

	/* code originally written by R.F. Carswell, IOA Cambridge */

	/* this will be number of times grid calls cloudy */
	optimize.nOptimiz = 0;

	/* variables with optimizer */
	for( i=0; i < LIMPAR; i++ )
	{
		optimize.OptIncrm[i] = 0.;
		optimize.varang[i][0] = -FLT_MAX;
		optimize.varang[i][1] = FLT_MAX;
		/* this should be overwritten by format of vary line */
		strcpy( optimize.chVarFmt[i], "error - no optimizer line image was set" );
	}

	/* necessary to do this to keep all lines in */
	prt.lgFaintOn = false;
	conv.LimFail = 1000;

	/* this initializes variables at the start of each simulation
	* in a grid, before the parser is called - this must set any values
	* that may be changed by the command parser */
	InitDefaultsPreparse();

	/* Read the isotope data and allocate the required space */
	LoadIsotopes();

	optimize.lgInitialParse = true;

	/* call READR the first time to scan off all variable options */
	/* this is just an initial parsing to get the number of iterations and
	 * the number of varied parameters.  The other Init* routines are not 
	 * called after this because this is all done again later for each grid point */
	ParseCommands();

	optimize.lgInitialParse = false;

	/* >>chng 00 aug 09, return memory allocated for grains, they are not used, PvH */
	gv.clear();

	optimize.nvary = optimize.nparm;

	/* option to change default increments; if zero then leave as is */
	for( i=0; i < LIMPAR; i++ )
	{
		if( optimize.OptIncrm[i] != 0. )
		{
			optimize.vincr[i] = optimize.OptIncrm[i];
		}
	}

	if( called.lgTalk )
	{
		/* check that at least 1 observed quantity was entered */
		unsigned long nObsQuant = optimize.xLineInt_Obs.size() + optimize.ContNFnu.size() +
			optimize.temp_obs.size() + optimize.ColDen_Obs.size();
		if( optimize.lgOptLum )
			nObsQuant++;
		if( optimize.lgOptDiam )
			nObsQuant++;
		if( nObsQuant == 0 && !grid.lgGrid )
		{
			fprintf( ioQQQ, " The input stream has vary commands, but\n" );
			fprintf( ioQQQ, " no observed quantities were entered.  Whats up?\n" );
			fprintf( ioQQQ, " Use the NO VARY command if you intended to disable optimization.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that the total number of parameters to vary is greater than 1 */
		if( optimize.nvary < 1 )
		{
			fprintf( ioQQQ, " No parameters to vary were entered. Whats up?\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( optimize.nvary > long(nObsQuant) && !grid.lgGrid )
		{
			fprintf( ioQQQ, " PROBLEM - More parameters are varied then there are observables.\n" );
			fprintf( ioQQQ, " PROBLEM - This run may not converge as a result.\n" );
			fprintf( ioQQQ, " PROBLEM - Please reduce the number of free parameters,"
				 " or add more observables.\n" );
		}

		if( strcmp(optimize.chOptRtn,"XSPE") == 0 && optimize.nRangeSet != optimize.nvary )
		{
			fprintf( ioQQQ, " Every parameter with a VARY option must have a GRID specified,\n" );
			fprintf( ioQQQ, " and the GRID must be specified after the VARY option.\n" );
			fprintf( ioQQQ, " These requirements were not satisfied for %ld parameter(s).\n",
				 abs(optimize.nvary - optimize.nRangeSet) );
			cdEXIT(EXIT_FAILURE);
		}

		/* lgTrOptm set with trace grid command */
		if( trace.lgTrOptm )
		{
			for( i=0; i < optimize.nvary; i++ )
			{
				/*print the command format as debugging aid */
				fprintf( ioQQQ, "%s\n", optimize.chVarFmt[i]);

				/* now generate the actual command with parameters */
				string chLine = MakeInputLine(i);

				/* print the resulting command line*/
				fprintf( ioQQQ, "%s\n", chLine.c_str() );
			}
		}

		/* say who we are */
		if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
			fprintf( ioQQQ, "%58cGrid  Driver\n", ' ' );
		else
			fprintf( ioQQQ, "%54cOptimization  Driver\n", ' ' );
		int indent = (int)((122 - t_version::Inst().chVersion.length())/2);
		fprintf( ioQQQ, "%*cCloudy %s\n\n",indent,' ',t_version::Inst().chVersion.c_str());
		fprintf( ioQQQ, "%23c**************************************%7.7s**************************************\n",
			 ' ', t_version::Inst().chDate.c_str() );
		fprintf( ioQQQ, "%23c*%81c*\n", ' ', ' ' );

		/* now echo initial input quantities with flag for vary */
		/* first loop steps over all command lines entered */
		for( i=0; i < long(input.crd.size()); i++ )
		{
			/* put space to start line, overwrite if vary found */
			strcpy( chNote, "       " );
			/* loop over all vary commands, see if this is one */
			for( j=0; j < optimize.nvary; j++ )
			{
				if( i == optimize.nvfpnt[j] )
				{
					/* this is a vary command, put keyword at start */
					strcpy( chNote, "VARY>>>" );
				}
			}

			// exclude lines from init files
			if( input.crd[i]->InclLevel == 0 )
				fprintf( ioQQQ, "%22.7s * %-80s*\n", chNote, input.crd[i]->chCardSav.c_str() );
		}
		fprintf( ioQQQ, "%23c*%81c*\n", ' ', ' ' );
		fprintf( ioQQQ, "%23c***********************************************************************************\n\n\n", ' ' );

		/* option to trace logical flow within this sub */
		if( optimize.lgOptimFlow )
		{
			for( j=0; j < optimize.nvary; j++ )
			{
				i = optimize.nvfpnt[j];
				fprintf( ioQQQ, " trace:%80.80s\n", input.crd[i]->chCardSav.c_str());
				fprintf( ioQQQ, "%80.80s\n", optimize.chVarFmt[j]);
				fprintf( ioQQQ, " number of variables on line:%4ld\n", 
					 optimize.nvarxt[j] );
				fprintf( ioQQQ, " Values:" );
				for( ii=1; ii <= optimize.nvarxt[j]; ii++ )
				{
					fprintf( ioQQQ, "%10.2e", optimize.vparm[ii-1][j] );
				}
				fprintf( ioQQQ, "\n" );
			}
		}

		if( strcmp(optimize.chOptRtn,"PHYM") == 0 )
		{
			fprintf( ioQQQ, " Up to %ld iterations will be performed,\n", 
				 optimize.nIterOptim );
			fprintf( ioQQQ, " and the final version of the input file will be written to the file %s\n", 
				 chOptimFileName.c_str() );

			fprintf( ioQQQ, " The Phymir method will be used" );
			if( optimize.lgParallel )
			{
				if( cpu.i().lgMPI() )
					fprintf( ioQQQ, " in MPI mode.\n" );
				else
					fprintf( ioQQQ, " in parallel mode.\n" );

				fprintf( ioQQQ, " The maximum no. of CPU's to be used is %ld.\n",
					 optimize.useCPU );
			}
			else
				fprintf( ioQQQ, " in sequential mode.\n" );
		}

		else if( strcmp(optimize.chOptRtn,"SUBP") == 0 )
		{
			fprintf( ioQQQ, " Up to %ld iterations will be performed,\n", 
				 optimize.nIterOptim );
			fprintf( ioQQQ, " and the final version of the input file will be written to the file %s\n", 
				 chOptimFileName.c_str() );

			fprintf( ioQQQ, " The Subplex method will be used.\n" );
		}

		else if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
		{
			fprintf( ioQQQ, " Producing grid output.\n" );
		}

		else
		{
			fprintf( ioQQQ, " I do not understand what method to use.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		fprintf( ioQQQ, "\n   %ld parameter(s) will be varied.  The first lines, and the increments are:\n", 
			 optimize.nvary );

		for( i=0; i < optimize.nvary; i++ )
		{
			optimize.varmax[i] = -FLT_MAX;
			optimize.varmin[i] = FLT_MAX;

			/* now generate the actual command with parameters */
			string chLine = MakeInputLine(i);

			fprintf( ioQQQ, "\n %s\n", chLine.c_str() );
			if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
			{
				if( grid.paramValuesFromList[i].size() != 0U )
				{
					fprintf( ioQQQ, " %li parameters read from list.\n", grid.numParamValues[i] );
				}				
				else
				{
					fprintf( ioQQQ, " %s increment is %.3g, the limits are %.3g to %.3g\n",
						 grid.lgLinearIncrements[i] ? "Linear" : "Log",
						 grid.paramIncrements[i], grid.paramLimits[i][0], grid.paramLimits[i][1] );
				}
			}
			else
				fprintf( ioQQQ, " Initial increment is %.3g, the limits are %.3g to %.3g\n", 
					 optimize.vincr[i], optimize.varang[i][0], optimize.varang[i][1] );
		}
	}

	if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
	{
		if( called.lgTalk )
		{
			if( cpu.i().lgMPI() )
				fprintf( ioQQQ, "\n Running in MPI grid mode on %ld CPUs. ", cpu.i().nCPU() );
			else if( grid.lgParallel )
				fprintf( ioQQQ, "\n Running in parallel grid mode on %d CPUs. ", grid.useCPU );
			else
				fprintf( ioQQQ, "\n Running in single-CPU grid mode. " );
			fprintf( ioQQQ, "I will now start to write the input files.\n\n" );
		}

		for( j=0; j < optimize.nvary; j++ )
			ptem[j] = grid.paramLimits[j][0]; 
		for( j=optimize.nvary; j < LIMPAR; j++ )
		{
			ptem[j] = 0.f; 
			grid.paramIncrements[j] = 0.f;
			grid.lgLinearIncrements[j] = false;
		}

		gridXspec(ptem,optimize.nvary);

		if( called.lgTalk )
		{
			fprintf( ioQQQ, " **************************************************\n" );
			fprintf( ioQQQ, " **************************************************\n" );
			fprintf( ioQQQ, " **************************************************\n" );
			fprintf( ioQQQ, "\n Writing input files has been completed.\n\n\n" );
		}
	}
	else
	{
		called.lgTalk = false;
		/* this flag is needed to turn print on to have effect */
		called.lgTalkIsOK = false;

		optimize_do();
	}
}
