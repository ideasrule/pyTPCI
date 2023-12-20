/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*optimize_do main driver for optimization runs*/
#include "cddefines.h"
#define	NPLXMX	(LIMPAR*(LIMPAR+6)+1)
#include "input.h"
#include "called.h"
#include "prt.h"
#include "optimize.h"

/* called by cdDrive */
void optimize_do()
{
	long int iflag, 
	  ii, 
	  iworke[NPLXMX], 
	  j, 
	  mode, 
	  need, 
	  nfe, 
	  np;
	chi2_type fret;
	realnum fx,
	  ptem[LIMPAR], 
	  delta[LIMPAR], 
	  toler,
	  worke[NPLXMX];

	DEBUG_ENTRY( "optimize_do()" );

	/* main driver for optimization runs
	 * Drives cloudy to optimize variables;*/

	/* code originally written by R.F. Carswell, IOA Cambridge */

	toler = (realnum)log10(1. + optimize.OptGlobalErr);

	if( strcmp(optimize.chOptRtn,"PHYM") == 0 )
	{
		/* Phymir method */
		for( j=0; j < optimize.nvary; j++ )
		{
			ptem[j] = optimize.vparm[0][j];
			delta[j] = optimize.vincr[j];
		}
		/* >>chng 06 jan 02, fix uninitialized var problem detected by valgrind/purify, PvH */
		for( j=optimize.nvary; j < LIMPAR; j++ )
		{
			ptem[j] = -FLT_MAX;
			delta[j] = -FLT_MAX;
		}
		optimize_phymir(ptem,delta,optimize.nvary,&fret,toler);
		for( j=0; j < optimize.nvary; j++ )
			optimize.vparm[0][j] = ptem[j];
	}

	else if( strcmp(optimize.chOptRtn,"SUBP") == 0 )
	{
		fprintf( ioQQQ, " Begin optimization with SUBPLEX\n" );
		need = 2*optimize.nvary + optimize.nvary*(optimize.nvary + 4) + 1;
		if( need > NPLXMX )
		{
			fprintf( ioQQQ, " Increase size of NPLXMX in parameter statements to handle this many variables.\n" );
			fprintf( ioQQQ, " I need at least %5ld\n", need );
			cdEXIT(EXIT_FAILURE);
		}
		for( j=0; j < optimize.nvary; j++ )
			ptem[j] = optimize.vparm[0][j];

		/* The subroutine SUBPLX input into cloudy 8/4/94.
		 * The program itself is very well commented.
		 * The mode must set to 0 for the default values.
		 * The switch iflag tells if the program terminated normally.   */
		mode = 0;

		/*  >>chng 97 dec 08, remove first arg, optimize_func, since not used in routines */
		optimize_subplex(
			/* the number of parameters to vary */
			optimize.nvary,
			/* the relative error, single number */
			toler,
			/* maximum number of function evaluations before giving up */
			optimize.nIterOptim,
			/* mode of operation, we simply set to zero */
			mode,
			/* the initial changes in the guessed best coefficients, typically 0.2 to 1  */
			optimize.vincr,
			/* a vector of nvary initial parameters that are the starting guesses for the parameters */
			ptem,
			/* a realnum, this is simply ignored */
			&fx,
			/* another parameter that is simply ignored, a long int */
			&nfe,
			/* a realnum that is NPLXMX long, used for working space by the routine */
			/* an array that is 20*26 + 1 elements long, used for working space */
			worke,
			/* a long int that is NPLXMX long, used for working space by the routine */
			/* an array that is 20*26 + 1 elements long, used for working space */
			iworke,
			/* a long int - says what happened, if -1 then exceeded nIterOptim iteration */
			&iflag);

		if( iflag == -1 )
		{
			fprintf( ioQQQ, " SUBPLEX exceeding maximum iterations.\n"
				 " This can be reset with the OPTIMZE ITERATIONS command.\n" );
		}

		for( j=0; j < optimize.nvary; j++ )
			optimize.vparm[0][j] = ptem[j];

		if( optimize.lgOptimFlow )
		{
			fprintf( ioQQQ, " trace return optimize_subplex:\n" );
			for( j=0; j < optimize.nvary; j++ )
			{
				fprintf( ioQQQ, " Values:" );
				for( ii=1; ii <= optimize.nvarxt[j]; ii++ )
					fprintf( ioQQQ, " %.2e", optimize.vparm[ii-1][j] );
				fprintf( ioQQQ, "\n" );
			}
		}
	}
	else
		TotalInsanity();

	// turn output back on for final model
	called.lgTalk = cpu.i().lgMPI_talk();
	called.lgTalkIsOK = cpu.i().lgMPI_talk();
	prt.lgFaintOn = true;

	if( called.lgTalk )
	{
		fprintf( ioQQQ, " **************************************************\n" );
		fprintf( ioQQQ, " **************************************************\n" );
		fprintf( ioQQQ, " **************************************************\n" );
		fprintf( ioQQQ, "\n Cloudy was called %4ld times.\n\n", optimize.nOptimiz );

		for( long i=0; i < optimize.nvary; i++ )
		{
			np = optimize.nvfpnt[i];

			/* now generate the actual command with parameters */
			input.crd[np]->chCardSav = MakeInputLine(i);

			fprintf( ioQQQ, " Optimal command: %s\n", input.crd[np]->chCardSav.c_str());
			fprintf( ioQQQ, "  Smallest value:%10.2e Largest value:%10.2e Allowed range %10.2e to %10.2e\n", 
				 optimize.varmin[i], optimize.varmax[i], optimize.varang[i][0], 
				 optimize.varang[i][1] );
		}
	}

	if( cpu.i().lgMaster() )
	{
		/* save optimal parameters */
		FILE *ioOptim = open_data( chOptimFileName, "w" );
		for( size_t i=0; i < input.crd.size(); i++ )
			if( input.crd[i]->InclLevel == 0 )
				fprintf( ioOptim, "%s\n", input.crd[i]->chCardSav.c_str());
		fclose( ioOptim );

		/* recalculate values in cloudy for the best fit, and print out
		 * all the information */
		fprintf( ioQQQ, "\f" );

		for( long i=0; i < optimize.nvary; i++ )
			ptem[i] = optimize.vparm[0][i];

		(void)optimize_func(ptem);
	}
}
