/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseGrid parse the grid command lines */
#include "cddefines.h"
#include "grid.h"
#include "input.h"
#include "parser.h"

/* ParseGrid - called from ParseCommands if GRID command found */
void ParseGrid(
	/* command line, which was changed to all caps in main parsing routine */
	Parser &p)
{
	DEBUG_ENTRY( "ParseGrid()" );

	/* RP fake optimizer to run a grid of calculations, also accepts
	 * keyword XSPEC */
	strcpy( optimize.chOptRtn, "XSPE" );
	grid.lgGrid = true;

	if( p.nMatch("REPE") )
	{
		/* just keep repeating, don't actually change the values in the grid.
		 * useful for debugging unintentional crosstalk */
		grid.lgStrictRepeat = true;
	}

	if( p.nMatch("SEPA") )
	{
		/* keep main output files separate, this cuts down the
		 * gathering time substantially for huge grid runs */
		grid.lgKeepMainOutputSeparate = true;
	}

	if( p.nMatch("SEQU") )
	{
		/* run grid in sequential mode */
		grid.lgParallel = false;
	}

	/* 06 aug 22, change to accept three parameters: lower and upper limit and number of points. */
	/* scan off range for the previously selected variable */
	if( optimize.nparm > 0 )
	{
		ASSERT( optimize.nparm <= LIMPAR );
	
		grid.lgLinearIncrements[optimize.nparm-1] = p.nMatch("LINE") ? true : false ;
	
		string chLabel;
		if( p.nMatch("LIST") )
		{
			if( p.GetQuote( chLabel ) )
				p.StringError();
	
			vector<double> tmpVector;
			input_readvector( chLabel, tmpVector );
			if( tmpVector.size() < 2u )
			{
				fprintf(ioQQQ, "PROBLEM there must be at least two grid steps but only one was found.\n");
				cdEXIT( EXIT_FAILURE );
			}
	
			for( auto it = tmpVector.begin(); it != tmpVector.end(); ++it )
				grid.paramValuesFromList[optimize.nparm-1].emplace_back( (realnum)*it );
			grid.numParamValues[optimize.nparm-1] = grid.paramValuesFromList[optimize.nparm-1].size();
			++optimize.nRangeSet;
	
			// Set original limits to end points.
			grid.paramLimits[optimize.nparm-1][0] = *min_element( grid.paramValuesFromList[optimize.nparm-1].begin(),
																  grid.paramValuesFromList[optimize.nparm-1].end() );
			grid.paramLimits[optimize.nparm-1][1] = *max_element( grid.paramValuesFromList[optimize.nparm-1].begin(),
																  grid.paramValuesFromList[optimize.nparm-1].end() );
			// Set increment to zero.
			grid.paramIncrements[optimize.nparm-1] = 0;
		}
		else
		{
			grid.paramLimits[optimize.nparm-1][0] = (realnum)p.FFmtRead();
			grid.paramLimits[optimize.nparm-1][1] = (realnum)p.FFmtRead();
			grid.paramIncrements[optimize.nparm-1] = (realnum)p.FFmtRead();
			if( grid.paramIncrements[optimize.nparm-1] < 0_r )
				grid.lgNegativeIncrements = true;

			/* the increase step should not be 0 */
			if( grid.paramIncrements[optimize.nparm-1] == 0_r )
			{
				fprintf( ioQQQ," The increment (third parameter) should not be zero.\n" );
				fprintf( ioQQQ," Sorry.\n" );
				cdEXIT( EXIT_FAILURE );
			}

			if( p.lgEOL() )
			{
				fprintf( ioQQQ," This command has changed since the definition given in Porter et al. 2006, PASP, 118, 920.\n" );
				fprintf( ioQQQ," The grid command now requires three parameters: lower limit, upper limit, and increment.\n" );
				fprintf( ioQQQ," The keywords RANGE and STEPS are no longer necessary.\n" );
				fprintf( ioQQQ," Sorry.\n" );
				cdEXIT( EXIT_FAILURE );
			}
			else
			{
				++optimize.nRangeSet;
			}

			realnum ratio =	(grid.paramLimits[optimize.nparm-1][1] - grid.paramLimits[optimize.nparm-1][0])/
			      grid.paramIncrements[optimize.nparm-1];

			/* Alert if the uplimit and lowlimit are wrong */
			if( ratio < 0_r )
			{
				fprintf( ioQQQ, "The increment (third parameter) has the wrong sign. \
					It doesn't take you from the initial to the final grid value (first and second parameter, resp.).\n" );
				fprintf( ioQQQ," Sorry.\n" );
				cdEXIT( EXIT_FAILURE );
			}

			// this takes care of the blowup in the error due to cancellation in limits[1]-limits[0]
			// it assumes that limits[1]-limits[0] is accurate within 3*eps*(limits[0]+limits[1])/2
			// which should be a very conservative estimate...
			realnum feps = safe_div(1.5_r*
				(grid.paramLimits[optimize.nparm-1][1] + grid.paramLimits[optimize.nparm-1][0]),
				(grid.paramLimits[optimize.nparm-1][1] - grid.paramLimits[optimize.nparm-1][0]));

			// this will blow for pathologically narrow grid ranges
			ASSERT( abs(feps) < realnum(INT_MAX) );
			int eps = max(int(abs(feps)),3);

			// take special care if step is integer fraction of max-min (which is nearly always the case)
			if( fp_equal( ratio, realnum(nint(ratio)), int(eps) ) )
				grid.numParamValues[optimize.nparm-1] = nint(ratio) + 1;
			else
				grid.numParamValues[optimize.nparm-1] = long(ratio) + 1;

			if( grid.numParamValues[optimize.nparm-1] < 2 )
			{
				fprintf( ioQQQ, " There must be at least two grid points in each dimension.\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT( EXIT_FAILURE );
			}

			grid.numParamValues[optimize.nparm-1] = MAX2( 2, grid.numParamValues[optimize.nparm-1] );
		}

		if( p.nMatch("CYCL") )
		{
			/* cycle through the grid multiple times, only used for testing */
			grid.nCycle = nint(p.FFmtRead());
			if( p.lgEOL() )
				grid.nCycle = 2;
			if( grid.nCycle < 2 )
			{
				fprintf( ioQQQ, " Invalid repetion number for cycle: %ld\n", grid.nCycle );
				fprintf( ioQQQ, " Usage: grid <p1> <p2> <p3> cycle [ <n> ] with n >= 2.\n" );
				fprintf( ioQQQ, "    or: grid list \"filename\" cycle [ <n> ] with n >= 2.\n" );
				cdEXIT( EXIT_FAILURE );
			}
		}

		if( p.nMatch("NCPU") )
		{
			/* set number of CPUs to be used for the grid */
			long dum = nint(p.FFmtRead());
			if( grid.lgParallel )
			{
				grid.useCPU = p.lgEOL() ? cpu.i().nCPU() : max(dum,0);
			}
			else
			{
				grid.useCPU = 1;
			}
			if( grid.useCPU == 0 )
			{
				fprintf( ioQQQ, " Invalid number of CPUs: %ld\n", dum );
				fprintf( ioQQQ, " Usage: grid <p1> <p2> <p3> ncpus [ <n> ] with n >= 1.\n" );
				fprintf( ioQQQ, "    or: grid list \"filename\" ncpus [ <n> ] with n >= 1.\n" );
				cdEXIT( EXIT_FAILURE );
			}
		}

		// Create some buffer area in the allowed range of parameter values to prevent
		// accidentally going over the limit due to roundoff error. The buffer is 1/10th
		// of a step, so should still guard against doing too many steps due to bugs.
		realnum safety = 0.001f*grid.paramIncrements[optimize.nparm-1];

		if( grid.lgLinearIncrements[optimize.nparm-1] )
		{
			if( grid.paramLimits[optimize.nparm-1][0]-safety<=0. )
			{
				fprintf(ioQQQ,"The current implementation of the grid command works with log parameter values even when you specify LINEAR.\n");
				fprintf(ioQQQ,"A non-positive value was entered.  The grid command cannot deal with this.\n");
				cdEXIT( EXIT_FAILURE );
			}
			optimize.varang[optimize.nparm-1][0] = log10(grid.paramLimits[optimize.nparm-1][0]-safety);
			optimize.varang[optimize.nparm-1][1] = log10(grid.paramLimits[optimize.nparm-1][1]+safety);
		}
		else
		{
			optimize.varang[optimize.nparm-1][0] = grid.paramLimits[optimize.nparm-1][0]-safety;
			optimize.varang[optimize.nparm-1][1] = grid.paramLimits[optimize.nparm-1][1]+safety;
		}
	}

	return;
}
