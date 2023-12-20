/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParsePowerlawContinuum parse the power law continuum command */
#include "cddefines.h"
#include "rfield.h"
#include "optimize.h"
#include "input.h"
#include "parser.h"

void ParsePowerlawContinuum(Parser &p)
{
	DEBUG_ENTRY( "ParsePowerlawContinuum()" );

	/* power law with cutoff and X-ray continuum */
	strcpy( rfield.chSpType[rfield.nShape], "POWER" );

	/* first parameter is slope of continuum, probably should be negative */
	rfield.slope[rfield.nShape] = p.FFmtRead();
	if( p.lgEOL() )
	{
		fprintf( ioQQQ, " There should have been a number on this line.   Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( rfield.slope[rfield.nShape] >= 0. )
	{
		fprintf( ioQQQ, " Is the slope of this power law correct?\n" );
	}

	/* second optional parameter is high energy cut off */
	rfield.cutoff[rfield.nShape][0] = p.FFmtRead();

	/* no cutoff if eof hit */
	if( p.lgEOL() )
	{
		/* no extra parameters at all, so put in extreme cutoffs */
		rfield.cutoff[rfield.nShape][0] = 1e4;
		rfield.cutoff[rfield.nShape][1] = 1e-4;
	}
	else
	{
		/* first cutoff was present, check for second */
		rfield.cutoff[rfield.nShape][1] = p.FFmtRead();
		if( p.lgEOL() )
			rfield.cutoff[rfield.nShape][1] = 1e-4;
	}

	/* check that energies were entered in the correct order */
	if( rfield.cutoff[rfield.nShape][1] > rfield.cutoff[rfield.nShape][0] )
	{
		fprintf( ioQQQ, " The optional cutoff energies do not appear to be in the correct order.  Check Hazy.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* check for keyword KELVIN to interpret cutoff energies as degrees */
	if( p.nMatch("KELV") )
	{
		/* temps are log if first le 10 */
		if( rfield.cutoff[rfield.nShape][0] <= 10. || p.nMatch(" LOG") )
			rfield.cutoff[rfield.nShape][0] = exp10(rfield.cutoff[rfield.nShape][0]);
		if( rfield.cutoff[rfield.nShape][1] <= 10. || p.nMatch(" LOG") )
			rfield.cutoff[rfield.nShape][1] = exp10(rfield.cutoff[rfield.nShape][1]);
		rfield.cutoff[rfield.nShape][0] /= TE1RYD;
		rfield.cutoff[rfield.nShape][1] /= TE1RYD;
	}

	if( rfield.cutoff[rfield.nShape][0] < 0. || 
		rfield.cutoff[rfield.nShape][1] < 0. )
	{
		fprintf( ioQQQ, " A negative cutoff energy is not physical.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( rfield.cutoff[rfield.nShape][1] == 0. && 
		rfield.slope[rfield.nShape] <= -1. )
	{
		fprintf( ioQQQ, " A power-law with this slope, and no low energy cutoff, may have an unphysically large\n brightness temperature in the radio.\n" );
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		if( p.nMatch("VARYB") )
		{
			/* this test is for key "varyb", meaning to vary second parameter
			 * the cutoff temperature
			 * this is the number of parameters to feed onto the input line */
			optimize.nvarxt[optimize.nparm] = 1;
			sprintf( optimize.chVarFmt[optimize.nparm], "POWER LAW %f KELVIN %%f %f LOG", 
				 rfield.slope[rfield.nShape],
				 log10(rfield.cutoff[rfield.nShape][1]*TE1RYD) );
			optimize.vparm[0][optimize.nparm] =
				(realnum)log10(rfield.cutoff[rfield.nShape][0]*TE1RYD);
			optimize.varang[optimize.nparm][0] = (realnum)log10(rfield.cutoff[rfield.nShape][1]*TE1RYD);
			optimize.varang[optimize.nparm][1] = FLT_MAX;
			optimize.vincr[optimize.nparm] = 0.2f;
		}
		else if( p.nMatch("VARYC") )
		{
			/* the keyword was "varyc"
			 * this is the number of parameters to feed onto the input line */
			optimize.nvarxt[optimize.nparm] = 1;
			sprintf( optimize.chVarFmt[optimize.nparm], "POWER LAW %f KELVIN %f %%f LOG", 
				 rfield.slope[rfield.nShape],
				 log10(rfield.cutoff[rfield.nShape][0]*TE1RYD) );
			optimize.vparm[0][optimize.nparm] =
				(realnum)log10(rfield.cutoff[rfield.nShape][1]*TE1RYD);
			optimize.varang[optimize.nparm][0] = -FLT_MAX;
			optimize.varang[optimize.nparm][1] = (realnum)log10(rfield.cutoff[rfield.nShape][0]*TE1RYD);
			optimize.vincr[optimize.nparm] = 0.2f;
		}
		else
		{
			/* vary the first parameter only, but still are two more
			 * this is the number of parameters to feed onto the input line */
			optimize.nvarxt[optimize.nparm] = 1;
			sprintf( optimize.chVarFmt[optimize.nparm], "POWER LAW %%f KELVIN %f %f LOG", 
				 log10(rfield.cutoff[rfield.nShape][0]*TE1RYD),
				 log10(rfield.cutoff[rfield.nShape][1]*TE1RYD) );
			optimize.vparm[0][optimize.nparm] = (realnum)rfield.slope[rfield.nShape];
			optimize.vincr[optimize.nparm] = 0.2f;
		}
		++optimize.nparm;
	}

	/*>>chng 06 nov 10, BUGFIX, nShape was incremented before previous branch
	 * and so crashed with log10 0 since nSpage was beyond set values
	 * caused fpe domain function error with log 0
	 * fpe caught by Pavel Abolmasov */
	++rfield.nShape;
	if( rfield.nShape >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}
