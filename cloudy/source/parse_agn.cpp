/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAgn parse parameters for the AGN continuum shape command */
#include "cddefines.h"
#include "rfield.h"
#include "input.h"
#include "optimize.h"
#include "parser.h"

void ParseAgn(Parser &p)
{
	double BigBump, 
	  Ratio, 
	  XRays, 
	  xnu;

	DEBUG_ENTRY( "ParseAgn()" );

	/* this radiation field will be something like an AGN */
	strcpy( rfield.chSpType[rfield.nShape], "AGN  " );

	/* there were no numbers on the line - this could be the Kirk option,
		* to use the numbers for the continuum in the atlas paper */
	if( p.nMatch("KIRK") )
	{
		/* million degree cutoff, but in rydbergs */
		rfield.slope[rfield.nShape] = 1e6 /  TE1RYD;

		/* cutoff is second parameter is really alpha ox */
		rfield.cutoff[rfield.nShape][0] = -1.40;

		/* bb slope is third parameter */
		rfield.cutoff[rfield.nShape][1] = -0.50;

		/* slope of X-Ray component is last parameter */
		rfield.cutoff[rfield.nShape][2] = -1.0;
	}
	else
	{
		/* first parameter is temperature of big bump
		* second parameter is alpha ox */
		/* slope is first parameter is really temperature of bump */
		rfield.slope[rfield.nShape] = p.FFmtRead();
		if( p.lgEOL() )
		{

			fprintf( ioQQQ, " The big bump temperature should have been on this line.   Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( rfield.slope[rfield.nShape] <= 0. )
		{
			fprintf( ioQQQ, " Non positive temperature not allowed.   Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* temps are log if first le 10 */
		if( rfield.slope[rfield.nShape] <= 10. )
			rfield.slope[rfield.nShape] = 
			exp10(rfield.slope[rfield.nShape]);

		/* want cutoff in ryd not kelvin */
		rfield.slope[rfield.nShape] /= TE1RYD;

		/* cutoff is second parameter is really alpha ox */
		rfield.cutoff[rfield.nShape][0] = p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " alpha ox should have been on this line.   Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( rfield.cutoff[rfield.nShape][0] > 3. || 
			rfield.cutoff[rfield.nShape][0] < -3. )
		{
			fprintf( ioQQQ, " An alpha ox of%10.2e looks funny to me.  Check Hazy to make sure its ok.\n", 
			rfield.cutoff[rfield.nShape][0] );
		}

		if( rfield.cutoff[rfield.nShape][0] >= 0. )
		{
			fprintf( ioQQQ, " The sign of alpha ox is almost certainly incorrect.   Check Hazy.\n" );
		}

		/* bb slope is third parameter */
		rfield.cutoff[rfield.nShape][1] = p.FFmtRead();
		if( p.lgEOL() )
			rfield.cutoff[rfield.nShape][1] = -0.5f;

		/* slope of X-Ray component is last parameter */
		rfield.cutoff[rfield.nShape][2] = p.FFmtRead();
		if( p.lgEOL() )
			rfield.cutoff[rfield.nShape][2] = -1.0f;
	}

	/* 403.3 is ratio of energies where alpha ox defined,
	 * assumed to be 2500A and 2keV */
	Ratio = pow(403.3,rfield.cutoff[rfield.nShape][0] - 1.);

	/* following code must be kept parallel with ffun1 */
	xnu = 0.3645;
	BigBump = pow(xnu,-1. + rfield.cutoff[rfield.nShape][1])*
	  sexp(xnu/rfield.slope[rfield.nShape]);
	xnu = 147.;

	/* XRays = xnu**(-2.) */
	XRays = pow(xnu,rfield.cutoff[rfield.nShape][2] - 1.);
	if( BigBump <= 0. )
	{
		fprintf( ioQQQ, " Big Bump had zero flux at .3645 Ryd.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	realnum SaveCutoff = (realnum)rfield.cutoff[rfield.nShape][0];
	rfield.cutoff[rfield.nShape][0] = (Ratio/(XRays/BigBump));

	/* vary option */
	if( optimize.lgVarOn )
	{
		/* AGN with its four parameters */
		strcpy( optimize.chVarFmt[optimize.nparm], "AGN LOG T=%f, a(ox)=%f a(uv)=%f a(x)=%f" );
		optimize.nvarxt[optimize.nparm] = 4;
		// specified in K but stored in Ryd, need log of temp stored here
		optimize.vparm[0][optimize.nparm] = (realnum)log10(rfield.slope[rfield.nShape]*TE1RYD);
		optimize.vparm[1][optimize.nparm] = SaveCutoff;
		optimize.vparm[2][optimize.nparm] = (realnum)rfield.cutoff[rfield.nShape][1];
		optimize.vparm[3][optimize.nparm] = (realnum)rfield.cutoff[rfield.nShape][2];

		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* the increment in the first steps away from the original value */
		optimize.vincr[optimize.nparm] = 0.5f;
		++optimize.nparm;
	}


	/* lastly increment number of spectra and check that we still have room in the vector */
	++rfield.nShape;
	if( rfield.nShape >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
