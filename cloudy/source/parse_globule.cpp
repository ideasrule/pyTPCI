/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseGlobule parse parameters off the globule command */
#include "cddefines.h"
#include "radius.h"
#include "dense.h"
#include "optimize.h"
#include "input.h"
#include "parser.h"

void ParseGlobule(Parser &p)
{
	DEBUG_ENTRY( "ParseGlobule()" );

	if( dense.gas_phase[ipHYDROGEN] > 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER More than one density command was entered.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* globule with density increasing inward
	 * parameters are outer density, radius of globule, and density power */
	radius.glbden = (realnum)p.FFmtRead();
	radius.glbden = p.lgEOL() ? 1.f : exp10(radius.glbden);
	dense.SetGasPhaseDensity( ipHYDROGEN, radius.glbden );

	if( dense.gas_phase[ipHYDROGEN] <= 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER Hydrogen density must be > 0.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	radius.glbrad = (realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		radius.glbrad = 3.086e18f;
	}
	else
	{
		radius.glbrad = exp10(radius.glbrad);
	}

	/* this is largest zone thickness, used to set first zone thickness */
	radius.sdrmax = radius.glbrad/25.;
	radius.lgSdrmaxRel = false;

	/* turn off min dr checking in NEXTDR */
	radius.lgDrMnOn = false;
	radius.glbpow = (realnum)p.FFmtRead();
	if( p.lgEOL() )
		radius.glbpow = 1.;
	strcpy( dense.chDenseLaw, "GLOB" );

	/* this is distance to globule */
	radius.glbdst = radius.glbrad;

	/* vary option */
	if( optimize.lgVarOn )
	{
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;

		/* this is the number of parameters to feed onto the input line */
		optimize.nvarxt[optimize.nparm] = 3;
		// the keyword LOG is not used above, but is checked elsewhere
		strcpy( optimize.chVarFmt[optimize.nparm], "GLOBULE %f LOG %f %f" );

		/* param is log of abundance by number relative to hydrogen */
		optimize.vparm[0][optimize.nparm] = (realnum)log10(radius.glbden);
		optimize.vparm[1][optimize.nparm] = (realnum)log10(radius.glbrad);
		optimize.vparm[2][optimize.nparm] = radius.glbpow;
		optimize.vincr[optimize.nparm] = 0.2f;
		++optimize.nparm;
	}
	return;
}
