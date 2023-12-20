/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCoronal parse parameters off coronal equilibrium command */
#include "cddefines.h"
#include "thermal.h"
#include "input.h"
#include "optimize.h"
#include "phycon.h"
#include "dynamics.h"
#include "parser.h"
#include "atmdat.h"

/*ParseCoronal parse parameters off coronal equilibrium command */
void ParseCoronal(Parser &p)
{
	double a;

	DEBUG_ENTRY( "ParseCoronal()" );

	/* use coronal command to establish initial conditions in a cooling
	 * time-varying cloud */
	if( p.nMatch( "INIT"  ) && p.nMatch( "TIME"  ) )
	{
		dynamics.lgTimeDependentStatic = true;
		dynamics.lgRecom = true;
		dynamics.lg_coronal_time_init = true;
		if( p.nMatch( "TRAC"  ) )
			dynamics.lgTracePrint = true;
	}

	/* coronal equilibrium; set constant temperature to number on line */
	thermal.lgTemperatureConstant = true;
	thermal.lgTemperatureConstantCommandParsed = true;

	// kinetic temperatures for a given ion are higher for coronal equilibrium
	// simulations - large Fe chianti/stout models are needed to get the full cooling
	if( !atmdat.lgChiantiLevelsSet )
	{
		atmdat.nChiantiMaxLevelsFe = atmdat.nDefaultCollLevelsFe;
		atmdat.nChiantiMaxLevels = atmdat.nDefaultCollLevels;
	}
	if( !atmdat.lgStoutLevelsSet )
	{
		atmdat.nStoutMaxLevelsFe = atmdat.nDefaultCollLevelsFe;
		atmdat.nStoutMaxLevels = atmdat.nDefaultCollLevels;
	}

	a = p.FFmtRead();
	if( p.lgEOL() )
	{
		fprintf( ioQQQ, " There should be a temperature on this line.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* numbers less than or equal to 10 are the log of the temperature */
	if( (a <= 10. && !p.nMatch("LINE")) || p.nMatch(" LOG") )
	{
		thermal.ConstTemp = (realnum)exp10(a);
	}
	else
	{
		thermal.ConstTemp = (realnum)a;
	}

	/* check temperature bounds */
	if( thermal.ConstTemp < phycon.TEMP_LIMIT_LOW )
	{
		thermal.ConstTemp = (realnum)(1.0001*phycon.TEMP_LIMIT_LOW);
		fprintf( ioQQQ, " PROBLEM Te too low, reset to %g K.\n",
			 thermal.ConstTemp );
	}
	if( thermal.ConstTemp > phycon.TEMP_LIMIT_HIGH )
	{
		thermal.ConstTemp = (realnum)(0.9999*phycon.TEMP_LIMIT_HIGH);
		fprintf( ioQQQ, " PROBLEM Te too high, reset to %g K.\n",
			 thermal.ConstTemp );
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		/*  no luminosity options on vary */
		optimize.nvarxt[optimize.nparm] = 1;
		strcpy( optimize.chVarFmt[optimize.nparm], "COROnal equilibrium %f LOG" );
		if( dynamics.lg_coronal_time_init )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME INIT" );

		/*  pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;

		/*  log of temp will be pointer */
		optimize.vparm[0][optimize.nparm] = (realnum)log10(thermal.ConstTemp);
		optimize.vincr[optimize.nparm] = 0.1f;
		optimize.varang[optimize.nparm][0] = (realnum)log10(1.00001*phycon.TEMP_LIMIT_LOW);
		optimize.varang[optimize.nparm][1] = (realnum)log10(0.99999*phycon.TEMP_LIMIT_HIGH);
		++optimize.nparm;
	}
	return;
}
