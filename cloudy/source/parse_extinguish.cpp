/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseExtinguish parse the extinguish command */
#include "cddefines.h"
#include "rfield.h"
#include "parser.h"

/*ParseExtinguish parse the extinguish command */
void ParseExtinguish( Parser &p )
{
	// extinguish ionizing continuum by absorbing column AFTER
	// setting luminosity or intensity.  First number is the column
	// density (log), second number is leakage (def=0%)
	// last number is lowest energy (ryd), last two may be omitted
	// from right to left 
	// 
	// extinction is actually done in extin, which is called by ContSetIntensity

	DEBUG_ENTRY( "ParseExtinguish()" );

	rfield.ExtinguishColumnDensity = (realnum)p.FFmtRead();
	if( p.lgEOL() )
		p.NoNumb("extinguishing column");

	// default is for the number to be the log of the column. 
	// there is a linear option for the column or optical depth,
	// if linear does not occur then log, and convert to linear */
	if( !p.nMatch("LINE"  ) )
	{
		if( rfield.ExtinguishColumnDensity>35. )
		{
			fprintf(ioQQQ,
				" The first parameter on this command line is the log of either the column density or optical depth.\n");
			fprintf(ioQQQ,
				" The value seems pretty big to me - please check it.\n");
			/* flush it since we will probably crash */
			fflush(ioQQQ);
		}
		rfield.ExtinguishColumnDensity = exp10(rfield.ExtinguishColumnDensity);
	}

	/* option to set leakage - default is 0. */
	rfield.ExtinguishLeakage = (realnum)p.FFmtRead();
	if( p.lgEOL() )
		rfield.ExtinguishLeakage = 0.;

	/* negative leaks are logs */
	if( rfield.ExtinguishLeakage < 0. )
		rfield.ExtinguishLeakage = exp10(rfield.ExtinguishLeakage);

	if( rfield.ExtinguishLeakage > 1. )
	{
		/* but leaks greater than 1 are not allowed */
		fprintf( ioQQQ, " A leakage of%9.0f%% was entered - this must be less than 100%%\n", 
			rfield.ExtinguishLeakage*100. );
		cdEXIT(EXIT_FAILURE);
	}
	// user input check that H-ionizing radiation is blocked if
	// table Draine used */
	rfield.lgBlockHIon = true;

	/* option to set lowest energy for absorber */
	rfield.ExtinguishLowEnergyLimit = (realnum)p.FFmtRead();
	if( p.lgEOL() )
		rfield.ExtinguishLowEnergyLimit = 0.99946f;
	else
	{
		if( rfield.ExtinguishLowEnergyLimit <= 0. )
			rfield.ExtinguishLowEnergyLimit = exp10(rfield.ExtinguishLowEnergyLimit);
		if( rfield.ExtinguishLowEnergyLimit < 0.99946 )
			fprintf( ioQQQ, " Energy less than 1 Ryd!!\n" );
	}

	/* specify optical depth at 1 Ryd rather than column density */
	if( p.nMatch("OPTI"  ) )
	{
		/* convert the optical depth into the proper column density */
		rfield.ExtinguishColumnDensity /= (realnum)(rfield.ExtinguishConvertColDen2OptDepth*
			pow(rfield.ExtinguishLowEnergyLimit,rfield.ExtinguishEnergyPowerLow) );
	}

	return;
}
