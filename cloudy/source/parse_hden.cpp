/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseHDEN parse the hden command */
#include "cddefines.h"
#include "input.h"
#include "dense.h"
#include "optimize.h"
#include "parser.h"

void ParseHDEN(Parser &p )
{
	DEBUG_ENTRY( "ParseHDEN()" );

	if( dense.gas_phase[ipHYDROGEN] > 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER More than one density command was entered.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* log of hydrogen density */
	dense.SetGasPhaseDensity( ipHYDROGEN, (realnum)p.FFmtRead() );
	if( p.lgEOL() )
	{
		fprintf( ioQQQ, " DISASTER The density MUST be entered with this command. STOP\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* check for further options */
	if( ! p.nMatch("LINE") )
	{
		/* check size of density - will we crash? */
		if( dense.gas_phase[ipHYDROGEN] > log10(FLT_MAX) ||
		    dense.gas_phase[ipHYDROGEN] < log10(FLT_MIN) )	
		{
			fprintf(ioQQQ,
				" DISASTER - the log of the entered hydrogen density is %.3f - too extreme for this processor.\n",
				dense.gas_phase[ipHYDROGEN]);
			if( dense.gas_phase[ipHYDROGEN] > 0. )
				fprintf(ioQQQ,
					" DISASTER - the log of the largest hydrogen density this processor can do is %.3f.\n",
					log10(FLT_MAX) );
			else
				fprintf(ioQQQ,
					" DISASTER - the log of the smallest hydrogen density this processor can do is %.3f.\n",
					log10(FLT_MIN) );
			fprintf(ioQQQ," Sorry.\n" );
			
			cdEXIT(EXIT_FAILURE);
		}

		dense.SetGasPhaseDensity( ipHYDROGEN, exp10(dense.gas_phase[ipHYDROGEN]) );
	}

	if( dense.gas_phase[ipHYDROGEN] > MAX_DENSITY )
	{
		fprintf( ioQQQ, "This density is too high.  This version of Cloudy does not permit densities greater than %e cm-3.\n", MAX_DENSITY );
		cdEXIT(EXIT_FAILURE);
	}

	if( dense.gas_phase[ipHYDROGEN] <= 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER Hydrogen density must be > 0.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* this is the linear initial density */
	dense.den0 = dense.gas_phase[ipHYDROGEN];

	/* check if exponent entered */
	dense.DensityPower = (realnum)p.FFmtRead();
	/* this branch when exponent was entered - do something with it */
	if( !p.lgEOL() )
	{
		/*  not constant density
		 *  some sort of power law density distribution */
		if( p.nMatch("COLU") )
		{
			/* density will depend on column density to a power
			 * number entered is col den, convert to scale radius
			 * at this point HDEN is LOG10 of hydrogen density */
			dense.rscale = (realnum)exp10(p.FFmtRead());
			if( p.lgEOL() )
			{
				fprintf( ioQQQ, " The column density MUST be set if the col den option is to be used.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			strcpy( dense.chDenseLaw, "POWC" );
		}
		else if( p.nMatch("DEPT") )
		{
			/* depth option, sets scale radius, log cm */
			dense.rscale = (realnum)exp10(p.FFmtRead());
			if( p.lgEOL() )
			{
				fprintf( ioQQQ, " The scale depth MUST be set if the depth option is to be used.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			strcpy( dense.chDenseLaw, "POWD" );
		}
		else
		{
			/* radius option, will be relative to inner radius */
			strcpy( dense.chDenseLaw, "POWR" );
		}
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		/*  pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)log10(dense.gas_phase[ipHYDROGEN]);
		optimize.vincr[optimize.nparm] = 1.;

		/* these are varios options for density laws, 
		 * first is constant density or pressre*/
		if( strcmp(dense.chDenseLaw ,"CDEN") == 0 || 
			 strcmp(dense.chDenseLaw ,"CPRE") == 0 ||
			 strcmp(dense.chDenseLaw ,"WIND") == 0 ||
			 strcmp(dense.chDenseLaw ,"DYNA") == 0 
			)
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "HDEN=%f LOG" );
			optimize.nvarxt[optimize.nparm] = 1;
		}

		/* power law density distrution */
		else if( strcmp(dense.chDenseLaw,"POWR") == 0 )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "HDEN=%f LOG, power=%f" );
			optimize.nvarxt[optimize.nparm] = 2;
			optimize.vparm[1][optimize.nparm] = dense.DensityPower;
		}

		/* power law with density scale depending on column density */
		else if( strcmp(dense.chDenseLaw,"POWC") == 0 )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "HDEN=%f LOG, power=%f, column=%f" );
			optimize.nvarxt[optimize.nparm] = 3;
			optimize.vparm[1][optimize.nparm] = dense.DensityPower;
			optimize.vparm[2][optimize.nparm] = (realnum)log10(dense.rscale);
		}

		/* power law with density scale depending on depth */
		else if( strcmp(dense.chDenseLaw,"POWD") == 0 )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "HDEN=%f LOG, power=%f, depth=%f" );
			optimize.nvarxt[optimize.nparm] = 3;
			optimize.vparm[1][optimize.nparm] = dense.DensityPower;
			optimize.vparm[2][optimize.nparm] = (realnum)log10(dense.rscale);
		}

		/* could not identify an option */
		else
		{
			fprintf( ioQQQ, " Internal error in HDEN\n" );
			cdEXIT(EXIT_FAILURE);
		}
		++optimize.nparm;
	}
	return;
}
