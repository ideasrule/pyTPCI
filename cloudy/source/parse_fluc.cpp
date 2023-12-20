/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseFluc parse the fluctuations command, which affects either density or abundances */
#include "cddefines.h"
#include "dense.h"
#include "parser.h"
#include "physconst.h"

void ParseFluc(Parser &p )
{
	double flmax, 
	  flmin, 
	  period, 
	  temp;

	DEBUG_ENTRY( "ParseFluc()" );

	/* rapid density fluctuations
	 * first parameter is log of period, 2 is log den max, 3 log Nmin */
	if( p.nMatch("ABUN") )
	{
		/* abundances varied, not density */
		dense.lgDenFlucOn = false;
	}
	else
	{
		/* density is varied */
		dense.lgDenFlucOn = true;
	}

	/* optional keyword COLUMN makes sin over column density rather than radius */
	if( p.nMatch("COLU") )
	{
		/* found key, not fluc over radius, over col den instead */
		dense.lgDenFlucRadius = false;
	}
	else
	{
		/* no key, use default of radius */
		dense.lgDenFlucRadius = true;
	}

	/* 1st number log of period in centimeters */
	period = exp10(p.FFmtRead());
	dense.flong = (realnum)(PI2/period);
	temp = p.FFmtRead();

	/* check size of density - will we crash? */
	if( temp > log10(FLT_MAX) || temp < log10(FLT_MIN) )
	{
		fprintf(ioQQQ,
			" DISASTER - the log of the entered max hydrogen density is %.3f - too extreme for this processor.\n",
			temp);
		if( temp > 0. )
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

	/* 2nd number log of max hydrogen density */
	flmax = exp10(temp);

	temp = p.FFmtRead();

	/* check size of density - will we crash? */
	if( temp > log10(FLT_MAX) || temp < log10(FLT_MIN) )
	{
		fprintf(ioQQQ,
			" DISASTER - the log of the entered min hydrogen density is %.3f - too extreme for this processor.\n",
			temp);
		if( temp > 0. )
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

	/* 3rd number log of min hydrogen density */
	flmin = exp10(temp);

	if( flmax/flmin > 100. )
	{
		fprintf( ioQQQ, "This range of density probably will not work.\n" );
	}
	if( flmax > 1e15 )
	{
		fprintf( ioQQQ, "These parameters look funny to me.  Please check Hazy.\n" );
	}
	if( p.lgEOL() || (flmin > flmax) )
	{
		fprintf( ioQQQ, "There MUST be three numbers on this line.\n" );
		fprintf( ioQQQ, "These must be the period(cm), max, min densities, all logs, in that order.\n" );
		if( flmin > flmax )
			fprintf( ioQQQ, "The max density must be greater or equal than the min density.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* this is optional phase shift for the command */
	dense.flcPhase = (realnum)p.FFmtRead();

	/* FacAbunSav = (cfirst * COS( depth*flong+flcPhase ) + csecnd) */
	dense.cfirst = (realnum)((flmax - flmin)/2.);
	dense.csecnd = (realnum)((flmax + flmin)/2.);
	/* these will be added together with the first mult by sin - which goes to
	 * -1 - must not have a negative density */
	ASSERT( dense.cfirst < dense.csecnd );
	/* >>chng 96 jul 13 moved depset to SetAbundances fac
	 * if( lgDenFlucOn ) then
	 * this is a pressure law
	 * chCPres = 'SINE'
	 * else
	 * this is the metallicity of the gas
	 * do i=3,limelm
	 * depset(i) = flmax
	 * end do
	 * endif
	 *
	 * now get density if this is density option (not abundances) */
	if( dense.lgDenFlucOn )
	{
		strcpy( dense.chDenseLaw, "SINE" );

		if( dense.gas_phase[ipHYDROGEN] > 0. )
		{
			fprintf( ioQQQ, " PROBLEM DISASTER More than one density command was entered.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* depth is zero for first zone */
		dense.SetGasPhaseDensity( ipHYDROGEN, dense.cfirst*(realnum)cos(dense.flcPhase) + dense.csecnd );

		if( dense.gas_phase[ipHYDROGEN] <= 0. )
		{
			fprintf( ioQQQ, " PROBLEM DISASTER Hydrogen density must be > 0.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	return;
}
