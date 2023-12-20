/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseRatio derive continuum luminosity of this continuum relative to previous */
#include "cddefines.h"
#include "optimize.h"
#include "input.h"
#include "rfield.h"
#include "parser.h"

void ParseRatio(Parser &p)
{
	bool lgAoxOn;
	double aox;

	DEBUG_ENTRY( "ParseRatio()" );

	/* enter a continuum luminosity as a ratio of
	 * nuFnu for this continuum relative to a previous continuum
	 * format; first number is ratio of second to first continuum
	 * second number is energy for this ratio
	 * if third numbewr on line, then 2nd number is energy of
	 * first continuum, while 3rd number is energy of second continuum */

	if( p.m_nqh == 0 )
	{
		fprintf( ioQQQ, " Can\'t form ratio since this is first continuum.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* per sq cm or luminosity is really irrelevant */
	strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
	strcpy( rfield.chSpNorm[p.m_nqh], "RATI" );

	/* this option is to specify alpha ox */
	if( p.nMatch("ALPH") )
	{
		/* lgAoxOn is flag saying that we will spicify alpha ox */
		lgAoxOn = true;
		/* only one parameter to be recognized, alpha ox */
		aox = p.FFmtRead();

		/* 403.3 is ratio of energies where alpha ox defined,
		 * assumed to be 2500A and 2keV */
		rfield.totpow[p.m_nqh] = pow(403.3,aox);
		rfield.range[p.m_nqh][0] = 0.3645;
		rfield.range[p.m_nqh][1] = 147.;
	}

	else
	{
		/* set flag saying that alpha ox will not be specified */
		lgAoxOn = false;
		/* set this to impossible number since not used, but lint needs a value */
		aox = -DBL_MAX;
		/* specify ratio, two energies */
		rfield.totpow[p.m_nqh] = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("continuum ratio");

		/* assumed to be a log if negative or zero */
		if( rfield.totpow[p.m_nqh] <= 0. || p.nMatch(" LOG" ) )
		{
			rfield.totpow[p.m_nqh] = exp10(rfield.totpow[p.m_nqh]);
		}

		rfield.range[p.m_nqh][0] = p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There must be at least 2 numbers on this line.\n" );
			fprintf( ioQQQ, " The ratio, and one or two energies\n" );
			cdEXIT(EXIT_FAILURE);
		}

		rfield.range[p.m_nqh][1] = p.FFmtRead();
		/* if only one number then assume same for each */
		if( p.lgEOL() )
			rfield.range[p.m_nqh][1] = rfield.range[p.m_nqh][0];

		if( rfield.range[p.m_nqh][0] < rfield.emm() || 
		    rfield.range[p.m_nqh][1] < rfield.emm() )
		{
			fprintf( ioQQQ, " One of the energies is too low, outside the range of the code.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		if( rfield.range[p.m_nqh][0] > rfield.egamry() ||
		    rfield.range[p.m_nqh][1] > rfield.egamry() )
		{
			fprintf( ioQQQ, " One of the energies is too high, outside the range of the code.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME"  ) )
		rfield.lgTimeVary[p.m_nqh] = true;

	/* vary option */
	if( optimize.lgVarOn )
	{
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		if( lgAoxOn )
		{
			/* this is the number of parameters to feed onto the input line */
			optimize.nvarxt[optimize.nparm] = 1;
			/* vary alpha ox */
			strcpy( optimize.chVarFmt[optimize.nparm], "RATIO alphox %f" );
			optimize.lgOptimizeAsLinear[optimize.nparm] = true;
			/* param is linear scale factor */
			optimize.vparm[0][optimize.nparm] = (realnum)aox;
			optimize.vincr[optimize.nparm] = 0.2f;
		}
		else
		{
			/* this is the number of parameters to feed onto the input line */
			optimize.nvarxt[optimize.nparm] = 3;
			strcpy( optimize.chVarFmt[optimize.nparm], "RATIO LOG %f %f %f" );
			/* param is log of abundance by number relative to hydrogen */
			optimize.vparm[0][optimize.nparm] = (realnum)log10(rfield.totpow[p.m_nqh]);
			optimize.vparm[1][optimize.nparm] = (realnum)rfield.range[p.m_nqh][0];
			optimize.vparm[2][optimize.nparm] = (realnum)rfield.range[p.m_nqh][1];
			optimize.vincr[optimize.nparm] = 0.2f;
		}
		if( rfield.lgTimeVary[p.m_nqh] )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		++optimize.nparm;
	}

	++p.m_nqh;
	if( p.m_nqh >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}
