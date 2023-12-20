/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseIonPar parse the ionization parameter command */
#include "cddefines.h"
#include "radius.h"
#include "optimize.h"
#include "rfield.h"
#include "input.h"
#include "parser.h"

void ParseIonParI(Parser &p)
{
	ParseIonPar(p,'I');
}
void ParseIonParX(Parser &p)
{
	ParseIonPar(p,'X');
}

void ParseIonPar(Parser &p,
					  char chType)
{

	DEBUG_ENTRY( "ParseIonPar()" );

	/*  check not too many continua */
	if( p.m_nqh >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* this is counter for where to start scanning number on line - different
	 * for XI than for IONIZ */
	/* say that continuum is per unit area, ionization parameter*/
	strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
	if( chType == 'I' )
	{
		/* this is the usual ionization parameter, U */
		strcpy( rfield.chSpNorm[p.m_nqh], "IONI" );
	}
	else if( chType == 'X' )
	{
		/* the X-Ray ionization parameter, xi, defined by Tarter+69,
		 * The explicit energy bounds, 1-1000 Ryd, are quoted in Kallman & Bautista 2001 */
		strcpy( rfield.chSpNorm[p.m_nqh], "IONX" );
		rfield.range[(p.m_nqh)][0] = 1.;
		rfield.range[(p.m_nqh)][1] = 1000.;
	}
	else
	{
		fprintf(ioQQQ," ParseIonPar hit chCard insanity.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* get the ionization parameter*/
	rfield.totpow[p.m_nqh] = p.FFmtRead();
	if( p.lgEOL() )
		p.NoNumb("ionization parameter");

	/* check for linear option, if present take log since rfield.totpow[p.m_nqh]
	 * being log ionization parameter is the default */
	if( p.nMatch( "LINE"  ) )
		rfield.totpow[p.m_nqh] = log10(rfield.totpow[p.m_nqh]);

	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME"  ) )
		rfield.lgTimeVary[p.m_nqh] = true;

	/* vary option */
	if( optimize.lgVarOn )
	{
		if( chType == 'I' )
		{
			/* this is the usual ionization parameter, U */
			strcpy( optimize.chVarFmt[optimize.nparm], "IONIZATION PARAMETER= %f LOG" );
		}
		else if( chType == 'X' )
		{
			/* the X-Ray ionization parameter, xi */
			strcpy( optimize.chVarFmt[optimize.nparm], "XI= %f LOG" );
		}
		else
		{
			fprintf( ioQQQ, " Insanity in detecting which ionization parameter.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		if( rfield.lgTimeVary[p.m_nqh] )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)rfield.totpow[p.m_nqh];
		optimize.vincr[optimize.nparm] = 0.5;
		optimize.nvarxt[optimize.nparm] = 1;
		++optimize.nparm;
	}

	/* increment nmber of specifications of continuum intensities, */
	++p.m_nqh;
	return;
}
