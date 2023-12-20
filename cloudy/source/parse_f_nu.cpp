/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseF_nu parse intensity command parameters */
#include "cddefines.h"
#include "rfield.h"
#include "radius.h"
#include "parser.h"

void ParseF_nu(
  Parser &p, 
  const char *chType, 
  bool lgNU2)
{
	double fr;

	DEBUG_ENTRY( "ParseF_nu()" );

	/* flux density of this continuum source, at optional frequency */

	strcpy( rfield.chRSpec[p.m_nqh], chType );

	rfield.totpow[p.m_nqh] = p.FFmtRead();

	/* large luminosity but per sq cm */
	if( rfield.totpow[p.m_nqh] > 37. && 
		strcmp(rfield.chRSpec[p.m_nqh],"SQCM") == 0 )
	{
		fprintf( ioQQQ, " This intensity is VERY large.  Problems?  Was luminosity intended??\n" );
	}

	if( p.lgEOL() )
	{
		p.NoNumb("flux density");
	}

	strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

	/* this is optional frequency in Rydbergs */
	rfield.range[p.m_nqh][0] = p.FFmtRead();

	/* >>chng 96 dec 18, was 1 changed to H mass nuc ryd
	 * if( lgEOL ) range(nqh,1) = 1. */
	if( p.lgEOL() )
	{
		rfield.range[p.m_nqh][0] = HIONPOT;
	}

	if( rfield.range[p.m_nqh][0] <= 0. )
	{
		rfield.range[p.m_nqh][0] = exp10(rfield.range[p.m_nqh][0]);
	}

	if( lgNU2 )
	{
		/* range is now freq in ryd, totpow is log of product nu*f_nu */
		fr = log10(rfield.range[p.m_nqh][0]*FR1RYD);
		rfield.totpow[p.m_nqh] -= fr;
	}

	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME"  ) )
		rfield.lgTimeVary[p.m_nqh] = true;

	++p.m_nqh;
	if( p.m_nqh >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
