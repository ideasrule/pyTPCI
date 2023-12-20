/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAbsMag parse the absolute magnitude command */
#include "cddefines.h"
#include "rfield.h"
#include "called.h"
#include "parser.h"

void ParseAbsMag(Parser &p)
{
	DEBUG_ENTRY( "ParseAbsMag()" );

	/* enter luminosity in absolute magnitudes */
	strcpy( rfield.chRSpec[p.m_nqh], "4 PI" );
	rfield.totpow[p.m_nqh] = p.FFmtRead();
	if( p.lgEOL() )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " There should have been a number on this line.   Sorry.\n" );
		}
		cdEXIT(EXIT_FAILURE);
	}
	if( p.nMatch("BOLO") )
	{
		strcpy( rfield.chSpNorm[p.m_nqh], "LUMI" );
		rfield.range[p.m_nqh][0] = rfield.emm();
		rfield.range[p.m_nqh][1] = rfield.egamry();
		/* page 197 allen 76 */
		rfield.totpow[p.m_nqh] = ((4.75 - rfield.totpow[p.m_nqh])/
		  2.5 + 33.5827);
	}
	else if( p.nMatch("VISU") )
	{
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );
		/* this is 5550A, the center of the V filter */
		rfield.range[p.m_nqh][0] = 0.164f;
		/* totpow(nqh) = (4.79 - totpow(nqh)) / 2.5  + 18.758
		 * page 197, allen 76, 3rd line from bottom */
		rfield.totpow[p.m_nqh] = (-rfield.totpow[p.m_nqh]/2.5 + 
		  20.65296);
	}
	else
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " Keyword BOLOmetric or VISUal must appear.\n" );
		}
		cdEXIT(EXIT_FAILURE);
	}

	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME"  ) )
		rfield.lgTimeVary[p.m_nqh] = true;

	++p.m_nqh;
	if( p.m_nqh >= LIMSPC )
	{
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
