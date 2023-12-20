/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseMap parse map command to produce map of heating and cooling,
 * map is produced by calling punt(" map") */
#include "cddefines.h"
#include "hcmap.h"
#include "parser.h"

void ParseMap(Parser &p )
{
	DEBUG_ENTRY( "ParseMap()" );

	/* say output goes to stdout */
	ioMAP = ( ioQQQ == NULL ) ? stdout : ioQQQ;

	/* do cooling space map for specified zones
	 * if no number, or <0, do map and punch out without doing first zone */
	hcmap.MapZone = (long)p.FFmtRead();
	if( p.lgEOL() )
	{
		hcmap.MapZone = 0;
		return;
	}

	if( p.nMatch("RANG") )
	{
		bool lgLogOn;
		hcmap.RangeMap[0] = (realnum)p.FFmtRead();
		if( hcmap.RangeMap[0] <= 10. )
		{
			hcmap.RangeMap[0] = exp10(hcmap.RangeMap[0]);
			lgLogOn = true;
		}
		else
		{
			lgLogOn = false;
		}
		hcmap.RangeMap[1] = (realnum)p.FFmtRead();
		if( lgLogOn )
			hcmap.RangeMap[1] = exp10(hcmap.RangeMap[1]);

		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There must be a zone number, followed by two temperatures, on this line.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		if( hcmap.RangeMap[1] <= hcmap.RangeMap[0] )
		{
			fprintf( ioQQQ, " The upper temperature limit must be larger than the lower: "
				 "found lower=%g upper=%g.\n", hcmap.RangeMap[0], hcmap.RangeMap[1] );
			cdEXIT(EXIT_FAILURE);
		}
	}
}
