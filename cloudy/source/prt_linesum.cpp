/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtLineSum parse print line sum command to enter set of lines into sum  */
#include "cddefines.h"
#include "cddrive.h"
#include "radius.h"
#include "lines.h"
#include "parser.h"
#include "prt.h"

static vector<LineID> lineids;
static vector<long> ipLine;

void ParseLineList(Parser &p, vector<LineID>& lines)
{
	DEBUG_ENTRY( "ParseLineList()" );

	/* now read in lines */
	lines.clear();
	while( true )
	{
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
		if( p.hasCommand("END") )
			break;

		lines.emplace_back(p.getLineID());
		if( !p.lgReachedEnd() )
		{
			fprintf( ioQQQ, "ParsePrtLineSum: found junk at end of input line:\n" );
			p.showLocation();
			cdEXIT(EXIT_FAILURE);
		}
	}
}

void ParsePrtLineSum(Parser &p)
{
	DEBUG_ENTRY( "ParsePrtLineSum()" );

	ParseLineList(p, lineids);
	ipLine.resize(lineids.size());
}

double PrtLineSum()
{
	DEBUG_ENTRY( "PrtLineSum()" );

	double sum = 0.;
	/* this can be called during setup mode, in which case we do nothing */
	if( LineSave.ipass < 0 )
		return sum;
	
	if( LineSave.ipass == 0 )
	{
		bool lgFail = false;
		for( size_t i=0; i < lineids.size(); i++ )
		{
			/* save the array index for each line */
			if( (ipLine[i] = LineSave.findline(lineids[i])) <= 0 )
			{
				fprintf( ioQQQ, " PrtLineSum could not find line " );
				prt_line_err( ioQQQ, lineids[i] );
				lgFail = true;
			}
		}
		if( lgFail )
			cdEXIT(EXIT_FAILURE);
		return sum;
	}
	
	/* now sum the line */
	for( size_t i=0; i < lineids.size(); i++ )
	{
		/* this version of chLine uses index, does not search*/
		double absint, relint;
		cdLine_ip(ipLine[i], &relint, &absint);
		absint /= radius.Conv2PrtInten;
		sum += absint;
	}
	return sum;
}

