/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAge parse parameters off the age command */
#include "cddefines.h"
#include "timesc.h"
#include "parser.h"

#define NUMBEROF(a) (sizeof(a)/sizeof((a)[0]))

namespace Time {
	const double YEAR=3.15569e7,
		MILLENIUM=YEAR*1000.,
		CENTURY=YEAR*100.,
		MONTH=YEAR/12.,
		FORTNIGHT=(24.*3600.*14.),
		WEEK=(24.*3600.*7.),
		DAY=(24.*3600.),
		HOUR=3600.,
		MINUTE=60.,
		SECOND=1.;

	KeyAction<UnitConverter> TimeUnits[] =
	{ 
		MakeKeyAction("MILL", UnitConverter(MILLENIUM)),
		MakeKeyAction("CENT", UnitConverter(CENTURY)),
		MakeKeyAction("YEAR", UnitConverter(YEAR)),
		MakeKeyAction("MONT", UnitConverter(MONTH)),
		MakeKeyAction("FORT", UnitConverter(FORTNIGHT)),
		MakeKeyAction("WEEK", UnitConverter(WEEK)),
		MakeKeyAction("DAY ", UnitConverter(DAY)),
		MakeKeyAction("HOUR", UnitConverter(HOUR)),
		MakeKeyAction("MINU", UnitConverter(MINUTE)),
		MakeKeyAction("SECO", UnitConverter(SECOND)),		
	};

}

void ParseAge( Parser &p )
{
	DEBUG_ENTRY( "ParseAge()" );

	/* set age for the cloud
	 * various timescales will be checked in AgeCheck, called in comment */

	realnum value = (realnum)p.FFmtRead();

	/* key " off" turns age off */
	if( p.lgEOL() && (!p.nWord(" OFF")) )
	{
		fprintf( ioQQQ, " The age must be on this line.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* check if log of age */
	if( p.nWord(" LOG") )
	{
		value = exp10(value);
	}

	parserProcess(p, Time::TimeUnits, NUMBEROF(Time::TimeUnits), &value);

	timesc.CloudAgeSet = value;

	return;
}
