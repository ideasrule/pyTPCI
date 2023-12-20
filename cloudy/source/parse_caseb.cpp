/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCaseB - parse the Case A, B, or C command */
#include "cddefines.h"
#include "opacity.h"
#include "parser.h"

/*ParseCaseB - parse the Case A, B, or C command */
void ParseCaseB( Parser &p )
{
	DEBUG_ENTRY( "ParseCaseB()" );

	/* which case do we do? */
	if( p.nMatch(" A ") )
	{
		opac.lgCaseB = false;
//		rfield.lgInducProcess = false;
	}
	else if( p.nMatch(" B ") )
	{
		opac.lgCaseB = true;
//		rfield.lgInducProcess = false;
	}
	else if( p.nMatch(" C ") )
		opac.lgCaseB = false;

	/* difference between Case A and Case C is in induced processes -
	 * do not include them in Case A or Case B */
	/* optional optical depth in Ly a */
	opac.tlamin = (realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		/* optical depths not specified - use defaults for cases A and B */ 
		if( opac.lgCaseB )
			/* set default to 1e5 to get more realistic conditions in H+ region.
			 * Very large tau caused extreme Lya trapping, photoionization rates, 
			 * & radiation pressure*/
			opac.tlamin = 1e5f;
		else
			/* Case A or Case C - Lyman lines optically thin */
			opac.tlamin = 1e-5f;
	}
	else
		opac.tlamin = exp10(opac.tlamin);

	/* Hummer and Storey case B; no collisions from n=1, 2 (usually in) */
	if( p.nMatch("HUMM") )
		opac.lgCaseB_HummerStorey = true;

	/* the NO PHOTOIONIZATION option, turns off excited state photoionization */
	if( p.nMatch("NO PH") )
		opac.lgCaseB_no_photo = true;

	/* the NO PDEST option, turns off line destruction by background opacities */
	if( p.nMatch("NO PDE") )
		opac.lgCaseB_no_pdest = true;

	return;
}
