/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseSphere parse parameters on sphere command */
#include "cddefines.h"
#include "geometry.h"
#include "opacity.h"
#include "parser.h"

void ParseSphere(Parser &p )
{
	double fac;

	DEBUG_ENTRY( "ParseSphere()" );

	/* compute a spherical model, diffuse field from other side included */
	geometry.lgSphere = true;

	/* turn off electron scattering opacity */
	opac.lgScatON = false;

	/* if "STATIC" is specified then set to case B for H lines */
	if( p.nMatch("STAT") )
	{
		geometry.lgStatic = true;
		//opac.tlamin = 1e5;
		/* this is option to not check on iterations, used for debugging H atom */
		if( p.nMatch("(OK)") )
		{
			geometry.lgStaticNoIt = true;
		}
	}

	/* set the covering facto to full coverage */
	geometry.covgeo = 1.;
	geometry.covrt = 1.;

	/* check for optional covering factor, which is no longer parsed here */
	fac = p.FFmtRead();

	if( !p.lgEOL() )
	{
		/* >>chng 01 jul 16, remove covering factor, only on covering factor command */
		fprintf(ioQQQ," The number %g appeared on the SPHERE command.\n", fac);
		fprintf(ioQQQ," The covering factor can no longer be set with the SPHERE command.\n");
		fprintf(ioQQQ," The number has been ignored.\n");
	}

	/* if the "BEAM" or "SLIT" option is specified then only part 
	 * of the sphere is observed, and intensities
	 * should not be increased by r^2.  There are two limiting cases, SLIT in which
	 * the slit is longer than the diameter of the nebula and the contibution to the
	 * detected luminosity goes as r^1, and BEAM when the contribution is r^0, 
	 * or same as plane parallel */
	if( p.nMatch("SLIT") || p.nMatch("BEAM") )
	{
		/* >>chng 01 jul 16, remove options, put on APERTURE command */
		fprintf(ioQQQ," The SLIT and BEAM options are now part of the APERTURE command.\n");
		fprintf(ioQQQ," The syntax is the same.\n");
		fprintf(ioQQQ," This option has been ignored.\n");
	}
	return;
}
