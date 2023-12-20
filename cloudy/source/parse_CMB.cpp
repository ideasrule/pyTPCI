/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCMB parse parameters from fireball command */
#include "cddefines.h"
#include "parser.h"
#include "cosmology.h"
#include "radius.h"
#include "rfield.h"

void ParseCMB(double z, 
  long int *nqh)
{
	double a, 
	  rlogl;

	DEBUG_ENTRY( "ParseCMB()" );

	/* check that stack of shape and luminosity specifications
	 * is parallel, stop if not - this happens is background comes
	 * BETWEEN another set of shape and luminosity commands */
	if( rfield.nShape != *nqh )
	{
		fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* put in a black body */
	strcpy( rfield.chSpType[rfield.nShape], "BLACK" );
	/* >>chng 03 may 23, CMB temp from 2.756 to 2.725 */
	rfield.slope[rfield.nShape] = (CMB_TEMP*(1. + z));
	rfield.cutoff[rfield.nShape][0] = 0.;
	rfield.cutoff[rfield.nShape][1] = 0.;
	strcpy( rfield.chSpNorm[*nqh], "LUMI" );
	a = log10(rfield.slope[rfield.nShape]);
	rlogl = log10(4.*STEFAN_BOLTZ) + 4.*a;
	strcpy( rfield.chRSpec[*nqh], "SQCM" );
	rfield.range[*nqh][0] = rfield.emm();
	rfield.range[*nqh][1] = rfield.egamry();
	rfield.totpow[*nqh] = rlogl;
	/* this is an isotropic radiation field */
	rfield.lgBeamed[*nqh] = false;
	rfield.Illumination[*nqh] = Illuminate::ISOTROPIC;

	++rfield.nShape;
	++*nqh;
	if( *nqh >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* this flag says that CMB has been set */
	rfield.lgCMB_set = true;
	return;
}
