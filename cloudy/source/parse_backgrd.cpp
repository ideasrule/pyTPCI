/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseBackgrd parse options for the BACKGROUND command - this actually enters two continua*/
#include "cddefines.h"
#include "cosmology.h"
#include "radius.h"
#include "rfield.h"
#include "parser.h"

void ParseBackgrd(Parser &p)
{
	double a, 
	  fac, 
	  rlogl, 
	  z;

	DEBUG_ENTRY( "ParseBackgrd()" );

	/* check that stack of shape and luminosity specifications
	 * is parallel, stop if not - this happens is background comes
	 * BETWEEN another set of shape and luminosity commands */
	if( rfield.nShape != p.m_nqh )
	{
		fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* diffuse x-ray background from 
	 * >>refer	cosmic	background	Ostriker and Ikeuchi ApJL 268, L63.
	 * parameter on cmnd is redshift, followed by optional J21 
	 * >>refer	cosmic	background	Ikeuchi, S.; Ostriker, J. P., 1986, ApJ 301, 522 */

	/* >>chng 02 aug 12, do not try to enter bare powerlaw, use table power law default instead */
	{
		// There must be a neater way than this to extract the
		// functionality from ParseTable...
		p.set_point(0);
		string cmdSave = p.getRawTail();
		p.setline("TABLE POWERLAW  ");
		ParseTable( p );
		p.setline(cmdSave.c_str());
		p.set_point(4);
	}

	/* now generate equivalent of luminosity command
	 * enter phi(h), the number of h-ionizing photons/cm2 */

	strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
	strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );
	/* this is an isotropic radiation field */
	rfield.lgBeamed[p.m_nqh] = false;
	rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;

	/* this will be the redshift */
	z = p.FFmtRead();
	if( p.lgEOL() )
		z = 0.;

	/* optional scale factor */
	fac = p.FFmtRead();
	if( p.lgEOL() )
		fac = 1.0;

	/* this equation from Ostriker and Ikeuchi Ap.J.L 268, L63.
	 * this should be J_nu into 4\pi sr
	 * >>chng 96 jul 23, from ostriker to vedel, et al mnras 271, 743. */
	rfield.totpow[p.m_nqh] = (log10(PI4*fac*1.e-21/
		(1.+powi(5./(1.+z),4))));

	/* this is at 1 ryd for H
	 * range(nqh,1) = 1.
	 *>>chgn 96 dec 18, to H ioniz pot from 1 ryd */
	rfield.range[p.m_nqh][0] = HIONPOT;

	++p.m_nqh;
	if( p.m_nqh >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* add fireball unless NO FIREBALL or NO CMP in */
	if( !(p.nMatch("NO FI") || p.nMatch("NO CM") ) )
	{

		/* put in a black body */
		strcpy( rfield.chSpType[rfield.nShape], "BLACK" );
		/* >>chng 03 may 23, CMB temp from 2.756 to 2.725 */
		rfield.slope[rfield.nShape] = (CMB_TEMP*(1. + z));
		rfield.cutoff[rfield.nShape][0] = 0.;
		rfield.cutoff[rfield.nShape][1] = 0.;
		strcpy( rfield.chSpNorm[p.m_nqh], "LUMI" );
		a = log10(rfield.slope[rfield.nShape]);
		rlogl = log10(4.*STEFAN_BOLTZ) + 4.*a;
		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		rfield.range[p.m_nqh][0] = rfield.emm();
		rfield.range[p.m_nqh][1] = rfield.egamry();
		rfield.totpow[p.m_nqh] = rlogl;

		/* this flag says that CMB has been set */
		rfield.lgCMB_set = true;

		++rfield.nShape;
		++p.m_nqh;
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
}
