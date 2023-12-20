/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "emission.h"
#include "opacity.h"

/*EmLineJunk set all elements of transition struc to dangerous values */
/*EmLineZero set all elements of transition struc to zero */

/*EmLineJunk set all elements of transition struc to dangerous values */
void EmLineJunk( EmissionList::reference t )
{

	DEBUG_ENTRY( "EmLineJunk()" );

	/* optical depth in continuum to ill face */
	t.TauCon() = -FLT_MAX;

	/* inward and total line optical depths */
	t.TauIn() = -FLT_MAX;
	t.TauInSpecific() = -FLT_MAX;
	t.TauTot() = -FLT_MAX;

	/* type of redistribution function, */
	t.iRedisFun() = INT_MIN;

	/* array offset for line within fine opacity */
	t.ipFine() = -10000;

	/* inward fraction */
	t.FracInwd() = -FLT_MAX;

	/* continuum pumping rate */
	t.pump() = -FLT_MAX;

	/* line intensity */
	t.xIntensity() = -FLT_MAX;
	t.xObsIntensity() = -FLT_MAX;

	/* gf value */
	t.gf() = -FLT_MAX;

	/* escape and destruction probs */
	t.Pesc() = -FLT_MAX;
	t.Pdest() = -FLT_MAX;
	t.Pelec_esc() = -FLT_MAX;

	/* damping constant */
	t.damp() = -FLT_MAX;

	/* ratio of collisional to radiative excitation*/
	t.ColOvTot() = -FLT_MAX;

	/* auto-ionization fraction */
	t.AutoIonizFrac() = -FLT_MAX;

	/* line opacity */
	t.opacity() = -FLT_MAX;
	t.mult_opac() = -FLT_MAX;

	t.PopOpc() = -FLT_MAX;

	t.VoigtLineCen() = -FLT_MAX;

	/* transition prob, Einstein A upper to lower */
	t.Aul() = -FLT_MAX;

	/* ots rate */
	t.ots() = -FLT_MAX;
	return;
}

/*TauZero set initial values of inward, outward, and local-to-continuum-source optical depths */
void TauZero( EmissionList::reference t )
{
	DEBUG_ENTRY( "TauZero()" );

	/* total optical depth in all overlapping lines to illuminated face,
	 * used for pumping */
	t.TauCon() = opac.taumin;

	/* inward and total line optical depths */
	/* >>chng 03 feb 14, from 0 to opac.taumin */
	t.TauIn() = opac.taumin;
	t.TauInSpecific() = opac.taumin;

	/* total optical depths */
	t.TauTot() = 1e20f;
	
	return;
}

/*EmLineZero zeros out the emission line structure */
void EmLineZero( EmissionList::reference t )
{
	DEBUG_ENTRY( "EmLineZero()" );

	/* inward fraction */
	/* >>chng 03 feb 14, from 0 to 1 */
	t.FracInwd() = 1.;

	/* continuum pumping rate */
	t.pump() = 0.;

	/* line intensity */
	t.xIntensity() = 0.;
	t.xObsIntensity() = 0.;

	/* escape and destruction probs */
	/* >>chng 03 feb 14, change from 0 to 1 */
	t.Pesc() = 1.;
	t.Pdest() = 0.;
	t.Pelec_esc() = 0.;

	/* ratio of collisional to radiative excitation*/
	t.ColOvTot() = 1.;

	t.VoigtLineCen() = 1.;

	/* pop that enters net opacity */
	t.PopOpc() = 0.;
	t.mult_opac() = 0.;

	/* ots rate */
	t.ots() = 0.;
	return;
}
