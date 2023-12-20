/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_collis fill in collisional ionization rates, and resulting cooling */
#include "cddefines.h"
#include "phycon.h"
#include "rfield.h"
#include "heavy.h"
#include "atmdat_adfa.h"
#include "ionbal.h"
#include "dense.h"

void ion_collis(
	/* element number on c scale, H is 0 */
	long int nelem )
{
	long int ion,
		limit;
	double DimaRate, crate;

	DEBUG_ENTRY( "ion_collis()" );

	/* compute collisional ionization rate */

	/* CollidRate[nelem][ion][0] is collisional ionization rate, s-1
	 * CollidRate[nelem][ion][1] is collisional ionization cooling, erg/s
	 */

	/* zero out rates below lowest ionization stage we will consider */
	for( ion=0; ion < (dense.IonLow[nelem] - 1); ion++ )
	{
		ionbal.CollIonRate_Ground[nelem][ion][0] = 0.;
		ionbal.CollIonRate_Ground[nelem][ion][1] = 0.;
	}

	/* chng logic changed to be precisely same as ion_solver */
	/* >>chng 02 nov 08, change 2 to NISO */
	/*limit = MIN2(nelem-2,dense.IonHigh[nelem]-1);*/
	limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);
	ASSERT( limit < LIMELM );

	for( ion=dense.IonLow[nelem]; ion <= limit; ion++ )
	{
		//Get the rate coefficients using either Dima or Hybrid
		DimaRate = t_ADfA::Inst().coll_ion_wrapper( nelem, ion , phycon.te );

		crate = DimaRate*dense.EdenHCorr;

		/* total collisional ionization rate 
		 * with only thermal suprathermal electrons */
		ionbal.CollIonRate_Ground[nelem][ion][0] = crate;

		/* cooling due to collisional ionization, which only includes thermal */
		ionbal.CollIonRate_Ground[nelem][ion][1] = (crate*
			rfield.anu(Heavy.ipHeavy[nelem][ion]-1)* EN1RYD);
	}

	for( ion=dense.IonHigh[nelem]; ion <= nelem; ion++ )
	{
		ionbal.CollIonRate_Ground[nelem][ion][0] = 0.;
		ionbal.CollIonRate_Ground[nelem][ion][1] = 0.;
	}

	/* check not rates are negative - in release mode this loop will optimize out */
	for( ion=0; ion <= nelem; ion++ )
	{
		/* there can be no negative rates */
		ASSERT( ionbal.CollIonRate_Ground[nelem][ion][0] >= 0. );
	}
	return;
}
