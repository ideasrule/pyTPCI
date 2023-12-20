/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_setup convert level 1 and level 2 line parameters and pointers into internal form 
 * used by code, line data were read in by atmdat_readin */
#include "cddefines.h"
#include "lines.h"
#include "taulines.h"
#include "opacity.h"
#include "lines_service.h"

void lines_setup(void)
{
	static bool lgFirst = true;

	DEBUG_ENTRY( "lines_setup()" );

	/* only do this one time, and only if number of atom_level2 lines is positive */
	if( lgFirst && nWindLine>0)
	{ 
		lgFirst = false;
		/* these are the massive set of op lines, with g-bar approx cs
		 * confirm that input data are valid */

		for( long i=0; i < nWindLine; i++ )
		{
			/* this information was read in in createdata */
			ASSERT( (*TauLine2[i].Hi()).nelem() > 0 );
			ASSERT( (*TauLine2[i].Hi()).nelem() <= (int)LIMELM );

			ASSERT( (*TauLine2[i].Hi()).IonStg() > 0 );
			ASSERT( (*TauLine2[i].Hi()).IonStg() <= (int)LIMELM );

			ASSERT( (*TauLine2[i].Lo()).g() >0. );

			ASSERT( (*TauLine2[i].Hi()).g() > 0. );

			/* check that energy is positive*/
			ASSERT( TauLine2[i].EnergyWN() > 0 );

			/* TauLine2[i].Emis().gf() this is gf if positive, A if negative */
			/* test whether a or gf entered, convert A to gf */
			if( TauLine2[i].Emis().gf() < 0. )
			{
				/* convert A (=-gf) into real gf */
				TauLine2[i].Emis().gf() *= (realnum)((-(*TauLine2[i].Hi()).g())/TRANS_PROB_CONST/POW2(TauLine2[i].EnergyWN()));
			}

			/*now put into standard format */
			TauLine2[i].WLAng() = wn2ang( double( TauLine2[i].EnergyWN() ) );
			(*TauLine2[i].Lo()).Pop() = 0.;
			(*TauLine2[i].Hi()).Pop() = 0.;
			TauLine2[i].Emis().iRedisFun() = ipPRD;

			/* these are line optical depth arrays
			 * inward optical depth */
			TauLine2[i].Emis().TauIn() = opac.taumin;
			TauLine2[i].Emis().TauCon() = opac.taumin;
			TauLine2[i].Emis().ColOvTot() = 0.;
			/* outward optical depth */
			TauLine2[i].Emis().TauTot() = 1e20f;
			/* escape probability */
			TauLine2[i].Emis().Pesc() = 1.;
			/* inward part of line */
			TauLine2[i].Emis().FracInwd() = 1.;
			/* destruction probability */
			TauLine2[i].Emis().Pdest() = 0.;
			TauLine2[i].Emis().Pelec_esc() = 0.;
			/* line pumping rate */
			TauLine2[i].Emis().pump() = 0.;
			/* population of lower level */
			(*TauLine2[i].Lo()).Pop() = 0.;
			/* population of upper level */
			(*TauLine2[i].Hi()).Pop() = 0.;
			/* population of lower level with correction for stim emission */
			TauLine2[i].Emis().PopOpc() = 0.;
			/* following two heat exchange excitation, deexcitation */
			TauLine2[i].Coll().cool() = 0.;
			TauLine2[i].Coll().heat() = 0.;
			/* intensity of line */
			TauLine2[i].Emis().xIntensity() = 0.;
			TauLine2[i].Emis().xObsIntensity() = 0.;
			/* ots rate */
			TauLine2[i].Emis().ots() = 0.;
		}
	}

	for( size_t i=0; i < UTALines.size(); i++ )
	{
		/* this information was read in in createdata */
		ASSERT( (*UTALines[i].Hi()).nelem() > 0 );
		ASSERT( (*UTALines[i].Hi()).nelem() <= (int)LIMELM );

		ASSERT( (*UTALines[i].Hi()).IonStg() > 0 );
		ASSERT( (*UTALines[i].Hi()).IonStg() <= (int)LIMELM );

		ASSERT( (*UTALines[i].Lo()).g() > 0. );

		ASSERT( (*UTALines[i].Hi()).g() > 0. );

		/* check that energy is positive*/
		ASSERT( UTALines[i].EnergyWN() > 0 );

		(*UTALines[i].Lo()).Pop() = 0.;
		(*UTALines[i].Hi()).Pop() = 0.;
		UTALines[i].Emis().iRedisFun() = ipPRD;

		/* these are line optical depth arrays
		 * inward optical depth */
		UTALines[i].Emis().TauIn() = opac.taumin;
		UTALines[i].Emis().TauCon() = opac.taumin;
		UTALines[i].Emis().ColOvTot() = 0.;
		/* outward optical depth */
		UTALines[i].Emis().TauTot() = 1e20f;
		/* escape probability */
		UTALines[i].Emis().Pesc() = 1.;
		/* inward part of line */
		UTALines[i].Emis().FracInwd() = 1.;
		/* destruction probability */
		UTALines[i].Emis().Pdest() = 0.;
		UTALines[i].Emis().Pelec_esc() = 0.;
		/* line pumping rate */
		UTALines[i].Emis().pump() = 0.;
		/* population of lower level */
		(*UTALines[i].Lo()).Pop() = 0.;
		/* population of upper level */
		(*UTALines[i].Hi()).Pop() = 0.;
		/* population of lower level with correction for stim emission */
		UTALines[i].Emis().PopOpc() = 0.;
		/* following two heat exchange excitation, deexcitation */
		UTALines[i].Coll().cool() = 0.;
		/* heat is the net heat per pump and was set when data read in
		 * this is different from other lines with this structure 
		UTALines[i].Coll().heat() = 0.;*/
		/* intensity of line */
		UTALines[i].Emis().xIntensity() = 0.;
		UTALines[i].Emis().xObsIntensity() = 0.;
		UTALines[i].Emis().ots() = 0.;
	}

	for( size_t i=0; i < HFLines.size(); i++ )
	{
		/* this information was read in in createdata */
		ASSERT( (*HFLines[i].Hi()).nelem() > 0 );
		ASSERT( (*HFLines[i].Hi()).nelem() <= (int)LIMELM );

		ASSERT( (*HFLines[i].Hi()).IonStg() > 0 );
		ASSERT( (*HFLines[i].Hi()).IonStg() <= (int)LIMELM );

		ASSERT( (*HFLines[i].Lo()).g() > 0. );

		ASSERT( (*HFLines[i].Hi()).g() > 0. );

		/* check that energy is positive*/
		ASSERT( HFLines[i].EnergyWN() > 0 );
		ASSERT( HFLines[i].Emis().Aul()>0 );
		ASSERT( HFLines[i].Emis().damp()>0 );

		/* HFLines[i].Emis->gf() this is gf if positive, A if negative */
		/* test whether a or gf entered, convert A to gf */
		if( HFLines[i].Emis().gf() < 0. )
		{
			/* convert A (=-gf) into real gf */
			HFLines[i].Emis().gf() *= (realnum)(-(*HFLines[i].Hi()).g()/TRANS_PROB_CONST/POW2(HFLines[i].EnergyWN()));
		}

		/*now put into standard format */
		HFLines[i].WLAng() = 1.e8f/HFLines[i].EnergyWN();
		(*HFLines[i].Lo()).Pop() = 0.;
		(*HFLines[i].Hi()).Pop() = 0.;
		/* change from partial to complete redistribution */
		HFLines[i].Emis().iRedisFun() = ipCRD;

		/* these are line optical depth arrays
		 * inward optical depth */
		HFLines[i].Emis().TauIn() = opac.taumin;
		HFLines[i].Emis().TauCon() = opac.taumin;
		HFLines[i].Emis().ColOvTot()=0;
		/* outward optical depth */
		HFLines[i].Emis().TauTot() = 1e20f;
		/* escape probability */
		HFLines[i].Emis().Pesc() = 1.;
		/* inward part of line */
		HFLines[i].Emis().FracInwd() = 1.;
		/* destruction probability */
		HFLines[i].Emis().Pdest() = 0.;
		HFLines[i].Emis().Pelec_esc() = 0.;
		/* line pumping rate */
		HFLines[i].Emis().pump() = 0.;
		/* population of lower level */
		(*HFLines[i].Lo()).Pop() = 0.;
		/* population of upper level */
		(*HFLines[i].Hi()).Pop() = 0.;
		/* population of lower level with correction for stim emission */
		HFLines[i].Emis().PopOpc() = 0.;
		/* following two heat exchange excitation, deexcitation */
		HFLines[i].Coll().cool() = 0.;
		HFLines[i].Coll().heat() = 0.;
		/* intensity of line */
		HFLines[i].Emis().xIntensity() = 0.;
		HFLines[i].Emis().xObsIntensity() = 0.;
		HFLines[i].Emis().ots() = 0.;
	}
	return;
}
