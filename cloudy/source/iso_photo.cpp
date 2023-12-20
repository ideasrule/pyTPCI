/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_photo do photoionization rates for element nelem on the ipISO isoelectronic sequence */
#include "cddefines.h"
#include "hydrogenic.h"
#include "rfield.h"
#include "opacity.h"
#include "trace.h"
#include "ionbal.h"
#include "thermal.h"
#include "iso.h"
#include "gammas.h"
#include "freebound.h"
#include "cosmology.h"

void iso_photo(
	/* iso sequence, hydrogen or helium for now */
	long ipISO , 
	/* the chemical element, 0 for hydrogen */
	long int nelem)
{
	long int limit ,
		n;

	t_phoHeat photoHeat;

	DEBUG_ENTRY( "iso_photo()" );

	/* check that we were called with valid charge */
	ASSERT( nelem >= 0 && nelem < LIMELM );
	ASSERT( ipISO < NISO );

	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/* do photoionization rates */
	/* induced recombination; FINDUC is integral of
	 * pho rate times EXP(-hn/kt) for induc rec
	 * CIND is this times hnu-hnu0 to get ind rec cooling
	 * ionbal.lgPhotoIoniz_On is 1, set to 0 with 'no photoionization' command
	 * ipSecIon points to 7.353 Ryd, lowest energy where secondary ioniz
	 * of hydrogen is possible */

	// photoionization of ground, this is different from excited states because 
	// there will eventually be more than one shell, (when li-like done)
	sp->fb[0].gamnc = GammaBn(sp->fb[0].ipIsoLevNIonCon,
		rfield.nflux,
		sp->fb[0].ipOpac,
		sp->fb[0].xIsoLevNIonRyd,
		&sp->fb[0].RecomInducRate,
		&sp->fb[0].RecomInducCool_Coef,
		&photoHeat)*
		ionbal.lgPhotoIoniz_On;

	if( cosmology.lgDo )
	{
		photoHeat.HeatNet = 0.;
		photoHeat.HeatLowEnr = 0.;
		photoHeat.HeatHiEnr = 0.;
		sp->fb[0].gamnc = 0.;
		sp->fb[0].RecomInducRate = 0.;
		sp->fb[0].RecomInducCool_Coef = 0.;
	}
	
	/* heating due to photo of ground */
	sp->fb[0].PhotoHeat = photoHeat.HeatNet*ionbal.lgPhotoIoniz_On;

	/* save these rates into ground photo rate vector */
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][0] = sp->fb[ipH1s].gamnc;
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][1] = photoHeat.HeatLowEnr*ionbal.lgPhotoIoniz_On;
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][2] = photoHeat.HeatHiEnr*ionbal.lgPhotoIoniz_On;

	/* CompRecoilIonRate is direct photioniz rate due to 
	 * bound compton scattering of very hard x-rays+Compton scat */
	/* now add in compton recoil, to ground state, save heating as high energy heat */
	ASSERT( ionbal.CompRecoilIonRate[nelem][nelem-ipISO]>=0. &&
		ionbal.CompRecoilHeatRate[nelem][nelem-ipISO]>= 0. );
	sp->fb[0].gamnc += ionbal.CompRecoilIonRate[nelem][nelem-ipISO];
	sp->fb[0].PhotoHeat += ionbal.CompRecoilHeatRate[nelem][nelem-ipISO];

	/* now add bound compton scattering to ground term photoionization rate */
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][0] += ionbal.CompRecoilIonRate[nelem][nelem-ipISO];
	/* add heat to high energy heating term */
	ionbal.PhotoRate_Shell[nelem][nelem-ipISO][0][2] += ionbal.CompRecoilHeatRate[nelem][nelem-ipISO];

	/* option to print ground state photoionization rates */
	if( trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO]) )
	{
		GammaPrt(sp->fb[0].ipIsoLevNIonCon,
			rfield.nflux,
			sp->fb[0].ipOpac,
			ioQQQ,
			sp->fb[0].gamnc,
			sp->fb[0].gamnc*0.05);
	}

	limit = rfield.nflux; 
	/* >>chng 06 aug 17, to numLevels_local instead of _max. */
	for( n=1; n < sp->numLevels_local; n++ )
	{
		/* continuously update rates for n <=3, but only update
		 * rates for higher levels when redoing static opacities */
		if( 0 && sp->st[n].n()>4 && !opac.lgRedoStatic && !sp->lgMustReeval )
			break;

		/** \todo	2	- hydro.lgHInducImp should depend on iso and nelem,
		 * even better - just call one gamnc and within that code
		 * check to see whether induced is important by looking
		 * at occnum near threshold */
		if( hydro.lgHInducImp )
		{
			sp->fb[n].gamnc = 
				GammaBn(
				sp->fb[n].ipIsoLevNIonCon,
				limit,
				sp->fb[n].ipOpac,
				sp->fb[n].xIsoLevNIonRyd,
				&sp->fb[n].RecomInducRate,
				&sp->fb[n].RecomInducCool_Coef,
				&photoHeat)*
				ionbal.lgPhotoIoniz_On;
		}
		else
		{
			sp->fb[n].gamnc = 
				GammaK(sp->fb[n].ipIsoLevNIonCon,
				limit,
				sp->fb[n].ipOpac,1.,
				&photoHeat)*
				ionbal.lgPhotoIoniz_On;

			/* these are zero */
			sp->fb[n].RecomInducRate = 0.;
			sp->fb[n].RecomInducCool_Coef = 0.;
		}
		sp->fb[n].PhotoHeat = photoHeat.HeatNet*ionbal.lgPhotoIoniz_On;

		ASSERT( sp->fb[n].gamnc>= 0. );
		ASSERT( sp->fb[n].PhotoHeat>= 0. );
		/* this loop only has excited states */
	}

	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			if( nelem==ipHYDROGEN )
			{
				fprintf(ioQQQ," buggbugg hphotodebugg%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
					nzone,
				  sp->fb[0].gamnc,
				  sp->fb[1].gamnc,
				  sp->fb[3].gamnc,
				  sp->fb[4].gamnc,
				  sp->fb[5].gamnc,
				  sp->fb[6].gamnc);
			}
		}
	}

	/* >>chng 02 jan 19, kill excited state photoionization with case b no photo */
	/* option for case b conditions, kill all excited state photoionizations */
	if( opac.lgCaseB_no_photo )
	{
		for( n=1; n < sp->numLevels_max; n++ )
		{
			sp->fb[n].gamnc = 0.;
			sp->fb[n].PhotoHeat = 0.;
			sp->fb[n].RecomInducRate = 0.;
			sp->fb[n].RecomInducCool_Coef = 0.;
		}
	}
	{
		/* this block turns off induced recom for some element */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && ipISO==1 && nelem==5)
		{
			/* this debugging block is normally not active */
			for( n=0; n < sp->numLevels_max; n++ )
			{
				sp->fb[n].RecomInducRate = 0.;
			}
		}
	}

	if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
	{
		fprintf( ioQQQ, "     iso_photo, ipISO%2ld nelem%2ld low, hi=",ipISO,nelem);
		fprintf( ioQQQ,PrintEfmt("%9.2e", sp->fb[ipH1s].gamnc));
		ASSERT(nelem>=ipISO);
		fprintf( ioQQQ,PrintEfmt("%9.2e", ionbal.CompRecoilIonRate[nelem][nelem-ipISO]));
		fprintf( ioQQQ, " total=");
		fprintf( ioQQQ,PrintEfmt("%9.2e",sp->fb[ipH1s].gamnc ));
		fprintf( ioQQQ, "\n");
	}
	return;
}
