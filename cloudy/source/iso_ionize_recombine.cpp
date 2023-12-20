/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_ionize_recombine find state specific creation and destruction rates for iso sequences */
/*ChargeTransferUpdate update rate of ct ionization and recombination for H atoms */
#include "cddefines.h"
#include "ionbal.h"
#include "conv.h"
#include "atmdat.h"
#include "secondaries.h"
#include "iso.h"
#include "phycon.h"
#include "rfield.h"
#include "trace.h"
#include "freebound.h"
#include "dense.h"
/*lint -e778 constant expression evaluates to zero - in HMRATE macro */
/*iso_charge_transfer_update update rate of ct ionization and recombination for H atoms */
void iso_charge_transfer_update(long nelem1)
{
	DEBUG_ENTRY( "iso_charge_transfer_update()" );

	if( nelem1 >= t_atmdat::NCX )
		return;

	atmdat.CharExcIonTotal[nelem1] = 0.;
	atmdat.CharExcRecTotal[nelem1] = 0.;

	// find total charge transfer rate coefficients
	// charge transfer reactions of lighter species
	for ( long nelem=ipHYDROGEN; nelem < nelem1; ++nelem)
	{
		atmdat.CharExcIonTotal[nelem1] += atmdat.CharExcIonOf[nelem][nelem1][0]*dense.xIonDense[nelem][1];
		long ipISO=nelem;
		atmdat.CharExcRecTotal[nelem1] += atmdat.CharExcRecTo[nelem][nelem1][0]*iso_sp[ipISO][nelem].st[0].Pop();
	}

	// and of heavier species
	for( long nelem=nelem1+1; nelem<LIMELM; ++nelem)
	{
		for( long ion=0; ion<=nelem; ++ion )
		{
			// charge transfer ionization of nelem1, recombination for other species, cm-3 s-1
			atmdat.CharExcIonTotal[nelem1] += atmdat.CharExcRecTo[nelem1][nelem][ion]*dense.xIonDense[nelem][ion+1]; 
			// charge transfer recombination of nelem1, ionization for other species, cm-3 s-1
			atmdat.CharExcRecTotal[nelem1] += atmdat.CharExcIonOf[nelem1][nelem][ion]*dense.xIonDense[nelem][ion];
		}
	}

	return;
}

/*iso_ionize_recombine find state specific creation and destruction rates for iso sequences */
void iso_ionize_recombine(
	/* iso sequence, hydrogen or helium for now */
	long ipISO , 
	/* the chemical element, 0 for hydrogen */
	long int nelem )
{
	long int level;
	double Recom3Body;

	DEBUG_ENTRY( "iso_ionize_recombine()" );

	ASSERT( ipISO >= 0 && ipISO < NISO );
	ASSERT( nelem >= 0 && nelem < LIMELM );

	/* find recombination and ionization elements, will use to get simple estimate
	 * of ionization ratio below */
	/* >>chng 06 jul 20, level should go to numLevels_local instead of numLevels_max */

	fixit("must apply error to iso.gamnc");

	for( level=ipH1s; level< iso_sp[ipISO][nelem].numLevels_local; ++level)
	{
		// This term is intended to capture LTE in DR-dominated plasmas.
		// It is scaled by Occ num / Occ num at STE to ensure correct limits w.r.t. ambient radiation. 
		double DR_reverse = 0.;
		if( ipISO > ipH_LIKE && iso_sp[ipISO][nelem].fb[level].PopLTE > SMALLFLOAT )
		{ 
			long indexIP = iso_sp[ipISO][nelem].fb[level].ipIsoLevNIonCon-1;
			double OccNum =	rfield.OccNumbIncidCont[indexIP] + rfield.OccNumbContEmitOut[indexIP];
			double OccNumBB = 1./(exp( iso_sp[ipISO][nelem].fb[level].xIsoLevNIonRyd / phycon.te_ryd ) - 1. );
			double RelOccNum = MIN2( 1., OccNum / OccNumBB );
			//fprintf( ioQQQ, "DEBUGGG DR_reverse %li\t%li\t%e\t%e\t%e\t%e\t%e\t%e\n", level, indexIP, RelOccNum, OccNum, OccNumBB, 
			//	rfield.anu(indexIP), iso_sp[ipISO][nelem].fb[level].xIsoLevNIonRyd, rfield.anu(indexIP+1) );
			DR_reverse = iso_sp[ipISO][nelem].fb[level].DielecRecomb * RelOccNum /  
				iso_sp[ipISO][nelem].fb[level].PopLTE; 
		}
	
		/* all process moving level to continuum, units s-1 */
		iso_sp[ipISO][nelem].fb[level].RateLevel2Cont = iso_sp[ipISO][nelem].fb[level].gamnc + 
		  	iso_sp[ipISO][nelem].fb[level].ColIoniz* dense.EdenHCorr + 
			DR_reverse +
			secondaries.csupra[nelem][nelem-ipISO]*iso_ctrl.lgColl_ionize[ipISO];

		/* all processes from continuum to level n, units s-1 */
		iso_sp[ipISO][nelem].fb[level].RateCont2Level = (
			/* radiative recombination */
			iso_sp[ipISO][nelem].fb[level].RadRecomb[ipRecRad]*
			iso_sp[ipISO][nelem].fb[level].RadRecomb[ipRecNetEsc] + 

			/* dielectronic recombination */
			iso_sp[ipISO][nelem].fb[level].DielecRecomb +

			/* induced recombination */
			iso_sp[ipISO][nelem].fb[level].RecomInducRate*iso_sp[ipISO][nelem].fb[level].PopLTE + 

			/* collisional or three body recombination */
			/* PopLTE(level,nelem) is only LTE pop when mult by n_e n_H */
			iso_sp[ipISO][nelem].fb[level].ColIoniz*dense.EdenHCorr*iso_sp[ipISO][nelem].fb[level].PopLTE
			) * dense.eden;
	
		ASSERT( iso_sp[ipISO][nelem].fb[level].RadRecomb[ipRecRad] >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[level].RadRecomb[ipRecNetEsc] >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[level].DielecRecomb >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[level].RecomInducRate >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[level].PopLTE >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[level].ColIoniz >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[level].RateCont2Level >= 0. );

		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			iso_sp[ipISO][nelem].fb[level].RateCont2Level *=
				iso_sp[ipISO][nelem].ex[ iso_sp[ipISO][nelem].numLevels_max ][level].ErrorFactor[IPRAD];
		}
	}

	/* now create sums of recombination and ionization rates for ISO species */
	ionbal.RateRecomIso[nelem][ipISO] = 0.;
	ionbal.RR_rate_coef_used[nelem][nelem-ipISO] = 0.;
	Recom3Body = 0.;
	/* >>chng 06 jul 20, level should go to numLevels_local instead of numLevels_max */
	for( level=0; level< iso_sp[ipISO][nelem].numLevels_local; ++level)
	{

		/* units of ionbal.RateRecomTot are s-1,
		 * equivalent ionization term is done after level populations are known */
		ionbal.RateRecomIso[nelem][ipISO] += iso_sp[ipISO][nelem].fb[level].RateCont2Level;

		/* just the radiative recombination rate coef, cm3 s-1 */
		ionbal.RR_rate_coef_used[nelem][nelem-ipISO] += iso_sp[ipISO][nelem].fb[level].RadRecomb[ipRecRad]*
			iso_sp[ipISO][nelem].fb[level].RadRecomb[ipRecNetEsc];

		/* >>chng 05 jul 11, from > to >=, some very opt thick sims did block escape to zero */
		ASSERT( ionbal.RR_rate_coef_used[nelem][nelem-ipISO]>= 0. );

		/* this is three-body recombination rate coef by itself - 
		 * need factor of eden to become rate */
		Recom3Body += iso_sp[ipISO][nelem].fb[level].ColIoniz*dense.EdenHCorr*iso_sp[ipISO][nelem].fb[level].PopLTE;
	}

	/* fraction of total recombinations due to three body - when most are due to this
	 * small changes in temperature can produce large changes in recombination coefficient,
	 * and so in ionization */
	iso_sp[ipISO][nelem].RecomCollisFrac = Recom3Body* dense.eden / ionbal.RateRecomIso[nelem][ipISO];

	/* very first pass through here rate RateIoniz not yet evaluated */
	if( conv.nTotalIoniz==0 )
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = iso_sp[ipISO][nelem].fb[0].RateLevel2Cont;

	/* get simple estimate of atom to ion ratio */
	if( ionbal.RateRecomIso[nelem][ipISO] > 0. )
	{
		iso_sp[ipISO][nelem].xIonSimple = ionbal.RateIonizTot(nelem,nelem-ipISO)/ionbal.RateRecomIso[nelem][ipISO];
	}
	else
	{
		iso_sp[ipISO][nelem].xIonSimple = 0.;
	}

	if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
	{
		fprintf( ioQQQ, "     iso_ionize_recombine iso=%2ld Z=%2ld Level2Cont[0] %10.2e RateRecomTot %10.2e xIonSimple %10.2e\n", 
			ipISO, nelem, iso_sp[ipISO][nelem].fb[0].RateLevel2Cont, ionbal.RateRecomIso[nelem][ipISO], iso_sp[ipISO][nelem].xIonSimple );
	}

	return;
}
/*lint +e778 constant expression evaluates to zero - in HMRATE macro */
