/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "ionbal.h"
#include "dense.h"
#include "warnings.h"

t_ionbal ionbal;

void t_ionbal::alloc()
{
	ipCompRecoil.reserve(LIMELM);
	RateRecomIso.reserve(LIMELM);
	RateIoniz.reserve(LIMELM);
	CollIonRate_Ground.reserve(LIMELM);
	PhotoRate_Shell.reserve(LIMELM);
	/* create arrays for ions */
	for( long nelem=0; nelem < LIMELM; ++nelem )
	{
		ipCompRecoil.reserve(nelem, nelem+1);
		RateRecomIso.reserve(nelem, NISO);
		RateIoniz.reserve(nelem, nelem+1);
		CollIonRate_Ground.reserve(nelem, nelem+1);
		PhotoRate_Shell.reserve(nelem, nelem+1);

		for( long ion=0; ion < nelem+1; ++ion )
		{
			RateIoniz.reserve(nelem, ion, nelem+2);
			CollIonRate_Ground.reserve(nelem, ion, 2);
			PhotoRate_Shell.reserve(nelem, ion, NSHELLS);
			for( long ns=0; ns < NSHELLS; ++ns )
				PhotoRate_Shell.reserve(nelem, ion, ns, 3);
		}
	}
	/* these will save bound electron recoil information data */
	ipCompRecoil.alloc();
	CompRecoilIonRate.alloc(ipCompRecoil.clone());
	CompRecoilIonRateSave.alloc(ipCompRecoil.clone());
	CompRecoilHeatRate.alloc(ipCompRecoil.clone());
	CompRecoilHeatRateSave.alloc(ipCompRecoil.clone());
	PhotoRate_Shell.alloc();
	CollIonRate_Ground.alloc();
	ExcitationGround.alloc(ipCompRecoil.clone());
	UTA_ionize_rate.alloc(ipCompRecoil.clone());
	UTA_heat_rate.alloc(ipCompRecoil.clone());
	/* space for ionization recombination arrays */
	RateIoniz.alloc();
	RateRecomTot.alloc(ipCompRecoil.clone());
	RateRecomIso.alloc();
	RR_rate_coef_used.alloc(ipCompRecoil.clone());
	RR_Verner_rate_coef.alloc(ipCompRecoil.clone());
	/* rate coefficients [cm3 s-1] for Badnell DR recombination */
	DR_Badnell_rate_coef.alloc(ipCompRecoil.clone());
	RR_Badnell_rate_coef.alloc(ipCompRecoil.clone());
	CX_recomb_rate_used.alloc(ipCompRecoil.clone());
	DR_Badnell_suppress_fact.alloc(ipCompRecoil.clone());

	RateIoniz = 0.;
	RateRecomIso = 0.;

	/* >>chng 03 aug 09, set these to impossible values */
	RateRecomTot = -1.;
	UTA_ionize_rate = -1.;
	UTA_heat_rate = -1.;
	CompRecoilIonRate = -1.;
	CompRecoilIonRateSave = -1.;
	CompRecoilHeatRate = -1.;
	CompRecoilHeatRateSave = -1.;
			
	ipCompRecoil = -100000;
	DR_Badnell_rate_coef = 0.;
	RR_Badnell_rate_coef = 0.;

	RR_rate_coef_used.invalidate();
	RR_Verner_rate_coef.invalidate();
	CX_recomb_rate_used.invalidate();
}

void t_ionbal::zero()
{
	/* now zero out these arrays */
	CompRecoilHeatRate = 0.;
	CompRecoilIonRate = 0.;
	UTA_ionize_rate = 0.;
	UTA_heat_rate = 0.;
	CollIonRate_Ground = 0.;
	ExcitationGround = 0.;
	RateRecomTot = 0.;
	/* must be zero since ion routines use these when
	 * not yet defined */
	PhotoRate_Shell = 0.;

	/* limits for highest and lowest stages of ionization in TrimStage */
	lgTrimhiOn = true;
	trimhi = 1e-6;
	lgTrimloOn = true;
	trimlo = 1e-10;
	lgNewTrim = true;
	
	lgPhotoIoniz_On = true;
	lgCompRecoil = true;

	lgDRsup = true;

	lgNoCota = false;
	for( long nelem = 0; nelem < LIMELM; ++nelem )
	{
		CotaRate[nelem] = 0.;
	}
	ilt = 0;
	iltln = 0;
	ilthn = 0;
	ihthn = 0;
	ifail = 0;
	lgGrainIonRecom = true;

	/* option to print recombination coefficient then exit */
	lgRecom_Badnell_print = false;
	guess_noise = 0.;
}

void t_ionbal::comment(t_warnings& w)
{
	if( ionbal.CompHeating_Max > 0.05 )
	{
		ostringstream chLine;
		chLine <<  "  !Bound Compton heating reached " << fixed << setprecision(2);
		chLine << ionbal.CompHeating_Max*100. << "% of the local heating.";
		w.bangin(chLine.str());
	}
	else if( ionbal.CompHeating_Max > 0.01 )
	{
		ostringstream chLine;
		chLine <<  "   Bound Compton heating reached " << fixed << setprecision(2);
		chLine << ionbal.CompHeating_Max*100. << "% of the local heating.";
		w.notein(chLine.str());
	}

	/* say if 3-body recombination coefficient function outside range of validity
	 * tripped if te/z**2 < 100 or approx 10**13: => effect >50% of radiative
	 * other integers defined in source for da */
	if( ionbal.ifail > 0 && ionbal.ifail <= 10 )
	{
		ostringstream chLine;
		chLine << "   3 body recombination coefficient outside range " << ionbal.ifail;
		w.notein(chLine.str());
	}
	else if( ionbal.ifail > 10 )
	{
		ostringstream chLine;
		chLine << " C-3 body recombination coefficient outside range " << ionbal.ifail;
		w.caunin(chLine.str());
	}
}

double t_ionbal::RateIonizTot( long nelem, long ion ) const
{
	double sum = 0.;
	
	for( long ion_to=ion+1; ion_to<=dense.IonHigh[nelem]; ion_to++ )
		sum += RateIoniz[nelem][ion][ion_to];
	
	return sum;
}
