/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "grainvar.h"

GrainVar gv;

void ShellData::p_clear0()
{
	p.clear();
	y01.clear();
	AvNr.clear();
	Ener.clear();
	y01A.clear();
}

void ShellData::p_clear1()
{
	nelem = LONG_MIN;
	ns = LONG_MIN;
	ionPot = -DBL_MAX;
	ipLo = LONG_MIN;
	nData = 0;
}

void AEInfo::p_clear0()
{
	nData.clear();
	IonThres.clear();
	AvNumber.clear();
	Energy.clear();
}

void AEInfo::p_clear1()
{
	nSubShell = 0;
}

void ChargeBin::p_clear0()
{
	yhat.clear();
	yhat_primary.clear();
	ehat.clear();
	cs_pdt.clear();
	fac1.clear();
	fac2.clear();
}

void ChargeBin::p_clear1()
{
	DustZ = LONG_MIN;
	nfill = 0;
	FracPop = -DBL_MAX;
	tedust = 1.f;
}

void GrainBin::p_clear0()
{
	pool.clear();
	dstab1.clear();
	pure_sc1.clear();
	asym.clear();
	dstab1_x_anu.clear();
	sd.clear();
	y0b06.clear();
	inv_att_len.clear();
}

void GrainBin::p_clear1()
{
	pool.resize(NCHS);
	nDustFunc = DF_STANDARD;
	lgPAHsInIonizedRegion = false;
	avDGRatio = 0.;
	dstfactor = 1.f;
	dstAbund = -FLT_MAX;
	GrnDpth = 1.f;
	cnv_H_pGR = -DBL_MAX;
	cnv_H_pCM3 = -DBL_MAX;
	cnv_CM3_pGR = -DBL_MAX;
	cnv_CM3_pH = -DBL_MAX;
	cnv_GR_pH = -DBL_MAX;
	cnv_GR_pCM3 = -DBL_MAX;
	/* used to check that the energy grid resolution scale factor in
	 * grains opacity files is the same as current cloudy scale */
	RSFCheck = 0.;
	memset( dstems, 0, NDEMS*sizeof(dstems[0]) );
	memset( dstslp, 0, NDEMS*sizeof(dstslp[0]) );
	memset( dstslp2, 0, NDEMS*sizeof(dstslp2[0]) );
	lgTdustConverged = false;
	/* >>chng 00 jun 19, tedust has to be greater than zero
	 * to prevent division by zero in GrainElecEmis and GrainCollHeating, PvH */
	tedust = 1.f;
	TeGrainMax = FLT_MAX;
	avdust = 0.;
	LowestZg = LONG_MIN;
	nfill = 0;
	sd.reserve(15);
	AveDustZ = -DBL_MAX;
	dstpot = -DBL_MAX;
	RateUp = -DBL_MAX;
	RateDn = -DBL_MAX;
	StickElecNeg = -DBL_MAX;
	StickElecPos = -DBL_MAX;
	avdpot = 0.;
	le_thres = FLT_MAX;
	BolFlux = -DBL_MAX;
	GrainCoolTherm = -DBL_MAX;
	GasHeatPhotoEl = -DBL_MAX;
	GrainHeat = DBL_MAX/10.;
	GrainHeatColl = -DBL_MAX;
	GrainGasCool = DBL_MAX/10.;
	ChemEn = -DBL_MAX;
	ChemEnH2 = -DBL_MAX;
	thermionic = -DBL_MAX;
	lgQHeat = false;
	lgUseQHeat = false;
	lgEverQHeat = false;
	lgQHTooWide = false;
	QHeatFailures = 0;
	qnflux = LONG_MAX;
	qnflux2 = LONG_MAX;
	qtmin = -DBL_MAX;
	qtmin_zone1 = -DBL_MAX;
	HeatingRate1 = -DBL_MAX;
	memset( DustEnth, 0, NDEMS*sizeof(DustEnth[0]) );
	memset( EnthSlp, 0, NDEMS*sizeof(EnthSlp[0]) );
	memset( EnthSlp2, 0, NDEMS*sizeof(EnthSlp2[0]) );
	rate_h2_form_grains_HM79 = 0.;
	rate_h2_form_grains_CT02 = 0.;
	rate_h2_form_grains_ELRD = 0.;
	/* >>chng 04 feb 05, zero this rate in case "no molecules" is set, will.in, PvH */
	rate_h2_form_grains_used = 0.;
	DustDftVel = 1.e3f;
	avdft = 0.;
	/* NB - this number should not be larger than NCHU */
	nChrgOrg = gv.nChrgRequested;
	nChrg = nChrgOrg;
	for( int nz=0; nz < NCHS; nz++ )
		ichrg[nz] = nz;
}

void GrainVar::p_clear0()
{
	ReadRecord.clear();
	dstab.clear();
	dstsc.clear();
	for( int nelem=0; nelem < LIMELM; nelem++ )
		AugerData[nelem].clear();
	GrainEmission.clear();
	GraphiteEmission.clear();
	SilicateEmission.clear();
	bin.clear();
}

void GrainVar::p_clear1()
{
	lgAnyDustVary = false;
	TotalEden = 0.;
	dHeatdT = 0.;
	lgQHeatAll = false;
	/* lgGrainElectrons - should grain electron source/sink be included in overall electron sum?
	 * default is true, set false with no grain electrons command */
	lgGrainElectrons = true;
	lgQHeatOn = true;
	lgDHetOn = true;
	lgDColOn = true;
	GrainMetal = 1.;
	dstAbundThresholdNear = 1.e-6f;
	dstAbundThresholdFar = 1.e-3f;
	lgWD01 = false;
	nChrgRequested = NCHRG_DEFAULT;
	/* by default grains always reevaluated - command grains reevaluate off sets to false */
	lgReevaluate = true;
	/* flag saying neg grain drift vel found */
	lgNegGrnDrg = false;

	/* counts how many times GrainDrive has been called */
	nCalledGrainDrive = 0;

	/* this is sest true with "set PAH Bakes" command - must also turn off
	 * grain heating with "grains no heat" to only get their results */
	lgBakesPAH_heat = false;

	/* this is option to turn off all grain physics while leaving
	 * the opacity in, set false with no grain physics command */
	lgGrainPhysicsOn = true;

	/* scale factor set with SET GRAINS HEAT command to rescale grain photoelectric
	 * heating as per Allers et al. 2005 */
	GrainHeatScaleFactor = 1.f;

	/* the following entries define the physical behavior of each type of grains
	 * (entropy function, expression for Zmin and ionization potential, etc) */
	which_enth[MAT_CAR] = ENTH_CAR;
	which_zmin[MAT_CAR] = ZMIN_CAR;
	which_pot[MAT_CAR] = POT_CAR;
	which_ial[MAT_CAR] = IAL_CAR;
	which_pe[MAT_CAR] = PE_CAR;
	which_strg[MAT_CAR] = STRG_CAR;
	which_H2distr[MAT_CAR] = H2_CAR;

	which_enth[MAT_SIL] = ENTH_SIL;
	which_zmin[MAT_SIL] = ZMIN_SIL;
	which_pot[MAT_SIL] = POT_SIL;
	which_ial[MAT_SIL] = IAL_SIL;
	which_pe[MAT_SIL] = PE_SIL;
	which_strg[MAT_SIL] = STRG_SIL;
	which_H2distr[MAT_SIL] = H2_SIL;

	which_enth[MAT_PAH] = ENTH_PAH;
	which_zmin[MAT_PAH] = ZMIN_CAR;
	which_pot[MAT_PAH] = POT_CAR;
	which_ial[MAT_PAH] = IAL_CAR;
	which_pe[MAT_PAH] = PE_CAR;
	which_strg[MAT_PAH] = STRG_CAR;
	which_H2distr[MAT_PAH] = H2_CAR;

	which_enth[MAT_CAR2] = ENTH_CAR2;
	which_zmin[MAT_CAR2] = ZMIN_CAR;
	which_pot[MAT_CAR2] = POT_CAR;
	which_ial[MAT_CAR2] = IAL_CAR;
	which_pe[MAT_CAR2] = PE_CAR;
	which_strg[MAT_CAR2] = STRG_CAR;
	which_H2distr[MAT_CAR2] = H2_CAR;

	which_enth[MAT_SIL2] = ENTH_SIL2;
	which_zmin[MAT_SIL2] = ZMIN_SIL;
	which_pot[MAT_SIL2] = POT_SIL;
	which_ial[MAT_SIL2] = IAL_SIL;
	which_pe[MAT_SIL2] = PE_SIL;
	which_strg[MAT_SIL2] = STRG_SIL;
	which_H2distr[MAT_SIL2] = H2_SIL;

	which_enth[MAT_PAH2] = ENTH_PAH2;
	which_zmin[MAT_PAH2] = ZMIN_CAR;
	which_pot[MAT_PAH2] = POT_CAR;
	which_ial[MAT_PAH2] = IAL_CAR;
	which_pe[MAT_PAH2] = PE_CAR;
	which_strg[MAT_PAH2] = STRG_CAR;
	which_H2distr[MAT_PAH2] = H2_CAR;

	which_enth[MAT_SIC] = ENTH_SIC;
	which_zmin[MAT_SIC] = ZMIN_CAR;
	which_pot[MAT_SIC] = POT_CAR;
	which_ial[MAT_SIC] = IAL_CAR;
	which_pe[MAT_SIC] = PE_CAR;
	which_strg[MAT_SIC] = STRG_CAR;
	which_H2distr[MAT_SIC] = H2_CAR;

	for( int nelem=0; nelem < LIMELM; nelem++ )
	{
		for( int ion=0; ion <= nelem+1; ion++ )
		{
			for( int ion_to=0; ion_to <= nelem+1; ion_to++ )
			{
				GrainChTrRate[nelem][ion][ion_to] = 0.f;
			}
		}
	}

	/* this sets the default abundance dependence for PAHs,
	 * proportional to n(H0) / n(Htot) 
	 * changed with SET PAH command */
	chPAH_abundance = "H";
}
