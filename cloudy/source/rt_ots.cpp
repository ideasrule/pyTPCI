/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_OTS compute diffuse fields due to H, He atoms, ion, triplets, metal recombination,
 * called by ConvBase  */
/*RT_OTS_AddLine add local destruction of lines to ots field */
/*RT_OTS_AddCont add local destruction of continuum to ots field */
/*RT_OTS_Update sum flux, otscon, otslin, ConInterOut, outlin, to form SummeDif, SummedCon SummedOcc */
/*RT_OTS_Zero - zero out some vectors - 
 * this is only called when code initialized by ContSetIntensity */
/*RT_OTS_ChkSum sanity check confirms summed continua reflect contents of individuals */
#include "cddefines.h"
#include "taulines.h"
#include "opacity.h"
#include "dense.h"
#include "iso.h"
#include "hmi.h"
#include "h2.h"
#include "rfield.h"
#include "conv.h"
#include "rt.h"
#include "heavy.h"
#include "he.h"
#include "trace.h"
#include "mole.h"
#include "freebound.h"
#include "two_photon.h"
#include "cosmology.h"

/* this flag may be used for debugging ots rate changes */
static int nOTS_Line_type = -1;
static int nOTS1=-1 , nOTS2=-1;
/*add local destruction of continuum to ots field */
STATIC void RT_OTS_AddCont(
	/* the ots rate itself */
	realnum ots, 
	/* pointer to continuum cell for ots, on f scale */
	long int ip);

/* =================================================================== */

void RT_OTS(void)
{
	long int
		ipla,
		ipISO ,
		nelem,
		n;
	realnum 
		difflya,
		esc,
		ots;

	/* the Bowen HeII yield 
	 * >>chng 06 aug 08, from 0.3 to 0.4, update with netzer */
	const realnum BOWEN = 0.5f;
	long int ipHi, 
	  ipLo;

	double bwnfac;
	double ots660;
	realnum cont_phot_destroyed;
	double save_lya_dest,
	  save_he2lya_dest = 0.;

	double save_he2rec_dest = 0.;

	/*static long int nCall=0;
	++nCall;
	fprintf(ioQQQ,"debugggtos enter %li\n", nCall );*/

	DEBUG_ENTRY( "RT_OTS()" );

	for( long i=0; i < rfield.nflux; i++ )
	{
		rfield.otslin[i] = 0.;
		rfield.otscon[i] = 0.;
	}

	/**************************************************************************
	 *
	 * the bowen HeII - OIII fluorescence problem
	 *
	 **************************************************************************/
	nOTS_Line_type = 0;
	nelem = ipHELIUM;
	if( dense.lgElmtOn[nelem] )
	{
		/* conversion per unit atom to OIII, at end of sub we divide by it,
		 * to fix lines back to proper esc/dest probs */
		bwnfac = BOWEN * MAX2(0.f,1.f- iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Pesc_total());

		/* the last factor accounts for fact that two photons are produced,
		 * and the branching ratio */
		ots660 = iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Aul()*
		  iso_sp[ipH_LIKE][nelem].st[ipH2p].Pop()*
		  /*>>chng 06 aug 08, rm 0.8 factor from below, renorm aft discuss with Netzer */
		  iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Pdest() *BOWEN*2.0;

		/* now add this to the ots field */
		if( ots660 > SMALLFLOAT ) 
			RT_OTS_AddLine(ots660 , he.ip660 );

		/* decrease the destruction prob by the amount we will add elsewhere,
		 * ok since dest probs updated on every iteration*/
		iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Pdest() *= (realnum)bwnfac;
		ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Pdest() >= 0. );
		{
			/* debugging code for line oscillation problems 
			 * most often Lya OTS oscillations*/
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ,"DEBUG HeII Bowen nzone %li bwnfac:%.2e bwnfac esc:%.2e ots660 %.2e\n",
					nzone,
					bwnfac , 
					bwnfac/BOWEN , 
					ots660 );
			}
		}
	}

	else
	{
		bwnfac = 1.;
	}

	/* save Lya loss rates so we can reset at end */
	save_lya_dest = iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Pdest();

	/* this is option to kill Lya ots rates, 
	 * rfield.lgLyaOTS is usually true (1), and set false (0) with
	 * no lya ots command */
	iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Pdest() *= rfield.lgLyaOTS;

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* option to kill heii lya and rec continua ots */
		save_he2lya_dest = iso_sp[ipH_LIKE][ipHELIUM].trans(ipH2p,ipH1s).Emis().Pdest();
		iso_sp[ipH_LIKE][ipHELIUM].trans(ipH2p,ipH1s).Emis().Pdest() *= rfield.lgHeIIOTS;

		/* option to kill heii lya and rec continua ots */
		save_he2rec_dest = iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].RadRecomb[ipRecRad];
		iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].RadRecomb[ipRecRad] *= rfield.lgHeIIOTS;
	}

	nOTS_Line_type = 1;

	/* make ots fields due to lines and continua of species treated with unified 
	 * isoelectronic sequence */
	/* loop over all elements */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			nOTS1 = ipISO;
			nOTS2 = nelem;
			/* do if this stage exists */
			/** \todo	2	should also check whether IonLo is in bounds - in func_set_ion test
			 * He0 is set to zero, so this does not do anything.  as NISO grows this
			 * will become larger waste of time */
			if( (dense.IonHigh[nelem] >= nelem+1-ipISO )  )
			{
				/* generate line ots rates */
				/* now loop over all possible levels, but skip non-radiative
				* since there is no pointer to this continuum */
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
				{
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						/* this signifies a fake line */
						/* >>chng 03 may 19, DEST0 is the smallest possible
						 * dest prob - not a real number, don't add to ots field */
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() < 1 ||
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pdest()<= DEST0 )
							continue;

						/* ots rates, the destp prob was set in hydropesc */
						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ots() = 
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul()*
							iso_sp[ipISO][nelem].st[ipHi].Pop()*
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pdest();

						ASSERT( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ots() >= 0. );
						/* way to kill lyalpha ots rates
						if( nelem==ipHYDROGEN && ipHi==ipH2p && ipLo==ipH1s )
						{
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ots() = 0.;
						} */

						/* finally dump the ots rate into the stack */
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ots() > SMALLFLOAT ) 
							RT_OTS_AddLine(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ots(),
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() );
					}
				}
				{
					/* debugging code for line oscillation problems 
					 * most often Lya OTS oscillations*/
					/*@-redef@*/
					enum {DEBUG_LOC=false};
					/*@+redef@*/
					if( DEBUG_LOC )
					{
						long ip;
						if( ipISO==0 && nelem==0  && nzone>500  ) 
						{
							ipHi = 2;
							ipLo = 0;
							ip = iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont();
							fprintf(ioQQQ,"DEBUG hlyaots\t%.2f\tEdenTrue\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
								fnzone, 
								dense.EdenTrue ,
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ots(),
								opac.opacity_abs[ip-1],
								iso_sp[ipISO][nelem].st[ipHi].Pop(),
								dense.xIonDense[nelem][nelem+1-ipISO],
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pdest(),
								rfield.otslin[ip-1]);
						}
					}
				}

				/**************************************************************************
				*
				* ots recombination bound-free b-f continua continuum
				*
				**************************************************************************/

				/* put in OTS continuum */
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( n=0; n < iso_sp[ipISO][nelem].numLevels_local; n++ )
				{
					cont_phot_destroyed = (realnum)(iso_sp[ipISO][nelem].fb[n].RadRecomb[ipRecRad]*
						(1. - iso_sp[ipISO][nelem].fb[n].RadRecomb[ipRecEsc])*dense.eden*
						dense.xIonDense[nelem][nelem+1-ipISO]);
					ASSERT( cont_phot_destroyed >= 0. );

					/* continuum energy index used in this routine is decremented by one there */
					RT_OTS_AddCont(cont_phot_destroyed,iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon);
					/* debugging code for rec continua */
					{
						/*@-redef@*/
						enum {DEBUG_LOC=false};
						/*@+redef@*/
						if( DEBUG_LOC && nzone > 400 && nelem==0 && n==2 )
						{
							long ip = iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon-1;
							fprintf(ioQQQ,"hotsdebugg\t%.3f\t%li\th con ots\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
								fnzone, 
								n ,
								iso_sp[ipISO][nelem].st[n].Pop(),
								cont_phot_destroyed,
								cont_phot_destroyed/opac.opacity_abs[ip],
								rfield.otscon[ip] ,
								opac.opacity_abs[ip] ,
								findspecieslocal("H-")->den ,
								hmi.HMinus_photo_rate);
						}
					}
				}
			}
		}
	}
	/* more debugging code for rec continua */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			nelem = 0;
			fprintf(ioQQQ,"hotsdebugg %li \t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				nzone, 
				rfield.otscon[iso_sp[ipH_LIKE][nelem].fb[0].ipIsoLevNIonCon-1],
				rfield.otscon[iso_sp[ipH_LIKE][nelem].fb[1].ipIsoLevNIonCon-1],
				rfield.otscon[iso_sp[ipH_LIKE][nelem].fb[3].ipIsoLevNIonCon-1],
				rfield.otscon[iso_sp[ipH_LIKE][nelem].fb[4].ipIsoLevNIonCon-1],
				rfield.otscon[iso_sp[ipH_LIKE][nelem].fb[5].ipIsoLevNIonCon-1],
				rfield.otscon[iso_sp[ipH_LIKE][nelem].fb[6].ipIsoLevNIonCon-1],
				opac.opacity_abs[iso_sp[ipH_LIKE][nelem].fb[6].ipIsoLevNIonCon-1]);
		}
	}

	/* now reset Lya dest prob in case is was clobbered by rfield.lgHeIIOTS */
	iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Pdest() = (realnum)save_lya_dest;
	if( dense.lgElmtOn[ipHELIUM] )
	{
		iso_sp[ipH_LIKE][ipHELIUM].trans(ipH2p,ipH1s).Emis().Pdest() = (realnum)save_he2lya_dest;
		iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].RadRecomb[ipRecRad] = save_he2rec_dest;
		if( bwnfac > 0. )
		{
			/* increase the destruction prob by the amount we decreased it above */
			iso_sp[ipH_LIKE][ipHELIUM].trans(ipH2p,ipH1s).Emis().Pdest() /= (realnum)bwnfac;
		}
	}

	if( trace.lgTrace )
	{
		fprintf(ioQQQ,"     RT_OTS Pdest %.2e ots rate %.2e in otslin %.2e con opac %.2e\n",
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Pdest(), 
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().ots() , 
			rfield.otslin[iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont()-1]  ,
			opac.opacity_abs[iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont()-1]
			);
	}

	nOTS_Line_type = 2;
	/* recombination Lya for all elements not yet converted into std isoelectronc form */
	for( nelem=NISO; nelem < LIMELM; nelem++ )
	{
		long int ion;
		/* do not include species treated in iso-electronic fashion in the following,
		 * these were treated above */
		for( ion=0; ion < nelem+1-NISO; ion++ )
		{
			if( dense.xIonDense[nelem][ion+1] > 0. )
			{
				nOTS1 = nelem;
				nOTS2 = ion;
				/* now do the recombination Lya */
				ipla = Heavy.ipLyHeavy[nelem][ion];
				ASSERT( ipla>0 );
				esc = opac.ExpmTau[ipla-1];
				/* xLyaHeavy is set to a fraction of total rad rec in ion_recomb, includes eden */
				difflya = Heavy.xLyaHeavy[nelem][ion]*dense.xIonDense[nelem][ion+1];
				/* >>chng 00 dec 22, from MIN2 to MAX2, MIN2 had effect of always
				 * setting the ots rates to zero */
				ots = difflya*MAX2(0.f,1.f-esc);
				/*if( nelem==6 && ion==2 )
					fprintf(ioQQQ," debugggnly\t %.2f\t%.2e\n",fnzone, ots );*/
				ASSERT( ots >= 0.);
				/*if( iteration == 2 && nzone>290 && ipla== 2339 )
					fprintf(ioQQQ,"recdebugg1 %.2e %li %li %.2e %.2e \n", 
					ots, nelem, ion,
					esc , dense.xIonDense[nelem][ion+1]);*/
				if( ots > SMALLFLOAT ) 
					RT_OTS_AddLine(ots,ipla);

				/* now do the recombination balmer lines */
				ipla = Heavy.ipBalHeavy[nelem][ion];
				esc = opac.ExpmTau[ipla-1];
				/* xLyaHeavy is set to a fraction of total rad rec in ion_recomb, includes eden */
				difflya = Heavy.xLyaHeavy[nelem][ion]*dense.xIonDense[nelem][ion+1];
				/* >>chng 00 dec 22, from MIN2 to MAX2, MIN2 had effect of always
				 * setting the ots rates to zero */
				ots = difflya*MAX2(0.f,1.f-esc);
				ASSERT( ots >= 0.);
				/*if( iteration == 2 &&nzone==294 && ipla== 2339 )
					fprintf(ioQQQ,"recdebugg2 %.2e %li %li\n", ots, nelem, ion );*/
				if( ots > SMALLFLOAT ) 
					RT_OTS_AddLine(ots,ipla);
			}
		}
	}

	nOTS_Line_type = 5;
	/* do the level2 level 2 lines */
	for( nOTS1=0; nOTS1 < nWindLine; nOTS1++ )
	{
		if( (*TauLine2[nOTS1].Hi()).IonStg() < (*TauLine2[nOTS1].Hi()).nelem()+1-NISO )
		{
			TauLine2[nOTS1].Emis().ots() = (*TauLine2[nOTS1].Hi()).Pop() * TauLine2[nOTS1].Emis().Aul() * TauLine2[nOTS1].Emis().Pdest();
			if( TauLine2[nOTS1].Emis().ots() > SMALLFLOAT ) 
				RT_OTS_AddLine( TauLine2[nOTS1].Emis().ots() , TauLine2[nOTS1].ipCont());
		}
	}

	nOTS_Line_type = 6;
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( dBaseSpecies[ipSpecies].lgActive )
		{
			for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
				 tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local || (*tr).ipCont() <= 0)
					continue;
				(*tr).Emis().ots() = (*(*tr).Hi()).Pop() * (*tr).Emis().Aul() * (*tr).Emis().Pdest();
				RT_OTS_AddLine( (*tr).Emis().ots() , (*tr).ipCont());
			}
		}
	}

	nOTS_Line_type = 7;
	/* the large H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_RT_OTS();

	return;
}

/* =================================================================== */

void RT_OTS_AddLine(double ots, 
	/* pointer on the f scale */
  long int ip )
{

	DEBUG_ENTRY( "RT_OTS_AddLine()" );

	/* add ots due to line destruction to radiation field */

	/* return if outside bounds of this continuum source, ip > rfield.nflux
	 * first case ip==0 happens when called with dummy line */
	if( ip==0 || ip > rfield.nflux )
	{ 
		return;
	}

	/*the local ots rate must be non-negative */
	ASSERT( ots >= 0. );
	/* continuum pointer must be positive */
	ASSERT( ip > 0 );

	/* add locally destroyed flux of photons to line OTS array */
	/* check whether local gas opacity  (units cm-1) is positive, if so
	 * convert line destruction rate into ots rate by dividing by it */
	if( opac.opacity_abs[ip-1] > 0. )
	{
		rfield.otslin[ip-1] += (realnum)(ots/opac.opacity_abs[ip-1]);
	}
	/* first iteration is 1, second is two */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && /*iteration == 2 && nzone>294 &&*/ ip== 2363 /*&& ots > 1e16*/)
		{
			fprintf(ioQQQ,"DEBUG ots, opc, otsr %.3e\t%.3e\t%.3e\t",
				ots ,
				opac.opacity_abs[ip-1],
				ots/opac.opacity_abs[ip-1] );
			fprintf(ioQQQ,"iteration %li type %i %i %i \n", 
				iteration, 
				nOTS_Line_type,
				nOTS1,nOTS2 );
		}
	}
	return;
}

/* =================================================================== */

/*add local destruction of continuum to ots field */
STATIC void RT_OTS_AddCont(
					/* the ots rate itself */
					realnum ots, 
					/* pointer to continuum cell for ots, on f scale */
					long int ip)
{

	DEBUG_ENTRY( "RT_OTS_AddCont()" );

	/* 
	 * routine called to add ots due to continuum destruction to
	 * radiation field
	 */

	/* check if outside bounds of this continuum source */
	if( ip > rfield.nflux )
	{ 
		return;
	}

	ASSERT( ip > 0 );
	ASSERT( ots >= 0. );
	ASSERT( ip <= rfield.nflux_with_check );

	/* add locally destroyed flux of photons to continuum OTS array */
	/* check whether local gas opacity  (units cm-1) is positive, if so
	 * convert continuum destruction rate into ots rate by dividing by it */
	if( opac.opacity_abs[ip-1] > 0. )
	{
		rfield.otscon[ip-1] += (realnum)(ots/opac.opacity_abs[ip-1]);
	}
	return;
}

/* =================================================================== */

/*RT_OTS_Update update ots rates, called in ConvBase */
void RT_OTS_Update(double *SumOTS) /* summed ots rates */
{
	long int i;

	DEBUG_ENTRY( "RT_OTS_Update()" );

	/* option to kill ots rates with no ots lines command */
	if( rfield.lgKillOTSLine )
	{
		for( i=0; i < rfield.nflux; i++ )
		{
			rfield.otslin[i] = 0.;
		}
	}
	
	vzero(rfield.ConOTS_local_photons);
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* >>chng 01 sep 23, rewrote for iso sequences */
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{	
			if( dense.IonHigh[nelem] >= nelem+1-ipISO  )
			{
				t_iso_sp* sp = &iso_sp[ipISO][nelem];
				
				for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
				{
					CalcTwoPhotonEmission( *tnu, rfield.lgInducProcess && iso_ctrl.lgInd2nu_On );
					for( long nu=0; nu < tnu->ipTwoPhoE; nu++ )
					{
						rfield.ConOTS_local_photons[nu] += tnu->local_emis[nu] * (1.f - opac.ExpmTau[nu]);
					}
				}
			}
		}
	}

	/* remember largest change in ots rates */
	*SumOTS = 0.;
	/* now update new ots rates */
	for( i=0; i < rfield.nflux; ++i )
	{
		double CurrentInverseOpacity = 1./MAX2( SMALLDOUBLE , opac.opacity_abs[i] );

		/* this is local ots continuum created by destroyed diffuse continua, 
		 * currently only two-photon */
		rfield.ConOTS_local_OTS_rate[i] = (realnum)((double)rfield.ConOTS_local_photons[i]*CurrentInverseOpacity);		

		/* remember sum of ots rates for convergence criteria */
		*SumOTS += (rfield.otscon[i] + rfield.otslin[i])*opac.opacity_abs[i];

		if( !cosmology.lgDo )
		{
			rfield.SummedDif[i] = rfield.otscon[i] + rfield.otslin[i] + rfield.outlin_noplot[i]+ 
				rfield.ConInterOut[i]*rfield.lgOutOnly + rfield.outlin[0][i] +
				rfield.ConOTS_local_OTS_rate[i];
		}

		rfield.SummedCon[i] = rfield.flux[0][i] + rfield.SummedDif[i];
		rfield.SummedOcc[i] = rfield.SummedCon[i]*rfield.convoc[i];
	}

	rfield.setTrimming();

	/* sum of accumulated flux from particular frequency to infinity */
	rfield.flux_accum[rfield.nflux-1] = 0.;
	for( i=1; i < rfield.nflux; i++ )
	{
		rfield.flux_accum[rfield.nflux-i-1] = rfield.flux_accum[rfield.nflux-i] +
			rfield.SummedCon[rfield.nflux-i-1];
	}

	/* >>chng 02 jul 23, set to black body at local temp if in optically thick continuum,
	 * between plasma frequency and energy where brems is thin */
	ASSERT( rfield.ipPlasma > 0 );

	/* all radiation fields are zero below plasma frequency */
	for( i=0; i < rfield.ipPlasma-1; i++ )
	{
		rfield.otscon[i] = 0.;
		rfield.ConOTS_local_OTS_rate[i] = 0.;
		rfield.ConOTS_local_photons[i] = 0.;
		rfield.otslin[i] = 0.;
		rfield.SummedDif[i] = 0.;
		rfield.SummedCon[i] = 0.;
		rfield.SummedOcc[i] = 0.;
		rfield.ConInterOut[i] = 0.;
	}
	return;
}

/* =================================================================== */

/*RT_OTS_Zero zero out some vectors - this is only called when code 
 * initialized by ContSetIntensity */
void RT_OTS_Zero( void )
{
	long int i;

	DEBUG_ENTRY( "RT_OTS_Zero()" );

	/* this loop goes up to nflux itself (<=) since the highest cell
	 * will be used to pass unity through the code to verify integrations */
	for( i=0; i <= rfield.nflux; i++ )
	{
		rfield.SummedDif[i] = 0.;
		/* the main ots vectors */
		rfield.otscon[i] = 0.;
		rfield.otslin[i] = 0.;

		rfield.ConInterOut[i] = 0.;
		rfield.outlin[0][i] = 0.;
		rfield.outlin_noplot[i] = 0.;
		rfield.SummedDif[i] = 0.;
		/* "zero" for the summed con will be just the incident radiation field */
		rfield.SummedCon[i] = rfield.flux[0][i];
		rfield.SummedOcc[i] = rfield.SummedCon[i]*rfield.convoc[i];
		rfield.ConOTS_local_photons[i] = 0.;
		rfield.ConOTS_local_OTS_rate[i] = 0.;
	}

	rfield.setTrimming();

	rfield.resetCoarseTransCoef();
	return;
}

/* =================================================================== */

/*RT_OTS_ChkSum sanity check confirms summed continua reflect contents of individuals */
void RT_OTS_ChkSum(long int ipPnt)
{
	static long int nInsane=0;
	long int i;
	double phisig;
	const int LIM_INSAME_PRT = 30;

	DEBUG_ENTRY( "RT_OTS_ChkSum()" );

	/* this entire sub is a sanity check */
	/* >>chng 02 jul 23, lower bound from 0 to rfield.ipEnergyBremsThin - since now
	 * set radiation field to black body below this energy */
	for( i=rfield.ipEnergyBremsThin; i < rfield.nflux; i++ )
	{
		phisig = rfield.otscon[i] + rfield.otslin[i] + rfield.ConInterOut[i]*rfield.lgOutOnly + 
		  rfield.outlin[0][i]+
		  rfield.outlin_noplot[i]+
		  rfield.ConOTS_local_OTS_rate[i];
		/* >>chng 02 sep 19, add sec test on SummedDif since it can be zero whild 
		 * phisig is just above small float */
		if( phisig > 0. && rfield.SummedDif[i]> 0.)
		{
			if( fabs(rfield.SummedDif[i]/phisig-1.) > 1e-3 )
			{
				++nInsane;
				/* limit amount of printout - in many failures would print entire
				 * continuum array */
				if( nInsane < LIM_INSAME_PRT )
				{
					fprintf( ioQQQ, " PROBLEM RT_OTS_ChkSum insane SummedDif at energy %.5e error= %.2e i=%4ld\n", 
					rfield.anu(i), rfield.SummedDif[i]/phisig - 1., i );
					fprintf( ioQQQ, " SummedDif, sum are%11.4e%11.4e\n", 
					rfield.SummedDif[i], phisig );
					fprintf( ioQQQ, " otscon otslin ConInterOut outlin are%11.4e%11.4e%11.4e%11.4e\n", 
					rfield.otscon[i], rfield.otslin[i]+rfield.outlin_noplot[i], rfield.ConInterOut[i], 
					rfield.outlin[0][i]+rfield.outlin_noplot[i] );
					fprintf( ioQQQ, " line continuum here are %4.4s %4.4s\n", 
					rfield.chLineLabel[i].c_str(), rfield.chContLabel[i].c_str() );
				}
			}
		}

		phisig += rfield.flux[0][i];
		/* >>chng 02 sep 19, add sec test on SummedDif since it can be zero when 
		 * phisig is just above small float */
		if( phisig > 0. && rfield.SummedDif[i]> 0. )
		{
			if( fabs(rfield.SummedCon[i]/phisig-1.) > 1e-3 )
			{
				++nInsane;
				/* limit amount of printout - in many failures would print entire
				 * continuum array */
				if( nInsane < LIM_INSAME_PRT )
				{
					fprintf( ioQQQ, " PROBLEM RT_OTS_ChkSum %3ld, insane SummedCon at energy %.5e error=%.2e i=%ld\n", 
					ipPnt, rfield.anu(i), rfield.SummedCon[i]/phisig - 1., i );
					fprintf( ioQQQ, " SummedCon, sum are %.4e %.4e\n", 
					rfield.SummedCon[i], phisig );
					fprintf( ioQQQ, " otscon otslin ConInterOut outlin flux are%.4e %.4e %.4e %.4e %.4e\n", 
					rfield.otscon[i], rfield.otslin[i]+rfield.outlin_noplot[i], rfield.ConInterOut[i], 
					rfield.outlin[0][i]+rfield.outlin_noplot[i], rfield.flux[0][i] );
					fprintf( ioQQQ, " line continuum here are %s %s\n", 
					rfield.chLineLabel[i].c_str(), rfield.chContLabel[i].c_str()
					);
				}
			}
		}
	}

	if( nInsane > 0 )
	{
		fprintf( ioQQQ, " PROBLEM RT_OTS_ChkSum too much insanity to continue.\n");
		/* TotalInsanity exits after announcing a problem */
		TotalInsanity();
	}
	return;
}

/* =================================================================== */

/*RT_OTS_PrtRate print continuum and line ots rates when trace ots is on */
void RT_OTS_PrtRate(
	  /* weakest rate to print */
	  double weak ,
	  /* flag, 'c' continuum, 'l' line, 'b' both */
	  int chFlag )
{
	long int i;

	DEBUG_ENTRY( "RT_OTS_PrtRate()" );

	/* arg must be one of these three */
	ASSERT( chFlag=='l' || chFlag=='c' || chFlag=='b' );

	/* 
	 * both printouts have cell number (on C array scale)
	 * energy in ryd
	 * the actual value of the ots rate 
	 * the ots rate relative to the continuum at that energy 
	 * rate times opacity 
	 * all are only printed if greater than weak 
	 */

	/*===================================================================*/
	/* first print ots continua                                          */
	/*===================================================================*/
	if( chFlag == 'c' || chFlag == 'b' )
	{
		fprintf( ioQQQ, "     DEBUG OTSCON array, anu, otscon, opac, OTS*opac limit:%.2e zone:%.2f IonConv?%c\n",
					weak,fnzone ,TorF(conv.lgConvIoniz()) );

		for( i=0; i < rfield.nflux_with_check; i++ )
		{
			if( rfield.otscon[i]*opac.opacity_abs[i] > weak )
			{
				fprintf( ioQQQ, "     %4ld%12.4e%12.4e%12.4e%12.4e %s \n", 
				  i, 
				  rfield.anu(i), 
				  rfield.otscon[i], 
				  opac.opacity_abs[i], 
				  rfield.otscon[i]*opac.opacity_abs[i], 
				  rfield.chContLabel[i].c_str());

			}
		}
	}

	/*===================================================================*/
	/* second print ots line rates                                       */
	/*===================================================================*/
	if( chFlag == 'l' || chFlag == 'b' )
	{
		fprintf( ioQQQ, "DEBUG density He %.2e He+2 %.2e O+2 %.2e\n",
			dense.gas_phase[ipHELIUM] , dense.xIonDense[ipHELIUM][2],
			dense.xIonDense[ipOXYGEN][2] );
		fprintf( ioQQQ, "     DEBUG OTSLIN array, anu, otslin, opac, OTS*opac Lab nLine limit:%.2e zone:%.2f IonConv?%c\n",
					weak,fnzone,TorF(conv.lgConvIoniz()) );

		for( i=0; i < rfield.nflux_with_check; i++ )
		{
			if( rfield.otslin[i]*opac.opacity_abs[i] > weak )
			{
				fprintf( ioQQQ, "     %4ld%12.4e%12.4e%12.4e%12.4e %s %3li\n", 
				  i, 
				  rfield.anu(i), 
				  rfield.otslin[i], 
				  opac.opacity_abs[i],
				  rfield.otslin[i]*opac.opacity_abs[i],
				  rfield.chLineLabel[i].c_str() ,
				  rfield.line_count[i] );
			}
		}
	}
	return;
}
