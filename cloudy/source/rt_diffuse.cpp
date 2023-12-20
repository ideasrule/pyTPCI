/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_diffuse evaluate local diffuse emission for this zone,
 * fill in ConEmitLocal[depth][energy] with diffuse emission,
 * called by Cloudy, this routine adds energy to the outward beam 
 * OTS rates for this zone were set in RT_OTS - not here */
#include "cddefines.h"
#include "taulines.h"
#include "grains.h"
#include "grainvar.h"
#include "iso.h"
#include "dense.h"
#include "opacity.h"
#include "trace.h"
#include "coolheavy.h"
#include "rfield.h"
#include "phycon.h"
#include "hmi.h"
#include "radius.h"
#include "atmdat.h"
#include "heavy.h"
#include "h2.h"
#include "rt.h"
#include "freebound.h"
#include "two_photon.h"
#include "lines_service.h"
#include "atmdat_gaunt.h"
#include "vectorize.h"
#include "ipoint.h"

STATIC void RT_iso_integrate_RRC( );

void RT_diffuse(void)
{
	/* arrays used in this routine
	 * rfield.ConEmitLocal[depth][energy] local emission per unit vol 
	 * rfield.DiffuseEscape is the spectrum of diffuse emission that escapes this zone,
	 * at end of this routine part is thrown into the outward beam
	 * by adding to rfield.ConInterOut
	 * units are photons s-1 cm-3 
	 * one-time init done on first call */

	/* rfield.DiffuseEscape and rfield.ConEmitLocal are same except that 
	 * rfield.ConEmitLocal is local emission, would be source function if div by opac
	 * rfield.DiffuseEscape is part that escapes so has RT built into it 
	 * rfield.DiffuseEscape is used to define rfield.ConInterOut below as per this statement
	 * rfield.ConInterOut[nu] += rfield.DiffuseEscape[nu]*(realnum)radius.dVolOutwrd;
	 */
	/* \todo	0	define only rfield.ConEmitLocal as it is now done, 
	 * do not define rfield.DiffuseEscape at all
	 * at bottom of this routine use inward and outward optical depths to define
	 * local and escaping parts 
	 * this routine only defines 
	 * rfield.ConInterOut - set to rfield.DiffuseEscape times vol element 
	 * so this is only var that
	 * needs to be set 
	 */

	long int ip=-100000, 
	  ipla=-100000,
	  limit=-100000, 
	  nu=-10000;

	double EdenAbund, 
	  difflya, 
	  fac, 
	  factor, 
	  gamma, 
	  gion, 
	  gn, 
	  photon; 

	DEBUG_ENTRY( "RT_diffuse()" );

	/* many arrays were allocated to nupper, and we will add unit flux to [nflux] -
	 * this must be true to work */
	ASSERT( rfield.nflux < rfield.nflux_with_check );

	/* this routine evaluates the local diffuse fields
	 * it fills in all of the following vectors */

	realnum *conEmitZone = &rfield.ConEmitLocal[nzone][0];

	vzero(rfield.DiffuseEscape);
	memset(conEmitZone, 0, (unsigned)rfield.nflux_with_check*sizeof(conEmitZone[0]));
	vzero(rfield.TotDiff2Pht);
	vzero(rfield.DiffuseLineEmission);

	// calculate recombination spectra and cooling
	RT_iso_integrate_RRC();

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* >>chng 01 sep 23, rewrote for iso sequences */
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			// calculate recombination spectra and cooling
			// RT_iso_integrate_RRC( ipISO, nelem );
			/* the product of the densities of the parent ion and electrons */
			EdenAbund = dense.eden*dense.xIonDense[nelem][nelem+1-ipISO];

			/* recombination continua for all iso seq - 
			 * if this stage of ionization exists */
			if( dense.IonHigh[nelem] >= nelem+1-ipISO  )
			{
				t_iso_sp* sp = &iso_sp[ipISO][nelem];

				// add line emission from the model iso atoms
				for( long ipHi=1; ipHi < sp->numLevels_local; ipHi++ )
				{
					for( long ipLo=0; ipLo < ipHi; ipLo++ )
					{
						// skip non-radiative transitions
						if( sp->trans(ipHi,ipLo).ipCont() < 1 )
							continue;

						/* number of photons in the line has not been defined up until now,
						 * do so now.  this is redone in lines.  */
						sp->trans(ipHi,ipLo).Emis().xIntensity() = 
							sp->trans(ipHi,ipLo).Emis().Aul()*
							sp->st[ipHi].Pop()*
							sp->trans(ipHi,ipLo).Emis().Pesc_total() *
							sp->trans(ipHi,ipLo).EnergyErg();

						// Would be better to enable checks (and remove argument) --
						// present state is to ensure backwards compatibility with previous
						// unchecked code.
						// First argument is fraction of line not emitted by scattering --
						// would be better to do this on the basis of line physics rather than
						// fiat...
						const bool lgDoChecks = false;
						sp->trans(ipHi,ipLo).outline(1.0, lgDoChecks );
					}
				}

				/*Iso treatment of two photon emission.  */
				/* NISO could in the future be increased, but we want this assert to blow
				 * so that it is understood this may not be correct for other iso sequences,
				 * probably should break since will not be present */
				ASSERT( ipISO <= ipHE_LIKE );

				/* upper limit to 2-phot is energy of 2s to ground */
				for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
				{
					CalcTwoPhotonEmission( *tnu, rfield.lgInducProcess && iso_ctrl.lgInd2nu_On );

					for( nu=0; nu < tnu->ipTwoPhoE; nu++ )
					{
						/* information - only used in save output */
						rfield.TotDiff2Pht[nu] += tnu->local_emis[nu];

						/* total local diffuse emission */
						conEmitZone[nu] += tnu->local_emis[nu];

						/* this is escaping part of two-photon emission, 
						 * as determined from optical depth to illuminated face */
						rfield.DiffuseEscape[nu] += tnu->local_emis[nu] * opac.ExpmTau[nu];
					}
					enum {DEBUG_LOC=false};
					if( DEBUG_LOC )
					{
						fprintf( ioQQQ, "Two-photon emission coefficients - ipISO, nelem = %2li, %2li\n", ipISO, nelem );
						PrtTwoPhotonEmissCoef( *tnu, EdenAbund );
					}
				}
			}
		}
	}

	/* add recombination continua for elements heavier than those done with iso seq */
	for( long nelem=NISO; nelem < LIMELM; nelem++ )
	{
		// zero out all stages since dense.IonLow[nelem] may have been lower last time around
		for( long ion=0; ion < nelem-NISO+1; ion++ )
		{
			Heavy.RadRecCon[nelem][ion] = 0.;
		}

		/* do not include species with iso-sequence in following */
		/* >>chng 03 sep 09, upper bound was wrong, did not include NISO */
		for( long ion=dense.IonLow[nelem]; ion < nelem-NISO+1; ion++ )
		{
			if( dense.xIonDense[nelem][ion+1] > 0. )
			{
				long int ns, nshell,igRec , igIon,
					iplow , iphi , ipop;

				ip = Heavy.ipHeavy[nelem][ion]-1;
				ASSERT( ip >= 0 );

				/* nflux was reset upward in ConvInitSolution to encompass all
				 * possible line and continuum emission.  this test should not
				 * possibly fail.  It could if the ionization were to increase with depth
				 * although the continuum mesh is designed to deal with this.
				 * This test is important because the nflux cell in ConInterOut 
				 * is used to carry out the unit integration, and if it gets 
				 * clobbered by diffuse emission the code will declare 
				 * insanity in PrtComment */
				if( ip >= rfield.nflux )
					continue;

				/* get shell number, stat weights for this species */
				atmdat_outer_shell( nelem+1 , nelem+1-ion , &nshell, &igRec , &igIon );
				gn = (double)igRec;
				gion = (double)igIon;

				/* shell number */
				ns = Heavy.nsShells[nelem][ion]-1;
				ASSERT( ns == (nshell-1) );

				/* lower and upper energies, and offset for opacity stack */
				iplow = opac.ipElement[nelem][ion][ns][0]-1;
				iphi = opac.ipElement[nelem][ion][ns][1];
				iphi = MIN2( iphi , rfield.nflux );
				ipop = opac.ipElement[nelem][ion][ns][2];

				/* now convert ipop to the offset in the opacity stack from threshold */
				ipop = ipop - iplow;

				EdenAbund = dense.eden*dense.xIonDense[nelem][ion+1];
				gamma = 0.5*MILNE_CONST*gn/gion/phycon.te/phycon.sqrte;

				/* this is ground state continuum from stored opacities */
				if( rfield.ContBoltz[iplow] > SMALLFLOAT )
				{
					for( nu=iplow; nu < iphi; ++nu )
					{
						photon = gamma*rfield.ContBoltz[nu]/rfield.ContBoltz[iplow]*
							rfield.widflx(nu)*opac.OpacStack[nu+ipop]*rfield.anu2(nu);
						/* add heavy rec to ground in active beam,*/
						/** \todo	2	should use ConEmitLocal for all continua, but not followed
						 * by rfield.DiffuseEscape - put that at the end.  Once continua all
						 * bundled this way, it will be easy to save them as a function
						 * of depth and then do exact rt */
						conEmitZone[nu] += (realnum)photon*EdenAbund;
						rfield.DiffuseEscape[nu] += (realnum)photon*EdenAbund*opac.ExpmTau[nu];

						// escaping RRC
						Heavy.RadRecCon[nelem][ion] += rfield.anu(nu) *
						  	emergent_line( photon*EdenAbund/2. , photon*EdenAbund/2. , 
							// energy on fortran scale
							nu+1 );
					}
				}
				// units erg cm-3 s-1
				Heavy.RadRecCon[nelem][ion] *= EN1RYD;

				/* now do the recombination Lya */
				ipla = Heavy.ipLyHeavy[nelem][ion]-1;
				ASSERT( ipla >= 0 );
				/* xLyaHeavy is set to a fraction of the total rad rec in ion_recomb, includes eden */
				difflya = Heavy.xLyaHeavy[nelem][ion]*dense.xIonDense[nelem][ion+1];
				rfield.DiffuseLineEmission[ipla] += (realnum)difflya;

				/* >>chng 03 jul 10, here and below, use outlin_noplot */
				rfield.outlin_noplot[ipla] += (realnum)(difflya*radius.dVolOutwrd*opac.tmn[ipla]*opac.ExpmTau[ipla]);

				/* now do the recombination Balmer photons */
				ipla = Heavy.ipBalHeavy[nelem][ion]-1;
				ASSERT( ipla >= 0 );
				/* xLyaHeavy is set to fraction of total rad rec in ion_recomb, includes eden */
				difflya = Heavy.xLyaHeavy[nelem][ion]*dense.xIonDense[nelem][ion+1];
				rfield.outlin_noplot[ipla] += (realnum)(difflya*radius.dVolOutwrd*opac.tmn[ipla]*opac.ExpmTau[ipla]);
			}
		}
	}

	/* free-free free free brems emission for all ions */
	limit = MIN2( rfield.ipMaxBolt , rfield.nflux );

	if( CoolHeavy.lgFreeOn )
	{
		t_brems_den sum;
		t_gaunt::Inst().brems_sum_ions(sum);
		// there is no need to keep these separate here...
		sum.den_ion[1] += sum.den_Hp + sum.den_Hep;
		sum.den_ion[2] += sum.den_Hepp;

		double fac = dense.eden * FREE_FREE_EMIS / phycon.sqrte;
	
		vector<double> TotBrems( limit );
		/* First add H- brems.  Reaction is H(1s) + e -> H(1s) + e + hnu. */
		t_gaunt::Inst().brems_rt( -1, phycon.te, fac*sum.den_Hm, TotBrems );
		/* chng 02 may 16, by Ryan...do all brems for all ions in one fell swoop,
		 * using gaunt factors from t_gaunt.gff.	*/
		for( long ion=1; ion < LIMELM+1; ++ion )
			if( sum.den_ion[ion] > 0. )
				t_gaunt::Inst().brems_rt( ion, phycon.te, fac*sum.den_ion[ion], TotBrems );

		for( nu=0; nu < limit; nu++ )
		{
			/* >>chng 05 feb 20, move into this test on brems opacity - should not be
			 * needed since would use expmtau to limit outward beam */
			/* >>chng 01 jul 01, move thick brems back to ConEmitLocal but do not add
			* to outward beam - ConLocNoInter array removed as result
			* if problems develop with very dense BLR clouds, this may be reason */
			conEmitZone[nu] += realnum(TotBrems[nu]/rfield.anu(nu)) +
				pow2(dense.eden)*rfield.eeBremsDif[nu];

			/* do not add optically thick part to outward beam since self absorbed
			* >>chng 96 feb 27, put back into outward beam since do not integrate
			* over it anyway. */
			/* >>chng 99 may 28, take back out of beam since DO integrate over it
			* in very dense BLR clouds */
			/* >>chng 01 jul 10, add here, in only one loop, where optically thin */
			rfield.DiffuseEscape[nu] += realnum(TotBrems[nu]/rfield.anu(nu)) +
				pow2(dense.eden)*rfield.eeBremsDif[nu];
		}
	}

	/* grain dust emission */
	/* >>chng 01 nov 22, moved calculation of grain flux to qheat.c, PvH */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		/* this calculates diffuse emission from grains,
		 * and stores the result in gv.GrainEmission */
		GrainMakeDiffuse();

		for( nu=0; nu < rfield.nflux; nu++ )
		{
			conEmitZone[nu] += gv.GrainEmission[nu];
			rfield.DiffuseEscape[nu] += gv.GrainEmission[nu];
		}
	}

	/* hminus emission */
	fac = dense.eden*(double)dense.xIonDense[ipHYDROGEN][0];
	gn = 1.;
	gion = 2.;
	gamma = 0.5*MILNE_CONST*gn/gion/phycon.te/phycon.sqrte;
	/* >>chng 00 dec 15 change limit to -1 of H edge */
	limit = ipoint(atmdat.EIonPot[ipHELIUM][0]);

	if( rfield.ContBoltz[hmi.iphmin-1] > 0. )
	{
		for( nu=hmi.iphmin-1; nu < limit; nu++ )
		{
			/* H- flux photons cm-3 s-1 
			 * ContBoltz is ratio of Boltzmann factor for each freq */
			factor = gamma*rfield.ContBoltz[nu]/rfield.ContBoltz[hmi.iphmin-1]*rfield.widflx(nu)*
			  opac.OpacStack[nu-hmi.iphmin+opac.iphmop]*
			  rfield.anu2(nu)*fac;
			conEmitZone[nu] += (realnum)factor;
			rfield.DiffuseEscape[nu] += (realnum)factor;
		}
	}
	else
	{
		for( nu=hmi.iphmin-1; nu < limit; nu++ )
		{
			double arg = MAX2(0.,TE1RYD*(rfield.anu(nu)-HMINUSIONPOT)/phycon.te);
			/* this is the limit sexp normally uses */
			if( arg > SEXP_LIMIT ) 
				break;
			/* H- flux photons cm-3 s-1 
			 * flux is in photons per sec per ryd */
			factor = gamma*exp(-arg)*rfield.widflx(nu)*
				opac.OpacStack[nu-hmi.iphmin+opac.iphmop]*
			  rfield.anu2(nu)*fac;
			conEmitZone[nu] += (realnum)factor;
			rfield.DiffuseEscape[nu] += (realnum)factor;
		}
	}

	for( long ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] && iso_ctrl.lgDielRecom[ipISO] ) 
			{
				for( long i=0; i<iso_sp[ipISO][nelem].numLevels_local; i++ )
				{
					const TransitionList::iterator& tr = SatelliteLines[ipISO][nelem].begin()+ipSatelliteLines[ipISO][nelem][i];
					(*tr).Emis().xIntensity() = 
						(*tr).Emis().Aul()*
						(*(*tr).Hi()).Pop()*
						(*tr).Emis().Pesc_total()*
						(*tr).EnergyErg();

					(*tr).outline_resonance();
				}
			}
		}
	}

	/* outward level 2 line photons */
	for( long i=0; i < nWindLine; i++ )
	{
		/* must not also do lines that were already done as part
		 * of the isoelectronic sequences */
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			{
				enum {DEBUG_LOC=false};
				if( DEBUG_LOC /*&& nzone > 10*/ && i==4821 )
				{
					/* set up to dump the Fe 9 169A line */
					fprintf(ioQQQ,"DEBUG dump lev2 line %li\n", i );
					DumpLine( TauLine2[i] );/**/
					fprintf(ioQQQ,"DEBUG dump %.3e %.3e %.3e\n",
						rfield.outlin[0][TauLine2[i].ipCont()-1],
						phots( TauLine2[i] )*TauLine2[i].Emis().FracInwd()*radius.BeamInOut*opac.tmn[i]*TauLine2[i].Emis().ColOvTot(),
						phots( TauLine2[i] )*(1. - TauLine2[i].Emis().FracInwd())*radius.BeamOutOut* TauLine2[i].Emis().ColOvTot() );
				}
			}
			TauLine2[i].outline_resonance();
			/*if( i==2576 ) fprintf(ioQQQ,"DEBUG dump %.3e %.3e \n",
				rfield.outlin[0][TauLine2[i].ipCont()-1] , rfield.outlin_noplot[TauLine2[i].ipCont()-1]);*/
		}
	}

	/* outward hyperfine structure line photons */
	for( size_t i=0; i < HFLines.size(); i++ )
	{
		HFLines[i].outline_resonance();
	}

	/* external database lines */
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( dBaseSpecies[ipSpecies].lgActive )
		{
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{	
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local || (*tr).ipCont() <= 0)
					continue;
				(*tr).outline_resonance();
			}
		}
	}

	/* H2 emission */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_RT_diffuse();

	/** \todo	2	add fegrain to outward beams, but within main formalism by including grains
	 * in all x-ray processes */

	if( trace.lgTrace )
		fprintf( ioQQQ, " RT_diffuse returns.\n" );

	/* >>chng 02 jul 25, zero out all light below plasma freq */
	for( nu=0; nu < rfield.ipPlasma-1; nu++ )
	{
		rfield.flux_beam_const[nu] = 0.;
		rfield.flux_beam_time[nu] = 0.;
		rfield.flux_isotropic[nu] = 0.;
		rfield.flux[0][nu] = 0.;
		conEmitZone[nu] = 0.;
		rfield.otscon[nu] = 0.;
		rfield.otslin[nu] = 0.;
		rfield.outlin[0][nu] = 0.;
		rfield.outlin_noplot[nu] = 0.;
		rfield.reflin[0][nu] = 0.;
		rfield.TotDiff2Pht[nu] = 0.;
		rfield.ConInterOut[nu] = 0.;
	}

	/* find occupation number, also assert that no continua are negative */
	for( nu=0; nu < rfield.nflux; nu++ )
	{
		/* >>chng 00 oct 03, add diffuse continua */
		/* local diffuse continua */
		rfield.OccNumbDiffCont[nu] = 
			conEmitZone[nu]*rfield.convoc[nu];

		/* units are photons cell-1 cm-2 s-1  */
		rfield.ConSourceFcnLocal[nzone][nu] = 
			/* units photons cm-3 s-1 cell-1, */
			safe_div( conEmitZone[nu],
			/* units cm-1 */
			(realnum)opac.opacity_abs[nu],
			0_r );

		/* confirm that all are non-negative */
		ASSERT( rfield.flux_beam_const[nu] >= 0.);
		ASSERT( rfield.flux_beam_time[nu] >= 0.);
		ASSERT( rfield.flux_isotropic[nu] >= 0.);
		ASSERT( rfield.flux[0][nu] >= 0.);
		ASSERT( conEmitZone[nu] >= 0.);
		ASSERT( rfield.otscon[nu] >= 0.);
		ASSERT( rfield.otslin[nu] >= 0.);
		ASSERT( rfield.outlin[0][nu] >= 0.);
		ASSERT( rfield.outlin_noplot[nu] >= 0.);
		ASSERT( rfield.reflin[0][nu] >= 0.);
		ASSERT( rfield.TotDiff2Pht[nu] >= 0.);
		ASSERT( rfield.ConInterOut[nu] >= 0.);
	}

	/* option to kill outward lines with no outward lines command*/
	if( rfield.lgKillOutLine )
	{
		for( nu=0; nu < rfield.nflux; nu++ )
		{
			rfield.outlin[0][nu] = 0.;
			rfield.outlin_noplot[nu] = 0.;
		}
	}

	/* option to kill outward continua with no outward continua command*/
	if( rfield.lgKillOutCont )
	{
		for( nu=0; nu < rfield.nflux; nu++ )
		{
			rfield.ConInterOut[nu] = 0.;
		}
	}
	return;
}

STATIC void RT_iso_integrate_RRC()
{
	DEBUG_ENTRY( "RT_iso_integrate_RRC()" );

	realnum* conEmitZone = &rfield.ConEmitLocal[nzone][0];
	/* loop over iso-sequences of all elements 
	 * to add all recombination continua and lines*/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* >>chng 01 sep 23, rewrote for iso sequences */
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			ASSERT( nelem >= ipISO );
			ASSERT( nelem < LIMELM );
			
			// recombination continua for all iso seq - 
			// if this stage of ionization exists 
			if( dense.IonHigh[nelem] < nelem+1-ipISO  )
				continue;
			/* this will be the sum of recombinations to all excited levels */
			double SumCaseB = 0.;
			
			/* the product of the densities of the parent ion and electrons */
			double EdenAbund = dense.eden*dense.xIonDense[nelem][nelem+1-ipISO];
			
			t_iso_sp* sp = &iso_sp[ipISO][nelem];

			avx_ptr<double> arg(sp->numLevels_local), val(sp->numLevels_local);
			for( long n=0; n < sp->numLevels_local; n++ )
			{
				long ipLo = sp->fb[n].ipIsoLevNIonCon-1;
				double thresh = sp->fb[n].xIsoLevNIonRyd;
				double widflx = rfield.anumax(ipLo) - thresh;
				arg[n] = -widflx/phycon.te_ryd;
			}
			vexpm1( arg.ptr0(), val.ptr0(), 0, sp->numLevels_local );
			
			// loop over all levels to include recombination diffuse continua,
			// pick highest energy continuum point that opacities extend to 
			long ipHi = rfield.nflux;
			// >>chng 06 aug 17, should go to numLevels_local instead of _max. 
			for( long n=0; n < sp->numLevels_local; n++ )
			{
				double Sum1level = 0.;
				double RadRecCon = 0.;
				// the number is (2 pi me k/h^2) ^ -3/2 * 8 pi/c^2 / ge - it includes
				// the stat weight of the free electron in the demominator 
				double gamma = 0.5*MILNE_CONST*sp->st[n].g()/iso_ctrl.stat_ion[ipISO]/phycon.te/phycon.sqrte;
				
				// loop over all recombination continua 
				// escaping part of recombinations are added to rfield.ConEmitLocal 
				// added to ConInterOut at end of routine
				long ipLo = sp->fb[n].ipIsoLevNIonCon-1;
				double thresh = sp->fb[n].xIsoLevNIonRyd;
				long offset = -sp->fb[n].ipIsoLevNIonCon + sp->fb[n].ipOpac;
				double RadRecomb = sp->fb[n].RadRecomb[ipRecEsc];
				ASSERT( rfield.anumin(ipLo) <= thresh && thresh < rfield.anumax(ipLo) );
				double efac = 0.;

				fixit("need to include induced recombination in diffuse spectrum.");
				// Probably best just to do the whole thing here
				// then no need for induced terms in iso_photo.

				for( long nu=ipLo; nu < ipHi; nu++ )
				{
					double bfac;
					if( nu == ipLo )
					{
						double widflx = rfield.anumax(ipLo) - thresh;
						// fac is the fraction of the cell width that is above threshold
						double fac = widflx/rfield.widflx(ipLo);
						bfac = -fac*phycon.te_ryd*val[n]/widflx;
						efac = exp(-widflx/phycon.te_ryd);
					}
					else
					{
						// efac = exp(-(rfield.anumin(nu)-thresh)/phycon.te_ryd);
						bfac = efac*rfield.ContBoltzHelp1[nu];
						efac *= rfield.ContBoltzHelp2[nu];
					}
					/* photon is in photons cm^3 s^-1 per cell */
					double photon = gamma*bfac*rfield.widflx(nu)*
						opac.OpacStack[nu+offset] * rfield.anu2(nu);
					
					Sum1level += photon;

					// no point in wasting time on adding tiny numbers...
					if( photon/Sum1level < 1.e-20 )
						break;
					
					/* total local diffuse emission units photons cm-3 s-1 cell-1,*/
					conEmitZone[nu] += (realnum)(photon*EdenAbund);
					
					// sp->fb[n].RadRecomb[ipRecEsc] is escape probability
					// rfield.DiffuseEscape is local emission that escapes this zone 
					rfield.DiffuseEscape[nu] += (realnum)(photon*EdenAbund*RadRecomb);
					
					// total RRC radiative recombination continuum, pointer must be on fortran scale
					RadRecCon += rfield.anu(nu) *
						emergent_line( photon*EdenAbund/2., photon*EdenAbund/2., nu+1 );
				}
				
				// convert to erg cm-3 s-1
				sp->fb[n].RadRecCon = RadRecCon*EN1RYD;
				/* this will be used below to confirm case B sum */
				if( n > 0 )
				{
					/* SumCaseB will be sum to all excited */
					SumCaseB += Sum1level;
				}
			}
			
			/* this is check on self-consistency */
			sp->CaseBCheck = MAX2(sp->CaseBCheck, (realnum)(SumCaseB/sp->RadRec_caseB));
		}
	}
	
	return;
}
