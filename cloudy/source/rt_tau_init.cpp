/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_tau_init set initial outward optical depths at start of first iteration,
 * it is only called by cloudy one time per complete calculation, just after
 * continuum set up and before start of convergence attempts. 
 * 
 * RT_tau_reset after first iteration, updates the optical depths, mirroring this
 * routine but with the previous iteration's variables */
#include "cddefines.h"
#include "taulines.h"
#include "doppvel.h"
#include "iso.h"
#include "h2.h"
#include "rfield.h"
#include "dense.h"
#include "opacity.h"
#include "thermal.h"
#include "geometry.h"
#include "stopcalc.h"
#include "ipoint.h"
#include "conv.h"
#include "rt.h"
#include "trace.h"
#include "freebound.h"

static const realnum TAULIM = 1e8;

void RT_tau_init(void)
{
	long int nelem,
	  ipISO,
	  ipHi, 
	  ipLo,
	  nHi;

	bool lgBalmerTauOn;

	DEBUG_ENTRY( "RT_tau_init()" );

	ASSERT( dense.eden > 0. );

	/* Zero lines first */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				if( iso_ctrl.lgDielRecom[ipISO] )
				{
					// SatelliteLines are indexed by lower level
					for( ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_max; ipLo++ )
					{
						SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]].Zero();
					}
				}

				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Zero();
					}
				}
				for( ipHi=2; ipHi <iso_ctrl.nLyman[ipISO]; ipHi++ )
				{
					ExtraLymanLines[ipISO][nelem][ipExtraLymanLines[ipISO][nelem][ipHi]].Zero();
				}
			}
		}
	}

	/* this is a dummy optical depth array for non-existant lines 
	 * when this goes over to struc, make sure all are set to zero here since
	 * init in grid may depend on it */
	(*TauDummy).Zero();

	/* lines in cooling function with Mewe approximate collision strengths */
	for( long i=0; i < nWindLine; i++ )
	{
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			/* inward optical depth */
			TauLine2[i].Zero();
		}
	}

	/* inner shell lines */
	for( size_t i=0; i < UTALines.size(); i++ )
	{
		/* these are line optical depth arrays
		 * inward optical depth */
		/* heat is special for this array - it is heat per pump */
		double hsave = UTALines[i].Coll().heat();
		UTALines[i].Zero();
		UTALines[i].Coll().heat() = hsave;
	}

	/* hyperfine structure lines */
	for( size_t i=0; i < HFLines.size(); i++ )
	{
		HFLines[i].Zero();
	}

	/* external database lines */
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		// can't filter by lgActive, must do all to properly reset in grids
		//if( dBaseSpecies[ipSpecies].lgActive )
		{
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				(*tr).Zero();
			}
		}
	}

	/* initialize optical depths in H2 */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_LineZero();

	/* ==================================================================*/
	/* end setting lines to zero */

	/* >>chng 02 feb 08, moved to here from opacitycreateall, which was called later.
	 * bug inhibited convergence in some models.  Caught by PvH */
	/* this is option to stop at certain optical depth */
	if( StopCalc.taunu > 0. )
	{
		StopCalc.iptnu = ipoint(StopCalc.taunu);
		StopCalc.iptnu = MIN2(StopCalc.iptnu,rfield.nflux_with_check-1);
	}
	else
	{
		StopCalc.iptnu = rfield.nflux_with_check;
		/* >>chng 04 dec 21, remove from here and init to 1e30 in zero */
		/*StopCalc.tauend = 1e30f;*/
	}

	/* if taunu was set with a stop optical depth command, then iptnu must
	 * have also been set to a continuum index before this code is reaced */
	ASSERT( StopCalc.taunu == 0. || StopCalc.iptnu >= 0 );

	/* set initial and total optical depths in arrays
	 * TAUNU is set when lines read in, sets stopping radius */
	if( StopCalc.taunu > 0. )
	{
		/*  an optical depth has been specified */
		if( StopCalc.iptnu >= iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon )
		{
			/* at ionizing energies */
			for( long i=0; i < (iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - 1); i++ )
			{
				/* taumin can be reset with taumin command */
				opac.TauAbsGeo[1][i] = opac.taumin;
				opac.TauScatGeo[1][i] = opac.taumin;
				opac.TauTotalGeo[1][i] = opac.taumin;
			}

			for( long i=iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1; i < rfield.nflux_with_check; i++ )
			{
				/* TauAbsGeo(i,2) = tauend * (anu(i)/anu(iptnu))**(-2.43) */
				opac.TauAbsGeo[1][i] = StopCalc.tauend;
				opac.TauScatGeo[1][i] = opac.taumin;
				opac.TauTotalGeo[1][i] = opac.TauAbsGeo[1][i] + opac.TauScatGeo[1][i];
			}
		}

		else
		{
			/* not specified at ionizing energies */
			for( long i=0; i < (iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - 1); i++ )
			{
				opac.TauAbsGeo[1][i] = StopCalc.tauend;
				opac.TauScatGeo[1][i] = StopCalc.tauend;
				opac.TauTotalGeo[1][i] = StopCalc.tauend;
			}

			for( long i=iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1; i < rfield.nflux_with_check; i++ )
			{
				opac.TauAbsGeo[1][i] = (realnum)(TAULIM*pow(rfield.anu(i),-2.43));
				opac.TauScatGeo[1][i] = opac.taumin;
				opac.TauTotalGeo[1][i] = opac.TauAbsGeo[1][i] + opac.TauScatGeo[1][i];
			}

		}
	}

	else
	{
		/*  ending optical depth not specified, assume 1E8 at 1 Ryd */
		for( long i=0; i < (iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - 1); i++ )
		{
			opac.TauAbsGeo[1][i] = opac.taumin;
			opac.TauScatGeo[1][i] = opac.taumin;
			opac.TauTotalGeo[1][i] = opac.taumin;
		}

		for( long i=iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1; i < rfield.nflux_with_check; i++ )
		{
			opac.TauAbsGeo[1][i] = (realnum)(TAULIM*pow(rfield.anu(i),-2.43));
			opac.TauScatGeo[1][i] = opac.taumin;
			opac.TauTotalGeo[1][i] = opac.TauAbsGeo[1][i] + opac.TauScatGeo[1][i];
		}
	}

	/* if lgSphere then double outer, set inner to half
	 * assume will be optically thin at sub-ionizing energies */
	if( geometry.lgSphere || opac.lgCaseB )
	{
		for( long i=0; i < (iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - 1); i++ )
		{
			opac.TauAbsGeo[0][i] = opac.taumin;
			opac.TauAbsGeo[1][i] = opac.taumin*2.f;
			opac.TauScatGeo[0][i] = opac.taumin;
			opac.TauScatGeo[1][i] = opac.taumin*2.f;
			opac.TauTotalGeo[0][i] = 2.f*opac.taumin;
			opac.TauTotalGeo[1][i] = 4.f*opac.taumin;
		}

		for( long i=iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1; i < rfield.nflux_with_check; i++ )
		{
			opac.TauAbsGeo[0][i] = opac.TauAbsGeo[1][i];
			opac.TauAbsGeo[1][i] *= 2.;
			opac.TauScatGeo[0][i] = opac.TauScatGeo[1][i];
			opac.TauScatGeo[1][i] *= 2.;
			opac.TauTotalGeo[0][i] = opac.TauTotalGeo[1][i];
			opac.TauTotalGeo[1][i] *= 2.;
		}

		if( StopCalc.taunu > 0. )
		{
			/* ending optical depth specified, and lgSphere also */
			StopCalc.tauend *= 2.;
		}
	}

	/* fix escape prob for He, metals, first set log of Te, needed by RECEFF
	 * do not do this if temperature has been set by command */
	if( !thermal.lgTemperatureConstant )
	{
		double TeNew;
		/* this is a typical temperature for the H+ zone, and will use it is
		 * the line widths & opt depth in the following.  this is not meant to be the first
		 * temperature during the search phase. */
		TeNew = 2e4;
		/* >>chng 05 jul 19, in PDR models the guess of the optical depth in Lya could
		 * be very good since often limiting column density is set, but must use
		 * the temperature of H0 gas.  in a dense BLR this is roughly 7000K and
		 * closer to 100K for a low-density PDR - use these guesses */
		if( dense.gas_phase[ipHYDROGEN] >= 1e9 )
		{
			/* this is a typical BLR H0 temp */
			TeNew = 7000.;
		}
		else if( dense.gas_phase[ipHYDROGEN] <= 1e5 )
		{
			/* this is a typical PDR H0 temp */
			TeNew = 100.;
		}
		else
		{
			/* power law interpolation between them */
			TeNew = 0.5012 * pow( (double)dense.gas_phase[ipHYDROGEN], 0.46 );
		}

		/* propagate this temperature through the code */
		/* must not call tfidle at this stage since not all vectors have been allocated */
		TempChange( TeNew );
	}

	SumDensities();
	ASSERT( dense.xNucleiTotal > 0. );	
	/* set inward optical depths for hydrogenic ions to small number proportional to abundance */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* now get actual optical depths */
			double AbunRatio = dense.gas_phase[nelem]/dense.xNucleiTotal;
			for( ipLo=ipH1s; ipLo < (iso_sp[ipH_LIKE][nelem].numLevels_max - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
				{
					/* set all inward optical depths to taumin, regardless of abundance */
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() = opac.taumin;
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() = 2.f*opac.taumin;
				}
			}

			/* La may be case B, tlamin set to taumin and reset with Case B
			 * command to 1e5.  Case A and C set it to 1e-5 */
			iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn() = opac.tlamin;

			/* scale factor so that all other Lyman lines are appropriate for this Lya optical depth*/
			realnum f = opac.tlamin/iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().opacity();
			fixit("this appears to be redundant to code below.");

			for( nHi=3; nHi<=iso_sp[ipH_LIKE][nelem].n_HighestResolved_max; nHi++ )
			{
				ipHi = iso_sp[ipH_LIKE][nelem].QN2Index(nHi, 1, 2);
				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() = MAX2( opac.taumin, 
					f*iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().opacity() );
			}
			for( ipHi=iso_sp[ipH_LIKE][nelem].numLevels_max - iso_sp[ipH_LIKE][nelem].nCollapsed_max; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
			{
				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() = MAX2( opac.taumin, 
					f*iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().opacity() );
			}

			/* after this set of if's the total Lya optical depth will be known,
			 * common code later on will set rest of Lyman lines
			 * if case b then set total Lyman to twice inner */
			if( opac.lgCaseB )
			{
				/* force outer optical depth to twice inner if case B */
				iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() = 
					(realnum)(2.*iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn());
				/* force off Balmer et al optical depths */
				lgBalmerTauOn = false;
			}

			else
			{
				/* usual case for H LyA; try to guess outer optical depth */
				if( StopCalc.colnut < 6e29 )
				{
					/* neutral column is defined */
					iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() = (realnum)(StopCalc.colnut*
					  rt.DoubleTau*iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().opacity()/
					  GetDopplerWidth(dense.AtomicWeight[nelem])*AbunRatio);
					ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() >= 0. );
				}
				/* has optical depth at some energy where we want to stop been specified?
				 * taunu is energy where
				 * this has been specified - this checks on Lyman continuum, taunu = 1 */
				else if( StopCalc.taunu < 3. && StopCalc.taunu >= 0.99 )
				{
					/* Lyman continuum optical depth defined */
					iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() = (realnum)(StopCalc.tauend*
					  1.2e4*1.28e6/GetDopplerWidth(dense.AtomicWeight[nelem])*rt.DoubleTau*AbunRatio);
					ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() >= 0. );
				}
				else if( StopCalc.HColStop < 6e29 )
				{
					/* neutral column not defined, guess from total col and U */
					double coleff = StopCalc.HColStop - 
						MIN2(MIN2(rfield.qhtot/dense.eden,1e24)/2.6e-13,0.6*StopCalc.HColStop);

					iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() = (realnum)(coleff*
					  7.6e-14*AbunRatio);
					ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() >= 0. );
				}
				else
				{
					/* no way to estimate 912 optical depth, set to large number */
					iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() = (realnum)(1e20*
					  AbunRatio);
					ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() >= 0. );
				}
				/* allow Balmer et al. optical depths */
				lgBalmerTauOn = true;
			}

			realnum TAddHLya = 0.f;
			bool lgLyaContinuumCorrection = false;
			if (lgLyaContinuumCorrection)
			{
				/* Lya total optical depth now known, is it optically thick?*/
				if (iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() > 1. )
				{
					TAddHLya = (realnum)MIN2(1.,iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot()/
													 1e4);
					iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn() += TAddHLya;
				}
				else
				{
					TAddHLya = opac.tlamin;
				}
			}

			/* this scale factor is to set other lyman lines, given the Lya optical depth */
			f = iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot()/
				iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().opacity();

			ipISO = ipH_LIKE;
			ASSERT( ipISO<NISO && nelem < LIMELM );
			for( nHi=3; nHi<=iso_sp[ipISO][nelem].n_HighestResolved_max; nHi++ )
			{
				ipHi = iso_sp[ipISO][nelem].QN2Index(nHi, 1, 2);
				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauTot() = MAX2( opac.taumin, 
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().opacity() * f );

				/* increment inward optical depths by rt all lyman lines, inward
				 * optical depth was set above, this adds to it.  this was originally
				 * set some place above */
				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() += TAddHLya*
				  iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().opacity()/
				  iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().opacity();

				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() = MAX2( 
					opac.taumin, iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() );
			}
			for( ipHi=iso_sp[ipH_LIKE][nelem].numLevels_max - iso_sp[ipH_LIKE][nelem].nCollapsed_max; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
			{
				/* set total optical depth for higher lyman lines */
				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauTot() = MAX2( opac.taumin, 
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().opacity() * f );

				/* increment inward optical depths by rt all lyman lines, inward
				 * optical depth was set above, this adds to it.  this was originally
				 * set some place above */
				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() += TAddHLya*
				  iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().opacity()/
				  iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().opacity();

				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() = MAX2( 
					opac.taumin, iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() );
			}

			/* try to guess what Balmer cont optical guess,
			 * first branch is case where we will stop at a balmer continuum optical
			 * depth - taunu is energy where tauend was set */
			if( StopCalc.taunu > 0.24 && StopCalc.taunu < 0.7 )
			{
				iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauTot() = (realnum)(StopCalc.tauend*
				  3.7e4*lgBalmerTauOn*AbunRatio + opac.taumin);

				iso_sp[ipH_LIKE][nelem].trans(ipH3s,ipH2p).Emis().TauTot() = (realnum)(StopCalc.tauend*
				  3.7e4*lgBalmerTauOn*AbunRatio + opac.taumin);

				iso_sp[ipH_LIKE][nelem].trans(ipH3d,ipH2p).Emis().TauTot() = (realnum)(StopCalc.tauend*
				  3.7e4*lgBalmerTauOn*AbunRatio + opac.taumin);
			}

			else
			{
				/* this is a guess based on Ferland&Netzer 1979, but it gets very large */
				double balc = rfield.qhtot*2.1e-19*lgBalmerTauOn*AbunRatio + opac.taumin;

				iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauTot() = 
					(realnum)(MIN2(2e5, balc*3.7e4*lgBalmerTauOn+opac.taumin));

				iso_sp[ipH_LIKE][nelem].trans(ipH3s,ipH2p).Emis().TauTot() = 
					(realnum)(MIN2(2e5, balc*3.7e4*lgBalmerTauOn+opac.taumin));

				iso_sp[ipH_LIKE][nelem].trans(ipH3d,ipH2p).Emis().TauTot() = 
					(realnum)(MIN2(2e5, balc*3.7e4*lgBalmerTauOn+opac.taumin));

				ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauTot() >= 0.);
				ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH3s,ipH2p).Emis().TauTot() >= 0.);
				ASSERT( iso_sp[ipH_LIKE][nelem].trans(ipH3d,ipH2p).Emis().TauTot() >= 0.);
			}

			/* transitions down to 2s are special since 'alpha' (2s-2p) has
			 * small optical depth, so use 3 to 2p instead */
			f = iso_sp[ipH_LIKE][nelem].trans(ipH3s,ipH2p).Emis().TauTot()/
				iso_sp[ipH_LIKE][nelem].trans(ipH3s,ipH2p).Emis().opacity()* lgBalmerTauOn;

			ipLo = ipH2s;
			for( ipHi=ipLo + 1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
			{
				if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
					continue;

				iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() = opac.taumin +
					f* iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().opacity();
				ASSERT(iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() >= 0.);
			}

			/* this is for rest of lines, scale from opacity */
			for( ipLo=ipH2p; ipLo < (iso_sp[ipH_LIKE][nelem].numLevels_max - 1); ipLo++ )
			{
				long ipNS, ipNPlusOneP;
#if 0
				/* scale everything with same factor we use for (n+1)P -> nS */
				long ipISO = ipH_LIKE;
				long n = N_(ipLo);
				ipNS = iso_sp[ipH_LIKE][nelem].QN2Index(n, 0, 2);
				ipNPlusOneP = iso_sp[ipH_LIKE][nelem].QN2Index(n+1, 1, 2);
				if( ipNPlusOneP < 0 )
					ipNPlusOneP = ipNS + 1;
#else
				/* old way */
				ipNS = ipLo;
				ipNPlusOneP = ipNS + 1;
#endif

				if( iso_sp[ipH_LIKE][nelem].trans(ipNPlusOneP,ipNS).ipCont() <= 0 )
				{	
					f = SMALLFLOAT;
				}
				else
				{
					f = iso_sp[ipH_LIKE][nelem].trans(ipNPlusOneP,ipNS).Emis().TauTot()/
					  iso_sp[ipH_LIKE][nelem].trans(ipNPlusOneP,ipNS).Emis().opacity()*
					  lgBalmerTauOn;
				}

				for( ipHi=ipLo + 1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
				{
					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
						continue;

					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() = opac.taumin +
						f* iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().opacity();
					ASSERT(iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() >= 0.);
				}
			}

			/* this loop is over all possible H lines, do some final cleanup */
			for( ipLo=ipH1s; ipLo < (iso_sp[ipH_LIKE][nelem].numLevels_max - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
				{
					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
						continue;

					/* TauCon is line optical depth to inner face used for continuum pumping rate,
					 * not equal to TauIn in case of static sphere since TauCon will never
					 * count the far side line optical depth */
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauCon() = 
						iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn();

					/* make sure inward optical depth is not larger than half total */
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() = 
						MIN2(	iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() ,
						iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot()/2.f );
					ASSERT(iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() >= 0.f);

					/* this is fraction of line that goes inward, not known until second iteration*/
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().FracInwd() = 0.5;
				}
			}
		}
	}

	/* initialize all he-like optical depths */
	for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			for( ipLo=0; ipLo < (iso_sp[ipHE_LIKE][nelem].numLevels_max - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < iso_sp[ipHE_LIKE][nelem].numLevels_max; ipHi++ )
				{
					/* set all inward optical depths to taumin, regardless of abundance */
					iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() = opac.taumin;
					iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() = 2.f*opac.taumin;
				}
			}
		}
	}

	/* now do helium like sequence if case b */
	if( opac.lgCaseB )
	{
		for( nelem=1; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				/* La may be case B, tlamin set to 1e9 by default with case b command */
				iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauIn() = opac.tlamin;

				iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauCon() = 
					iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauIn();

				iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauTot() = 
					2.f*iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauIn();

				realnum ratio = opac.tlamin/iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().opacity();

				/* >>chng 02 jan 05, remove explicit list of lyman lines, use As to guess
				 * which are which - this will work for any number of levels */
				for( long i=1; i < iso_sp[ipHE_LIKE][nelem].numLevels_max; i++ )
				{
					if( iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).ipCont() <= 0 )
						continue;

					iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauIn() =
						ratio*iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().opacity();
					/* reset line optical depth to continuum source */
					iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauCon() =
						iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauIn();
					iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauTot() =
						2.f*iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauIn();
					/*fprintf(ioQQQ,
						"RT_tau_init optical depths  Z=%li i=%li  tot=%g in=%g \n",
						nelem, i, iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauTot(),
						iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauIn() );*/
				}
			}
		}
	}

	/* begin sanity check, total greater than 1/0.9 time inward */
	bool lgHit = false;
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			for( ipLo=ipH1s; ipLo < (iso_sp[ipH_LIKE][nelem].numLevels_max - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
				{
					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
						continue;

					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() < 
						0.9f*iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() )
					{
						fprintf(ioQQQ,
							"RT_tau_init insanity for h line, Z=%li lo=%li hi=%li tot=%g in=%g \n",
							nelem , ipLo, ipHi , iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauTot() , 
							iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() );
						lgHit = true;
					}
				}
			}
		}
	}
	if( lgHit )
	{
		fprintf( ioQQQ," stopping due to insanity in RT_tau_init\n");
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}
	/*end sanity check */

	/* fix offset for effective column density optical depth */
	rt.tauxry = opac.TauAbsGeo[0][rt.ipxry-1];

	/* this is flag detected by dest prob routines to see whether ots rates are
	 * oscillating - damp them out if so */
	conv.lgOscilOTS = false;

	/* now that optical depths have been incremented, do escape prob, this
	 * is located here instead on in cloudy.c (where it belongs) because
	 * information generated by RT_line_all is needed for the following printout */
	conv.lgFirstSweepThisZone = true;
	RT_line_all_escape( NULL );

	/* rest of routine is printout in case of trace */
	if( trace.lgTrace )
	{
		/* default is h-like */
		ipISO = ipH_LIKE;
		if( trace.lgIsoTraceFull[ipHE_LIKE] )
			ipISO = ipHE_LIKE;

		if( trace.lgIsoTraceFull[ipISO] )
		{
			t_iso_sp *sp= &iso_sp[ipISO][trace.ipIsoTrace[ipISO]];
			fprintf( ioQQQ, "\n\n   up TauTot array as set in RT_tau_init ipZTrace (nelem)= %ld\n",
				trace.ipIsoTrace[ipH_LIKE] );
			long upper_limit = iso_sp[ipISO][ trace.ipIsoTrace[ipISO] ].numLevels_local;
			for( ipHi=2; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ, " %3ld", ipHi );
				for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
				{
					if( sp->trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					fprintf( ioQQQ,PrintEfmt("%9.2e",
						sp->trans(ipHi,ipLo).Emis().TauTot() ));
				}
				fprintf( ioQQQ, "\n" );
			}

			fprintf( ioQQQ, "\n\n TauIn array\n" );
			for( ipHi=1; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ, " %3ld", ipHi );
				for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
				{
					if( sp->trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					fprintf( ioQQQ,PrintEfmt("%9.2e", 
						sp->trans(ipHi,ipLo).Emis().TauIn() ));
				}
				fprintf( ioQQQ, "\n" );
			}

			fprintf( ioQQQ, "\n\n Aul As array\n" );
			for( ipHi=1; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ, " %3ld", ipHi );
				for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
				{
					if( sp->trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					fprintf( ioQQQ,PrintEfmt("%9.2e", 
						sp->trans(ipHi,ipLo).Emis().Aul()) );
				}
				fprintf( ioQQQ, "\n" );
			}

			fprintf( ioQQQ, "\n\n Aul*esc array\n" );
			for( ipHi=1; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ, " %3ld", ipHi );
				for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
				{
					if( sp->trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					fprintf( ioQQQ,PrintEfmt("%9.2e", 
						sp->trans(ipHi,ipLo).Emis().Aul()*
						sp->trans(ipHi,ipLo).Emis().Ploss()));
				}
				fprintf( ioQQQ, "\n" );
			}

			fprintf( ioQQQ, "\n\n H opac array\n" );
			for( ipHi=1; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ, " %3ld", ipHi );
				for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
				{
					if( sp->trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					fprintf( ioQQQ,PrintEfmt("%9.2e", 
						sp->trans(ipHi,ipLo).Emis().opacity() ));
				}
				fprintf( ioQQQ, "\n" );
			}
		}

		else
		{
			fprintf( ioQQQ, " RT_tau_init called.\n" );
		}
	}

	ASSERT( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauIn()>= 0. );
	ASSERT( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().PopOpc()>= 0. );
	return;
}
