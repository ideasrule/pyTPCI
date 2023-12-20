/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* iso_solve main routine to call iso_level and determine iso level balances */
/* iso_renorm - renormalize iso sequences so that they agree with the ionization balance */
/* AGN_He1_CS routine to save table needed for AGN3 - collision strengths of HeI */
#include "cddefines.h"
#include "atmdat.h"
#include "conv.h"
#include "elementnames.h"
#include "helike_cs.h"
#include "hmi.h"
#include "hydrogenic.h"
#include "ionbal.h"
#include "iso.h"
#include "helike.h"
#include "phycon.h"
#include "rfield.h"
#include "secondaries.h"
#include "thermal.h"
#include "trace.h"
#include "mole.h"
#include "freebound.h"
#include "two_photon.h"
#include "dense.h"

void iso_collapsed_update( void )
{	
	for( long nelem=ipHYDROGEN; nelem<=ipHELIUM; nelem++ )
	{
		if (dense.lgElmtOn[nelem])
		{
			for( long ipISO=ipH_LIKE; ipISO<MIN2(NISO,nelem+1); ipISO++ )
			{
				if ((dense.IonHigh[nelem] >= nelem - ipISO &&
					  dense.IonLow[nelem] <= nelem - ipISO) || !conv.nTotalIoniz)
				{
					iso_cascade( ipISO, nelem );	
				}
			}
		}
	}
}

void iso_update_rates( void )
{
	for( long nelem=ipHYDROGEN; nelem<LIMELM; nelem++ )
	{
		if (!dense.lgElmtOn[nelem])
			continue;
		for( long ipISO=ipH_LIKE; ipISO<MIN2(NISO,nelem+1); ipISO++ )
		{
			if ((dense.IonHigh[nelem] >= nelem - ipISO &&
				  dense.IonLow[nelem] <= nelem - ipISO) || !conv.nTotalIoniz)
			{
				
				/* evaluate collisional rates */
				iso_collide( ipISO,  nelem );
				
				/* truncate atom if physical conditions limit the maximum principal quantum number of a
				 * bound electron to a number less than the allocated size */
				if( iso_ctrl.lgContinuumLoweringEnabled[ipISO] && !conv.nPres2Ioniz )
					iso_continuum_lower( ipISO, nelem );

				/* evaluate recombination rates -- needs to precede iso_photo because of topoff fix */
				iso_radiative_recomb( ipISO , nelem );

				/* evaluate photoionization rates */
				iso_photo( ipISO , nelem );

				/* Generate Gaussian errors if turned on. */
				if( iso_ctrl.lgRandErrGen[ipISO] && nzone==0 && !iso_sp[ipISO][nelem].lgErrGenDone )
				{
					iso_error_generation(ipISO, nelem );
				}

				iso_radiative_recomb_effective( ipISO, nelem );

				iso_ionize_recombine( ipISO , nelem );

				ionbal.RateRecomTot[nelem][nelem-ipISO] = ionbal.RateRecomIso[nelem][ipISO];
			}
			
			/** \todo 2 the indices for the two-photon rates must be changed for further iso sequences. */  
			ASSERT( ipISO <= ipHE_LIKE );
			// two-photon processes
			t_iso_sp* sp = &iso_sp[ipISO][nelem]; 
			for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
			{
				CalcTwoPhotonRates( *tnu, rfield.lgInducProcess && iso_ctrl.lgInd2nu_On );
			}
		}
	}
}

void iso_solve(long ipISO, long nelem, double &maxerr)
{
	DEBUG_ENTRY( "iso_solve()" );

	maxerr = 0.;
	/* do not consider elements that have been turned off */
	if( dense.lgElmtOn[nelem] )
	{
		/* note that nelem scale is totally on c not physical scale, so 0 is h */
		/* evaluate balance if ionization reaches this high */
		if( (dense.IonHigh[nelem] >= nelem - ipISO) &&
			(dense.IonLow[nelem] <= nelem - ipISO) )
		{
			/* solve for the level populations */
			double renorm;
			iso_level( ipISO, nelem, renorm );
			if (fabs(renorm-1.0) > maxerr)
				maxerr = fabs(renorm-1.0);

			/* this just contains a bunch of trace statements. */
			HydroLevel(ipISO, nelem);
		}
		else
		{
			/* zero it out since no population*/
			iso_sp[ipISO][nelem].st[0].Pop() = 0.;
			for( long ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
			{
				iso_sp[ipISO][nelem].st[ipHi].Pop() = 0.;
				for( long ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					/* population of lower level rel to ion, corrected for stim em */
					iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().PopOpc() =  0.;
					iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().mult_opac() = 0.;
				}
			}
		}

		ASSERT( (*iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Lo()).Pop() == iso_sp[ipISO][nelem].st[0].Pop() );

		if( 0 )
		{
			/*Print effective recombination coefficients to 2s level for H */
			if ( nelem == 0 && ipISO == 0 )
			{
				double alphaeff2s=0., pop, alphaeff2p = 0;
				long numLevs = iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max;
				fprintf( ioQQQ,"Effective recombination 2S/2P, ipISO=%li, nelem=%li, Te = %e, dens=%e\t%e, numlevels=%li\n", ipISO, nelem, phycon.te,
						dense.eden, dense.xIonDense[nelem][ipISO+1],numLevs );
				fprintf( ioQQQ, "N\tL\tS\tAlphaEffec\n" );

				for( long ipHi=3; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					//if (iso_sp[ipISO][nelem].st[ipHi].l()==1)
					pop = iso_sp[ipISO][nelem].st[ipHi].Pop();
					if (ipHi > numLevs )
						pop /=(2.*pow2((double)iso_sp[ipISO][nelem].st[ipHi].n()));

					alphaeff2s +=  pop*iso_sp[ipISO][nelem].trans(ipHi,1).Emis().Aul();

					alphaeff2p +=  pop*iso_sp[ipISO][nelem].trans(ipHi,2).Emis().Aul();

				}

				alphaeff2s /= (dense.eden*dense.xIonDense[nelem][ipISO+1]);
				alphaeff2s += iso_sp[ipISO][nelem].fb[1].RadRecomb[ipRecRad];

                                alphaeff2p /= (dense.eden*dense.xIonDense[nelem][ipISO+1]);
                                alphaeff2p += iso_sp[ipISO][nelem].fb[2].RadRecomb[ipRecRad];


				fprintf( ioQQQ, "%li\t%li\t%li\t%e\t\n", iso_sp[ipISO][nelem].st[1].n(),iso_sp[ipISO][nelem].st[1].l() ,
						iso_sp[ipISO][nelem].st[1].S(),alphaeff2s );
				fprintf( ioQQQ, "%li\t%li\t%li\t%e\t\n", iso_sp[ipISO][nelem].st[2].n(),iso_sp[ipISO][nelem].st[2].l() ,
                                                iso_sp[ipISO][nelem].st[2].S(),alphaeff2p );
				fprintf( ioQQQ, "\n" );
			}
		}


	}

	return;
}

void IonHydro( void )
{
	DEBUG_ENTRY( "IonHydro()" );

	/* ============================================================================== */
	/* rest is for hydrogen only */

	{
		/*@-redef@*/
		/* often the H- route is the most efficient formation mechanism for H2,
		 * will be through rate called ratach
		 * this debug print statement is to trace h2 oscillations */
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if(DEBUG_LOC )
		{
			fprintf(ioQQQ,"DEBUG  \t%.2f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
				fnzone,
				hmi.H2_total ,
				findspecieslocal("H2")->den,
				dense.xIonDense[ipHYDROGEN][0],
				dense.xIonDense[ipHYDROGEN][1],
				hmi.H2_Solomon_dissoc_rate_used_H2g,
				hmi.H2_Solomon_dissoc_rate_BD96_H2g,
				hmi.H2_Solomon_dissoc_rate_TH85_H2g);
		}
	}
#if 0
	/* >>chng 01 may 09, add option to force abundance, with element name ioniz cmmnd */
	if( dense.lgSetIoniz[ipHYDROGEN] )
	{
		realnum dense_ation = dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHYDROGEN][1];
		dense.xIonDense[ipHYDROGEN][1] = dense.SetIoniz[ipHYDROGEN][1]*dense_ation;
		dense.xIonDense[ipHYDROGEN][0] = dense.SetIoniz[ipHYDROGEN][0]*dense_ation;

		/* initialize ground state pop too */
		iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop() = dense.xIonDense[ipHYDROGEN][0];
	}
	else
	{
		/* 
		 * >> chng 03 jan 15 rjrw:- terms are now in ion_solver, to allow for
		 * molecular sources and sinks of H and H+.  ion_solver renormalizes
		 * to keep the total H abundance correct -- only the molecular
		 * network is allowed to change this. 
		 */
		ion_solver( ipHYDROGEN , false );
	}

	fixit(); /* this is called in HydroLevel above, is it needed in both places? */
	/* >>hcng 05 mar 24,
	 * renormalize the populations and emission of H atom to agree with chemistry */
	double renorm;
	iso_renorm( ipHYDROGEN, ipH_LIKE, renorm );
#endif		
	ion_solver( ipHYDROGEN , false );

	/* remember the ratio of pops of 2p to 1s for possible printout in prtComment
	 * and to obtain Lya excitation temps.  the pop of ground is not defined if
	 * NO ionization at all since these pops are relative to ion */
	/* >>chng 99 jun 03, added MAX2 to protect against totally neutral gas */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()/MAX2(SMALLDOUBLE,iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()) > 0.1 &&
		iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop() > SMALLDOUBLE )
	{
		hydro.lgHiPop2 = true;
		hydro.pop2mx = (realnum)MAX2(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()/iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop(),
		  hydro.pop2mx);
	}

	double gamtot = iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].gamnc + secondaries.csupra[ipHYDROGEN][0];

	double coltot = iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ColIoniz + 
		iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Coll().ColUL( colliders ) / dense.EdenHCorr *
		4. * iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Boltzmann();

	/* if ground state destruction rate is significant, recall different dest procceses */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].RateLevel2Cont > SMALLFLOAT )
	{
		hydro.H_ion_frac_photo = 
			(realnum)(iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].gamnc/iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].RateLevel2Cont );

		/* fraction of ionizations of H from ground, due to thermal collisions */
		hydro.H_ion_frac_collis = 
			(realnum)(iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ColIoniz*dense.eden/iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].RateLevel2Cont);

		/* this flag is used in ConvBase to decide whether we
		 * really need to converge the secondary ionization rates */
		secondaries.sec2total = 
			(realnum)(secondaries.csupra[ipHYDROGEN][0] / iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].RateLevel2Cont);

		/* frac of ionizations due to ct */
		atmdat.HIonFrac = atmdat.CharExcIonTotal[ipHYDROGEN] / iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].RateLevel2Cont;
	}
	else
	{
		hydro.H_ion_frac_collis = 0.;
		hydro.H_ion_frac_photo = 0.;
		secondaries.sec2total = 0.;
		atmdat.HIonFrac = 0.;
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "       Hydrogenic return %.2f ",fnzone);
		fprintf(ioQQQ,"H0:%.3e ", dense.xIonDense[ipHYDROGEN][0]);
		fprintf(ioQQQ,"H+:%.3e ", dense.xIonDense[ipHYDROGEN][1]);
		fprintf(ioQQQ,"H2:%.3e ", hmi.H2_total);
		fprintf(ioQQQ,"H-:%.3e ", findspecieslocal("H-")->den);
		fprintf(ioQQQ,"ne:%.3e ", dense.eden);
		fprintf( ioQQQ, " REC, COL, GAMT= ");
		/* recomb rate coef, cm^3 s-1 */
		fprintf(ioQQQ,"%.2e ", iso_sp[ipH_LIKE][ipHYDROGEN].RadRec_effec );
		fprintf(ioQQQ,"%.2e ", coltot);
		fprintf(ioQQQ,"%.2e ", gamtot);
		fprintf( ioQQQ, " CSUP=");
		PrintE82( ioQQQ, secondaries.csupra[ipHYDROGEN][0]);
		fprintf( ioQQQ, "\n");
	}

	return;
}

/* iso_renorm - renormalize iso sequences so that they agree with the ionization balance */
void iso_renorm( long nelem, long ipISO, double &renorm )
{
	DEBUG_ENTRY( "iso_renorm()" );

	const bool lgJustAssert = false, lgJustTest = true;
	double sum_atom_iso;

	renorm = 1.0;
	
	if (ipISO > nelem)
		return;

	/*>>chng 04 mar 23, add this renorm */
	/* renormalize the state specific populations, so that they
	 * add up to the results that came from ion_solver 
	 * units at first is sum div by H+ density, since Pop2Ion defn this way */
	sum_atom_iso = 0.;
	/* >> chng 06 aug 31, from numLevels_max to _local */
	for( long level=0; level < iso_sp[ipISO][nelem].numLevels_local; level++ )
	{
		/* cm-3 - total population in iso solved model */
		sum_atom_iso += iso_sp[ipISO][nelem].st[level].Pop();
	}

	// If total iso population is tiny, populate ground state (probably due to
	// e.g. ++dense.IonHigh[nelem] somewhere)
	if (sum_atom_iso <= SMALLFLOAT)
	{
		// Ensure this is different from the ion density, so it is signalled as non-convergence
		if (dense.xIonDense[nelem][nelem-ipISO] > 2.0*SMALLFLOAT)
			sum_atom_iso = 0.5*dense.xIonDense[nelem][nelem-ipISO];
		else
			sum_atom_iso = 1.0;
		if ( !lgJustTest )
			iso_sp[ipISO][nelem].st[0].Pop() = sum_atom_iso;
	}

	/* >>chng 04 may 25, humunculus sum_atom_iso is zero */
	renorm = dense.xIonDense[nelem][nelem-ipISO] / sum_atom_iso;

	if (lgJustTest)
		return;

	if (0)
		fprintf (ioQQQ, "Iso_Renorm %ld %ld %g %g %g[%g]\n",
					ipISO,nelem,renorm-1.,
					dense.xIonDense[nelem][nelem-ipISO], sum_atom_iso,
					iso_sp[ipISO][nelem].st[0].Pop());
	if (lgJustAssert) 
	{
		ASSERT(fp_equal(renorm,1.0));
		return;
	}
	if (conv.lgConvIoniz() && dense.xIonDense[nelem][nelem-ipISO] > SMALLFLOAT && 
		 fabs(renorm - 1.0) > conv.IonizErrorAllowed)
	{
		conv.setConvIonizFail( "Iso vs. ion", 
									  dense.xIonDense[nelem][nelem-ipISO], 
									  sum_atom_iso);	
		//fprintf(ioQQQ,"%g %g %g\n",renorm-1.0,sum_atom_iso,dense.xIonDense[nelem][nelem-ipISO]);
	}

	/*fprintf(ioQQQ,"DEBUG renorm\t%.2f\t%.3e\n",fnzone, renorm);*/

	//ASSERT( renorm < BIGFLOAT );
	/* renormalize populations from iso-model atom so that they agree with ion solver & chemistry */
	/*fprintf(ioQQQ,"DEBUG h \t%.3e hydrogenic renorm %.3e\n", 
		iso_sp[ipH_LIKE][nelem].st[ipH1s].Pop ,
		iso_sp[ipH_LIKE][nelem].st[ipH1s].Pop/renorm );*/

	for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
	{
		iso_sp[ipISO][nelem].st[ipHi].Pop() *= renorm;
	}

	/* >> chng 06 aug 31, from numLevels_max to _local */
	for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
	{
		for( long ipLo=0; ipLo < ipHi; ipLo++ )
		{
			if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
				continue;

			/* population of lower level rel to ion, corrected for stim em */
			iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().PopOpc() *= renorm;
		}
	}

	return;
}

void iso_departure_coefficients( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_departure_coefficients()" );
		
	for( long level=0; level < iso_sp[ipISO][nelem].numLevels_local; level++ )
	{
		double denom = dense.xIonDense[nelem][nelem+1-ipISO]*
			iso_sp[ipISO][nelem].fb[level].PopLTE*dense.eden;

		if( iso_sp[ipISO][nelem].fb[level].PopLTE > 0. && denom > SMALLFLOAT )
			iso_sp[ipISO][nelem].st[level].DepartCoef() = safe_div( 
			  iso_sp[ipISO][nelem].st[level].Pop(), denom );
		else
			iso_sp[ipISO][nelem].st[level].DepartCoef() = 0.;
	}

	// These are levels above the lowered continuum. They are in equilibrium with the electron bath.
	for( long level=iso_sp[ipISO][nelem].numLevels_local; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
		iso_sp[ipISO][nelem].st[level].DepartCoef() = 1.;

	return;
}

/*iso_prt_pops print out iso sequence populations or departure coefficients */
void iso_prt_pops( long ipISO, long nelem, bool lgPrtDeparCoef )
{
	long int in, il, is, i, ipLo, nResolved, ipFirstCollapsed=LONG_MIN;
	char chPrtType[2][12]={"populations","departure"};
	/* first dimension is multiplicity */
	char chSpin[3][9]= {"singlets", "doublets", "triplets"};

#define ITEM_TO_PRINT(A_)	( lgPrtDeparCoef ? iso_sp[ipISO][nelem].st[A_].DepartCoef() : iso_sp[ipISO][nelem].st[A_].Pop() )

	DEBUG_ENTRY( "iso_prt_pops()" );

	ASSERT( ipISO < NISO );

	for( is = 1; is <= 3; ++is )
	{
		if( ipISO == ipH_LIKE && is != 2 )
			continue;
		else if( ipISO == ipHE_LIKE && is != 1 && is != 3 )
			continue;

		ipFirstCollapsed= iso_sp[ipISO][nelem].numLevels_local-iso_sp[ipISO][nelem].nCollapsed_local;
		nResolved = iso_sp[ipISO][nelem].st[ipFirstCollapsed-1].n();
		ASSERT( nResolved == iso_sp[ipISO][nelem].n_HighestResolved_local );
		ASSERT(nResolved > 0 );

		/* give element number and spin */
		fprintf(ioQQQ," %s %s  %s %s\n",
			iso_ctrl.chISO[ipISO],
			elementnames.chElementSym[nelem],
			chSpin[is-1],
			chPrtType[lgPrtDeparCoef]);

		/* header with the l states */
		fprintf(ioQQQ," n\\l=>    ");
		for( i =0; i < nResolved; ++i)
		{
			fprintf(ioQQQ,"%2ld         ",i);
		}
		fprintf(ioQQQ,"\n");

		/* loop over print quantum numbers, one per line, with l across */
		for( in = 1; in <= nResolved; ++in)
		{
			if( is==3 && in==1 )
				continue;

			fprintf(ioQQQ," %2ld      ",in);
			double TotalNumber=0.;
			for( il = 0; il < in; ++il)
			{
				if( ipISO==ipHE_LIKE && (in==2) && (il==1) && (is==3) )
				{
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipHe2p3P0) );
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipHe2p3P1) );
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipHe2p3P2) );
				}
				else
				{
					ipLo = iso_sp[ipISO][nelem].QN2Index(in, il, is);
					fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(ipLo) );
					TotalNumber += (2*il+1)*ITEM_TO_PRINT(ipLo);
				}
			}
			TotalNumber/=pow2((double)in);
			fprintf(ioQQQ,"\n Total = %9.3e",TotalNumber);
			fprintf(ioQQQ,"\n");
		}
	}
	/* above loop was over spin, now do collapsed levels, no spin or ang momen */
	for( il = ipFirstCollapsed; il < iso_sp[ipISO][nelem].numLevels_local; ++il)
	{
		in = iso_sp[ipISO][nelem].st[il].n();
		/* prin quan number of collapsed levels */
		fprintf(ioQQQ," %2ld      ",in);
		fprintf( ioQQQ, "%9.3e ", ITEM_TO_PRINT(il) );
		fprintf(ioQQQ,"\n");
	}
	return;
}

/* routine to save table needed for AGN3 - collision strengths of HeI */
void AGN_He1_CS( FILE *ioPun )
{

	long int i;

	/* list of temperatures where cs will be printed */
	const int NTE = 5;
	double TeList[NTE] = {6000.,10000.,15000.,20000.,25000.};
	double TempSave;

	DEBUG_ENTRY( "AGN_He1_CS()" );

	/* put on a header */
	fprintf(ioPun, "Te\t2 3s 33s\n");

	/* Restore the original temp when this routine done.	*/
	TempSave = phycon.te;

	for( i=0; i<NTE; ++i )
	{
		TempChange(TeList[i] , false);

		fprintf(ioPun , "%.0f\t", 
			TeList[i] );
		fprintf(ioPun , "%.2f\t", 
			HeCSInterp( 1 , ipHe3s3S , ipHe2s3S, ipELECTRON ) );
		fprintf(ioPun , "%.2f\t", 
			HeCSInterp( 1 , ipHe3p3P , ipHe2s3S, ipELECTRON ) );
		fprintf(ioPun , "%.2f\t", 
			HeCSInterp( 1 , ipHe3d3D , ipHe2s3S, ipELECTRON ) );
		fprintf(ioPun , "%.3f\t", 
			HeCSInterp( 1 , ipHe3d1D , ipHe2s3S, ipELECTRON ) );
		/*fprintf(ioPun , "%.1f\t%.1f\t%.1f\n", */
		fprintf(ioPun , "%.1f\n", 
			HeCSInterp( 1 , ipHe2p3P0 , ipHe2s3S, ipELECTRON ) +
			HeCSInterp( 1 , ipHe2p3P1 , ipHe2s3S, ipELECTRON ) +
			HeCSInterp( 1 , ipHe2p3P2 , ipHe2s3S, ipELECTRON ));
	}

	/* no need to force update since didn't do above	*/
	TempChange(TempSave , false);
	return;
}
