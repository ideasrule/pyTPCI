/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_solver solve the bi-diagonal matrix for ionization balance */
#include "cddefines.h"
#include "yield.h"
#include "prt.h"
#include "continuum.h"
#include "iso.h"
#include "dynamics.h"
#include "grainvar.h"
#include "newton_step.h"
#include "conv.h"
#include "secondaries.h"
#include "phycon.h"
#include "atmdat.h"
#include "heavy.h"
#include "elementnames.h"
#include "ionbal.h"
#include "trace.h"
#include "ion_trim.h"
#include "mole.h"
#include "freebound.h"
#include "dense.h"
#include "thirdparty.h"

STATIC bool lgTrivialSolution( long nelem, double abund_total );

STATIC void find_solution( long nelem, long ion_range, valarray<double> &xmat, valarray<double> &source);

STATIC void fill_array( long int nelem, long ion_range, valarray<double> &xmat, valarray<double> &auger );

STATIC void fill_ext_src_and_snk( long nelem, long ion_low, long ion_range, valarray<double> &xmat, valarray<double> &source );

STATIC double get_total_abundance_ions( long int nelem );

STATIC void HomogeneousSource( long nelem, long ion_low, long ion_range, valarray<double> &xmat, valarray<double> &source, double abund_total );

STATIC void renorm_solution( long nelem, long ion_range, long ion_low, double *source, double abund_total, bool *lgNegPop );

STATIC void store_new_densities( long nelem, long ion_range, long ion_low, double *source );

STATIC void clean_up( long nelem, double abund_total );

STATIC void PrintRates( long nelem, bool lgNegPop, double abund_total, valarray<double> &auger, bool lgPrintIt );

void solveions(double *ion, double *rec, double *snk, double *src,
			   long int nlev, long int nmax);

	/* this will be used to address the 2d arrays */
#	undef MAT
#	define MAT(M_,I_,J_)	((M_)[(I_)*(ion_range)+(J_)])

#	undef MAT1
#	define MAT1(M_,I_,J_)	((M_)[(I_)*(ion_range1)+(J_)])

#	undef MAT2
#	define MAT2(M_,I_,J_)	((M_)[(I_)*(ion_range2)+(J_)])

void ion_solver( long int nelem, bool lgPrintIt) 
{
	double abund_total = 0.0;
	long ion_low = dense.IonLow[nelem];
	bool lgNegPop = false;
	double error = 0.0;

	DEBUG_ENTRY( "ion_solver()" );

	iso_charge_transfer_update(nelem);
	
	long ion_range = dense.IonHigh[nelem]-dense.IonLow[nelem]+1;
	valarray<double> xmat(ion_range*ion_range);
	valarray<double> source(ion_range);
	valarray<double> auger(LIMELM+1);

	// V/W cycle would do ion solution *first*, then iso, then
	// ion again if iso made any changes
	// Seem to need to break out after iso at the moment.   Why?

	static ConvergenceCounter cctr=conv.register_("ION_SOLVES");
	++cctr;
	for (long it=0; it<4; ++it)
	{
		static ConvergenceCounter cctrl=conv.register_("ISO_LOOPS");
		++cctrl;

		for( long ipISO=ipH_LIKE; ipISO<MIN2(NISO,nelem+1); ipISO++ )
		{
			if( (dense.IonHigh[nelem] >= nelem - ipISO) &&
				  (dense.IonLow[nelem] <= nelem - ipISO)) 
			{
				iso_set_ion_rates(ipISO, nelem);
			}
		}
		
		//if ( it > 1 && error > 0.0)
		//	fprintf(ioQQQ,"ion nelem %ld loop %ld error %g nzone %ld\n",nelem,it,error,nzone);
				
		abund_total = get_total_abundance_ions( nelem );


		{
			/* this sets up a fake ionization balance problem, with a trivial solution,
			 * for debugging the ionization solver */
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC && nelem==ipCARBON )
			{
				broken();/* within optional debug print statement */
				dense.IonLow[nelem] = 0;
				dense.IonHigh[nelem] = 3;
				abund_total = 1.;
				long limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);
				/* make up ionization and recombination rates */
				for( long ion=dense.IonLow[nelem]; ion <= limit; ion++ )
				{
					double fac=1;
					if(ion)
						fac = 1e-10;
					ionbal.RateRecomTot[nelem][ion] = 100.;
					for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
					{
						/* direct photoionization of this shell */
						ionbal.PhotoRate_Shell[nelem][ion][ns][0] = fac;
					}
				}
			}
		}

		bool lgTrivial = lgTrivialSolution( nelem, abund_total );
		if( !lgTrivial )
		{
  			// fill xmat and source with appropriate terms
			fill_array( nelem, ion_range, xmat, auger );

			// Fill source and sink terms from chemistry.
			fill_ext_src_and_snk( nelem, ion_low, ion_range, xmat, source );
			
			// decide if matrix is homogeneous
			HomogeneousSource( nelem, ion_low, ion_range, xmat, source, abund_total );
			
			// Now find the solution
			find_solution( nelem, ion_range, xmat, source);
			
			renorm_solution( nelem, ion_range, ion_low, &source[0], abund_total, &lgNegPop );

			// save the results in the global density variables
			store_new_densities( nelem, ion_range, ion_low, &source[0] );

			clean_up( nelem, abund_total );

			ASSERT(ion_range >= dense.IonHigh[nelem]-dense.IonLow[nelem]+1);
			if (ion_range != dense.IonHigh[nelem]-dense.IonLow[nelem]+1)
			{
				ion_range = dense.IonHigh[nelem]-dense.IonLow[nelem]+1;
				xmat.resize(ion_range*ion_range);
				source.resize(ion_range);
			}			
		}

		error = 0.0;
		/* update and solve iso levels */
		for( long ipISO=ipH_LIKE; ipISO<MIN2(NISO,nelem+1); ipISO++ )
		{
			double thiserror;
			iso_solve( ipISO, nelem, thiserror );
			if (thiserror > error)
				error = thiserror;
		}

		if ( error < 1e-4 )
			break;
	}

	iso_satellite_update(nelem);

	for( long ipISO=ipH_LIKE; ipISO<MIN2(NISO,nelem+1); ipISO++ )
	{
		/* now evaluate departure coefficients */
		iso_departure_coefficients( ipISO, nelem );
	}

	if( prt.lgPrtArry[nelem] || lgPrintIt )
		PrintRates( nelem, lgNegPop, abund_total, auger, lgPrintIt );

	return;
}

STATIC bool lgTrivialSolution( long nelem, double abund_total )
{	
	double renorm;
	/* return if IonHigh is zero, since no ionization at all */
	if( dense.IonHigh[nelem] == dense.IonLow[nelem] )
	{
		/* set the atom to the total gas phase abundance */
		dense.xIonDense[nelem][dense.IonHigh[nelem]] = abund_total;
		for ( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
			iso_renorm(nelem, ipISO, renorm);
		return true;
	}
	else if( dense.lgSetIoniz[nelem] )
	{
		/* option to force ionization distribution with element name ioniz */
		for( long ion=0; ion<nelem+2; ++ion )
			dense.xIonDense[nelem][ion] = dense.SetIoniz[nelem][ion]*abund_total;
		for ( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
			iso_renorm(nelem, ipISO, renorm);
		return true;
	}
	else
		return false;
}

STATIC void find_solution( long nelem, long ion_range, valarray<double> &xmat, valarray<double> &source)
{
	int32 nerror;
	valarray<double> xmatsave(ion_range*ion_range);
	valarray<double> sourcesave(ion_range);
	valarray<int32> ipiv(ion_range);

	DEBUG_ENTRY( "find_solution()" );

	// save copy of xmat before continuing.
	for( long i=0; i< ion_range; ++i )
	{
		sourcesave[i] = source[i];
		for( long j=0; j< ion_range; ++j )
		{
			MAT( xmatsave, i, j ) = MAT( xmat, i, j );
		}
	}

	if (1)
	{
		nerror = solve_system(xmat,source,ion_range,NULL);
		
		if (nerror)
		{
			// solve_system doesn't return a valid solution, so just
			// provide initial values to keep things running
			fprintf(ioQQQ,"Error for nelem %ld, active ion range %ld--%ld\n",
					  nelem,dense.IonLow[nelem],dense.IonHigh[nelem]);
			fprintf(ioQQQ,"Initial ion abundances:");
			for( long j=0; j<nelem+2; ++j )
				fprintf(ioQQQ," %g",dense.xIonDense[nelem][j]);
			fprintf(ioQQQ,"\n");
			for( long j=0; j<ion_range; ++j )
				source[j] = dense.xIonDense[nelem][dense.IonLow[nelem]+j];
		}
	}
	else
	{
		/* this is the default solver - now get new solution */
		nerror = 0;
		/* Use general matrix solver */
		getrf_wrapper(ion_range, ion_range, &xmat[0], ion_range, &ipiv[0], &nerror);
		if( nerror != 0 )
		{
			fprintf( ioQQQ, 
					" DISASTER ion_solver: dgetrf finds singular or ill-conditioned matrix nelem=%li %s ion_range=%li, nConv %li IonLow %li IonHi %li\n",
					nelem , 
					elementnames.chElementSym[nelem],
					ion_range,
					conv.nTotalIoniz ,
                                dense.IonLow[nelem], dense.IonHigh[nelem]);
			fprintf( ioQQQ, " xmat follows\n");
			for( long i=0; i<ion_range; ++i )
			{
				for( long j=0;j<ion_range;j++ )
				{
					fprintf(ioQQQ,"%e\t",MAT(xmatsave,j,i));
				}
				fprintf(ioQQQ,"\n");
			}
			fprintf(ioQQQ,"source follows\n");
			for( long i=0; i<ion_range;i++ )
			{
				fprintf(ioQQQ,"%e\t",sourcesave[i]);
			}
			fprintf(ioQQQ,"\n");
			cdEXIT(EXIT_FAILURE);
		}
		getrs_wrapper('N', ion_range, 1, &xmat[0], ion_range, &ipiv[0], &source[0], ion_range, &nerror);
		if( nerror != 0 )
		{
			fprintf( ioQQQ, " DISASTER ion_solver: dgetrs finds singular or ill-conditioned matrix nelem=%li ionrange=%li\n",
					nelem , ion_range );
			cdEXIT(EXIT_FAILURE);
		}
	}
	
	{
		/* this is to debug following failed assert */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nelem == ipHYDROGEN )
		{
			fprintf(ioQQQ,"debuggg ion_solver1 %ld\t%.2f\t%.4e\t%.4e\tIon\t%.3e\tRec\t%.3e\n", 
					nelem,
					fnzone,
					phycon.te,
					dense.eden,
					ionbal.RateIonizTot(nelem,0) , 
					ionbal.RateRecomTot[nelem][0]);
			fprintf(ioQQQ," Msrc %.3e %.3e\n", mole.source[nelem][0], mole.source[nelem][1]);
			fprintf(ioQQQ," Msnk %.3e %.3e\n", mole.sink[nelem][0], mole.sink[nelem][1]);
			fprintf(ioQQQ," Poprat %.3e nomol %.3e\n",source[1]/source[0],
					ionbal.RateIonizTot(nelem,0)/ionbal.RateRecomTot[nelem][0]);
		}
	}
	
	for( long i=0; i<ion_range; i++ )
	{
		ASSERT( !isnan( source[i] ) );
		ASSERT( source[i] < MAX_DENSITY );
	}
	
	return;
}

STATIC void renorm_solution( long nelem, long ion_range, long ion_low, double *source, double abund_total, bool *lgNegPop ) 
{
	DEBUG_ENTRY( "renorm_solution()" );

	ASSERT( nelem >= 0 );
	ASSERT( nelem < LIMELM );
	ASSERT( ion_range <= nelem + 2 );
	ASSERT( ion_low >= 0 );
	ASSERT( ion_low <= nelem + 1 );

#if 0
#	define RJRW 0
	if( RJRW && 0 )
	{
		/* verify that the rates are sensible */
		double test;
		for(long i=0; i<ion_range; i++) {
			test = 0.;
			for(long j=0; j<ion_range; j++) {
				test = test+source[j]*MAT(xmatsave,j,i);
			}
			fprintf(ioQQQ,"%e\t",test);
		}
		fprintf(ioQQQ,"\n");

		test = 0.;
		fprintf(ioQQQ," ion %li abundance %.3e\n",nelem,abund_total);
		for( long ion=dense.IonLow[nelem]; ion < dense.IonHigh[nelem]; ion++ )
		{
			if( ionbal.RateRecomTot[nelem][ion] != 0 && source[ion-ion_low] != 0 )
				fprintf(ioQQQ," %li %.3e %.3e : %.3e\n",ion,source[ion-ion_low],
								source[ion-ion_low+1]/source[ion-ion_low],
								ionbal.RateIonizTot(nelem,ion)/ionbal.RateRecomTot[nelem][ion]);
			else
				fprintf(ioQQQ," %li %.3e [One ratio infinity]\n",ion,source[ion-ion_low]);
			test += source[ion-ion_low];
		}
	}
#endif

	/* 
	 * >> chng 03 jan 15 rjrw:- terms are now included for
	 * molecular sources and sinks.
	 *
	 * When the network is not in equilibrium, this will lead to a
	 * change in the derived abundance of coupled ions after the matrix
	 * solution.
	 *
	 * We therefore renormalize to keep the total abundance of the
	 * states treated by the ionization ladder constant -- only the
	 * molecular network is allowed to change this.
	 *
	 * The difference between `renorm' and 1. is a measure of the
	 * quality of the solution (it will be 1. if the rate of transfer
	 * into the ionization ladder species balances the rate of transfer
	 * out, for the consistent relative abundances)
	 *
	 */

	/* source[i] contains new solution for ionization populations
	 * save resulting abundances into main ionization density array,
	 * while checking whether any negative level populations occurred */
	*lgNegPop = false;
	for( long i=0; i < ion_range; i++ )
	{
		long ion = i+ion_low;

		if( source[i] < 0. )
		{
			/* >>chng 04 dec 04, put test on neg abund here, don't print unless value is very -ve */
			/* >>chng 06 feb 28, from -1e-10 to -1e-9, sim func_t10 had several negative
			 * during initial search, due to extremely high ionization */
			/* >>chng 06 mar 11, from 1e-9 to 2e-9 make many struc elements floats from double */
			if( source[i]<-2e-9 )
			{
				fprintf(ioQQQ,
				" PROBLEM negative ion population in ion_solver, nelem=%li, %s ion=%li val=%.3e Search?%c zone=%li iteration=%li\n",
				nelem ,
				elementnames.chElementSym[nelem],
				ion ,
				source[i] ,
				TorF(conv.lgSearch) ,
				nzone ,
				iteration );
				*lgNegPop = true;
				fixit("break PrintRates into one NegPop case and one trace? No auger defined here.");
				//PrintRates( nelem, *lgNegPop, abund_total, auger );
			}

			fixit("Kill this bit and force exit on negative populations.");
#if 1
			source[i] = 0.;
			/* if this is one of the iso seq model atoms then must also zero out pops */
			//if( ion == nelem+1-NISO ) //newmole had this should have been next line?
			if( ion > nelem-NISO && ion < nelem + 1 )
			{
				long int ipISO = nelem - ion;
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
					iso_sp[ipISO][nelem].st[level].Pop() = 0.;
			}
#endif
		}
	}

	double renormnew = 1.;
	{
		double dennew = 0.;

		/* find total population to renorm - also here check that negative pops did not occur */
		for( long i=0;i < ion_range; i++ )
		{
			dennew += source[i];
		}

		if(dennew > 0.)
		{
			renormnew = abund_total / dennew;
			/** \todo	2	renorm should == 1 when the molecules and
			 * ionization are in equilibrium.  Should monitor
			 * this figure of merit in calling routine.
			 * */
		}
		else
		{
			renormnew = 1.;
		}

		for( long i=0;i < ion_range; i++ )
		{
			source[i] *= renormnew;
			if( source[i] >= MAX_DENSITY )
			{
				long ion = i+ion_low;
				fprintf( ioQQQ, "PROBLEM DISASTER Huge density in ion_solver, nelem %ld ion %ld source %e renormnew %e\n",
					nelem, ion, source[i], renormnew );
			}
		}
	}
	/* check not negative, should be +ve, can be zero if species has become totally molecular.
	 * this happens for hydrogen if no cosmic rays, or cr ion set very low */
	if( renormnew < 0)
	{
		fprintf(ioQQQ,"PROBLEM impossible value of renorm \n");
	}
	ASSERT( renormnew>=0 );

	return;
}

STATIC void store_new_densities( long nelem, long ion_range, long ion_low, double *source )
{
	DEBUG_ENTRY( "store_new_densities()" );

	for( long i=0; i < ion_range; i++ )
	{
		long ion = i+ion_low;
		dense.xIonDense[nelem][ion] = source[i];
		ASSERT( dense.xIonDense[nelem][ion] < MAX_DENSITY );
	}

	return;
}

STATIC void clean_up( long nelem, double abund_total )
{
	DEBUG_ENTRY( "clean_up()" );

	fixit("this should only be done if trimming is not disabled?");

	/* Zero levels with abundances < 1e-25 which which will suffer numerical noise */
	ion_trim_small(nelem, abund_total);

	for ( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		double renorm;
		iso_renorm(nelem, ipISO, renorm);
	}

	// sanity check, either offset stages of low and high ionization
	ASSERT( (dense.IonLow[nelem] <= dense.IonHigh[nelem]) ||
		// both totally neutral
		(dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 ) ||
		// both fully stripped
		( dense.IonLow[nelem]==nelem+1 && dense.IonHigh[nelem]==nelem+1 ) );

	return;
}

STATIC double get_total_abundance_ions( long int nelem )
{
	DEBUG_ENTRY( "get_total_abundance_ions()" );

	ASSERT( nelem >= 0 );
	ASSERT( nelem < LIMELM );

	ionbal.elecsnk[nelem] = 0.;
	ionbal.elecsrc[nelem] = 0.;
	
	double abund_total = 0.;
	for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
	{
		abund_total += dense.xIonDense[nelem][ion];
	}

	realnum tot1 = dense.gas_phase[nelem];
	realnum tot2 = (realnum)(dense.xMolecules(nelem)+abund_total);
	if (0)
	{
		if (fabs(tot2-tot1) > conv.GasPhaseAbundErrorAllowed*tot1)
			fprintf(ioQQQ,"%ld %13.8g %13.8g %13.8g %13.8g\n",nelem,tot1,tot2,abund_total,tot2/tot1 - 1.0);
	}

	ASSERT( fp_equal_tol(tot1, tot2, realnum(conv.GasPhaseAbundErrorAllowed*tot1 + 100.f*FLT_MIN)) );

	ASSERT( abund_total < MAX_DENSITY );

	return abund_total;
}

STATIC void fill_array( long int nelem, long ion_range, valarray<double> &xmat, valarray<double> &auger )
{
	long int limit, 
	  IonProduced;
	double rateone;
	long ion_low;

	valarray<double> sink(ion_range);
	valarray<int32> ipiv(ion_range);

	DEBUG_ENTRY( "fill_array()" );

	/* this is on the c scale, so H is 0 */
	ASSERT( nelem >= 0);
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonHigh[nelem] >= 0 );

	/* impossible for HIonFrac[nelem] to be zero if IonHigh(nelem)=nelem+1
	 * HIonFrac(nelem) is stripped to hydrogen */
	/* >>chng 01 oct 30, to assert */
	//ASSERT( (dense.IonHigh[nelem] < nelem + 1) || dense.xIonDense[nelem][nelem+1-ipH_LIKE] > 0. );

	/* zero out the ionization and recombination rates that we will modify here,
	 * but not the iso-electronic stages which are done elsewhere,
	 * the nelem stage of ionization is he-like,
	 * the nelem+1 stage of ionization is h-like */

	/* loop over stages of ionization that we solve for here, 
	 * up through and including one less than nelem-NISO,
	 * never actually do highest NISO stages of ionization since they
	 * come from the ionization ratio from the next lower stage */
	limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);

	/* the full range of ionization - this is number of ionization stages */
	ASSERT( ion_range <= nelem+2 );

	ion_low = dense.IonLow[nelem];

	/* zero-out loop comes before main loop since there are off-diagonal
	 * elements in the main ionization loop, due to multi-electron processes,
	 * TotIonizRate and TotRecom were already set in h-like and he-like solvers 
	 * other recombination rates were already set by routines responsible for them */
	for( long ion_from=0; ion_from <= limit; ion_from++ )
	{
		for( long ion_to=0; ion_to < nelem+2; ion_to++ )
		{
			ionbal.RateIoniz[nelem][ion_from][ion_to] = 0.;
		}
	}

	/* auger is used only for debug printout - it is special because with multi-electron
	 * Auger ejection, very high stages of ionization can be produced, well beyond
	 * the highest stage considered here.  so we allocate to the limit */
	for( long i=0; i< LIMELM+1; ++i )
	{
		auger[i] = 0.;
	}

	/* zero out xmat */
	for( long i=0; i< ion_range; ++i )
	{
		for( long j=0; j< ion_range; ++j )
		{
			MAT( xmat, i, j ) = 0.;
		}
	}

	bool lgGrainsOn = gv.lgDustOn() && ionbal.lgGrainIonRecom && gv.lgGrainPhysicsOn;

	/* Now put in all recombination and ionization terms from CO_mole() that 
	 * come from molecular reactions. this traces molecular process that 
	 * change ionization stages with this ladder - but do not remove from 
	 * the ladder */
	for( long ion_to=dense.IonLow[nelem]; ion_to <= dense.IonHigh[nelem]; ion_to++ )
	{
		for( long ion_from=dense.IonLow[nelem]; ion_from <= dense.IonHigh[nelem]; ++ion_from )
		{
			/* do not do ion onto itself */
			if( ion_to != ion_from )
			{
				/* CT with molecules */
				rateone = mole.xMoleChTrRate[nelem][ion_from][ion_to] * atmdat.lgCTOn;
				MAT( xmat, ion_from-ion_low, ion_from-ion_low ) -= rateone;
				MAT( xmat, ion_from-ion_low, ion_to-ion_low ) += rateone;
				/* CT with grains */
				rateone = gv.GrainChTrRate[nelem][ion_from][ion_to]*lgGrainsOn;
				MAT( xmat, ion_from-ion_low, ion_from-ion_low ) -= rateone;
				MAT( xmat, ion_from-ion_low, ion_to-ion_low ) += rateone;
			}
		}
	}

	for( long ion=dense.IonLow[nelem]; ion <= limit; ion++ )
	{
		/* thermal & secondary collisional ionization */
		ionbal.RateIoniz[nelem][ion][ion+1] +=
			ionbal.CollIonRate_Ground[nelem][ion][0] +
			secondaries.csupra[nelem][ion] +
			/* inner shell ionization by UTA lines */
			ionbal.UTA_ionize_rate[nelem][ion];

		/* loop over all atomic sub-shells to include photoionization */
		for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
		{
			/* this is the primary ionization rate - add to diagonal element,
			 * test on ion stage is so that we don't include ionization from the very highest
			 * ionization stage to even higher - since those even higher stages are not considered
			 * this would appear as a sink - but populations of this highest level is ensured to
			 * be nearly trivial and neglecting it production of even higher ionization OK */
			/* >>chng 04 nov 29 RJRW, include following in this branch so only
			 * evaluated when below ions done with iso-sequence */
			if( ion+1-ion_low < ion_range )
			{
				/* this will be redistributed into charge states in following loop */ 

				/* t_yield::Inst().nelec_eject(nelem,ion,ns) is total number of electrons that can
				 * possibly be freed 
				 * loop over nej, the number of electrons ejected including the primary,
				 * nej = 1 is primary, nej > 1 includes primary plus Auger 
				 * t_yield::Inst().elec_eject_frac is prob of nej electrons */
				for( long nej=1; nej <= t_yield::Inst().nelec_eject(nelem,ion,ns); nej++ )
				{
					/* this is the ion that is produced by this ejection,
					 * limited by highest possible stage of ionization -
					 * do not want to ignore ionization that go beyond this */
					IonProduced = MIN2(ion+nej,dense.IonHigh[nelem]);
					rateone = ionbal.PhotoRate_Shell[nelem][ion][ns][0]*
						t_yield::Inst().elec_eject_frac(nelem,ion,ns,nej-1);

					/* direct photoionization of this shell */
					ionbal.RateIoniz[nelem][ion][IonProduced] += rateone;

					/* only used for possible printout - multiple electron Auger rate  -
					 * do not count one-electron as Auger */
					if( nej>1 )
						auger[IonProduced-1] += rateone;
				}
			}
		}
	}

	for( long ion=dense.IonLow[nelem]; ion < dense.IonHigh[nelem]; ion++ )
	{
		/* this is charge transfer recombination to H and He, ionization of other species */
		double ction;
		long ipISO=nelem-ion;

		if ( nelem < t_atmdat::NCX && nelem == ipISO )
		{
			ction = atmdat.CharExcIonTotal[nelem] * iso_sp[ipISO][nelem].st[0].Pop() / SDIV(dense.xIonDense[nelem][nelem-ipISO]);
		}
		else
		{
			ction=0;
			for (long nelem1=0; nelem1 < t_atmdat::NCX; ++nelem1)
				ction += atmdat.CharExcIonOf[nelem1][nelem][ion]*dense.xIonDense[nelem1][1];
		}
		/* depopulation processes enter with negative sign */
		MAT( xmat, ion-ion_low, ion-ion_low ) -= ction;
		MAT( xmat, ion-ion_low, ion+1-ion_low ) += ction;
	}

	for( long ion=dense.IonLow[nelem]; ion < dense.IonHigh[nelem]; ion++ )
	{
		/* this is charge transfer ionization of H and He, recombination of other species */
		double ctrec;
		long ipISO=nelem-ion;

		// NCX number of ions with CX rates, so this is test on hydrogenic H, He
		if ( nelem < t_atmdat::NCX && nelem == ipISO )
		{
			ctrec = atmdat.CharExcRecTotal[nelem];
		}
		else
		{
			// the remainder, so a loop over He^0 CX
			ctrec = 0.;
			for (long nelem1=0; nelem1<t_atmdat::NCX; ++nelem1)
			{
				long ipISO = nelem1;
				ctrec +=
					atmdat.CharExcRecTo[nelem1][nelem][ion]*iso_sp[ipISO][nelem1].st[0].Pop();
			}
		}

		MAT( xmat, ion+1-ion_low, ion+1-ion_low ) -= ctrec;
		MAT( xmat, ion+1-ion_low, ion-ion_low ) += ctrec;
	}

	for( long ion_from=dense.IonLow[nelem]; ion_from < dense.IonHigh[nelem]; ion_from++ )
	{
		for( long ion_to=ion_from+1; ion_to <= dense.IonHigh[nelem]; ion_to++ )
		{
			ionbal.elecsrc[nelem] += ionbal.RateIoniz[nelem][ion_from][ion_to]*dense.xIonDense[nelem][ion_from]*
				(ion_to-ion_from);
			/* depopulation processes enter with negative sign */
			MAT( xmat, ion_from-ion_low, ion_from-ion_low ) -= ionbal.RateIoniz[nelem][ion_from][ion_to];
			MAT( xmat, ion_from-ion_low, ion_to-ion_low ) += ionbal.RateIoniz[nelem][ion_from][ion_to];
		}
	}

	for( long ion=dense.IonLow[nelem]; ion<dense.IonHigh[nelem]; ion++ )
	{
		/* loss of next higher ion due to recombination to this ion stage */
		MAT( xmat, ion+1-ion_low, ion+1-ion_low ) -= ionbal.RateRecomTot[nelem][ion];
		MAT( xmat, ion+1-ion_low, ion-ion_low ) += ionbal.RateRecomTot[nelem][ion];
		ionbal.elecsnk[nelem] += ionbal.RateRecomTot[nelem][ion]*dense.xIonDense[nelem][ion+1];
	}

	return;
}

STATIC void fill_ext_src_and_snk( long nelem, long ion_low, long ion_range, valarray<double> &xmat, valarray<double> &source )
{
	DEBUG_ENTRY( "fill_ext_src_and_snk()" );

	for( long i=0; i<ion_range;i++ )
	{
		source[i] = 0.;
	}

	for( long i=0; i<ion_range;i++ )
	{
		long ion = i+ion_low;

		/* these are the external source and sink terms */
		/* need negative sign to get positive pops */
		source[i] -= mole.source[nelem][ion];
		MAT( xmat, i, i ) -= mole.sink[nelem][ion];
	}

	return;
}

STATIC void HomogeneousSource( long nelem, long ion_low, long ion_range, valarray<double> &xmat, valarray<double> &source, double abund_total )
{
	bool lgHomogeneous = true;
	double totsrc,
		totscl,
		maxdiag;

	DEBUG_ENTRY( "lgHomogeneousSource()" );

	/* this will be sum of source and sink terms, will be used to decide if
	 * matrix is singular */
	totsrc = totscl = maxdiag = 0.;
	for( long i=0; i<ion_range;i++ )
	{
		long ion = i + ion_low;

		totsrc -= source[i];
		fixit("need better way to determine total without molecular sinks");
		// In old newmole, totscl was calculated before
		// mole.sink terms are added into xmat, but those terms are now already done.
		// the above hack takes care of that, but a cleaner solution is needed.
		double diag = -(MAT( xmat, i, i )+mole.sink[nelem][ion]);
		if (diag > maxdiag)
			maxdiag = diag;
		totscl += diag*dense.xIonDense[nelem][ion];
	}

	// Largest condition number before we decide that conservation will get
	// a better answer than the raw linear solver
	const double CONDITION_SCALE = 1e8;
	double conserve_scale = maxdiag/CONDITION_SCALE;
	ASSERT( conserve_scale > 0.0 );

	/* matrix is not homogeneous if source is non-zero */
	if( totsrc > 1e-10*totscl )
		lgHomogeneous = false;

	fixit("dynamics rates need to be moved into fill_array?");
	/* chng 03 jan 13 rjrw, add in dynamics if required here,
	 * - only do advection if we have not overrun the radius scale */
	if( dynamics.doNonEquilibriumSolve( iteration ) )
	{
		for( long i=0; i<ion_range;i++ )
		{			 
			long ion = i+ion_low;
			MAT( xmat, i, i ) -= dynamics.Rate;
			source[i] -= dynamics.Source[nelem][ion];
			/* fprintf(ioQQQ," %li %li %.3e (%.3e %.3e)\n",i,i,MAT(*xmat,i,i),
			  dynamics.Rate, dynamics.Source[nelem][ion]);*/
		}
		lgHomogeneous = false;
	}

	/* >>chng 06 nov 21, for very high ionization parameter sims the H molecular
	* fraction can become so small that atom = atom + molecule.  In this case we
	* will not count system as an inhomogeneous case since linear algebra package
	* will fail.  If molecules are very small, count as homogeneous.
	* Note that we have already added sink terms to the main matrix and they
	* will not be removed, a possible source of error, but they must not
	* have been significant, given that the molecular fraction is so small */
	if( !lgHomogeneous && ion_range==2 )
	{
		/* solve algebraically */
		double a = MAT( xmat, 0, 0 ), 
			b = MAT( xmat, 1, 0 ) ,
			c = MAT( xmat, 0, 1 ) ,
			d = MAT( xmat, 1, 1 );
		b = SDIV(b);
		d = SDIV(d);
		double ratio1 = a/b , ratio2 = c/d , fratio1=fabs(a/b),fratio2=fabs(c/d);
		if( fabs(ratio1-ratio2) <= DBL_EPSILON*1e4*MAX2(fratio1,fratio2) )
		{
			lgHomogeneous = true;
		}
	}

#if 1
	if(
		fabs(dense.xMolecules(nelem)) <DBL_EPSILON*100.*dense.gas_phase[nelem] )
	{
		lgHomogeneous = true;
	}
#endif

	/* this is true if no source terms 
	 * we will use total population and species conservation to replace one
	 * set of balance equations since overdetermined */
	if( 1 || lgHomogeneous  )
	{
		double rate_ioniz=1., rate_recomb /*, scale = 0.*/;
		/* Simple estimate of most abundant ion */
		long jmax = 0;
		for( long i=0; i<ion_range-1;i++)
		{ 
			long ion = i+ion_low;
			rate_ioniz *= ionbal.RateIonizTot(nelem,ion);
			rate_recomb = ionbal.RateRecomTot[nelem][ion];
			/* trips if ion rate zero, so all the gas will be more neutral than this */
			if( rate_ioniz == 0)
				break;
			/* rec rate is zero */
			if( rate_recomb <= 0.) 
				break;

			rate_ioniz /= rate_recomb;
			if( rate_ioniz > 1.) 
			{
				/* this is peak ionization stage */
				jmax = i;
				rate_ioniz = 1.;
			}
		}

		/* replace its matrix elements with population sum */
		for( long i=0; i<ion_range;i++ )
		{
			MAT(xmat,i,jmax) -= conserve_scale;
		}
		source[jmax] -= abund_total*conserve_scale;
	}

#if 0
	if( false && nelem == ipHYDROGEN && (dynamics.lgAdvection || dynamics.lgTimeDependentStatic) && iteration>1 )
	{
		fprintf(ioQQQ,"DEBUGG Rate %.2f %.3e \n",fnzone,dynamics.Rate);
		fprintf(ioQQQ," %.3e %.3e\n", ionbal.RateIonizTot(nelem,0), ionbal.RateIonizTot(nelem,1) );
		fprintf(ioQQQ," %.3e %.3e\n", ionbal.RateRecomTot[nelem][0], ionbal.RateRecomTot[nelem][1]);
		fprintf(ioQQQ," %.3e %.3e %.3e\n\n", dynamics.Source[nelem][0], dynamics.Source[nelem][1], dynamics.Source[nelem][2]);
	}

	{
		/* option to print matrix */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone==1 && lgPrintIt )
		{
			fprintf( ioQQQ, 
				" DEBUG ion_solver: nelem=%li ion_range=%li, limit=%li, nConv %li xmat follows\n",
				nelem , ion_range,limit  , conv.nTotalIoniz );
			if( lgHomogeneous )
				fprintf(ioQQQ , "Homogeneous \n");
			for( long i=0; i<ion_range; ++i )
			{
				for( long j=0;j<ion_range;j++ )
				{
					fprintf(ioQQQ,"%e\t",MAT(xmat,i,j));
				}
				fprintf(ioQQQ,"\n");
			}
			fprintf(ioQQQ,"source follows\n");
			for( long i=0; i<ion_range;i++ )
			{
				fprintf(ioQQQ,"%e\t",source[i]);
			}
			fprintf(ioQQQ,"\n");
			fprintf(ioQQQ,"totsrc/totscl %g %g\n",totsrc,totscl);
		}
	}
#endif

}

STATIC void PrintRates( long nelem, bool lgNegPop, double abund_total, valarray<double> &auger, bool lgPrintIt )
{
	DEBUG_ENTRY( "PrintRates()" );

	long ion;

	/* this should not happen */
	if( lgNegPop )
	{
		fprintf( ioQQQ, " PROBLEM Negative population found for abundance of ionization stage of element %4.4s, ZONE=%4ld\n", 
		  elementnames.chElementNameShort[nelem], nzone );

		fprintf( ioQQQ, " Populations were" );
		for( ion=1; ion <= dense.IonHigh[nelem]+1; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", dense.xIonDense[nelem][ion-1] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " destroy vector =" );
		for( ion=1; ion <= dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", ionbal.RateIonizTot(nelem,ion-1) );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " CTHeavy  vector =" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", atmdat.CharExcIonOf[ipHELIUM][nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " CharExcIonOf[ipHYDROGEN] vtr=" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " CollidRate  vtr=" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", ionbal.CollIonRate_Ground[nelem][ion][0] );
		}
		fprintf( ioQQQ, "\n" );

		/* photo rates per subshell */
		fprintf( ioQQQ, " photo rates per subshell, ion\n" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%3ld", ion );
			for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
			{
				fprintf( ioQQQ, "%9.1e", ionbal.PhotoRate_Shell[nelem][ion][ns][0] );
			}
			fprintf( ioQQQ, "\n" );
		}

		/* now check out creation vector */
		fprintf( ioQQQ, " create  vector =" );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, "%9.1e", ionbal.RateRecomTot[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );
	}

	/* option to print ionization and recombination arrays
	 * prt flag set with print array print arrays command */
	if( lgPrintIt || prt.lgPrtArry[nelem] || lgNegPop )
	{
		const char* format = " %9.2e";
		/* say who we are, what we are doing .... */
		fprintf( ioQQQ, 
			"\n %s ion_solver ion/rec rt [s-1] %s nz%.2f Te%.4e ne%.4e Tot abun:%.3e ion abun%.2e mole%.2e Search? %c\n",
			elementnames.chElementSym[nelem],
			elementnames.chElementName[nelem],
			fnzone,
			phycon.te , 
			dense.eden,
			dense.gas_phase[nelem],
			abund_total ,
			dense.xMolecules(nelem) , TorF(conv.lgSearch)  );
		/* total ionization rate, all processes */
		fprintf( ioQQQ, " %s Ioniz total " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, ionbal.RateIonizTot(nelem,ion) );
		}
		fprintf( ioQQQ, "\n" );

		/* sum of all creation processes */
		fprintf( ioQQQ, " %s Recom total " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, ionbal.RateRecomTot[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* collisional ionization */
		fprintf( ioQQQ, " %s Coll ioniz  " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double ColIoniz = ionbal.CollIonRate_Ground[nelem][ion][0];
			if( ion > nelem - NISO )
			{
				long ipISO = nelem-ion;
				ASSERT( ipISO >=0 && ipISO < NISO );
				if( dense.xIonDense[nelem][nelem-ipISO] > 0. )
				{
					ColIoniz *= iso_sp[ipISO][nelem].st[0].Pop();
					for( long ipLevel=1; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
					{
						ColIoniz += iso_sp[ipISO][nelem].fb[ipLevel].ColIoniz * dense.EdenHCorr * iso_sp[ipISO][nelem].st[ipLevel].Pop();
					}
					ColIoniz /= dense.xIonDense[nelem][nelem-ipISO];
				}
			}
			fprintf( ioQQQ, format, ColIoniz );
		}
		fprintf( ioQQQ, "\n" );		

		/* UTA ionization */
		fprintf( ioQQQ, " %s UTA ioniz   " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, ionbal.UTA_ionize_rate[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* photo ionization */
		fprintf( ioQQQ, " %s Photoion snk" ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double PhotIoniz = 0.;
			for( long ipShell = 0; ipShell < Heavy.nsShells[nelem][ion]; ipShell++ )
				PhotIoniz += ionbal.PhotoRate_Shell[nelem][ion][ipShell][0];
			
			// still don't have the total if stage is in one of the iso-sequences
			if( ion > nelem - NISO )
			{
				long ipISO = nelem-ion;
				ASSERT( ipISO >=0 && ipISO < NISO );
				if( dense.xIonDense[nelem][nelem-ipISO]>0  )
				{
					PhotIoniz *= iso_sp[ipISO][nelem].st[0].Pop();
					for( long ipLevel=1; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
					{
						PhotIoniz += iso_sp[ipISO][nelem].fb[ipLevel].gamnc * iso_sp[ipISO][nelem].st[ipLevel].Pop();
					}
					PhotIoniz /= dense.xIonDense[nelem][nelem-ipISO];
				}
			}
			fprintf( ioQQQ, format, PhotIoniz );
		}
		fprintf( ioQQQ, "\n" );		

		/* photoionization source (source of this stage due to ionization of next stage) */
		fprintf( ioQQQ, " %s Photoion src" ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double source = 0.;
			if( ion>0 )
				source = ionbal.RateIoniz[nelem][ion-1][ion];
				
			fprintf( ioQQQ, format, source );
		}
		fprintf( ioQQQ, "\n" );
				
		/* auger production (of this stage due to auger ionization of another) */
		fprintf( ioQQQ, " %s Auger src   " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double source = 0.;
			if( ion>0 )
				source = auger[ion-1];
			
			fprintf( ioQQQ, format, source );
		}
		fprintf( ioQQQ, "\n" );

		/* secondary ionization */
		fprintf( ioQQQ, " %s Secon ioniz " ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, 
				secondaries.csupra[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* grain ionization - not total rate but should be dominant process */
		fprintf( ioQQQ, " %s ion on grn  "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, gv.GrainChTrRate[nelem][ion][ion+1] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge transfer ionization from chemistry */
		fprintf( ioQQQ, " %s ion xfr mol "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, mole.xMoleChTrRate[nelem][ion][ion+1] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge exchange ionization */
		fprintf( ioQQQ, " %s chr trn ion " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			/* sum has units s-1 */
			double sum;
			long ipISO=nelem-ion;

			if( nelem < t_atmdat::NCX && nelem == ipISO )
			{
				sum = atmdat.CharExcIonTotal[nelem] * iso_sp[ipISO][nelem].st[0].Pop() / SDIV(dense.xIonDense[nelem][nelem-ipISO]);
			}
			else
			{
				sum = 0.0;
				for (long nelem1=0; nelem1 < t_atmdat::NCX; ++nelem1)
					sum += atmdat.CharExcIonOf[nelem1][nelem][ion] * dense.xIonDense[nelem1][1];
			}
			
			fprintf( ioQQQ, format, sum );
		}
		fprintf( ioQQQ, "\n" );

		/* radiative recombination */
		fprintf( ioQQQ, " %s radiati rec "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, dense.eden*ionbal.RR_rate_coef_used[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* Badnell DR recombination */
		fprintf( ioQQQ, " %s drBadnel rec"  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < min(nelem-1,dense.IonHigh[nelem]); ion++ )
		{
			fprintf( ioQQQ, format, dense.eden*ionbal.DR_Badnell_rate_coef[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* Cota rate */
		fprintf( ioQQQ, " %s CotaRate rec"  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, dense.eden*ionbal.CotaRate[ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* grain recombination - not total but from next higher ion, should
		 * be dominant */
		fprintf( ioQQQ, " %s rec on grn  "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, gv.GrainChTrRate[nelem][ion+1][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge transfer recombination from chemistry */
		fprintf( ioQQQ, " %s rec xfr mol "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, mole.xMoleChTrRate[nelem][ion+1][ion] );
		}
		fprintf( ioQQQ, "\n" );

		/* charge exchange recombination */
		fprintf( ioQQQ, " %s chr trn rec "  ,elementnames.chElementSym[nelem]);
		for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
		{
			double sum;
			
			if( nelem==ipHELIUM && ion==0 )
			{
				sum = atmdat.CharExcRecTotal[nelem];
			}
			else if( nelem==ipHYDROGEN && ion==0 )
			{
				sum = atmdat.CharExcRecTotal[nelem];
			}
			else
			{
				sum = 0.0;
				for (long nelem1=0; nelem1<t_atmdat::NCX; ++nelem1)
				{
					long ipISO=nelem1;
					sum += atmdat.CharExcRecTo[nelem1][nelem][ion] * iso_sp[ipISO][nelem1].st[0].Pop();
				}
			}
			fprintf( ioQQQ, format, sum );
		}
		fprintf( ioQQQ, "\n" );

		{
			/* sources from the chemistry network */
			fprintf(ioQQQ," %s molecule src",
					elementnames.chElementSym[nelem]);
			for( ion=0; ion <= dense.IonHigh[nelem]; ion++ )
			{
				fprintf(ioQQQ,format, mole.source[nelem][ion]/SDIV(dense.xIonDense[nelem][ion]) );
			}
			fprintf( ioQQQ, "\n" );

			/* sinks from the chemistry network */
			fprintf(ioQQQ," %s molecule snk",
					  elementnames.chElementSym[nelem]);
			for( ion=0; ion <= dense.IonHigh[nelem]; ion++ )
			{
				fprintf(ioQQQ,format, mole.sink[nelem][ion] );
			}
			fprintf( ioQQQ, "\n" );
						
			if( dynamics.lgAdvection || dynamics.lgTimeDependentStatic )
			{
				/* source from the dynamcs */
				fprintf(ioQQQ," %s dynamics src",
						elementnames.chElementSym[nelem]);
				for( ion=0; ion <= dense.IonHigh[nelem]; ion++ )
				{
					fprintf(ioQQQ,format, dynamics.Source[nelem][ion]/SDIV(dense.xIonDense[nelem][ion]) );
				}
				fprintf( ioQQQ, "\n" );

				/* sinks from the dynamics */
				fprintf(ioQQQ," %s dynamics snk",
						  elementnames.chElementSym[nelem]);
				for( ion=0; ion <= dense.IonHigh[nelem]; ion++ )
				{
					fprintf(ioQQQ,format, dynamics.Rate );
				}
				fprintf( ioQQQ, "\n" );
			}
		}
		
		/* the "new" abundances the resulted from the previous ratio */
		fprintf( ioQQQ, " %s Abun [cm-3] " ,elementnames.chElementSym[nelem] );
		for( ion=0; ion <= dense.IonHigh[nelem]; ion++ )
		{
			fprintf( ioQQQ, format, dense.xIonDense[nelem][ion] );
		}
		fprintf( ioQQQ, "\n" );
	}

	if( lgNegPop )
	{
		ContNegative();
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	return;
}

/* 

	 Solve an ionization level system with specified ionization and
	 recombination rates between neighboring ions, and additional sink
	 and source terms.  The sink array is overwritten, and the results
	 appear in the source array.

	 Written in matrix form, the algorithm is equivalent to the
	 tridiagonal algorithm in Numerical Recipes applied to:

	 / i_0+a_0     -r_0          .           .    .  \ / x_0 \   / s_0 \
	 |  -i_0    i_1+a_1+r_0    -r_1          .    .  | | x_1 |   | s_1 |
	 |    .        -i_1      i_2+a_2+r_1   -r_2   .  | | x_2 |   | s_2 |
	 |    .          .       (etc....)               | | ... | = | ... |
	 \    .          .          .                    / \     /   \     /

	 where i, r are the ionization and recombination rates, s is the
	 source rate and a is the sink rate.

	 This matrix is diagonally dominant only when the sink terms are
	 large -- the alternative method coded here prevents rounding error
	 in the diagonal terms disturbing the solution when this is not the
	 case.

*/
	
/* solveions tridiagonal solver but optimized for structure of balance matrix */
void solveions(double *ion, double *rec, double *snk, double *src,
	       long int nlev, long int nmax)
{
	double kap, bet;
	long int i;

	DEBUG_ENTRY( "solveions()" );

	if(nmax != -1) 
	{
		/* Singular case */
		src[nmax] = 1.;
		for(i=nmax;i<nlev-1;i++)
			src[i+1] = src[i]*ion[i]/rec[i];
		for(i=nmax-1;i>=0;i--)
			src[i] = src[i+1]*rec[i]/ion[i];
	} 
	else 
	{
		kap = snk[0];    
		for(i=0;i<nlev-1;i++) 
		{
			bet = ion[i]+kap;
			if(bet == 0.)
			{
				fprintf(ioQQQ,"Ionization solver error\n");
				cdEXIT(EXIT_FAILURE);
			}
			bet = 1./bet;
			src[i] *= bet;
			src[i+1] += ion[i]*src[i];
			snk[i] = bet*rec[i];
			kap = kap*snk[i]+snk[i+1];
		}
		bet = kap;
		if(bet == 0.)
		{
			fprintf(ioQQQ,"Ionization solver error\n");
			cdEXIT(EXIT_FAILURE);
		}
		src[i] /= bet;

		for(i=nlev-2;i>=0;i--)
		{
			src[i] += snk[i]*src[i+1];
		}
	}
}

void ion_wrapper( long nelem )
{

	DEBUG_ENTRY( "ion_wrapper()" );

	ASSERT( nelem >= ipHYDROGEN );
	ASSERT( nelem < LIMELM );

	switch( nelem )
	{
	case ipHYDROGEN:
		IonHydro();
		break;
	case ipHELIUM:
		IonHelium();
		break;
	default:
		IonNelem(false,nelem);
		break;
	}
	
	if( trace.lgTrace && dense.lgElmtOn[nelem] && trace.lgHeavyBug )
	{
		fprintf( ioQQQ, "     ion_wrapper returns; %s fracs = ", elementnames.chElementSym[nelem] );
		for( long ion = 0; ion<=nelem+1; ion++ )
			fprintf( ioQQQ,"%10.3e ", dense.xIonDense[nelem][ion]/dense.gas_phase[nelem] );
		fprintf( ioQQQ, "\n" );
	}	

	ASSERT(lgElemsConserved());

	return;
}
