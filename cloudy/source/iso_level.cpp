/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_level solve for iso-sequence level populations */
#include "cddefines.h"
#include "atmdat.h"
#include "continuum.h"
#include "conv.h"
#include "doppvel.h"
#include "dynamics.h"
#include "elementnames.h"
#include "grainvar.h"
#include "he.h"
#include "ionbal.h"
#include "iso.h"
#include "opacity.h"
#include "phycon.h"
#include "rfield.h"
#include "trace.h"
#include "mole.h"
#include "freebound.h"
#include "two_photon.h"
#include "dense.h"
#include "vectorize.h"
#include "prt.h"
#include "save.h"
#include "iterations.h"
#include "species.h"


STATIC void iso_multiplet_opacities_one(
	const long int ipISO, const long int nelem);

/* solve for level populations  */
void iso_level( const long int ipISO, const long int nelem, double &renorm )
{
	long int ipHi,
		ipLo,
		i,
		level,
		level_error;

 	t_iso_sp* sp = &iso_sp[ipISO][nelem]; 
	const long int numlevels_local = sp->numLevels_local;

	double BigError;

	int32 nerror;
	bool lgNegPop=false;
	static vector<int32> ipiv;
	ipiv.resize(numlevels_local); 
	/* this block of variables will be obtained and freed here */
	double source=0., sink=0.;
	static vector<double> PopPerN;
	PopPerN.resize(sp->n_HighestResolved_local+1);

	DEBUG_ENTRY( "iso_level()" );
	renorm = 1.;

	/* check that we were called with valid charge */
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );

	/* these two collision rates must be the same or we are in big trouble,
	 * since used interchangeably */
	ASSERT( ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0]< SMALLFLOAT ||
		fabs( (sp->fb[0].ColIoniz* dense.EdenHCorr) /
		SDIV(ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0] ) - 1.) < 0.001 );

	if( (dense.IonHigh[nelem] >= nelem - ipISO) &&
		 (dense.IonLow[nelem] <= nelem - ipISO) &&
		 ionbal.RateRecomIso[nelem][ipISO] > 0. )
	{
		/* get simple estimate of atom to ion ratio */
		sp->xIonSimple = ionbal.RateIonizTot(nelem,nelem-ipISO)/ionbal.RateRecomIso[nelem][ipISO];
	}
	else
	{
		sp->xIonSimple = 0.;
	}

	/* which case atom to solve??? */
	if( dense.xIonDense[nelem][nelem+1-ipISO] < SMALLFLOAT )
	{
		/* don't bother if no ionizing radiation */
		strcpy( sp->chTypeAtomUsed, "zero " );
		if( trace.lgTrace && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_level iso=%2ld nelem=%2ld simple II/I=%10.2e so not doing equilibrium, doing %s.\n", 
				ipISO, nelem, sp->xIonSimple, sp->chTypeAtomUsed );
		}

		/* total ionization will just be the ground state */
		lgNegPop = false;
		sp->st[0].Pop() = dense.xIonDense[nelem][nelem-ipISO];
		for( long n=1; n < numlevels_local; n++ )
		{
			sp->st[n].Pop() =  0.;
		}
		sp->qTot2S = 0.;
	}
	else
	{
		/* fill in recombination vector - values were set in iso_ionize_recombine.cpp */
		static vector<double> creation, error, work, Save_creation;
		creation.resize(numlevels_local);
		error.resize(numlevels_local);
		work.resize(numlevels_local);
		Save_creation.resize(numlevels_local);

		ASSERT( dense.xIonDense[nelem][nelem+1-ipISO] >= 0.f );
		for( level=0; level < numlevels_local; level++ )
		{
			/* total recombination from once more ionized [cm-3 s-1] */
			creation[level] = sp->fb[level].RateCont2Level * dense.xIonDense[nelem][nelem+1-ipISO];
		}
		
		double ctsource=0.0, ctsink=0.0, ctrec=0.0;
		/* now charge transfer - all into/from ground, two cases, H and not H */
		if( nelem==ipHYDROGEN )
		{
			/* charge transfer, hydrogen onto everything else */
			/* charge exchange recombination */
			ctrec += atmdat.CharExcRecTotal[nelem];
			ctsink += atmdat.CharExcIonTotal[nelem];
		}
		else if( nelem==ipHELIUM && ipISO==ipHE_LIKE )
		{
			/* this is recom of He due to ct with all other gas constituents */
			ctrec += atmdat.CharExcRecTotal[nelem];
			ctsink += atmdat.CharExcIonTotal[nelem];
		}
		else
		{
			for (long nelem1=ipHYDROGEN; nelem1 < t_atmdat::NCX; ++nelem1)
			{
				long ipISO=nelem1;
				ctrec += atmdat.CharExcRecTo[nelem1][nelem][nelem-ipISO]*iso_sp[ipISO][nelem1].st[0].Pop();
				ctsink += atmdat.CharExcIonOf[nelem1][nelem][nelem-ipISO]*dense.xIonDense[nelem1][1];
			}					
		}
		ctsource += ctrec*dense.xIonDense[nelem][nelem+1-ipISO];
		
		if ( nelem > ipISO )
		{
			double ction=0.0;
			if( nelem==ipHELIUM )
			{
				ctsink += atmdat.CharExcRecTotal[nelem];
				ction += atmdat.CharExcIonTotal[nelem];
			}
			else
			{
				for (long nelem1=ipHYDROGEN; nelem1 < t_atmdat::NCX; ++nelem1)
				{
					long ipISO1=nelem1;
					ctsink += atmdat.CharExcRecTo[nelem1][nelem][nelem-(ipISO+1)]*iso_sp[ipISO1][nelem1].st[0].Pop();
					ction += atmdat.CharExcIonOf[nelem1][nelem][nelem-(ipISO+1)]*dense.xIonDense[nelem1][1];
				}
			}
			ctsource += ction * dense.xIonDense[nelem][nelem-(ipISO+1)];
		}
		
		/* now do the 2D array */
		multi_arr<double,2,C_TYPE> z, SaveZ;
		z.alloc(numlevels_local,numlevels_local);
		z.zero();

		/* this branch is main solver, full level populations 
		 * assert since this code must change if NISO ever increased */
		ASSERT( NISO == 2 );

		strcpy( sp->chTypeAtomUsed, "Popul" );
		if( trace.lgTrace && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_level iso=%2ld nelem=%2ld doing regular matrix inversion, %s\n", 
				ipISO, nelem, sp->chTypeAtomUsed );
		}

		//qList::const_iterator StElm = StatesElemNEW[nelem][nelem-ipISO].begin();
		qList::const_iterator StElm = sp->st.begin();
		static vector<double> Boltzmann_overg;
		Boltzmann_overg.resize(numlevels_local-1);
		for (unsigned lev = 0; lev < Boltzmann_overg.size(); ++lev)
			Boltzmann_overg[lev] = 1.0/(double)StElm[lev].g();
		
		/* >>chng 05 dec 21, rm eden to make into rate coefficient */
		sp->qTot2S = sp->fb[1].ColIoniz;

		static vector<double> coll_down, RadDecay, pump;
		coll_down.resize(numlevels_local);
		RadDecay.resize(numlevels_local);
		pump.resize(numlevels_local);
		avx_ptr<double> arg(1,numlevels_local), bstep(1,numlevels_local);

		for( level=1; level < numlevels_local; level++ )
		{
			fixit("Wouldn't need to mask this out if levels were in order");
			arg[level] = -max((StElm[level].energy().K()-StElm[level-1].energy().K()) / phycon.te, -38.);
		}
		vexp( arg.ptr0(), bstep.ptr0(), 1, numlevels_local );

		ColliderDensities colld(colliders);

		enum { DEBUG_RATES = false };

		if( DEBUG_RATES && iterations.lgLastIt )
			fprintf( stdout, "# ipISO\tnelem\tlevel\tcollDown\tcollIonz\tradRecom\tRadDecay\n" );

		/* master balance equation, use when significant population */
		for( level=0; level < numlevels_local; level++ )
		{
			double coll_down_total = 0.;

			/* all process depopulating level and placing into the continuum
			 * this does NOT include grain charge transfer ionization, added below */
			z[level][level] += sp->fb[level].RateLevel2Cont;

			if (level != 0)
			{
				for ( ipLo = 0; ipLo < level; ++ipLo )
					Boltzmann_overg[ipLo] *= bstep[level];
			}

			/* all processes populating level from below */
			for( ipLo=0; ipLo < level; ipLo++ )
			{
				coll_down[ipLo] = sp->trans(level,ipLo).Coll().ColUL( colld );
				if( DEBUG_RATES )
				{
					coll_down_total += coll_down[ipLo];
				}

				if ( rfield.lgPlasNu && sp->trans(level,ipLo).EnergyRyd()<rfield.plsfrq  )
				{
					RadDecay[ipLo] = iso_ctrl.SmallA;
					pump[ipLo] = iso_ctrl.SmallA;
				}
				else
				{
					RadDecay[ipLo] = MAX2( iso_ctrl.SmallA, sp->trans(level,ipLo).Emis().Aul()*
							       		  (sp->trans(level,ipLo).Emis().Ploss()) );
					pump[ipLo] = MAX2( iso_ctrl.SmallA, sp->trans(level,ipLo).Emis().pump() );
				}
			}

			if( iso_ctrl.lgRandErrGen[ipISO] )
			{
				for( ipLo=0; ipLo < level; ipLo++ )
				{
					coll_down[ipLo] *= sp->ex[level][ipLo].ErrorFactor[IPCOLLIS];
					RadDecay[ipLo] *= sp->ex[level][ipLo].ErrorFactor[IPRAD];
					pump[ipLo] *= sp->ex[level][ipLo].ErrorFactor[IPRAD];
				}
			}

			double glev = (double)StElm[level].g(), rglev = 1.0/glev;
			for( ipLo=0; ipLo < level; ipLo++ )
			{
				double coll_up = coll_down[ipLo] * glev *
					Boltzmann_overg[ipLo];

				z[ipLo][ipLo] += coll_up + pump[ipLo] ;
				z[ipLo][level] -= coll_up + pump[ipLo] ;

				double pump_down = pump[ipLo] *
					(double)StElm[ipLo].g() * rglev;

				z[level][level] += RadDecay[ipLo] + pump_down + coll_down[ipLo];
				z[level][ipLo] -= RadDecay[ipLo] + pump_down + coll_down[ipLo];

				/* find total collisions out of 2^3S to singlets. */
				if( (level == 1) && (ipLo==0) )
				{
					sp->qTot2S += coll_down[ipLo]/dense.EdenHCorr;
				}
				if( (ipLo == 1) && (ipISO==ipH_LIKE || StElm[level].S()==1) )
				{
					sp->qTot2S += coll_up/dense.EdenHCorr;
				}
			}

			if( DEBUG_RATES && iterations.lgLastIt )
			{
				fprintf( stdout,
					"%2ld\t%2ld\t%2ld\t"
					"%.4e\t"
					"%.4e\t"
					"%.4e\t"
					"%.4e\n",
					ipISO, nelem, level,
					coll_down_total,
					iso_sp[ipISO][nelem].fb[level].ColIoniz * dense.eden,
					iso_sp[ipISO][nelem].fb[level].RadRecomb[0] * dense.eden,
					1. / iso_sp[ipISO][nelem].st[level].lifetime() );
			}
		}

		if( DEBUG_RATES && iterations.lgLastIt )
			fprintf( stdout, "\n\n" );

		/** \todo 2 the indices for the two-photon rates must be changed for further iso sequences. */  
		ASSERT( ipISO <= ipHE_LIKE );
		for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
		{
			// induced two photon emission - special because upward and downward are
			// not related by ratio of statistical weights 
			// iso.lgInd2nu_On is controlled with SET IND2 ON/OFF command 

			fixit("need Pesc or the equivalent to multiply AulTotal?");
			// downward rate
			z[tnu->ipHi][tnu->ipLo] -= tnu->AulTotal;
			z[tnu->ipHi][tnu->ipHi] += tnu->AulTotal; 

		}
		
		if (iso_ctrl.lgInd2nu_On)
		{
			for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
			{
				// induced two photon emission - special because upward and downward are
				// not related by ratio of statistical weights 
				// iso.lgInd2nu_On is controlled with SET IND2 ON/OFF command 
				
				z[tnu->ipHi][tnu->ipLo] -= tnu->induc_dn;
				z[tnu->ipHi][tnu->ipHi] += tnu->induc_dn;
				
				// upward rate
				z[tnu->ipLo][tnu->ipHi] -= tnu->induc_up;
				z[tnu->ipLo][tnu->ipLo] += tnu->induc_up;
			}
		}
		
		/* grain charge transfer recombination and ionization to ALL other stages */
		for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
		{
			if( ion!=nelem-ipISO )
			{
				source += gv.GrainChTrRate[nelem][ion][nelem-ipISO] *
					dense.xIonDense[nelem][ion];
				sink += gv.GrainChTrRate[nelem][nelem-ipISO][ion];

				source += mole.xMoleChTrRate[nelem][ion][nelem-ipISO] * 
					dense.xIonDense[nelem][ion] * atmdat.lgCTOn;
				sink += mole.xMoleChTrRate[nelem][nelem-ipISO][ion] * atmdat.lgCTOn;

			}
		}
		
		/* add in source and sink terms from molecular network. */
		source += mole.source[nelem][nelem-ipISO];
		sink += mole.sink[nelem][nelem-ipISO];

#if	1
		/* >>chng 02 Sep 06 rjrw -- all elements have these terms */
		/*>>>chng 02 oct 01, only include if lgAdvection or lgTimeDependentStatic is set */
		if( dynamics.doNonEquilibriumSolve( iteration ) && dynamics.lgISO[ipISO] )
		{
			/* add in advection - these terms normally zero */
			source += dynamics.Source[nelem][nelem-ipISO];
			/* >>chng 02 Sep 06 rjrw -- advective term not recombination */
			sink += dynamics.Rate;
		}
#else
		/*>>>chng 02 oct 01, only include if lgAdvection or lgTimeDependentStatic is set */
		if( dynamics.doNonEquilibriumSolve( iteration ) && dynamics.lgISO[ipISO] )
		{
			for( level=0; level < numlevels_local; level++ )
			{
				creation[level] += dynamics.StatesElem[nelem][nelem-ipISO][level];
				z[level][level] += dynamics.Rate;
			}
		}
#endif

		/* ionization from/recombination from lower ionization stages */
		for(long ion_from=dense.IonLow[nelem]; ion_from < MIN2( dense.IonHigh[nelem], nelem-ipISO ) ; ion_from++ )
		{
			if( ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] >= 0. )
			{			  
				/* ionization from lower ionization stages, cm-3 s-1 */
				source += ionbal.RateIoniz[nelem][ion_from][nelem-ipISO] * dense.xIonDense[nelem][ion_from];
			}
			/* recombination to next lower ionization stage, s-1 */
			if( ion_from == nelem-1-ipISO )
				sink += ionbal.RateRecomTot[nelem][ion_from];
		}

		ASSERT( source >= 0.f );
		if (0)
		{
			/*
			 * Collisional ionization and photoionization can only be to the ground state in H iso-sequences
			 * to conserve energy.
			 */
			creation[0] += source;
			/*for( level=0; level < numlevels_local; level++ )
			{
				z[level][level] += sink;
			}*/
			/*recombination is only done from the ground state */
			z[0][0] += sink;
		}
		else
		{
			// Try Boltzmann weighting to capture LTE limit correctly
			t_iso_sp* sp = &iso_sp[ipISO][nelem];
			double partfun=0.0;
			for ( level = 0; level < numlevels_local; level++ )
			{
				partfun += sp->st[level].Boltzmann()*sp->st[level].g();
			}
			source /= partfun;
			for( level=0; level < numlevels_local; level++ )
				{
				creation[level] += source*
						sp->st[level].Boltzmann()*sp->st[level].g();
				z[level][level] += sink;
				}
		}
		creation[0] += ctsource;
		z[0][0] += ctsink;

		/* >>chng 04 nov 30, atom XX-like collisions off turns this off too */
		if( sp->trans(iso_ctrl.nLyaLevel[ipISO],0).Coll().rate_lu_nontherm() * iso_ctrl.lgColl_excite[ipISO] > 0. )
		{
			/* now add on supra thermal excitation */
			for( level=1; level < numlevels_local; level++ )
			{
				double RateUp , RateDown;

				RateUp = sp->trans(level,0).Coll().rate_lu_nontherm();
				RateDown = RateUp * (double)sp->st[ipH1s].g() /
					(double)sp->st[level].g();

				/* total rate out of lower level */
				z[ipH1s][ipH1s] += RateUp;

				/* rate from the upper level to ground */
				z[level][ipH1s] -= RateDown;

				/* rate from ground to upper level */
				z[ipH1s][level] -= RateUp;

				z[level][level] += RateDown;  
			}
		}

		/* =================================================================== 
		 *
		 * at this point all matrix elements have been established 
		 *
		 * ==================================================================== */
		/* save matrix, this allocates SaveZ */
		SaveZ = z;

		for( ipLo=0; ipLo < numlevels_local; ipLo++ )
			Save_creation[ipLo] = creation[ipLo];

		if( trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			const long numlevels_print = numlevels_local;
			fprintf( ioQQQ, "  pop level     others => (iso_level)\n" );
			for( long n=0; n < numlevels_print; n++ )
			{
				fprintf( ioQQQ, "  %s %s %2ld", iso_ctrl.chISO[ipISO], elementnames.chElementNameShort[nelem], n );
				for( long j=0; j < numlevels_print; j++ )
				{
					fprintf( ioQQQ,"\t%.9e", z[j][n] );
				}
				fprintf( ioQQQ, "\n" );
			}
			fprintf(ioQQQ," recomb ct %.2e co %.2e hectr %.2e hctr %.2e\n",
					  atmdat.CharExcRecTotal[ipHELIUM],
					  findspecieslocal("CO")->den ,
					  atmdat.CharExcRecTo[ipHELIUM][nelem][nelem-ipISO]*iso_sp[ipHE_LIKE][ipHELIUM].st[0].Pop() ,
					  atmdat.CharExcRecTo[ipHYDROGEN][nelem][nelem-ipISO]*iso_sp[ipH_LIKE][ipHYDROGEN].st[0].Pop() );
			fprintf(ioQQQ," recomb          ");
			for( long n=0; n < numlevels_print; n++ )
			{
				fprintf( ioQQQ,"\t%.9e", creation[n] );
			}
			fprintf( ioQQQ, "\n" );
		}


		if( sp->lgPrtMatrix )
		{
			valarray<double> c( get_ptr(creation), creation.size() );
			prt.matrix.prtRates( numlevels_local, z, c );
		}

		if( sp->lgImgMatrix && save.img_matrix.matchIteration( iteration ) &&
				save.img_matrix.matchZone( nzone ) )
		{
			valarray<double> SaveC( get_ptr(Save_creation), Save_creation.size() );
			save.img_matrix.createImage( iteration, nzone,
							numlevels_local, SaveZ, SaveC );
		}


		nerror = 0;

		getrf_wrapper(numlevels_local,numlevels_local,
						  z.data(),numlevels_local,&ipiv[0],&nerror);
		
		getrs_wrapper('N',numlevels_local,1,z.data(),numlevels_local,&ipiv[0],
						  &creation[0],numlevels_local,&nerror);

		if( nerror != 0 )
		{
			fprintf( ioQQQ, " iso_level: dgetrs finds singular or ill-conditioned matrix\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( sp->lgImgMatrix && save.img_matrix.matchIteration( iteration ) &&
				save.img_matrix.matchZone( nzone ) )
		{
			valarray<double> SavePops( get_ptr(creation), creation.size() );
			save.img_matrix.addImagePop_FITS( iteration, nzone,
								numlevels_local, SavePops );
		}


		/* check whether solution is valid */
		/* >>chng 06 aug 28, both of these from numLevels_max to _local. */
		for( level=ipH1s; level < numlevels_local; level++ )
		{
			double qn = 0., qx = 0.;
			error[level] = 0.;
			for( long n=ipH1s; n < numlevels_local; n++ )
			{
				double q = SaveZ[n][level]*creation[n];
				
				/* remember the largest size of element in sum to div by below */
				if ( q > qx )
					qx = q;
				else if (q < qn)
					qn = q;

				error[level] += q;
			}
			
			if (-qn > qx)
				qx = -qn;

			if( qx > 0. )
			{
				error[level] = (error[level] - Save_creation[level])/qx;
			}
			else
			{
				error[level] = 0.;
			}
		}

		/* remember largest residual in matrix inversion */
		BigError = -1.;
		level_error = -1;
		/* >>chng 06 aug 28, from numLevels_max to _local. */
		for( level=ipH1s; level < numlevels_local; level++ )
		{
			double abserror;
			abserror = fabs( error[level]);
			/* this will be the largest residual in the matrix inversion */
			if( abserror > BigError )
			{
				BigError = abserror;
				level_error = level;
			}
		}

		/* matrix inversion should be nearly as good as the accuracy of a double,
		 * but demand that it is better than epsilon for a float */
		if( BigError > FLT_EPSILON ) 
		{
			if( !conv.lgSearch )
				fprintf(ioQQQ," PROBLEM" );

			fprintf(ioQQQ,
				" iso_level zone %.2f - largest residual in iso=%li %s nelem=%li matrix inversion is %g "
				"level was %li Search?%c \n", 
				fnzone,
				ipISO,
				elementnames.chElementName[nelem],
				nelem , 
				BigError , 
				level_error,
				TorF(conv.lgSearch) );
		}
		
		// Force level balance to LTE
		if ( iso_ctrl.lgLTE_levels[ipISO] )
		{
			t_iso_sp* sp = &iso_sp[ipISO][nelem];
			double partfun=0.0;
			for ( level = 0; level < numlevels_local; level++ )
			{
				partfun += sp->st[level].Boltzmann()*sp->st[level].g();
			}
			double scale = dense.xIonDense[nelem][nelem-ipISO]/partfun;
			for ( level = 0; level < numlevels_local; level++ )
			{
				creation[level] = sp->st[level].Boltzmann()*sp->st[level].g()*scale;
			}
		}

		for( level = numlevels_local-1; level > 0; --level )
		{
			/* check for negative populations */
			if( creation[level] < 0. )
				lgNegPop = true;
		}

		if( lgNegPop && dense.lgSetIoniz[nelem] )
		{
			//  simulation can become unphysical if ionization is fixed.
			//  in this case, just put everything in ground.
			//  It's really the best we can do.
			for( level = 1; level < numlevels_local; ++level )
				creation[level] = 0.;
			creation[0] = dense.xIonDense[nelem][nelem-ipISO];
			// flip flag back
			lgNegPop = false;
		}

		/* put level populations into master array */
		for( level=0; level < numlevels_local; level++ )
		{
			ASSERT( creation[level] >= 0. );
			sp->st[level].Pop() = creation[level];

			if( sp->st[level].Pop() <= 0 && !conv.lgSearch )
			{
				fprintf(ioQQQ,
					"PROBLEM non-positive level pop for iso = %li, nelem = "
					"%li = %s, level=%li val=%.3e nTotalIoniz %li\n", 
					ipISO,
					nelem , 
					elementnames.chElementSym[nelem],
					level,
					sp->st[level].Pop() ,
					conv.nTotalIoniz);
			}
		}

		/* zero populations of unused levels. */
		for( level=numlevels_local; level < sp->numLevels_max; level++ )
		{
			sp->st[level].Pop() = 0.;
			/* >>chng 06 jul 25, no need to zero this out, fix limit to 3-body heating elsewhere. */
			/* sp->st[level].PopLTE = 0.; */
		}

		/* renormalize the populations to agree with ion solver */
		iso_renorm( nelem, ipISO, renorm );

		/* TotalPopExcited is sum of excited level pops */
		double TotalPopExcited = 0.;
		/* create sum of populations */
		for( level=1; level < numlevels_local; level++ )
			TotalPopExcited += sp->st[level].Pop();
		ASSERT( TotalPopExcited >= 0. );
		double TotalPop = TotalPopExcited + sp->st[0].Pop();

		/* option to force ionization */
		if( dense.lgSetIoniz[nelem] )
		{
			if( !fp_equal( TotalPop, (double)dense.xIonDense[nelem][nelem-ipISO] ) )
			{
				if( TotalPopExcited >= dense.xIonDense[nelem][nelem-ipISO] )
				{
					for( level=0; level < numlevels_local; level++ )
						sp->st[level].Pop() *=
							dense.xIonDense[nelem][nelem-ipISO] / TotalPop;
				}
				else
				{
					sp->st[0].Pop() = 
						MAX2( 1e-30 * dense.xIonDense[nelem][nelem-ipISO], 
						dense.xIonDense[nelem][nelem-ipISO] - TotalPopExcited );
				}
				sp->lgPopsRescaled = true;
			}	
			ASSERT( sp->st[0].Pop() >= 0. );
		}

		if( lgNegPop && !iso_ctrl.lgLTE_levels[ipISO] )
		{
			string species = makeChemical( nelem, nelem-ipISO );
			valarray<double> SaveC( get_ptr(Save_creation), Save_creation.size() );
			save.img_matrix.createImage( species, iteration, nzone,
							numlevels_local, SaveZ, SaveC, true );

			valarray<double> SavePops( get_ptr(creation), creation.size() );
			save.img_matrix.addImagePop_FITS( species, iteration, nzone,
							numlevels_local, SavePops, true );
		}
	}
	/* all solvers end up here */

	/* check on the sum of the populations */
	if( lgNegPop || dense.xIonDense[nelem][nelem-ipISO] < 0. )
	{
		fprintf( ioQQQ, 
			" DISASTER iso_level finds negative ion fraction for iso=%2ld nelem=%2ld "
			"%s using solver %s, IonFrac=%.3e simple=%.3e TE=%.3e ZONE=%4ld\n", 
			ipISO,
			nelem,
			elementnames.chElementSym[nelem],
			sp->chTypeAtomUsed,
			dense.xIonDense[nelem][nelem+1-ipISO]/SDIV(dense.xIonDense[nelem][nelem-ipISO]),
			sp->xIonSimple,
			phycon.te,
			nzone );

		fprintf( ioQQQ, " level pop are:\n" );
		for( i=0; i < numlevels_local; i++ )
		{
			fprintf( ioQQQ,PrintEfmt("%8.1e", sp->st[i].Pop() ));
			if( (i!=0) && !(i%10) ) fprintf( ioQQQ,"\n" );
		}
		fprintf( ioQQQ, "\n" );
		ContNegative();
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	for( ipHi=1; ipHi<numlevels_local; ++ipHi )
	{
		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			if( sp->trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
				continue;

			/* population of lower level, corrected for stimulated emission */
			sp->trans(ipHi,ipLo).Emis().PopOpc() = 
				sp->st[ipLo].Pop() - sp->st[ipHi].Pop()*
				sp->st[ipLo].g()/sp->st[ipHi].g();

			// don't allow masers from collapsed levels or in Case B only when MASER OFF has been specified
			if( iso_ctrl.lgNoMaser[ipISO][nelem] && ( N_(ipHi) > sp->n_HighestResolved_local || opac.lgCaseB ))
				sp->trans(ipHi,ipLo).Emis().PopOpc() = MAX2( 0., sp->trans(ipHi,ipLo).Emis().PopOpc() );
		}
	}

	// Zero PopOpc of inactive transitions.
	for( ipHi=numlevels_local; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( ipLo=0; ipLo<ipHi; ++ipLo )
		{
			sp->trans(ipHi,ipLo).Emis().PopOpc() = 0.; 
		}
	}
	return;
}

/** update multiplet opacities */
void iso_multiplet_opacities( void )
{
	for (long nelem = ipHYDROGEN; nelem < LIMELM; ++nelem)
	{
		if( ! dense.lgElmtOn[nelem] )
			continue;
		for ( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			if( (dense.IonHigh[nelem] >= nelem - ipISO) &&
				 (dense.IonLow[nelem] <= nelem - ipISO) )
			{
				iso_multiplet_opacities_one(ipISO, nelem);
			}
		}
	}
}

STATIC void iso_multiplet_opacities_one(
	const long int ipISO, const long int nelem)
{
	DEBUG_ENTRY( "iso_multiplet_opacities_one()" );
 	t_iso_sp* sp = &iso_sp[ipISO][nelem]; 
	// Allocate array for calculating n-n' multiplet opacities.
	multi_arr<double,2> MultOpacs;
	long nMax = sp->n_HighestResolved_max + sp->nCollapsed_max;
	MultOpacs.reserve( nMax+1 );
	for( long n=2; n <= nMax; ++n )
		MultOpacs.reserve( n, n+1 );
	MultOpacs.alloc();
	MultOpacs.zero();
	
	double rDopplerWidth = 1.0/GetDopplerWidth(dense.AtomicWeight[nelem]);
	
	// Sum n-n' multiplet opacities.	
	for( long ipHi=1; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( long ipLo=0; ipLo < ipHi; ++ipLo )
		{
			const TransitionProxy& tr = sp->trans(ipHi,ipLo);
			MultOpacs[ sp->st[ipHi].n() ][ sp->st[ipLo].n() ] += tr.Emis().PopOpc() *
				tr.Emis().opacity() * rDopplerWidth;
		}
	}

	// Now store n-n' multiplet opacities.
	for( long ipHi=1; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( long ipLo=0; ipLo < ipHi; ++ipLo )
		{
			const TransitionProxy& tr = sp->trans(ipHi,ipLo);
			tr.Emis().mult_opac() = MultOpacs[ sp->st[ipHi].n() ][ sp->st[ipLo].n() ];
		}
	}
}

void iso_set_ion_rates( long ipISO, long nelem)
{
	DEBUG_ENTRY( "iso_set_ion_rates()" );
 	t_iso_sp* sp = &iso_sp[ipISO][nelem]; 
	const long int numlevels_local = sp->numLevels_local;
	/* this is total ionization rate, s-1, of this species referenced to
	 * the total abundance */
	double TotalPop = 0.;
	ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = 0.;
	for( long level=0; level < numlevels_local; level++ )
	{
		/* sum of all ionization processes from this atom to ion, cm-3 s-1 now,
		 * but is divided by TotalPop below to become s-1 */
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] += 
			sp->st[level].Pop() * sp->fb[level].RateLevel2Cont;
		TotalPop += sp->st[level].Pop();
	}
	
	if( ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] > BIGDOUBLE )
	{
		fprintf( ioQQQ, "DISASTER RateIonizTot for Z=%li, ion %li is larger than BIGDOUBLE.  This is a big problem.",
					nelem+1, nelem-ipISO);
		cdEXIT(EXIT_FAILURE);
	}

	if (TotalPop <= SMALLFLOAT)
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] = sp->fb[0].RateLevel2Cont;
	else
		ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] /= TotalPop;

	if( ionbal.RateRecomIso[nelem][ipISO] > 0. )
		sp->xIonSimple = ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1]/ionbal.RateRecomIso[nelem][ipISO];
	else
		sp->xIonSimple = 0.;

	ASSERT( ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] >= 0. );

	if( ipISO == ipHE_LIKE && nelem==ipHELIUM && nzone>0 )
	{
		/* find fraction of He0 destructions due to photoionization of 2^3S */
		double ratio;
		double rateOutOf2TripS = sp->st[ipHe2s3S].Pop() * sp->fb[ipHe2s3S].RateLevel2Cont;
		if( rateOutOf2TripS > SMALLFLOAT )
		{
			ratio = rateOutOf2TripS /
				( rateOutOf2TripS + sp->st[ipHe1s1S].Pop()*ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] );
		}
		else
		{
			ratio = 0.;
		}
		if( ratio > he.frac_he0dest_23S )
		{
			/* remember zone where this happended and fraction, and frac due to photoionization */
			he.nzone = nzone;
			he.frac_he0dest_23S = ratio;
			rateOutOf2TripS = sp->st[ipHe2s3S].Pop() *sp->fb[ipHe2s3S].gamnc;
			if( rateOutOf2TripS > SMALLFLOAT )
			{
				he.frac_he0dest_23S_photo = rateOutOf2TripS /
					( rateOutOf2TripS + sp->st[ipHe1s1S].Pop()*ionbal.RateIoniz[nelem][nelem-ipISO][nelem-ipISO+1] );
			}
			else
			{
				he.frac_he0dest_23S_photo = 0.;
			}
		}
	}
}
