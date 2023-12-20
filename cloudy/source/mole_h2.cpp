/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_ContPoint set the ipCont struc element for the H2 molecule, called by ContCreatePointers */
/*H2_Accel radiative acceleration due to H2 */
/*H2_RadPress rad pressure due to h2 lines called in PresTotCurrent */
/*H2_InterEnergy internal energy of H2 called in PresTotCurrent */
/*H2_RT_diffuse do emission from H2 - called from RT_diffuse */
/*H2_itrzn - average number of H2 pop evaluations per zone */
/*H2_RTMake do RT for H2 - called from RT_line_all */
/*H2_RT_tau_inc increment optical depth for the H2 molecule, called from RT_tau_inc */
/*H2_LineZero initialize optical depths in H2, called from RT_tau_init */
/*H2_RT_tau_reset the large H2 molecule, called from RT_tau_reset */
/*H2_Colden maintain H2 column densities within X */
/*H2_LevelPops do level populations for H2, called by Hydrogenic */
/*H2_Level_low_matrix evaluate CO rotation cooling */
/*H2_cooling evaluate cooling and heating due to H2 molecule */
/*H2_X_coll_rate_evaluate find collisional rates within X */
/*cdH2_colden return column density in H2, negative -1 if cannot find state,
 * header is cddrive */
/*H2_DR choose next zone thickness based on H2 big molecule */
#include "cddefines.h" 
#include "cddrive.h" 
#include "atoms.h" 
#include "conv.h" 
#include "secondaries.h" 
#include "pressure.h" 
#include "trace.h" 
#include "hmi.h" 
#include "rt.h" 
#include "radius.h" 
#include "ipoint.h" 
#include "phycon.h" 
#include "thermal.h" 
#include "dense.h" 
#include "h2.h"
#include "mole.h"
#include "rfield.h"
#include "doppvel.h"
#include "lines_service.h"

/* turn this flag on to do minimal debug print of pops */
static const bool PRT_POPS = false;
/* this is limit to number of loops over H2 pops before giving up */
static const int LIM_H2_POP_LOOP = 10;
/* this is a typical dissociation cross section (cm2) for H2 + Hnu -> 2H + ke */
/* >>chng 05 may 11, had been 2.5e-19 */
static const realnum H2_DISS_ALLISON_DALGARNO = 6e-19_r;

static realnum collider_density[N_X_COLLIDER];
static realnum collider_density_total_not_H2;

void diatomics::H2_X_sink_and_source( void )
{
	DEBUG_ENTRY( "diatomics::H2_X_sink_and_source()" );

	/* this is total density of all colliders, is only used for collisional dissociation
	 * rates for H2 are not included here, will be added separately*/
	collider_density_total_not_H2 = collider_density[0] + 
		collider_density[1] + collider_density[4] + 
		dense.eden;

	for( long ipHi=0; ipHi<nLevels_per_elec[0]; ++ipHi )
	{
		H2_X_source[ipHi] = 0.;
		H2_X_sink[ipHi] = 0.;
	}

	double source_so_far = 0.;
	double source_so_far_s = 0.;
	double sink_so_far = 0.;
	double pop_tot = 0.;
	double sink_so_far_s = spon_diss_tot * H2_den_s;
	double pop_tot_s = 0.;

	for( long ipHi=0; ipHi < nLevels_per_elec[0]; ++ipHi )
	{
		/* array of energy sorted indices within X */
		long iVibHi = ipVib_H2_energy_sort[ipHi];
		long iRotHi = ipRot_H2_energy_sort[ipHi];

		/* count formation from grains and H- as a collisional formation process */
		/* cm-3 s-1, evaluated in mole_H2_form */
		H2_X_source[ipHi] += H2_X_formation[iVibHi][iRotHi];
		
		/*>>chng 05 sep 18, GS, H2 + e = H- + H*, H2_X_Hmin_back has units s-1 */
		H2_X_sink[ipHi] += H2_X_Hmin_back[iVibHi][iRotHi];
		
		/* this represents collisional dissociation into continuum of X,
		 * rates are just guesses */
		H2_X_sink[ipHi] += collider_density_total_not_H2 *
			H2_coll_dissoc_rate_coef[iVibHi][iRotHi] * lgColl_deexec_Calc;
		
		/*>>chng 05 jul 20, GS, collisional dissociation with H2g and H2s are added here*/
		H2_X_sink[ipHi] +=  hmi.H2_total*
			H2_coll_dissoc_rate_coef_H2[iVibHi][iRotHi] * lgColl_deexec_Calc;

		/* rate (s-1) out of this level */
		if( mole_global.lgStancil )
		{
			H2_X_sink[ipHi] += Cont_Dissoc_Rate[0][iVibHi][iRotHi];
		}
		else
			H2_X_sink[ipHi] += rfield.flux_accum[H2_ipPhoto[iVibHi][iRotHi]-1]*H2_DISS_ALLISON_DALGARNO;

		if ( states[ipHi].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
		{
			source_so_far_s += H2_X_source[ipHi];
			sink_so_far_s += H2_X_sink[ipHi] * states[ipHi].Pop();
			pop_tot_s += states[ipHi].Pop();
		}
		else
		{
			source_so_far += H2_X_source[ipHi];
			sink_so_far += H2_X_sink[ipHi] * states[ipHi].Pop();
			pop_tot += states[ipHi].Pop();
		}			
	}

	// cm-3 s-1
	double sink_tot = mole.sink_rate_tot(sp) * pop_tot;
	// cm-3 s-1
	double sink_left = sink_tot - sink_so_far;
	// divide by population in X to get units s-1
	ASSERT( pop_tot > 1e-10 * (*dense_total) );
	sink_left /= pop_tot;
	if( sink_left >= 0. )
	{
		for( long ipHi=0; ipHi < nLevels_per_elec[0]; ++ipHi )
		{
			if( states[ipHi].energy().WN() <= ENERGY_H2_STAR || !hmi.lgLeiden_Keep_ipMH2s )
			{
				H2_X_sink[ipHi] += sink_left;
			}
		}
	}
	
	// cm-3 s-1
	fixit("kill the second term (sp_star) when H2* is killed in chemistry");
	double sink_tot_s =  mole.sink_rate_tot(sp_star) * pop_tot_s;
	// cm-3 s-1
	double sink_left_s = sink_tot_s - sink_so_far_s;
	// divide by population in X to get units s-1
	if( pop_tot_s > 1e-30 * (*dense_total) )
		sink_left_s /= pop_tot_s;
	else
		sink_left_s = 0.;
	if( sink_left_s >= 0. )
	{
		for( long ipHi=0; ipHi < nLevels_per_elec[0]; ++ipHi )
		{
			if( states[ipHi].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
				H2_X_sink[ipHi] += sink_left_s;
		}
	}

	fixit("kill the second term (sp_star) when H2* is killed in chemistry");
	double source_tot = mole.source_rate_tot(sp);
	double source_left = source_tot - source_so_far;
	double source_tot_s = mole.source_rate_tot(sp_star);
	double source_left_s = source_tot_s - source_so_far_s;
	if( source_left+source_left_s >= 0. )
	{
		double rpop_lte = 1.;
		double rpop_lte_s = 0.;
		if (hmi.lgLeiden_Keep_ipMH2s)
		{
			double pop_lte = 0.;
			double pop_lte_s = 0.;
			for( long ipHi=0; ipHi < nLevels_per_elec[0]; ++ipHi )
			{
				long iElec = states[ipHi].n();
				long iVib = states[ipHi].v();
				long iRot = states[ipHi].J();
				if( states[ipHi].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
					pop_lte_s += H2_populations_LTE[iElec][iVib][iRot];
				else
					pop_lte += H2_populations_LTE[iElec][iVib][iRot];
			}
			rpop_lte = 1./SDIV(pop_lte);
			rpop_lte_s = 1./SDIV(pop_lte_s);
		}
		for( long ipHi=0; ipHi < nLevels_per_elec[0]; ++ipHi )
		{
			long iElec = states[ipHi].n();
			long iVib = states[ipHi].v();
			long iRot = states[ipHi].J();
			if( states[ipHi].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
				H2_X_source[ipHi] += source_left_s * H2_populations_LTE[iElec][iVib][iRot]*rpop_lte_s;
			else
				H2_X_source[ipHi] += source_left * H2_populations_LTE[iElec][iVib][iRot]*rpop_lte;
		}
	}

	return;
}

/*H2_X_coll_rate_evaluate find collisional rates within X - 
 * this is one time upon entry into H2_LevelPops */
void diatomics::H2_X_coll_rate_evaluate( void )
{
	DEBUG_ENTRY( "diatomics::H2_X_coll_rate_evaluate()" );

	/* set collider density 
	 * the colliders are:
	 * [0] = H
	 * [1], [5] = He (old and new cs data)
	 * [2] = H2 ortho
	 * [3] = H2 para
	 * [4] = H+ + H3+ */
	/* atomic hydrogen */
	collider_density[0] = dense.xIonDense[ipHYDROGEN][0];
	/* all ortho h2 */
	/* He - H2 */
	collider_density[1] = dense.xIonDense[ipHELIUM][0];
	/* H2 - H2(ortho) */
	collider_density[2] = h2.ortho_density_f;
	/* all para H2 */
	collider_density[3] = h2.para_density_f;
	/* protons - ionized hydrogen */
	collider_density[4] = dense.xIonDense[ipHYDROGEN][1];
	/* H3+ - assume that H3+ has same rates as proton */
	collider_density[4] += (realnum)findspecieslocal("H3+")->den;

	ASSERT( fp_equal( hmi.H2_total_f ,collider_density[2]+collider_density[3]) );

	if( nTRACE >= n_trace_full )
	{
		fprintf(ioQQQ," Collider densities are:");
		for( long nColl=0; nColl<N_X_COLLIDER; ++nColl )
		{
			fprintf(ioQQQ,"\t%.3e", collider_density[nColl]);
		}
		fprintf(ioQQQ,"\n");
	}

	H2_X_coll_rate.zero();

	for( long ipHi=0; ipHi < nLevels_per_elec[0]; ++ipHi )
	{
		if( lgColl_deexec_Calc )
		{
			/* excitation within X due to thermal particles */
			for( long ipLo=0; ipLo<ipHi; ++ipLo )
			{
				/* collisional interactions with upper levels within X */
				double colldown = 0.;
				mr3ci CollRate = CollRateCoeff.begin(ipHi, ipLo);
				for( long nColl=0; nColl<N_X_COLLIDER; ++nColl )
				{
					/* downward collision rate, units s-1 */
					colldown += CollRate[nColl]*collider_density[nColl];
					ASSERT( CollRate[nColl]*collider_density[nColl] >= 0. );
				}
				/* rate in from upper level, units cm-3 s-1 */
				H2_X_coll_rate[ipHi][ipLo] += colldown;
			}/* end loop over ipLo */
		}
	} 

	return;
}

/*H2_itrzn - average number of H2 pop evaluations per zone */
double diatomics::H2_itrzn( void )
{
	if( lgEnabled && nH2_zone>0 )
	{
		return( (double)nH2_pops / (double)nH2_zone );
	}
	else
	{
		return 0.;
	}
}

/* set the ipCont struc element for the H2 molecule, called by ContCreatePointers */
void diatomics::H2_ContPoint( void )
{
	if( !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_ContPoint()" );

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		ASSERT( (*tr).Emis().Aul() > 0. );
		(*tr).ipCont() = ipLineEnergy( (*tr).EnergyRyd(), label.c_str(), 0 );
		(*tr).Emis().ipFine() = ipFineCont( (*tr).EnergyRyd());
	}
	return;
}

/* ===================================================================== */
/* radiative acceleration due to H2 called in rt_line_driving */
double diatomics::H2_Accel(void)
{
	/* >>chng 05 jan 26, pops now set to LTE for small abundance case, so do this */
	if( !lgEnabled /*|| !nCall_this_zone*/ )
		return 0.;

	DEBUG_ENTRY( "H2_Accel()" );

	/* this routine computes the line driven radiative acceleration */

	double drive = 0.;
	
	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		ASSERT( (*tr).ipCont() > 0 );
		drive += (*tr).Emis().pump() * (*tr).Emis().PopOpc() * (*tr).EnergyErg();
	}

	return drive;
}

/* ===================================================================== */
/* rad pressure due to H2 lines called in PresTotCurrent */
double diatomics::H2_RadPress(void)
{
	/* will be used to check on size of opacity, was capped at this value */
	realnum smallfloat=SMALLFLOAT*10.f;

	/* radiation pressure sum is expensive - do not evaluate if we did not
	 * bother evaluating large molecule */
	if( !lgEnabled || !nCall_this_zone )
		return 0.;

	DEBUG_ENTRY( "H2_RadPress()" );

	realnum doppler_width =	GetDopplerWidth( mass_amu );
	double press = 0.;
	
	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		ASSERT( (*tr).ipCont() > 0 );
		if( (*(*tr).Hi()).Pop() > smallfloat && (*tr).Emis().PopOpc() > smallfloat )
		{
			press +=  PressureRadiationLine( *tr, doppler_width );
		}
	}

	if(nTRACE >= n_trace_full) 
		fprintf(ioQQQ,
		"  H2_RadPress returns, radiation pressure is %.2e\n", 
		press );
	return press;
}

#if 0
/* ===================== */
/* internal energy of H2 */
double diatomics::H2_InterEnergy(void)
{
	/* >>chng 05 jan 26, pops now set to LTE for small abundance case, so do this */
	if( !lgEnabled /*|| !nCall_this_zone*/ )
		return 0.;

	DEBUG_ENTRY( "H2_InterEnergy()" );

	double energy = 0.;
	for( qList::iterator st = trans.states.begin(); st != trans.states.end(); ++st )
		energy += st->Pop() * st->energy();
	
	return energy;
}
#endif

/*H2_RT_diffuse do emission from H2 - called from RT_diffuse */
void diatomics::H2_RT_diffuse(void)
{
	if( !lgEnabled || !nCall_this_zone )
		return;

	DEBUG_ENTRY( "H2_RT_diffuse()" );

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		qList::iterator Hi = (*tr).Hi();
		if( (*Hi).n() > 0 )
			continue;	
		(*tr).outline_resonance();
	}

	return;
}

/* RT for H2 lines */
void diatomics::H2_RTMake( linefunc line_one )
{
	if( !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_RTMake()" );

	realnum doppler_width =	GetDopplerWidth( mass_amu );
	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		/* >>chng 03 jun 18, added 4th parameter in call to this routine - says to not
		 * include self-shielding of line across this zone.  This introduces a dr dependent
		 * variation in the line pumping rate, which made H2 abundance fluctuate due to
		 * Solomon process having slight dr-caused mole. */
		line_one( *tr, false, 0.f, doppler_width ); 
	}

	return;
}

/* increment optical depth for the H2 molecule, called from RT_tau_inc which is called  by cloudy,
 * one time per zone */
void diatomics::H2_RT_tau_inc(void)
{
	/* >>chng 05 jan 26, now use LTE populations for small H2 abundance case, since electronic
	 * lines become self-shielding surprisingly quickly */
	if( !lgEnabled /*|| !nCall_this_zone*/ )
		return;

	DEBUG_ENTRY( "H2_RT_tau_inc()" );

	/* remember largest and smallest chemistry renormalization factor -
	 * if both networks are parallel will be unity,
	 * but only do this after we have stable solution */
	if( nzone > 0 && nCall_this_iteration>2 )
	{
		renorm_max = MAX2( H2_renorm_chemistry , renorm_max );
		renorm_min = MIN2( H2_renorm_chemistry , renorm_min );
	}

	realnum doppler_width =	GetDopplerWidth( mass_amu );
	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		ASSERT( (*tr).ipCont() > 0 );
		RT_line_one_tauinc( *tr,-9, -9, -9, -9, doppler_width );
	}

	return;
}


/* initialize optical depths in H2, called from RT_tau_init */
void diatomics::H2_LineZero( void )
{
	if( !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_LineZero()" );

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		(*tr).Zero();
	}

	return;
}

/* the large H2 molecule, called from RT_tau_reset */
void diatomics::H2_RT_tau_reset( void )
{
	if( !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_RT_tau_reset()" );

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		RT_line_one_tau_reset( *tr );
	}

	return;
}

/*H2_Level_low_matrix evaluate lower populations within X */
void diatomics::H2_Level_low_matrix(
	/* total abundance within matrix */
	realnum abundance )
{
	bool lgDoAs;
	int nNegPop;
	bool lgDeBug,
		lgZeroPop;
	double rot_cooling , dCoolDT;

	DEBUG_ENTRY( "H2_Level_low_matrix()" );

	/* option to not use the matrix */
	if( nXLevelsMatrix <= 1 )
	{
		return;
	}

	/* will need to allocate space for these but only on first call */
	if( lgFirst )
	{
		/* check that not more levels than there are in X */
		if( nXLevelsMatrix > nLevels_per_elec[0] )
		{
			/* number is greater than number of levels within X */
			fprintf( ioQQQ, 
				" The total number of levels used in the matrix solver must be <= %li, the number of levels within X.\n Sorry.\n",
				nLevels_per_elec[0]);
			cdEXIT(EXIT_FAILURE);
		}
		/* will never do this again */
		lgFirst = false;
		/* remember how much space we allocated in case ever called with more needed */
		/* >>chng 05 jan 19, allocate max number of levels
		ndim_allocated = nXLevelsMatrix;*/
		ndim_allocated = nLevels_per_elec[0];
		/* allocate the 1D arrays*/
		excit.resize( ndim_allocated );
		stat_levn.resize( ndim_allocated );
		pops.resize( ndim_allocated );
		create.resize( ndim_allocated );
		destroy.resize( ndim_allocated );
		depart.resize( ndim_allocated );
		/* create space for the 2D arrays */
		AulPump.alloc(ndim_allocated, ndim_allocated);
		CollRate_levn.alloc(ndim_allocated, ndim_allocated);
		AulDest.alloc(ndim_allocated, ndim_allocated);
		AulEscp.alloc(ndim_allocated, ndim_allocated);

		/* the statistical weights of the levels
		 * and excitation potentials of each level relative to ground */
		for( long j=0; j < ndim_allocated; j++ )
		{
			/* statistical weights for each level */
			stat_levn[j] = states[j].g();
			/* excitation energy of each level relative to ground, in K */
			excit[j] = states[j].energy().K();
		}
	}
	/* end allocate space and creating constant terms */

	/* this is test for call with too many rotation levels to handle - 
	 * logic needs for largest model atom to be called first */
	if( nXLevelsMatrix > ndim_allocated )
	{
		fprintf(ioQQQ," H2_Level_low_matrix has been called with the number of rotor levels greater than space allocated.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* all elements are used, and must be set to zero */
	for( long i=0; i < nXLevelsMatrix; i++ )
	{
		pops[i] = 0.;
		depart[i] = 0;
	}

	/* do we need to reevaluate radiative quantities?  only do this one time per zone */
	if( nzone!=nzoneAsEval || iteration!=iterationAsEval || nXLevelsMatrix!=levelAsEval)
	{
		lgDoAs = true;
		nzoneAsEval = nzone;
		iterationAsEval = iteration;
		levelAsEval = nXLevelsMatrix;
		ASSERT( levelAsEval <= ndim_allocated );
	}
	else
	{
		lgDoAs = false;
	}

	/* all elements are used, and must be set to zero */
	if( lgDoAs )
	{
		for( long i=0; i < nXLevelsMatrix; i++ )
		{
			pops[i] = 0.;
			depart[i] = 0;
			for( long j=0; j < nXLevelsMatrix; j++ )
			{
				AulEscp[j][i] = 0.;
				AulDest[j][i] = 0.;
				AulPump[j][i] = 0.;
				CollRate_levn[j][i] = 0.;
			}
		}
	}

	/* find all radiative interactions within matrix, and between
	 * matrix and upper X and excited electronic states */
	for( long ilo=0; ilo < nXLevelsMatrix; ilo++ )
	{
		long iRot = ipRot_H2_energy_sort[ilo];
		long iVib = ipVib_H2_energy_sort[ilo];

		/* H2_X_sink[ilo] includes all processes that destroy H2 in one step, 
		 * these include cosmic ray ionization and dissociation, photodissociation,
		 * BUT NOT THE SOLOMON process, which, directly, only goes to excited
		 * electronic states */
		destroy[ilo] = H2_X_sink[ilo];

		/* rates H2 is created from grains and H- units cm-3 s-1, evaluated in mole_H2_form */
		create[ilo] = H2_X_source[ilo];

		/* this loop does radiative decays from upper states inside matrix, 
		 * and upward pumps within matrix region into this lower level */
		if( lgDoAs )
		{
			for( long ihi=ilo+1; ihi<nXLevelsMatrix; ++ihi )
			{
				ASSERT( states[ihi].energy().WN() <= states[nXLevelsMatrix-1].energy().WN() );
				/* general case - but line may not actually exist */
				if( lgH2_radiative[ihi][ilo] )
				{
					const TransitionList::iterator&tr = trans.begin()+ ipTransitionSort[ihi][ilo] ;
					ASSERT( (*tr).ipCont() > 0 );

					/* NB - the destruction probability is included in 
					 * the total and the destruction is set to zero
					 * since we want to only count one ots rate, in 
					 * main calling routine, and do not want matrix 
					 * solver below to include it */
					AulEscp[ihi][ilo] = (*tr).Emis().Aul()*(
						(*tr).Emis().Pesc_total() );
					AulDest[ihi][ilo] = (*tr).Emis().Aul()*(*tr).Emis().Pdest();
					AulPump[ilo][ihi] = (*tr).Emis().pump();
					AulPump[ihi][ilo] = AulPump[ilo][ihi]*stat_levn[ilo]/stat_levn[ihi];
				}
			}
		}

		double rateout = 0.;
		double ratein = 0.;
		/* now do all levels within X, which are above nXLevelsMatrix,
		 * the highest level inside the matrix */
		for( long ihi=nXLevelsMatrix; ihi<nLevels_per_elec[0]; ++ihi )
		{
			if( lgH2_radiative[ihi][ilo] )
			{
				const TransitionList::iterator&tr = trans.begin()+ ipTransitionSort[ihi][ilo] ;
				ASSERT( (*tr).ipCont() > 0 );
				
				/* these will enter as net creation terms in creation vector, with
				 * units cm-3 s-1
				 * radiative transitions from above the matrix within X */
				ratein +=
					(*(*tr).Hi()).Pop() * (
					(*tr).Emis().Aul()*( (*tr).Emis().Ploss() ) +
					(*tr).Emis().pump() * (*(*tr).Lo()).g() / (*(*tr).Hi()).g() );
				/* rate out has units s-1 - destroys current lower level */
				rateout += (*tr).Emis().pump();
			}
		}

		/* all states above the matrix but within X */
		create[ilo] += ratein;

		/* rates out of matrix into levels in X but above matrix */
		destroy[ilo] += rateout;

		/* Solomon process, this sum dos all pump and decays from all electronic excited states */
		/* radiative rates [cm-3 s-1] from electronic excited states into X only vibration and rot */
		create[ilo] += H2_X_rate_from_elec_excited[iVib][iRot];

		/* radiative & cosmic ray rates [s-1] to electronic excited states from X only vibration and rot */
		destroy[ilo] += H2_X_rate_to_elec_excited[iVib][iRot];
	}

	/* this flag set with atom H2 trace matrix */
	if( nTRACE >= n_trace_matrix )
		lgDeBug = true;
	else
		lgDeBug = false;

	/* now evaluate the rates for all transitions within matrix */
	for( long ilo=0; ilo < nXLevelsMatrix; ilo++ )
	{
		double EloK = states[ilo].energy().K();

		if(lgDeBug)fprintf(ioQQQ,"DEBUG H2_Level_low_matrix, ilo=%li",ilo);
		for( long ihi=ilo+1; ihi < nXLevelsMatrix; ihi++ )
		{
			/* >>chng 05 may 31, replace with simple expresion */
			CollRate_levn[ihi][ilo] = H2_X_coll_rate[ihi][ilo];

			if(lgDeBug)fprintf(ioQQQ,"\t%.1e",CollRate_levn[ihi][ilo]);

			/* now get upward excitation rate - units s-1 */
			// use dsexp() explicitly rather than divide Boltzmann() factors
			// to avoid problems with underflow at low Te (see ticket #284)
			CollRate_levn[ilo][ihi] = CollRate_levn[ihi][ilo]*
				dsexp( (states[ihi].energy().K() - EloK) / phycon.te )*
				states[ihi].g() / states[ilo].g();
		}
		if(lgDeBug)fprintf(ioQQQ,"\n");

		/* now do all collisions for levels within X, which are above nXLevelsMatrix,
		 * the highest level inside the matrix */
		for( long ihi=nXLevelsMatrix; ihi<nLevels_per_elec[0]; ++ihi )
		{
			/* first do downward deexcitation rate */
			/* >>chng 04 sep 14, do all levels */
			/* >>chng 05 may 31, use summed rate */
			double ratein = H2_X_coll_rate[ihi][ilo];
			if(lgDeBug)fprintf(ioQQQ,"\t%.1e",ratein);

			/* now get upward excitation rate */
			// use dsexp() explicitly rather than divide Boltzmann() factors
			// to avoid problems with underflow at low Te (see ticket #284)
			double rateout = ratein *
				dsexp( (states[ihi].energy().K() - EloK) / phycon.te )*
				states[ihi].g()/states[ilo].g();

			/* these are general entries and exits going into vector */
			create[ilo] += ratein * states[ihi].Pop();
			destroy[ilo] += rateout;
		}
		if(lgDeBug)fprintf(ioQQQ,"\n");
	}

	/* H2 grain interactions */
	{
		for( long ihi=2; ihi < nXLevelsMatrix; ihi++ )
		{
			long iVibHi = ipVib_H2_energy_sort[ihi];
			long iRotHi = ipRot_H2_energy_sort[ihi];

			/* collisions with grains goes to either J=1 or J=0 depending on 
			* spin of upper level - this conserves op ratio - following
			* var is 1 if ortho, 0 if para, so this conserves op ratio
			* units are s-1 */
			CollRate_levn[ihi][H2_lgOrtho[0][iVibHi][iRotHi]] += rate_grain_op_conserve;
		}

		/* H2 ortho - para conversion on grain surface,
		 * rate (s-1) all v,J levels go to 0 or 1 */
		CollRate_levn[1][0] += 
			(realnum)(rate_grain_J1_to_J0);
	}

	/* now all levels in X above the matrix */
	for( long ihi=nXLevelsMatrix; ihi<nLevels_per_elec[0]; ++ihi )
	{
		long iVibHi = ipVib_H2_energy_sort[ihi];
		long iRotHi = ipRot_H2_energy_sort[ihi];

		/* these collisions all go into 0 or 1 depending on whether upper level was ortho or para 
		 * units are cm-3 s-1 - rate new molecules appear in matrix */
		create[H2_lgOrtho[0][iVibHi][iRotHi]] += states[ihi].Pop() * rate_grain_op_conserve;
	}

	/* debug print individual contributors to matrix elements */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC || lgDeBug)
		{
			fprintf(ioQQQ,"DEBUG H2 matexcit");
			for(long ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for(long ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",excit[ihi] );
			}
			fprintf(ioQQQ,"\n");
			for(long ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",stat_levn[ihi] );
			}
			fprintf(ioQQQ,"\n");

			fprintf(ioQQQ,"AulEscp[n][]\\[][n] = Aul*Pesc\n");
			for(long ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for(long ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"%li", ihi);
				for(long ilo=0; ilo<nXLevelsMatrix; ++ilo )
				{
					fprintf(ioQQQ,"\t%.2e",AulEscp[ilo][ihi] );
				}
				fprintf(ioQQQ,"\n");
			}

			fprintf(ioQQQ,"AulPump [n][]\\[][n]\n");
			for(long ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for(long ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"%li", ihi);
				for(long ilo=0; ilo<nXLevelsMatrix; ++ilo )
				{
					fprintf(ioQQQ,"\t%.2e",AulPump[ihi][ilo] );
				}
				fprintf(ioQQQ,"\n");
			}

			fprintf(ioQQQ,"CollRate_levn [n][]\\[][n]\n");
			for( long ilo=0; ilo<nXLevelsMatrix; ++ilo )
			{
				fprintf(ioQQQ,"\t%li",ilo );
			}
			fprintf(ioQQQ,"\n");
			for( long ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"%li", ihi);
				for( long ilo=0; ilo<nXLevelsMatrix; ++ilo )
				{
					fprintf(ioQQQ,"\t%.2e",CollRate_levn[ihi][ilo] );
				}
				fprintf(ioQQQ,"\n");
			}
			fprintf(ioQQQ,"SOURCE");
			for( long ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",create[ihi]);
			}
			fprintf(ioQQQ,"\nSINK");
			for( long ihi=0; ihi<nXLevelsMatrix;++ihi)
			{
				fprintf(ioQQQ,"\t%.2e",destroy[ihi]);
			}
			fprintf(ioQQQ,"\n");
		}
	}

	static Atom_LevelN atom_levelN;
	atom_levelN(
		nXLevelsMatrix,
		abundance,
		stat_levn,
		excit,
		'K',
		pops,
		depart,
		/* net transition rate, A * escape prob, s-1, indices are [upper][lower] */
		AulEscp, 
		AulDest,
		AulPump,
		CollRate_levn,
		create,
		destroy,
		&rot_cooling,
		&dCoolDT,
		" H2 ",
		lgPrtMatrix,
		lgImgMatrix,
		/* nNegPop positive if negative pops occurred, negative if too cold */
		&nNegPop,
		&lgZeroPop,
		lgDeBug ); /* option to print stuff - set to true for debug printout */

	if( nNegPop > 0 )
	{
		// Recovery procedure -- assume that sources have overshot
		fprintf(ioQQQ," H2_Level_low_matrix called atom_levelN which returned negative populations.\n");
		//ConvFail( "pops" , "H2" );
		double totpop = 0.;
		for( long i=0; i< nXLevelsMatrix; ++i )
		{
			if (pops[i] < 0.0)
				pops[i] = 0.0;
			totpop += pops[i];
		}
		totpop = abundance/totpop;
		for( long i=0; i< nXLevelsMatrix; ++i )
		{		
			pops[i] *= totpop;
		}
	}

	for( long i=0; i< nXLevelsMatrix; ++i )
	{
		states[i].Pop() = pops[i];
		fixit("need to also store DepartCoef above X");
		states[i].DepartCoef() = depart[i];
	}

	if( 0 && nTRACE >= n_trace_full) 
	{
		/* print pops that came out of matrix */
		fprintf(ioQQQ,"\n DEBUG H2_Level_lowJ dense_total: %.3e matrix rel pops\n", *dense_total);
		fprintf(ioQQQ,"v\tJ\tpop\n");
 		for( long i=0; i<nXLevelsMatrix; ++i )
		{
			long iRot = ipRot_H2_energy_sort[i];
			long iVib = ipVib_H2_energy_sort[i];
			fprintf(ioQQQ,"%3li\t%3li\t%.3e\t%.3e\t%.3e\n",
					  iVib , iRot , states[i].Pop(), create[i] , destroy[i]);
		}
	}
	
	/* nNegPop positive if negative pops occurred, negative if too cold */
	return;
}
/* do level populations for H2, called by Hydrogenic after ionization and H chemistry
 * has been recomputed */
void diatomics::H2_LevelPops( bool &lgPopsConverged, double &old_val, double &new_val )
{
	DEBUG_ENTRY( "H2_LevelPops()" );
	
	string convLabel;
	if( this==&h2 )
		convLabel = "H2_LOOPS";
	else if( this==&hd )
		convLabel = "HD_LOOPS";
	else
		TotalInsanity();
	static ConvergenceCounter cctr_diatom_l =conv.register_(convLabel);

	/* H2 not on, so space not allocated and return */
	if( !lgEnabled )
	{
		// need to do this even if not doing big model 
		CalcPhotoionizationRate();
		return;
	}

	double old_solomon_rate=-1.;
	long int n_pop_oscil = 0;
	int kase=0;
	bool lgConv_h2_soln,
		lgPopsConv_total,
		lgPopsConv_relative,
		lgHeatConv,
		lgSolomonConv,
		lgOrthoParaRatioConv;
	double quant_old=-1.,
		quant_new=-1.;

	bool lgH2_pops_oscil=false,
		lgH2_pops_ever_oscil=false;

	/* keep track of changes in population */
	double PopChgMax_relative=0. , PopChgMaxOld_relative=0., PopChgMax_total=0., PopChgMaxOld_total=0.;
	long int iRotMaxChng_relative , iVibMaxChng_relative,
		iRotMaxChng_total , iVibMaxChng_total,
		nXLevelsMatrix_save;
	double popold_relative , popnew_relative , popold_total , popnew_total;
	/* reason not converged */
	char chReason[100];

	/* these are convergence criteria - will be increased during search phase */
	double converge_pops_relative=conv.IonizErrorAllowed/3.,
		converge_pops_total=1e-3, 
		converge_ortho_para=1e-2;

	double dens_rel_to_lim_react = mole.species[sp->index].xFracLim;

	if(nTRACE >= n_trace_full ) 
	{
		fprintf(ioQQQ,
			"\n***************H2_LevelPops %s call %li this iteration, zone is %.2f, H2/H:%.e Te:%e ne:%e\n", 
			label.c_str(),
			nCall_this_iteration,
			fnzone,
			dens_rel_to_lim_react,
			phycon.te,
			dense.eden
			);
	}
	else if( nTRACE >= n_trace_final )
	{
		static long int nzone_prt=-1;
		if( nzone!=nzone_prt )
		{
			nzone_prt = nzone;
			fprintf(ioQQQ,"DEBUG zone %li species %s rel_to_lim:%.3e Te:%.3e *ne:%.3e n(%s):%.3e\n",
				nzone,
				label.c_str(),
				dens_rel_to_lim_react,
				phycon.te,
				dense.eden,
				label.c_str(),
				*dense_total );
		}
	}
	
	CalcPhotoionizationRate();

	/* evaluate Boltzmann factors and LTE unit population - for trivial abundances
	 * LTE populations are used in place of full solution */
	mole_H2_LTE();

	/* zero out populations and cooling, and return, if H2 fraction is small
	 * but, if H2 has ever been done, redo irregardless of abundance -
	 * if large H2 is ever evaluated then mole.H2_to_H_limit is ignored */
	if( (!lgEvaluated && dens_rel_to_lim_react < H2_to_H_limit )
		|| *dense_total < 1e-20 )
	{
		/* will not compute the molecule */
		if( nTRACE >= n_trace_full ) 
			fprintf(ioQQQ,
			"  H2_LevelPops %s pops too small, not computing, set to LTE and return, H2/H is %.2e and H2_to_H_limit is %.2e.\n",
			label.c_str(),
			dens_rel_to_lim_react,
			H2_to_H_limit);
		H2_zero_pops_too_low();
		fixit("set lgEvaluated = false here?");
		/* end of zero abundance branch */
		return;
	}

	/* check whether we need to update the H2_Boltzmann factors, LTE level populations,
	 * and partition function.  LTE level pops normalized by partition function,
	 * so sum of pops is unity */

	/* say that H2 has been computed, ignore previous limit to abund
	 * in future - this is to prevent oscillations as model is engaged */
	lgEvaluated = true;
	/* end loop setting H2_Boltzmann factors, partition function, and LTE populations */

	/* >>chng 05 jun 21,
	 * during search phase we want to use full matrix - save number of levels so that
	 * we can restore it */
	nXLevelsMatrix_save = nXLevelsMatrix;
	fixit("this does not appear to be necessary and may be counterproductive");
	if( conv.lgSearch )
	{
		nXLevelsMatrix = nLevels_per_elec[0];
	}

	/* 05 oct 27, had only reevaluated collision rates when 5% change in temperature
	 * caused temp failures in large G0 sims - 
	 * do not check whether we need to update the collision rates but
	 * reevaluate them always  
	 * >>chng 05 nov 04, above caused a 25% increase in the exec time for constant-T sims
	 * in test suite- original code had reevaluated if > 0.05 change in T - was too much
	 * change to 10x smaller, change > 0.005 */
	if( !fp_equal(phycon.te,TeUsedColl) )
	{
		H2_CollidRateEvalAll();
		TeUsedColl = phycon.te;
	}

	/* set the populations when this is the first call to this routine on 
	 * current iteration- will use LTE populations - populations were set by
	 * call to 	mole_H2_LTE before above block */
	if( nCall_this_iteration==0 || lgLTE )
	{
		/* very first call so use LTE populations */
		if(nTRACE >= n_trace_full ) 
			fprintf(ioQQQ,"%s 1st call - using LTE level pops\n", label.c_str() );

		H2_den_s = 0.;
		H2_den_g = 0.;
		for( qList::iterator st = states.begin(); st != states.end(); ++st )
		{
			long iElec = (*st).n();
			long iVib = (*st).v();
			long iRot = (*st).J();
			/* LTE populations are for unit H2 density, so need to multiply
			 * by total H2 density */
			double pop = H2_populations_LTE[iElec][iVib][iRot] * (*dense_total);
			H2_old_populations[iElec][iVib][iRot] = pop;
			(*st).Pop() = pop;
			/* find current population in H2s and H2g */
			if( (*st).energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
			{
				H2_den_s += pop;
			}
			else
			{
				H2_den_g += pop;
			}
		}

		/* first guess at ortho and para densities */
		ortho_density = 0.75* (*dense_total);
		para_density = 0.25* (*dense_total);
		{
			ortho_density_f = (realnum)ortho_density;
			para_density_f = (realnum)para_density;
		}
		ortho_para_current = ortho_density / SDIV( para_density );
		ortho_para_older = ortho_para_current;
		ortho_para_old = ortho_para_current;
		/* this is the fraction of the H2 pops that are within the levels done with a matrix */
		frac_matrix = 1.;
	}

	// make some population sums, normalize total to value handed from chemistry
	{
		pops_per_vib.zero();
		fill_n( pops_per_elec, N_ELEC, 0. );
		double pop_total = 0.;
		for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
		{
			long iElec = (*st).n();
			long iVib = (*st).v();
			
			pop_total += (*st).Pop();	
			pops_per_elec[iElec] += (*st).Pop();
			pops_per_vib[iElec][iVib] += (*st).Pop();
		}
		ASSERT( pops_per_elec[0]>SMALLFLOAT );
		// Now renorm the old populations to the correct current H2 density. 
		H2_renorm_chemistry = *dense_total/ SDIV(pop_total);
	}

	if(nTRACE >= n_trace_full) 
		fprintf(ioQQQ,
			"%s H2_renorm_chemistry is %.4e, *dense_total is %.4e pops_per_elec[0] is %.4e\n",
			label.c_str(),
			H2_renorm_chemistry ,
			*dense_total,
			pops_per_elec[0]);

	/* renormalize all level populations for the current chemical solution */
	for( qList::iterator st = states.begin(); st != states.end(); ++st )
	{
		long iElec = (*st).n();
		long iVib = (*st).v();
		long iRot = (*st).J();

		(*st).Pop() *= H2_renorm_chemistry;
		H2_old_populations[iElec][iVib][iRot] = (*st).Pop();
	}

	if(nTRACE >= n_trace_full )
		fprintf(ioQQQ,
		" H2 entry, old pops sumed to %.3e, renorm to htwo den of %.3e\n",
		pops_per_elec[0],
		*dense_total);

	/* >>chng 05 feb 10, reset convergence criteria if we are in search phase */
	fixit("I suspect this is counterproductive.  Test without it -- Ryan");
	if( conv.lgSearch )
	{
		converge_pops_relative *= 2.; /*def is 0.1 */
		converge_pops_total *= 3.;    /*def is 1e-3*/
		converge_ortho_para *= 3.;    /*def is 1e-2*/
	}

	if( !conv.nTotalIoniz )
		mole_update_species_cache();

	/* update state specific rates in H2_X_formation (cm-3 s-1) that H2 forms from grains and H- */
	mole_H2_form();

	/* evaluate total collision rates */
	H2_X_coll_rate_evaluate();

	/* this flag will say whether H2 populations have converged,
	 * by comparing old and new values */
	lgConv_h2_soln = false;
	/* this will count number of passes around following loop */
	long loop_h2_pops = 0;
	{
		if( nzone != nzoneEval )
		{
			nzoneEval = nzone;
			/* this is number of zones with H2 solution in this iteration */
			++nH2_zone;
		}
	}

	if( lgLTE )
		lgConv_h2_soln = true;

	/* begin - start level population solution
	 * first do electronic excited states, Lyman, Werner, etc
	 * using old solution for X
	 * then do matrix if used, then solve for pops of rest of X 
	 * >>chng 04 apr 06, subtract number of oscillations from limit - don't waste loops 
	 * if solution is unstable */
	while( loop_h2_pops < LIM_H2_POP_LOOP-n_pop_oscil && !lgConv_h2_soln && !lgLTE )
	{
		++cctr_diatom_l;

		/* this is number of trips around loop this time */
		++loop_h2_pops;
		/* this is number of times through this loop in entire iteration */
		++nH2_pops;

		/* radiative rates [cm-3 s-1] from electronic excited states into X vibration and rot */
		H2_X_rate_from_elec_excited.zero();
		/* radiative & cosmic ray rates [s-1] to electronic excited states from X */
		H2_X_rate_to_elec_excited.zero();
		H2_rad_rate_out.zero();
		H2_rad_rate_in.zero();
		pops_per_vib.zero();
		fill_n( pops_per_elec, N_ELEC, 0. );
		
		SolveExcitedElectronicLevels();

		/* evaluate rates that destroy or create ground electronic state */
		H2_X_sink_and_source();

		/* above set pops of excited electronic levels and found rates between them and X - 
		 * now solve highly excited levels within the X state by back-substitution */
		SolveSomeGroundElectronicLevels();

		/* now do lowest levels populations with matrix, 
		 * these should be collisionally dominated */
		if( nXLevelsMatrix )
		{
			H2_Level_low_matrix(
				/* the total abundance - frac_matrix is fraction of pop that was in these
				 * levels the last time this was done */
				(*dense_total) * (realnum)frac_matrix );
		}
		if(nTRACE >= n_trace_full) 
		{
			long iElecHi = 0;
			fprintf(ioQQQ," Rel pop(e=%li)" ,iElecHi);
		}

		/* find ortho and para densites, sum of pops in each vibration */
		/* this will become total pop is X, which will be renormed to equal *dense_total */
		{
			pops_per_elec[0] = 0.;
			for( md2i it = pops_per_vib.begin(0); it != pops_per_vib.end(0); ++it )
				*it = 0.;

			for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
			{
				long iElec = (*st).n();
				if( iElec > 0 ) continue;
				long iVib = (*st).v();
				pops_per_elec[iElec] += (*st).Pop();
				pops_per_vib[iElec][iVib] += (*st).Pop();
			}

			/* print sum of populations in each vibration if trace on */
			if(nTRACE >= n_trace_full) 
				for( md2ci it = pops_per_vib.begin(0); it != pops_per_vib.end(0); ++it )
					fprintf(ioQQQ,"\t%.2e", *it/(*dense_total));
	
			ASSERT( pops_per_elec[0] > SMALLFLOAT );
		}
		/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
		if(nTRACE >= n_trace_full) 
		{
			fprintf(ioQQQ,"\n");
			/* print the ground vibration state */
			fprintf(ioQQQ," Rel pop(0,J)");
			for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
			{
				long iElec = (*st).n();
				if( iElec > 0 ) continue;
				long iVib = (*st).v();
				if( iVib > 0 ) continue;
				fprintf(ioQQQ,"\t%.2e", (*st).Pop()/(*dense_total) );
			}
			fprintf(ioQQQ,"\n");
		}

		{
			double pop_total = 0.;
			for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
				pop_total += (*st).Pop();

			// The ratio of H2 that came out of the chemistry network to what we just obtained.
			double H2_renorm_conserve = *dense_total/ SDIV(pop_total);

			if (0)
				fprintf(ioQQQ,"DEBUG H2 %ld %ld %g %g %g %g %g\n", 
						  nzone, loop_h2_pops, H2_renorm_conserve,frac_matrix,pop_total,states[0].Pop(),states[1].Pop());

			/* renormalize populations  - were updated by renorm when routine entered, 
			 * before pops determined - but population determinations above do not have a sum rule on total
			 * population - this renorm is to preserve total population */
			pops_per_vib.zero();
			fill_n( pops_per_elec, N_ELEC, 0. );
			for( qList::iterator st = states.begin(); st != states.end(); ++st )
			{
				(*st).Pop() *= H2_renorm_conserve;
				long iElec = (*st).n();
				long iVib = (*st).v();
				pops_per_elec[iElec] += (*st).Pop();
				pops_per_vib[iElec][iVib] += (*st).Pop();
			}
		}

		/* now find population in states done with matrix - this is only used to pass
		 * to matrix solver */
		{
			double sum_pops_matrix = 0.;
			for( long i=0; i<nXLevelsMatrix; ++i )
			{
				sum_pops_matrix += states[i].Pop();
			}
			/* this is self consistent since pops_per_elec[0] came from current soln,
			 * as did the matrix.  pops will be renormalized by results from the chemistry
			 * a few lines down */
			frac_matrix = sum_pops_matrix / SDIV(*dense_total);
		}

		/* these will do convergence check */
		PopChgMaxOld_relative = PopChgMax_relative;
		PopChgMaxOld_total = PopChgMax_total;
		PopChgMax_relative = 0.;
		PopChgMax_total = 0.;
		iRotMaxChng_relative =-1;
		iVibMaxChng_relative = -1;
		iRotMaxChng_total =-1;
		iVibMaxChng_total = -1;
		popold_relative = 0.;
		popnew_relative = 0.;
		popold_total = 0.;
		popnew_total = 0.;

		// *****************************************
		// *****************************************
		// *****************************************
		// *****************************************
		//  Should be able to extract this loop!!!
		// *****************************************
		// *****************************************
		// *****************************************
		// *****************************************
		{
			/* this loop first checks for largest changes in populations, to determine whether
			 * we have converged, then updates the population array with a new value,
			 * which may be a mean of old and new
			 * update populations check convergence converged */
			double sumold = 0.;

			for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
			{
				long iElec = (*st).n();
				long iVib = (*st).v();
				long iRot = (*st).J();
				double pop = states[ ipEnergySort[iElec][iVib][iRot] ].Pop();
				/* keep track of largest relative change in populations to
				 * determines convergence */
				if(
					// >> chng 13 sep 21 -- pop test makes more sense as pre-condition
					pop>1e-6 * (*dense_total) && 
					/* >>chng 03 jul 19, this had simply been pop > SMALLFLOAT,
					 * change to relative pops > 1e-15, spent too much time converging
					 * levels at pops = 1e-37 */
					/* >>chng 03 dec 27, from rel pop 1e-15 to 1e-6 since converging heating will
					 * be main convergence criteria check convergence */
					fabs( pop - H2_old_populations[iElec][iVib][iRot])
					/* on first call some very high J states can have zero pop ,
					 * hence the SDIV, will retain sign for checks on oscilations,
					 * hence the fabs */
					> fabs(PopChgMax_relative)*pop
					  )
				{
					PopChgMax_relative = 
						(pop - H2_old_populations[iElec][iVib][iRot])/SDIV(pop);
					iRotMaxChng_relative = iRot;
					iVibMaxChng_relative = iVib;
					popold_relative = H2_old_populations[iElec][iVib][iRot];
					popnew_relative = pop;
				}
				/* >>chng 05 feb 08, add largest rel change in total, this will be converged
				 * down to higher accuracy than above 
				 * keep track of largest change in populations relative to total H2 to
				 * determine convergence check convergence */
				const double rel_change = (pop - H2_old_populations[iElec][iVib][iRot])/SDIV(*dense_total);
				/*  retain sign for checks on oscillations hence the fabs */
				if( fabs(rel_change) > fabs(PopChgMax_total) )
				{
					PopChgMax_total = rel_change;
					iRotMaxChng_total = iRot;
					iVibMaxChng_total = iVib;
					popold_total = H2_old_populations[iElec][iVib][iRot];
					popnew_total = pop;
				}

				kase = -1;
				/* update populations - we used the old populations to update the
				 * current new populations - will do another iteration if they changed
				 * by much.  here old populations are updated for next sweep through molecule */
				/* pop oscillations have occurred - use small changes */
				/* >>chng 04 may 10, turn this back on - now with min on how small frac new
				 * can become */
				const double abs_change = fabs( H2_old_populations[iElec][iVib][iRot] - pop );

				/* this branch very large changes, use mean of logs but onlly if both are positive*/
				if( abs_change > 3.*pop && H2_old_populations[iElec][iVib][iRot] * pop > 0 )
				{
					/* large changes or oscillations - take average in the log */
					H2_old_populations[iElec][iVib][iRot] = exp10(  
						log10(H2_old_populations[iElec][iVib][iRot])/2. +
						log10(pop)/2. );
					kase = 2;
				}

				/* modest change, use means of old and new */
				else if( abs_change> 0.1*pop )
				{
					realnum frac_old=0.25f;
					/* large changes or oscillations - take average */
					H2_old_populations[iElec][iVib][iRot] = 
						frac_old*H2_old_populations[iElec][iVib][iRot] +
						(1.f-frac_old)*pop;
					kase = 3;
				}
				else
				{
					/* small changes, use new value */
					H2_old_populations[iElec][iVib][iRot] = pop;
					kase = 4;
				}
				sumold += H2_old_populations[iElec][iVib][iRot];
			}

			/* will renormalize so that total population is correct */
			double H2_renorm_conserve_init = *dense_total/sumold;

			/* renormalize populations  - were updated by renorm when routine entered, 
			 * before pops determined - but population determinations above do not have a sum rule on total
			 * population - this renorm is to preserve total population */
			for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
			{
				long iElec = (*st).n();
				long iVib = (*st).v();
				long iRot = (*st).J();
				H2_old_populations[iElec][iVib][iRot] *= H2_renorm_conserve_init;
			}
		}

		/* get current ortho-para ratio, will be used as test on convergence */
		{
			ortho_density = 0.;
			para_density = 0.;
			H2_den_s = 0.;
			H2_den_g = 0.;

			for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
			{
				long iElec = (*st).n();
				long iVib = (*st).v();
				long iRot = (*st).J();
				const double& pop = (*st).Pop();
				/* find current population in H2s and H2g */
				if( (*st).energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
				{
					H2_den_s += pop;
				}
				else
				{
					H2_den_g += pop;
				}
				if( H2_lgOrtho[iElec][iVib][iRot] )
				{
					ortho_density += pop;
				}
				else
				{
					para_density += pop;
				}
			}
			ASSERT( fp_equal_tol( H2_den_s + H2_den_g, *dense_total, 1e-5 * (*dense_total) ) ); 
		}
		
		/* these will be used to determine whether solution has converged */
		ortho_para_older = ortho_para_old;
		ortho_para_old = ortho_para_current;
		ortho_para_current = ortho_density / SDIV( para_density );

		/* this will be evaluated in call to routine that follows - will check
		 * whether this has converged */
		old_solomon_rate = Solomon_dissoc_rate_g;

		/* >>chng 05 jul 24, break code out into separate routine for clarify
		 * located in mole_h2_etc.c - true says to only do Solomon rate */
		H2_Solomon_rate();

		/* are changes too large? must decide whether population shave converged,
		 * will check whether populations themselves have changed by much,
		 * but also change in heating by collisional deexcitation is stable */
		HeatChangeOld = HeatChange;
		HeatChange = HeatDexc_old - HeatDexc;
		{
			/* check whether pops are oscillating, as evidenced by change in
			 * heating changing sign */
			if( loop_h2_pops>2 && (
				(HeatChangeOld*HeatChange<0. ) ||
				(PopChgMax_relative*PopChgMaxOld_relative<0. ) ) )
			{
				lgH2_pops_oscil = true;
				if( loop_h2_pops > 6 )
				{
					loop_h2_oscil = loop_h2_pops;
					lgH2_pops_ever_oscil = true;
					++n_pop_oscil;
				}
			}
			else
			{
				lgH2_pops_oscil = false;
				/* turn off flag if no oscillations for a while */
				if( loop_h2_pops -  loop_h2_oscil > 4 )
				{
					lgH2_pops_ever_oscil = false;
				}
			}
		}

		/* reevaluate heating - cooling */
		HeatDexc_old = HeatDexc;
		H2_Cooling();

		/* begin check on whether solution is converged */
		lgConv_h2_soln = true;
		lgPopsConv_total = true;
		lgPopsConv_relative = true;
		lgHeatConv = true;
		lgSolomonConv = true;
		lgOrthoParaRatioConv = true;

		/* these are all the convergence tests 
		 * check convergence converged */
		if( fabs(PopChgMax_relative)>converge_pops_relative )
		{
			/*lgPopsConv = (fabs(PopChgMax_relative)<=0.1);*/
			lgConv_h2_soln = false;
			lgPopsConv_relative = false;
			/* >>chng 04 sep 08, set quant_new to new chng max gs */
			/*quant_old = PopChgMax_relative;*/
			quant_old = PopChgMaxOld_relative;
			/*quant_new = 0.;*/
			quant_new = PopChgMax_relative;

			strcpy( chReason , "rel pops changed" );
		}

		/* check largest change in a level population relative to total h2 
		 * population convergence converged check */
		else if( fabs(PopChgMax_total)>converge_pops_total)
		{
			lgConv_h2_soln = false;
			lgPopsConv_total = false;
			/* >>chng 04 sep 08, set quant_new to new chng max gs */
			/*quant_old = PopChgMax_relative;*/
			quant_old = PopChgMaxOld_total;
			/*quant_new = 0.;*/
			quant_new = PopChgMax_total;

			strcpy( chReason , "tot pops changed" );
		}

		/* >>chng 04 apr 30, look at change in ortho-para ratio, also that is not
		 * oscillating */
		/* >>chng 04 dec 15, only look at change, and don't make allowed change so tiny -
		 * these were attempts at fixing problems that were due to shielding not thin*/
		else if( fabs(ortho_para_current-ortho_para_old) / SDIV(ortho_para_current)> converge_ortho_para )
		/* else if( fabs(ortho_para_current-ortho_para_old) / SDIV(ortho_para_current)> 1e-3 
			&& (ortho_para_current-ortho_para_old)*(ortho_para_old-ortho_para_older)>0. )*/
		{
			lgConv_h2_soln = false;
			lgOrthoParaRatioConv = false;
			quant_old = ortho_para_old;
			quant_new = ortho_para_current;
			strcpy( chReason , "ortho/para ratio changed" );
		}
		/* >>chng 04 dec 16, reduce error allowed fm /5 to /2, to be similar to 
		 * logic in conv_base */
		else if( !thermal.lgTemperatureConstant &&
			fabs(HeatDexc-HeatDexc_old)/MAX2(thermal.ctot,thermal.htot) > 
			conv.HeatCoolRelErrorAllowed/10.
			/* >>chng 04 may 09, do not check on error in heating if constant temperature */
			/*&& !(thermal.lgTemperatureConstant || phycon.te <= phycon.TEMP_LIMIT_LOW  )*/ )
		{
			/* default on HeatCoolRelErrorAllowed is 0.02 */
			/*lgHeatConv = (fabs(HeatDexc-HeatDexc_old)/thermal.ctot <=
			 * conv.HeatCoolRelErrorAllowed/5.);*/
			lgConv_h2_soln = false;
			lgHeatConv = false;
			quant_old = HeatDexc_old/MAX2(thermal.ctot,thermal.htot);
			quant_new = HeatDexc/MAX2(thermal.ctot,thermal.htot);
			strcpy( chReason , "heating changed" );
			/*fprintf(ioQQQ,"DEBUG old new trip \t%.4e \t %.4e\n",
				HeatDexc_old,
				HeatDexc);*/
		}

		/* check on Solomon rate,
		 * >>chng 04 aug 28, do not do this check if induced processes are disabled,
		 * since Solomon process is then irrelevant */
		/* >>chng 04 sep 21, GS*/
		else if( rfield.lgInducProcess && 
			/* this is check that H2 abundance has not been set - if it has been
			 * then we don't care what the Solomon rate is doing */ 
			 hmi.H2_frac_abund_set==0 &&
			 /*>>chng 05 feb 10, rather than checking change in Solomon relative to Solomon,
			  * check it relative to total h2 destruction rate */
			fabs( Solomon_dissoc_rate_g - old_solomon_rate)/SDIV(hmi.H2_rate_destroy) > 
			conv.EdenErrorAllowed/5.)
		{
			lgConv_h2_soln = false;
			lgSolomonConv = false;
			quant_old = old_solomon_rate;
			quant_new = Solomon_dissoc_rate_g;
			strcpy( chReason , "Solomon rate changed" );
		}

		/* did we pass all the convergence test */
		if( !lgConv_h2_soln )
		{
			/* this branch H2 populations within X are not converged,
			 * print diagnostic */

			if( PRT_POPS || nTRACE >=n_trace_iterations )
			{
				/*fprintf(ioQQQ,"temppp\tnew\t%.4e\tnew\t%.4e\t%.4e\n",
					HeatDexc,
					HeatDexc_old,
					fabs(HeatDexc-HeatDexc_old)/thermal.ctot );*/
				fprintf(ioQQQ,"    %s loop %3li no conv oscl?%c why:%s ",
					label.c_str(),
					loop_h2_pops,
					TorF(lgH2_pops_ever_oscil),
					chReason );
				if( !lgPopsConv_relative )
					fprintf(ioQQQ," PopChgMax_relative:%.4e v:%li J:%li old:%.4e new:%.4e",
					PopChgMax_relative,
					iVibMaxChng_relative,
					iRotMaxChng_relative ,
					popold_relative ,
					popnew_relative );
				else if( !lgPopsConv_total )
					fprintf(ioQQQ," PopChgMax_total:%.4e v:%li J:%li old:%.4e new:%.4e",
					PopChgMax_total,
					iVibMaxChng_total,
					iRotMaxChng_total ,
					popold_total ,
					popnew_total );
				else if( !lgHeatConv )
					fprintf(ioQQQ," heat:%.4e old:%.4e new:%.4e",
					(HeatDexc-HeatDexc_old)/MAX2(thermal.ctot,thermal.htot), 
					quant_old , 
					quant_new);
				/* Solomon rate changed */ 
				else if( !lgSolomonConv )
					fprintf(ioQQQ," d(sol rate)/tot dest\t%2e",(old_solomon_rate - Solomon_dissoc_rate_g)/SDIV(hmi.H2_rate_destroy));
				else if( !lgOrthoParaRatioConv )
					fprintf(ioQQQ," current, old, older ratios are %.4e %.4e %.4e",
					ortho_para_current , ortho_para_old, ortho_para_older );
				else
					TotalInsanity();
				fprintf(ioQQQ,"\n");
			}
		}
		/* end convergence criteria */

		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, 
				"     H2 5lev %li Conv?%c",
				loop_h2_pops ,
				TorF(lgConv_h2_soln) );

			if( fabs(PopChgMax_relative)>0.1 )
				fprintf(ioQQQ," pops, rel chng %.3e",PopChgMax_relative);
			else
				fprintf(ioQQQ," rel heat %.3e rel chng %.3e H2 heat/cool %.2e",
					HeatDexc/thermal.ctot ,
					fabs(HeatDexc-HeatDexc_old)/thermal.ctot ,
					HeatDexc/thermal.ctot);

			fprintf( ioQQQ, 
				" Oscil?%c Ever Oscil?%c",
				TorF(lgH2_pops_oscil) ,
				TorF(lgH2_pops_ever_oscil) );
			fprintf(ioQQQ,"\n");
		}

		if( nTRACE >= n_trace_full ) 
		{
			fprintf(ioQQQ,
			"H2 loop\t%li\tkase pop chng\t%i\tchem renorm fac\t%.4e\tortho/para ratio:\t%.3e\tfrac of pop in matrix: %.3f\n",
			loop_h2_pops,
			kase,
			H2_renorm_chemistry,
			ortho_density / para_density ,
			frac_matrix);

			/* =======================INSIDE POPULATIONS CONVERGE LOOP =====================*/
			if( iVibMaxChng_relative>=0 && iRotMaxChng_relative>=0 && PopChgMax_relative>1e-10 )
				fprintf(ioQQQ,
					"end loop %li H2 max rel chng=%.3e from %.3e to %.3e at v=%li J=%li\n\n",
					loop_h2_pops,
					PopChgMax_relative , 
					H2_old_populations[0][iVibMaxChng_relative][iRotMaxChng_relative],
					states[ ipEnergySort[0][iVibMaxChng_relative][iRotMaxChng_relative] ].Pop(),
					iVibMaxChng_relative , iRotMaxChng_relative
					);
		}
	}
	/* =======================END POPULATIONS CONVERGE LOOP =====================*/

	/* >>chng 05 feb 08, do not print if we are in search phase */
	if( !lgConv_h2_soln && !conv.lgSearch )
	{
		conv.lgConvPops = false;
		lgPopsConverged = false;
		old_val = quant_old;
		new_val = quant_new;
	}

	for( qList::iterator st = states.begin(); st != states.end(); ++st )
	{
		ASSERT( (*st).Pop() >= 0. );
	}

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		/* following two heat exchange excitation, deexcitation */
		(*tr).Coll().cool() = 0.;
		(*tr).Coll().heat() = 0.;

		(*tr).Emis().PopOpc() = (*(*tr).Lo()).Pop() - (*(*tr).Hi()).Pop() * (*(*tr).Lo()).g() / (*(*tr).Hi()).g(); 

		/* number of photons in the line
		 * and line intensity */
		set_xIntensity( *tr );
	}
		
	average_energy_g = 0.;
	average_energy_s = 0.;
	/* determine average energy in ground and star */
	for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
	{
		double popTimesE = (*st).Pop() * (*st).energy().WN();
		if( (*st).energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
			average_energy_s += popTimesE;
		else
			average_energy_g += popTimesE;
	}
	/* average energy in ground and star */
	average_energy_g /= H2_den_g;
	if( H2_den_s > 1e-30 * (*dense_total) )
		average_energy_s /= H2_den_s;
	else
		average_energy_s = 0.;

	/* add up H2 + hnu => 2H, continuum photodissociation,
	 * this is not the Solomon process, true continuum */
	/* >>chng 05 jun 16, GS, add dissociation to triplet states*/
	photodissoc_BigH2_H2s = 0.;
	photodissoc_BigH2_H2g = 0.;
	/* >>chng 05 jul 20, GS, add dissociation by H2 g and H2s*/
	Average_collH_dissoc_g = 0.;
	Average_collH_dissoc_s = 0.;
	Average_collH2_dissoc_g = 0.;
	Average_collH2_dissoc_s = 0.;

	rel_pop_LTE_g =0.;
	rel_pop_LTE_s = 0.;
	
	double exp_disoc =  sexp(H2_DissocEnergies[0]/phycon.te_wn);

	/* >>chng 05 sep 12, TE, define a cutoff wavelength of 800 Angstrom 
	 * this is chosen as the cross sections given by 
	 *>>refer	H2	photo cs	Allison, A.C. & Dalgarno, A. 1969, Atomic Data, 1, 91 
	 * show a sharp decline in the cross section*/
	{
		static long ip_cut_off = -1;
		if( ip_cut_off < 0 )
		{
			/* one-time initialization of this pointer */
			ip_cut_off = ipoint( 1.14 );
		}

		/* >>chng 05 sep 12, TE, assume all H2s is at 2.5 eV
		 * the dissociation threshold is at 1.07896 Rydberg*/
		double flux_accum_photodissoc_BigH2_H2s = 0;
		fixit("this 2.5 seems like a pretty bad (and unnecessary) approximation.  Needs to be generalized at any rate.");
		long ip_H2_level = ipoint( 1.07896 - 2.5 / EVRYD);
		for( long i= ip_H2_level; i < ip_cut_off; ++i )
		{
			flux_accum_photodissoc_BigH2_H2s += ( rfield.flux[0][i-1] + rfield.ConInterOut[i-1]+ 
				rfield.outlin[0][i-1]+ rfield.outlin_noplot[i-1]  );
		}

		/* sum over all levels to obtain s and g populations and dissociation rates */
		for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
		{
			long iElec = (*st).n();
			if( iElec > 0 ) continue;
			long iVib = (*st).v();
			long iRot = (*st).J();
			const double &pop = (*st).Pop();
			fixit("generalize this factor (present value is (2m_e/m_H)^1.5/(2*2).  See Robin's Feb 7, 2009 email.");
			const double mass_stat_factor = 3.634e-5/(2*2);
			
			/* >>chng 05 mar 22, TE, this should be for H2* rather than total */
			/* this is the total rate of direct photo-dissociation of excited electronic states into 
			 * the X continuum - this is continuum photodissociation, not the Solomon process */
			/* >>chng 03 sep 03, make sum of pops of excited states */
			if( (*st).energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
			{
				double arg_ratio;
				photodissoc_BigH2_H2s += pop * flux_accum_photodissoc_BigH2_H2s;

				/* >>chng 05 july 20, GS, collisional dissociation, unit s-1*/
				Average_collH_dissoc_s += pop * H2_coll_dissoc_rate_coef[iVib][iRot];
				Average_collH2_dissoc_s += pop * H2_coll_dissoc_rate_coef_H2[iVib][iRot];

				/* >>chng 05 oct 17, GS, LTE populations of H2s*/
				arg_ratio = safe_div(exp_disoc,(*st).Boltzmann(),0.0);
				if( arg_ratio > 0. )
				{
					/* >>chng 05 oct 21, GS, only add ratio if Boltzmann factor > 0 */
					rel_pop_LTE_s += SAHA/SDIV(phycon.te32*arg_ratio)*
						(*st).g() * mass_stat_factor;
				}
			}
			else
			{
				double arg_ratio;
				/* >>chng 05 sep 12, TE, for H2g do the sum explicitly for every level*/
				double flux_accum_photodissoc_BigH2_H2g = 0;
				/* this is the dissociation energy needed for the level*/
				ip_H2_level = ipoint( 1.07896 - (*st).energy().Ryd() );

				for( long i= ip_H2_level; i < ip_cut_off; ++i )
				{
					flux_accum_photodissoc_BigH2_H2g += ( rfield.flux[0][i-1] + rfield.ConInterOut[i-1]+ 
						rfield.outlin[0][i-1]+ rfield.outlin_noplot[i-1] );
				}

				photodissoc_BigH2_H2g += pop * flux_accum_photodissoc_BigH2_H2g;

				/* >>chng 05 jun 28, TE, determine average energy level in H2g */
				average_energy_g += (pop * (*st).energy().WN() );
					
				/* >>chng 05 july 20, GS, collisional dissociation, unit s-1*/
				Average_collH_dissoc_g += pop * H2_coll_dissoc_rate_coef[iVib][iRot];
				Average_collH2_dissoc_g += pop * H2_coll_dissoc_rate_coef_H2[iVib][iRot];

				/* >>chng 05 oct 17, GS, LTE populations of H2g*/
				arg_ratio = safe_div(exp_disoc,(*st).Boltzmann(),0.0);
				if( arg_ratio > 0. )
				{
					rel_pop_LTE_g += SAHA/SDIV(phycon.te32*arg_ratio)*
						(*st).g() * mass_stat_factor;
				}
			}
		}
	}
	
	/* above sum was rate per unit vol since mult by H2 density, now div by H2* density to get rate s-1 */
	/* 0.25e-18 is wild guess of typical photodissociation cross section, from 
	 * >>refer	H2	dissoc	Allison, A.C. & Dalgarno, A. 1969, Atomic Data, 1, 91 
	 * this is based on an average of the highest v values they gave.  unfortunately, we want
	 * the highest J values - 
	 * final units are s-1*/ 
	
	if( H2_den_g > SMALLFLOAT )
	{
		Average_collH_dissoc_g /= SDIV(H2_den_g);/* unit cm3s-1*/
		Average_collH2_dissoc_g /= SDIV(H2_den_g);/* unit cm3s-1*/
		photodissoc_BigH2_H2g *= H2_DISS_ALLISON_DALGARNO / SDIV(H2_den_g);
	}
	else
	{
		Average_collH_dissoc_g = 0.;
		Average_collH2_dissoc_g = 0.;
		photodissoc_BigH2_H2g = 0.;
	}
	if( H2_den_s > SMALLFLOAT )
	{
		Average_collH_dissoc_s /= SDIV(H2_den_s);/* unit cm3s-1*/
		Average_collH2_dissoc_s /= SDIV(H2_den_s);/* unit cm3s-1*/
		photodissoc_BigH2_H2s *= H2_DISS_ALLISON_DALGARNO / SDIV(H2_den_s);
	}
	else
	{
		Average_collH_dissoc_s = 0.;
		Average_collH2_dissoc_s = 0.;
		photodissoc_BigH2_H2s = 0.;
	}


	// calculate some average rates from H2* to H2g 
	H2_Calc_Average_Rates();

	if( nTRACE )
	{
		fprintf(ioQQQ,"  H2_LevelPops exit1 %8.2f loops:%3li H2/H:%.3e Sol dis old %.3e new %.3e Sol dis star %.3e g-to-s %.3e photodiss star %.3e\n",
			fnzone ,
			loop_h2_pops ,
			dens_rel_to_lim_react,
			old_solomon_rate,
			Solomon_dissoc_rate_g,
			Solomon_dissoc_rate_s,
			gs_rate(),
			photodissoc_BigH2_H2s );
	}

	/* >>chng 03 sep 01, add this population - before had just used H2star from chem network */
	/* if big H2 molecule is turned on and used for this zone, use its
	 * value of H2* (pops of all states with v > 0 ) rather than simple network */

	/* update number of times we have been called */
	++nCall_this_iteration;

	/* this will say how many times the large H2 molecule has been called in this zone -
	 * if not called (due to low H2 abundance) then not need to update its line arrays */
	++nCall_this_zone;

	/* >>chng 05 jun 21,
	 * during search phase we want to use full matrix - save number of levels so that
	 * we can restore it */
	nXLevelsMatrix = nXLevelsMatrix_save;

	/* >>chng 05 jan 19, check how many levels should be in the matrix if first call on
	 * new zone, and we have a solution */
	/* end loop setting very first LTE populations */
	if( nCall_this_iteration && nzone != nzone_nlevel_set )
	{
		/* this is fraction of populations to include in matrix */
		const double FRAC = 0.99999;
		/* this loop is over increasing energy */
		double sum_pop = 0.;
		long nEner = 0;
		long iElec = 0;
		const bool PRT = false;
		if( PRT ) fprintf(ioQQQ,"DEBUG pops ");
		while( nEner < nLevels_per_elec[0] && sum_pop/(*dense_total) < FRAC )
		{
			/* array of energy sorted indices within X */
			ASSERT( iElec == ipElec_H2_energy_sort[nEner] );
			long iVib = ipVib_H2_energy_sort[nEner];
			long iRot = ipRot_H2_energy_sort[nEner];
			sum_pop += H2_old_populations[iElec][iVib][iRot];
			if( PRT ) fprintf(ioQQQ,"\t%.3e ", H2_old_populations[iElec][iVib][iRot]);
			++nEner;
		}
		if( PRT ) fprintf(ioQQQ,"\n");
		nzone_nlevel_set = nzone;
		/*fprintf(ioQQQ,"DEBUG zone %.2f old nmatrix %li proposed nmatrix %li sum_pop %.4e H2_total %.4e\n", 
			fnzone , nXLevelsMatrix ,nEner , sum_pop, *dense_total);
		nXLevelsMatrix = nEner;*/
	}

	return;
}
/*lint -e802 possible bad pointer */

void diatomics::SolveExcitedElectronicLevels( void )
{
	DEBUG_ENTRY( "diatomics::SolveExcitedElectronicLevels()" );

	multi_arr<double,3> rate_in;
	rate_in.alloc( H2_rad_rate_out.clone() );
	rate_in.zero();
	spon_diss_tot = 0.;
	double CosmicRayHILyaExcitationRate = ( hmi.lgLeidenCRHack ) ? secondaries.x12tot : 0.;

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		qList::iterator Lo = (*tr).Lo();	
		long iElecLo = (*Lo).n();
		long iVibLo = (*Lo).v();
		long iRotLo = (*Lo).J();
		qList::iterator Hi = (*tr).Hi();	
		long iElecHi = (*Hi).n();
		if( iElecHi < 1 ) continue;
		long iVibHi = (*Hi).v();
		long iRotHi = (*Hi).J();
		/* solve electronic excited state, 
		 * rate lower level in X goes to electronic excited state, s-1 
		 * first term is direct pump, second is cosmic ray excitation */
		/* collisional excitation of singlets by non-thermal electrons 
		 * this is stored ratio of electronic transition relative 
		 * cross section relative to the HI Lya cross section  */
		double rate_up = (*tr).Emis().pump() + CosmicRayHILyaExcitationRate * (*tr).Coll().col_str();
		double rate_down = 
			(*tr).Emis().Aul() * ( (*tr).Emis().Ploss() ) +
			rate_up * (*(*tr).Lo()).g() / (*(*tr).Hi()).g();

		/* this is a permitted electronic transition, must preserve nuclear spin */
		ASSERT( H2_lgOrtho[iElecHi][iVibHi][iRotHi] == H2_lgOrtho[iElecLo][iVibLo][iRotLo] );

		/* this is the rate [cm-3 s-1] electrons move into the upper level from X */
		rate_in[iElecHi][iVibHi][iRotHi] += H2_old_populations[iElecLo][iVibLo][iRotLo]*rate_up;

		/* rate [s-1] from levels within X to electronic excited states,
		 * includes photoexcitation and cosmic ray excitation */
		if( iElecLo==0 )
			H2_X_rate_to_elec_excited[iVibLo][iRotLo] += rate_up;
		H2_rad_rate_out[iElecLo][iVibLo][iRotLo] += rate_up;

		/* this is the rate [s-1] electrons leave the excited electronic upper level
		 * and decay into X - will be used to get pops of electronic excited states */
		H2_rad_rate_out[iElecHi][iVibHi][iRotHi] += rate_down; 
		ASSERT( rate_up >= 0. && rate_down >= 0. );
	}

	for( qList::iterator st = states.begin(); st != states.end(); ++st )
	{
		if( (*st).n() < 1 )
			continue;
	
		long iElec = (*st).n();
		long iVib = (*st).v();
		long iRot = (*st).J();

		H2_rad_rate_out[iElec][iVib][iRot] += H2_dissprob[iElec][iVib][iRot];
	
		/* update population [cm-3] of the electronic excited state this only includes 
		 * radiative processes between X and excited electronic states, and cosmic rays - 
		 * thermal collisions are neglected
		 * X is done below and includes all processes */
		double pop = rate_in[iElec][iVib][iRot] / SDIV( H2_rad_rate_out[iElec][iVib][iRot] );
		(*st).Pop() = pop;
		spon_diss_tot += pop * H2_dissprob[iElec][iVib][iRot];
		if( H2_old_populations[iElec][iVib][iRot]==0. )
			H2_old_populations[iElec][iVib][iRot] = pop;
		/* this is total pop in this vibration state */
		pops_per_vib[iElec][iVib] += pop;
		/* total pop in each electronic state */
		pops_per_elec[iElec] += pop;
	}

	fixit("uncomment and test");
	if( H2_den_s > 1e-30 * (*dense_total) )
		spon_diss_tot /= H2_den_s;
	else
		spon_diss_tot = 0.;

	if(nTRACE >= n_trace_full) 
	{
		for( long iElec=1; iElec<n_elec_states; ++iElec )
		{
			fprintf(ioQQQ," Pop(e=%li):",iElec);
			for( md2i it = pops_per_vib.begin(iElec); it != pops_per_vib.end(iElec); ++it )
				fprintf( ioQQQ,"\t%.2e", *it/(*dense_total) );
			fprintf(ioQQQ,"\n");
		}
	}

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		qList::iterator Lo = (*tr).Lo();	
		if( (*Lo).n() != 0 ) continue;
		qList::iterator Hi = (*tr).Hi();	
		if( (*Hi).n() < 1 ) continue;

		/* radiative rates [cm-3 s-1] from electronic excited states to X  */
		double rate = (*Hi).Pop() *
			((*tr).Emis().Aul() * ( (*tr).Emis().Ploss() ) +
			 (*tr).Emis().pump() * (*Lo).g() / (*Hi).g());
		H2_X_rate_from_elec_excited[(*Lo).v()][(*Lo).J()] += rate; 
		H2_rad_rate_in[(*Lo).v()][(*Lo).J()] += rate;
	}

	return;
}

void diatomics::SolveSomeGroundElectronicLevels( void )
{
	DEBUG_ENTRY( "diatomics::SolveSomeGroundElectronicLevels()" );

	/* these will be total rates into and out of the level */
	H2_col_rate_out.zero();
	H2_col_rate_in.zero();

	/* now evaluate total rates for all levels within X */
	for( long ipHi=0; ipHi<nLevels_per_elec[0]; ++ipHi)
	{
		/* array of energy sorted indices within X */
		ASSERT( ipElec_H2_energy_sort[ipHi] == 0 );
		long iVibHi = ipVib_H2_energy_sort[ipHi];
		long iRotHi = ipRot_H2_energy_sort[ipHi];

		realnum H2stat = states[ipHi].g();
		double EhiK = states[ipHi].energy().K();

		for( long ipLo=0; ipLo<ipHi; ++ipLo )
		{
			long iVibLo = ipVib_H2_energy_sort[ipLo];
			long iRotLo = ipRot_H2_energy_sort[ipLo];
	
			/* collision de-excitation [s-1] */
			realnum colldn = H2_X_coll_rate[ipHi][ipLo];
			/* inverse, rate up, [cm-3 s-1] */
			// use dsexp() explicitly rather than divide Boltzmann() factors
			// to avoid problems with underflow at low Te (see ticket #284)
			realnum collup = colldn *
				H2stat / states[ipLo].g() *
				dsexp( (EhiK - states[ipLo].energy().K()) / phycon.te );
			
			H2_col_rate_out[iVibHi][iRotHi] += colldn;
			H2_col_rate_in[iVibLo][iRotLo]  += colldn * H2_old_populations[0][iVibHi][iRotHi];

			H2_col_rate_out[iVibLo][iRotLo] += collup;
			H2_col_rate_in[iVibHi][iRotHi]  += collup * H2_old_populations[0][iVibLo][iRotLo];
		}
	}
		
	/* begin solving for X by back-substitution
	 * this is the main loop that determines populations within X 
	 * units of all rates in are cm-3 s-1, all rates out are s-1  
	 * nLevels_per_elec is number of levels within electronic 0 - so nEner is one
	 * beyond end of array here - but will be decremented at start of loop 
	 * this starts at the highest energy wihtin X and moves down to lower energies */
	long nEner = nLevels_per_elec[0];
	while( (--nEner) >= nXLevelsMatrix )
	{
		/* array of energy sorted indices within X - we are moving down
		 * starting from highest level within X */
		long iElec = ipElec_H2_energy_sort[nEner];
		ASSERT( iElec == 0 );
		long iVib = ipVib_H2_energy_sort[nEner];
		long iRot = ipRot_H2_energy_sort[nEner];

		if( nEner+1 < nLevels_per_elec[0] )
			ASSERT( states[nEner].energy().WN() < states[nEner+1].energy().WN() ||
				fp_equal( states[nEner].energy().WN(), states[nEner+1].energy().WN() ) );

		/* >>chng 05 apr 30,GS, Instead of *dense_total, the specific populations are used because high levels have much less
		 * populations than ground levels which consists most of the H2 population.
		 * only do this if working level is not v=0, J=0, 1 */ 
		if( nEner >1 )
		{
			H2_col_rate_out[iVib][iRot] += 
				/* H2 grain interactions
				 * rate (s-1) all v,J levels go to 0 or 1 preserving spin */
				(realnum)(rate_grain_op_conserve);

			/* this goes into v=0, and J=0 or 1 depending on whether initial
			 * state is ortho or para */
			H2_col_rate_in[0][H2_lgOrtho[0][iVib][iRot]] += 
				/* H2 grain interactions
				 * rate (cm-3 s-1) all v,J levels go to 0 or 1 preserving spin,
				 * in above lgOrtho says whether should go to 0 or 1 */
				(realnum)(rate_grain_op_conserve*H2_old_populations[0][iVib][iRot]);
		}
		else if( nEner == 1 )
		{
			/* this is special J=1 to J=0 collision, which is only fast at
			 * very low grain temperatures */
			H2_col_rate_out[0][1] += 
				/* H2 grain interactions
				 * H2 ortho - para conversion on grain surface,
				 * rate (s-1) all v,J levels go to 0 or 1, preserving nuclear spin */
				(realnum)(rate_grain_J1_to_J0);

			H2_col_rate_in[0][0] += 
				/* H2 grain interactions
				 * H2 ortho - para conversion on grain surface,
				 * rate (s-1) all v,J levels go to 0 or 1, preserving nuclear spin */
				(realnum)(rate_grain_J1_to_J0 *H2_old_populations[0][0][1]);
		}

		double pump_from_below = 0.;
		for( long ipLo = 0; ipLo<nEner; ++ipLo )
		{
			long iElecLo = ipElec_H2_energy_sort[ipLo];
			ASSERT( iElecLo == 0 );
			long iVibLo = ipVib_H2_energy_sort[ipLo];
			long iRotLo = ipRot_H2_energy_sort[ipLo];
			const TransitionList::iterator&tr = trans.begin() +ipTransitionSort[nEner][ipLo] ;

			/* the test on vibration is needed - the energies are ok but the space does not exist */
			if( ( abs(iRotLo-iRot) == 2 || iRotLo == iRot )  && (iVibLo <= iVib) && (*tr).ipCont() > 0 ) 
			{
				double rateone = (*tr).Emis().Aul() * ( (*tr).Emis().Ploss() );	
				// Pumping and cosmic-ray excitation from these levels up to higher (than X) levels is already included in H2_rad_rate_out before reaching here.
				// The following lines take care of that process within X
				double CosmicRayHILyaExcitationRate = ( hmi.lgLeidenCRHack ) ? secondaries.x12tot : 0.;
				double pump_up = (*tr).Emis().pump() + CosmicRayHILyaExcitationRate * (*tr).Coll().col_str();
				pump_from_below += pump_up * H2_old_populations[iElecLo][iVibLo][iRotLo];
				rateone += pump_up * (*(*tr).Lo()).g() / (*(*tr).Hi()).g();
				ASSERT( rateone >=0 );
				H2_rad_rate_out[0][iVib][iRot] += rateone;
				H2_rad_rate_in[iVibLo][iRotLo] += rateone * H2_old_populations[iElec][iVib][iRot];
			}
		}

		/* we now have the total rates into and out of this level, get its population 
		 * units cm-3 */
		states[nEner].Pop() =
			(H2_col_rate_in[iVib][iRot]+ H2_rad_rate_in[iVib][iRot]+H2_X_source[nEner]+pump_from_below) / 
			SDIV(H2_col_rate_out[iVib][iRot]+H2_rad_rate_out[0][iVib][iRot]+H2_X_sink[nEner]);

		ASSERT( states[nEner].Pop() >= 0. );
	}

	return;
}

/*H2_cooling evaluate cooling and heating due to H2 molecule */
#if defined(__ICC) && defined(__i386)
#pragma optimization_level 1
#endif
void diatomics::H2_Cooling( void )
{
	DEBUG_ENTRY( "H2_Cooling()" );

	/* nCall_this_iteration is not incremented until after the level
	 * populations have converged the first time.  so for the first n calls
	 * this will return zero, a good idea since populations will be wildly
	 * incorrect during search for first valid pops */
	if( !lgEnabled || !nCall_this_iteration )
	{
		HeatDexc = 0.;
		HeatDiss = 0.;
		HeatDexc_deriv = 0.;
		return;
	}
	
	HeatDiss = 0.;
	/* heating due to dissociation of electronic excited states */
	for( qList::iterator st = states.begin(); st != states.end(); ++st )
	{
		long iElec = (*st).n();
		long iVib = (*st).v();
		long iRot = (*st).J();
		HeatDiss += (*st).Pop() * H2_dissprob[iElec][iVib][iRot] * H2_disske[iElec][iVib][iRot];
	}
	
	/* dissociation heating was in eV - convert to ergs */
	HeatDiss *= EN1EV;

	/* now work on collisional heating due to bound-bound
	 * collisional transitions within X */
	HeatDexc = 0.;
	/* these are the colliders that will be considered as depopulating agents */
	/* the colliders are H, He, H2 ortho, H2 para, H+ */
	/* atomic hydrogen */

	/* this will be derivative */
	HeatDexc_deriv = 0.;
	long iElecHi = 0;
	for( long ipHi=1; ipHi<nLevels_per_elec[iElecHi]; ++ipHi )
	{
		realnum H2statHi = states[ipHi].g();
		double H2boltzHi = states[ipHi].Boltzmann();
		double H2popHi = states[ipHi].Pop();
		double ewnHi = states[ipHi].energy().WN();

		for( long ipLo=0; ipLo<ipHi; ++ipLo )
		{
			double rate_dn_heat = 0.;

			/* this sum is total downward heating summed over all colliders */
			mr3ci H2cr = CollRateCoeff.begin(ipHi, ipLo);
			for( long nColl=0; nColl<N_X_COLLIDER; ++nColl )
				/* downward collision rate */
				rate_dn_heat += H2cr[nColl]*collider_density[nColl];

			/* now get upward collisional cooling by detailed balance */
			double rate_up_cool = rate_dn_heat * states[ipLo].Pop() *
				/* rest converts into upward collision rate */
				H2statHi / states[ipLo].g() *
				safe_div(H2boltzHi, states[ipLo].Boltzmann(), 0.0 );

			rate_dn_heat *= H2popHi;

			/* net heating due to collisions within X - 
			 * positive if heating, negative is cooling
			 * this will usually be heating if X is photo pumped
			 * in printout and in save heating this is called "H2cX" */
			double conversion = (ewnHi - states[ipLo].energy().WN() ) * ERG1CM;
			double heatone = rate_dn_heat * conversion;
			double coolone = rate_up_cool * conversion;
			/* this is net heating, negative if cooling */
			double oneline = heatone - coolone;
			HeatDexc += oneline;

			/* derivative wrt temperature - assume exp wrt ground - 
			 * this needs to be divided by square of temperature in wn - 
			 * done at end of loop */
			HeatDexc_deriv +=  (realnum)(oneline * ewnHi);

			/* this would be a major logical error */
			ASSERT( 
				(rate_up_cool==0 && rate_dn_heat==0) || 
				(states[ipHi].energy().WN() > states[ipLo].energy().WN()) ); 
		}/* end loop over lower levels, all collisions within X */
	}/* end loop over upper levels, all collisions within X */

	/* this is inside h2 cooling, and is called extra times when H2 heating is important */
	if( PRT_POPS ) 
		fprintf(ioQQQ,
		"  DEBUG H2 heat fnzone\t%.2f\trenorm\t%.3e\tte\t%.4e\tdexc\t%.3e\theat/tot\t%.3e\n",
		fnzone , 
		H2_renorm_chemistry , 
		phycon.te , 
		HeatDexc,
		HeatDexc/thermal.ctot);

	/* this is derivative of collisional heating wrt temperature - needs 
	 * to be divided by square of temperature in wn */
	HeatDexc_deriv /=  (realnum)POW2(phycon.te_wn);

	if( nTRACE >= n_trace_full ) 
		fprintf(ioQQQ,
		" H2_Cooling Ctot\t%.4e\t HeatDiss \t%.4e\t HeatDexc \t%.4e\n" ,
		thermal.ctot , 
		HeatDiss , 
		HeatDexc );

	/* when we are very far from solution, during search phase, collisions within
	 * X can be overwhelmingly large heating and cooling terms, which nearly 
	 * cancel out.  Some dense cosmic ray heated clouds could not find correct
	 * initial solution due to noise introduced by large net heating which was
	 * the very noisy tiny difference between very large heating and cooling
	 * terms.  Do not include collisions with x as heat/cool during the
	 * initial search phase */
	if( conv.lgSearch )
	{
		HeatDexc = 0.;
		HeatDexc_deriv = 0.;
	}
	return;
}

/*cdH2_colden return column density in H2, negative -1 if cannot find state,
 * header is cdDrive */
double cdH2_colden( long iVib , long iRot )
{
	diatomics& diatom = h2;

	/*if iVib is negative, return
	 * total column density - iRot=0
	 * ortho column density - iRot 1
	 * para column density - iRot 2 
	 * else return column density in iVib, iRot */
	if( iVib < 0 )
	{
		if( iRot==0 )
		{
			/* return total column density */
			return( diatom.ortho_colden + diatom.para_colden );
		}
		else if( iRot==1 )
		{
			/* return ortho H2 column density */
			return diatom.ortho_colden;
		}
		else if( iRot==2 )
		{
			/* return para H2 column density */
			return diatom.para_colden;
		}
		else
		{
			fprintf(ioQQQ," iRot must be 0 (total), 1 (ortho), or 2 (para), returning -1.\n");
			return -1.;
		}
	}
	else if( diatom.lgEnabled )
	{
		return diatom.GetXColden( iVib, iRot );
	}
	/* error condition - no valid parameter */
	else
		return -1;
}

realnum diatomics::GetXColden( long iVib, long iRot )
{
	DEBUG_ENTRY( "diatomics::GetXColden()" );
	
	/* this branch want state specific column density, which can only result from
	 * evaluation of big molecule */
	int iElec = 0;
	if( iRot <0 || iVib >nVib_hi[iElec] || iRot > nRot_hi[iElec][iVib])
	{
		fprintf(ioQQQ," iVib and iRot must lie within X, returning -2.\n");
		fprintf(ioQQQ," iVib must be <= %li and iRot must be <= %li.\n",
			nVib_hi[iElec],nRot_hi[iElec][iVib]);
		return -2.;
	}
	else
		return H2_X_colden[iVib][iRot];
}

/*H2_Colden maintain H2 column densities within X */
void diatomics::H2_Colden( const char *chLabel )
{
	/* >>chng 05 jan 26, pops now set to LTE for small abundance case, so do this */
	if( !lgEnabled /*|| !nCall_this_zone*/ )
		return;

	DEBUG_ENTRY( "H2_Colden()" );

	if( strcmp(chLabel,"ZERO") == 0 )
	{
		/* zero out formation rates and column densites */
		H2_X_colden.zero();
		H2_X_colden_LTE.zero();
	}

	else if( strcmp(chLabel,"ADD ") == 0 )
	{
		/*  add together column densities */
		for( qList::iterator st = states.begin(); st != states.end(); ++st )
		{
			long iElec = (*st).n();
			if( iElec > 0 ) continue;
			long iVib = (*st).v();
			long iRot = (*st).J();
			/* state specific H2 column density */
			H2_X_colden[iVib][iRot] += (realnum)( (*st).Pop() * radius.drad_x_fillfac);
			/* LTE state specific H2 column density - H2_populations_LTE is normed to unity
			 * so must be multiplied by total H2 density */
			H2_X_colden_LTE[iVib][iRot] += (realnum)(H2_populations_LTE[0][iVib][iRot]*
				(*dense_total)*radius.drad_x_fillfac);
		}
	}

	/* we will not print column densities so skip that - if not print then we have a problem */
	else if( strcmp(chLabel,"PRIN") != 0 )
	{
		fprintf( ioQQQ, " H2_Colden does not understand the label %s\n", 
		  chLabel );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}

/*H2_DR choose next zone thickness based on H2 big molecule */
double diatomics::H2_DR(void)
{
	return BIGFLOAT;
}

/*H2_RT_OTS - add H2 ots fields */
void diatomics::H2_RT_OTS( void )
{
	/* do not compute if H2 not turned on, or not computed for these conditions */
	if( !lgEnabled || !nCall_this_zone )
		return;

	DEBUG_ENTRY( "H2_RT_OTS()" );

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		qList::iterator Hi = (*tr).Hi();	
		if( (*Hi).n() > 0 )
			continue;
								
		/* ots destruction rate */
		(*tr).Emis().ots() = (*(*tr).Hi()).Pop() * (*tr).Emis().Aul() * (*tr).Emis().Pdest();
								
		/* dump the ots rate into the stack - but only for ground electronic state*/
		RT_OTS_AddLine( (*tr).Emis().ots(), (*tr).ipCont() );
	}

	return;
}

void diatomics::H2_Calc_Average_Rates( void )
{
	/* >>chng 05 jul 09, GS*/ 
	/*  average Einstein value for H2* to H2g, GS*/
	double sumpop1 = 0.;
	double sumpopA1 = 0.;
	double sumpopcollH2O_deexcit = 0.;
	double sumpopcollH2p_deexcit = 0.;
	double sumpopcollH_deexcit = 0.;
	double popH2s = 0.;
	double sumpopcollH2O_excit = 0.;
	double sumpopcollH2p_excit = 0.;
	double sumpopcollH_excit = 0.;
	double popH2g = 0.;

	for( qList::const_iterator stHi = states.begin(); stHi != states.end(); ++stHi )
	{
		long iElecHi = (*stHi).n();
		if( iElecHi > 0 ) continue;
		long iVibHi = (*stHi).v();
		long iRotHi = (*stHi).J();
		double ewnHi = (*stHi).energy().WN();
		for( qList::const_iterator stLo = states.begin(); stLo != stHi; ++stLo )
		{
			long iVibLo = (*stLo).v();
			long iRotLo = (*stLo).J();
			double ewnLo2 = (*stLo).energy().WN();
			if( ewnHi > ENERGY_H2_STAR && ewnLo2 < ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
			{
				/* >>chng 05 jul 10, GS*/ 
				/*  average collisional rate for H2* to H2g, GS*/
				if( H2_lgOrtho[0][iVibHi][iRotHi] == H2_lgOrtho[0][iVibLo][iRotLo] )
				{ 
					long ihi = ipEnergySort[0][iVibHi][iRotHi];
					long ilo = ipEnergySort[0][iVibLo][iRotLo];
					const TransitionList::iterator&tr = 
						trans.begin()+ ipTransitionSort[ihi][ilo];
					double popHi = (*(*tr).Hi()).Pop();
					double popLo = (*(*tr).Lo()).Pop();

					/* sums of populations */
					popH2s += popHi;
					popH2g += popLo;

					/* sums of deexcitation rates - H2* to H2g */
					sumpopcollH_deexcit += popHi * CollRateCoeff[ihi][ilo][0];
					sumpopcollH2O_deexcit += popHi * CollRateCoeff[ihi][ilo][2];
					sumpopcollH2p_deexcit += popHi * CollRateCoeff[ihi][ilo][3];
						
					double temp = popLo * 
						(*stHi).g() / (*stLo).g() *
						safe_div((*stHi).Boltzmann(), (*stLo).Boltzmann(), 0.0 );

					/* sums of excitation rates - H2g to H2* */
					sumpopcollH_excit += temp * CollRateCoeff[ihi][ilo][0];
					sumpopcollH2O_excit += temp * CollRateCoeff[ihi][ilo][2];
					sumpopcollH2p_excit += temp * CollRateCoeff[ihi][ilo][3];

					if( lgH2_radiative[ihi][ilo] )
					{
						sumpop1 += popHi;
						sumpopA1 += popHi * (*tr).Emis().Aul();
					}
				}
			}
		}
	}
	Average_A = sumpopA1/SDIV(sumpop1);

	/* collisional excitation and deexcitation of H2g and H2s */
	Average_collH2_deexcit = (sumpopcollH2O_deexcit+sumpopcollH2p_deexcit)/SDIV(popH2s);
	Average_collH2_excit = (sumpopcollH2O_excit+sumpopcollH2p_excit)/SDIV(popH2g);
	Average_collH_excit = sumpopcollH_excit/SDIV(popH2g);
	Average_collH_deexcit = sumpopcollH_deexcit/SDIV(popH2s);

	/*fprintf(ioQQQ,
		"DEBUG Average_collH_excit sumpop = %.2e %.2e %.2e %.2e %.2e %.2e \n", 
		popH2g,popH2s,sumpopcollH_deexcit ,sumpopcollH_excit ,
		sumpopcollH_deexcit/SDIV(popH2s) ,sumpopcollH_excit/SDIV(popH2g));*/
	/*fprintf(ioQQQ,"sumpop = %le sumpopA = %le  Av= %le\n", 
	sumpop1,sumpopA1 , Average_A );*/

	return;
}

double diatomics::GetExcitedElecDensity( void )
{
	double H2_sum_excit_elec_den = 0.;
	for( long iElecHi=0; iElecHi<n_elec_states; ++iElecHi )
	{
		if( iElecHi > 0 )
			H2_sum_excit_elec_den += pops_per_elec[iElecHi];
	}

	return H2_sum_excit_elec_den;
}
/*lint +e802 possible bad pointer */

