/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvPresTempEdenIoniz solve for current pressure, calls PressureChange, ConvTempEdenIonize,
 * called by cloudy */
/*ConvFail handle convergence failure */
#include "cddefines.h"

#include "phycon.h"
#include "rt.h"
#include "dense.h"
#include "pressure.h"
#include "trace.h"
#include "conv.h"
#include "pressure_change.h"
#include "thermal.h"
#include "geometry.h" 
#include "grainvar.h"
#include "grains.h"

/*ConvPresTempEdenIoniz solve for current pressure, calls PressureChange, ConvTempEdenIoniz,
 * called by cloudy */
void ConvPresTempEdenIoniz()
{
	DEBUG_ENTRY( "ConvPresTempEdenIoniz()" );

	/* this will count number of times we call ConvBase in this zone,
	 * counter is incremented there 
	 * zero indicates first pass through solvers on this zone */
	conv.nPres2Ioniz = 0;
	conv.lgFirstSweepThisZone = true;
	conv.lgLastSweepThisZone = false;
	/* this will save the  history of density - pressure relationship
	 * for the current zone */
	if( nzone != conv.hist_pres_nzone )
	{
		/* first time in this zone - reset history */
		conv.hist_pres_nzone = nzone;
		conv.hist_pres_density.clear();
		conv.hist_pres_current.clear();
		conv.hist_pres_error.clear();
	}

	/* should still be true at end */
	conv.lgConvPops = true;

	/*rel_slope = 0.;*/

	if( trace.nTrConvg>=1  )
	{
		fprintf( ioQQQ, 
			" ConvPresTempEdenIoniz1 entered, will call ConvIoniz to initialize\n");
	}

	/* converge the ionization first, so that we know where we are, and have
	 * a valid foundation to begin the search */
	/* the true electron density dense.EdenTrue is set in eden_sum called by ConvBase */

	/* this evaluates current pressure, and returns whether or not 
	 * it is within tolerance of correct pressure */
	conv.lgConvPres = false; 

	/* convergence trace at this level */
	if( trace.nTrConvg>=1  )
	{
		fprintf( ioQQQ, 
			" ConvPresTempEdenIoniz1 ConvIoniz found following converged: Pres:%c, Temp:%c, Eden:%c, Ion:%c, Pops:%c\n", 
			TorF(conv.lgConvPres), 
			TorF(conv.lgConvTemp),
			TorF(conv.lgConvEden),
			TorF(conv.lgConvIoniz()),
			TorF(conv.lgConvPops));
	}

	/* trace convergence print at this level */
	if( trace.nTrConvg>=1  )
	{
		fprintf( ioQQQ, 
			"\n ConvPresTempEdenIoniz1 entering main pressure loop.\n");
	}

	AbundChange();

	if ( strcmp(dense.chDenseLaw,"CPRE") == 0 || 
		  strcmp(dense.chDenseLaw,"DYNA") == 0 ) 
	{
		/* set the initial temperature to the current value, so we will know
		 * if we are trying to jump over a thermal front */
		double TemperatureInitial = phycon.te;

		/* chng 02 dec 11 rjrw -- ConvIoniz => ConvTempEdenIoniz() here for consistency with inside loop */
		/* ConvIoniz; */
		ConvTempEdenIoniz();

		PresMode presmode;
		presmode.set();
		solverState st;
		st.press = pressureZone(presmode);

		/* this will be flag to check for pressure oscillations */
		bool lgPresOscil = false;
		bool lgStable = true;
		/* this is loop where it happened */
		long nloop_pres_oscil = 0;
		/* we will use these to check whether hden oscillating - would need to decrease step size */
		double hden_chng = 0.;
		/* this is dP_chng_factor, cut in half when pressure error changes sign */
		double dP_chng_factor = 1.;

		/* the limit to the number of loops */
		const int LOOPMAX = 100;
		/* this will be the limit, which we will increase if no oscillations occur */
		long LoopMax = LOOPMAX;
		long loop = 0;

		/* if start of calculation and we are solving for set pressure,
		 * allow a lot more iterations */
		if( nzone <= 1 && pressure.lgPressureInitialSpecified )
			LoopMax = 10*LOOPMAX;
		
		while( loop < LoopMax && !conv.lgConvPres )
		{
			/* there can be a pressure or density oscillation early in the search - if not persistent
			 * ok to clear flag */
			/* >>chng 01 aug 24, if large change in temperature allow lots more loops */
			if( fabs( TemperatureInitial - phycon.te )/phycon.te > 0.3 )
				LoopMax = 2*LOOPMAX;
			
			/* change current densities of all constituents if necessary, 
			 * PressureChange evaluates lgPresOK, true if pressure is now ok
			 * sets CurrentPressure and CorrectPressure */
			double hden_old = dense.gas_phase[ipHYDROGEN];
			/*pres_old = pressure.PresTotlCurr;*/
			
			/* this will evaluate current pressure, update the densities, 
			 * determine the wind velocity, and set conv.lgConvPres,
			 * return value is true if density was changed, false if no changes were necessary 
			 * if density changes then we must redo the temperature and ionization 
			 * PressureChange contains the logic that determines how to change the density to get
			 * the right pressure */

			static ConvergenceCounter cctr=conv.register_("PRES_CHANGES");
			++cctr;
			PressureChange(dP_chng_factor, presmode, st, lgStable);
			
			/* if product of these two is negative then hden is oscillating */
			double hden_chng_old = hden_chng;
			hden_chng = dense.gas_phase[ipHYDROGEN] - hden_old;
			if( fabs(hden_chng) < SMALLFLOAT )
				hden_chng = sign( (double)SMALLFLOAT, hden_chng );
			
			{
				/*@-redef@*/
				enum{DEBUG_LOC=false};
				/*@+redef@*/
				if( DEBUG_LOC && nzone > 150 && iteration > 1 )
				{
					fprintf(ioQQQ,"%li\t%.2e\t%.2e\n", 
							  nzone,
							  pressure.PresTotlCurr, 
							  pressure.PresTotlError*100.
						);
				}
			}
			
			/* check whether pressure is oscillating */
			if( ( ( hden_chng*hden_chng_old < 0. ) ) && loop > 1 )
			{
				/*fprintf(ioQQQ,"DEBUG\t%.2f\t%.2e\t%.2e\t%.2e\t%.2e\n",
				  fnzone,pres_chng,pres_chng_old ,hden_chng,hden_chng_old);*/
				lgPresOscil = true;
				nloop_pres_oscil = loop;
				/* >>chng 04 dec 09, go to this factor becoming smaller every time oscillation occurs */
				dP_chng_factor = MAX2( 0.125, dP_chng_factor * 0.5 );
				/* dP_chng_factor is how pressure changes with density - pass this to
				 * changing routine if it is stable */
				
				/*fprintf(ioQQQ,"oscilll %li %.2e %.2e %.2e %.2e dP_chng_factor %.2e\n", 
				  loop ,
				  pres_chng, 
				  pres_chng_old,
				  hden_chng , 
				  hden_chng_old ,
				  rel_slope);*/
			}
			
			/* convergence trace at this level */
			if( trace.nTrConvg>=1  )
			{
				fprintf( ioQQQ, 
							" ConvPresTempEdenIoniz1 %ld l:%3li nH:%.4e ne:%.4e PCurnt:%.4e err:%6.3f%% dPf:%.2e Te:%.4e Osc:%c Stb:%c\n", 
							nzone,
							loop, 
							dense.gas_phase[ipHYDROGEN], 
							dense.eden,
							pressure.PresTotlCurr, 
							/* this is percentage error */
							100.*pressure.PresTotlError,
							dP_chng_factor ,
							phycon.te,
							TorF(lgPresOscil),
							TorF(lgStable)  
					);
			}
			
			/* increment loop counter */
			++loop;
			
			/* there can be a pressure or density oscillation early in the search - if not persistent
			 * ok to clear flag */
			if( loop - nloop_pres_oscil > 4 )
				lgPresOscil = false;
			
			/* if we hit limit of loop, but no oscillations have happened, then we are
			 * making progress, and can keep going */
			if( loop == LoopMax && !lgPresOscil )
			{
				LoopMax = MIN2( LOOPMAX, LoopMax*2 );
			}
		}
	}
	else
	{
		double targetDensity = zoneDensity();
		double startingDensity = scalingDensity();
		double pdelta = conv.MaxFractionalDensityStepPerIteration;

		double logRatio = log(targetDensity/startingDensity);
		long nstep = (long) ceil(fabs(logRatio)/pdelta);
		// Ensure at least one update per zone
		if (nstep == 0)
			nstep = 1;
		double density_change_factor = exp(logRatio/nstep);

		for (long i=0; i<nstep; i++)
		{
			if (i == nstep - 1)
			{
				// Try to ensure exact answer at last iteration
				density_change_factor = targetDensity/scalingDensity();
			}

			/* this will evaluate current pressure, update the densities, 
			 * determine the wind velocity, and set conv.lgConvPres,
			 * return value is true if density was changed, false if no changes were necessary 
			 * if density changes then we must redo the temperature and ionization  */

			/* >>chng 04 feb 11, add option to remember current density and pressure */
			conv.hist_pres_density.push_back( dense.gas_phase[ipHYDROGEN] );
			conv.hist_pres_current.push_back( pressure.PresTotlCurr );
			conv.hist_pres_error.push_back( pressure.PresTotlError );
			
			if( trace.lgTrace )
			{
				fprintf( ioQQQ, 
							" DensityUpdate called, changing HDEN from %10.3e to %10.3e Set fill fac to %10.3e\n", 
							dense.gas_phase[ipHYDROGEN], density_change_factor*dense.gas_phase[ipHYDROGEN], 
							geometry.FillFac );
			}
			
			ScaleAllDensities((realnum) density_change_factor);
				
			/* must call TempChange since ionization has changed, there are some
			 * terms that affect collision rates (H0 term in electron collision) */
			TempChange(phycon.te , false);
			
			/* heating cooling balance while doing ionization, this is
			 * where the heavy lifting is done, this calls
			 * PresTotCurrent, which sets pressure.PresTotlCurr */
			ConvTempEdenIoniz();
			
			/* convergence trace at this level */
			if( trace.nTrConvg>=1  )
			{
				fprintf( ioQQQ, 
							" ConvPresTempEdenIoniz1 %.2f l:%3li nH:%.4e ne:%.4e PCurnt:%.4e err:%6.3f%% Te:%.4e Osc:%c\n", 
							fnzone,
							i, 
							dense.gas_phase[ipHYDROGEN], 
							dense.eden,
							pressure.PresTotlCurr, 
							/* this is percentage error */
							100.*pressure.PresTotlError,
							phycon.te,
							TorF(false)  );
			}
		}

		PresTotCurrent();
			
		conv.lgConvPres = true;
	}
	/* keep track of minimum and maximum temperature */
	thermal.thist = max((realnum)phycon.te,thermal.thist);
	thermal.tlowst = min((realnum)phycon.te,thermal.tlowst);

	/* >>chng 04 jan 31, now that all of the physics is converged, determine grain drift velocity */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
		GrainDrift();

	/* >>chng 01 mar 14, all ConvFail one more time, no matter how
	 * many failures occurred below.  Had been series of if, so multiple
	 * calls per failure possible. */
	/* >>chng 04 au 07, only announce pres fail here,
	 * we did not converge the pressure */
	if( !conv.lgConvIoniz() )
		ConvFail("ioni","");
	else if( !conv.lgConvEden )
		ConvFail("eden","");
	else if( !conv.lgConvTemp )
		ConvFail("temp","");
	else if( !conv.lgConvPres )
		ConvFail("pres","");

	/* this is only a sanity check that the summed continua accurately
	 * reflect all of the individual components.  Only include this
	 * when NDEBUG is not set, we are in not debug compile */
#	if !defined(NDEBUG)
	RT_OTS_ChkSum(0);
#	endif
}
