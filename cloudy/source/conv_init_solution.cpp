/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvInitSolution drive search for initial temperature, for illuminated face */
/*lgTempTooHigh returns true if temperature is too high */
/*lgCoolHeatCheckConverge return true if converged, false if change in heating cooling larger than set tolerance */
#include "cddefines.h"
#include "trace.h"
#include "struc.h"
#include "mole.h"
#include "dense.h"
#include "stopcalc.h"
#include "heavy.h"
#include "geometry.h"
#include "thermal.h"
#include "radius.h"
#include "phycon.h"
#include "pressure.h"
#include "conv.h"
#include "dynamics.h"
#include "rfield.h"

/* derivative of net cooling wrt temperature to check on sign oscillations */
static double dCoolNetDTOld = 0;

static double OxyInGrains , FracMoleMax;
/* determine whether chemistry is important - what is the largest
 * fraction of an atom in molecules?  Also is ice formation 
 * important
 * sets above two file static variables */

/*lgCoolHeatCheckConverge return true if converged, false if change in heating cooling larger than set tolerance */
STATIC bool lgCoolHeatCheckConverge( 
		double *CoolNet,
		bool lgReset )
{
	static double HeatOld=-1. , CoolOld=-1.;
	DEBUG_ENTRY( "lgCoolHeatCheckConverge()" );

	if( lgReset )
	{
		HeatOld = -1;
		CoolOld = -1;
	}

	double HeatChange = thermal.htot - HeatOld;
	double CoolChange = thermal.ctot - CoolOld;
	/* will use larger of heating or cooling in tests for relative convergence
	 * because in constant temperature models one may have trivially small value */
	double HeatCoolMax = MAX2( thermal.htot , thermal.ctot );

	/* is the heating / cooling converged?
	 * get max relative change, determines whether converged heating/cooling */
	double HeatRelChange = fabs(HeatChange)/SDIV( HeatCoolMax );
	double CoolRelChange = fabs(CoolChange)/SDIV( HeatCoolMax );
	bool lgConverged = true;
	if( MAX2( HeatRelChange , CoolRelChange ) > conv.HeatCoolRelErrorAllowed )
		lgConverged = false;

	*CoolNet = thermal.ctot - thermal.htot;

	HeatOld = thermal.htot;
	CoolOld = thermal.ctot;
	return lgConverged;
}

/* call lgCoolHeatCheckConverge until cooling and heating are converged */
STATIC bool lgCoolNetConverge( 
		double *CoolNet , 
		double *dCoolNetDT,
		bool lgReset )
{
	const long int LOOP_MAX=20;
	long int loop = 0;
	bool lgConvergeCoolHeat = false;

	DEBUG_ENTRY( "lgCoolNetConverge()" );

	while( loop < LOOP_MAX && !lgConvergeCoolHeat )
	{
		ConvEdenIoniz();
		lgConvergeCoolHeat = lgCoolHeatCheckConverge( CoolNet, lgReset && loop==0 );
		if( trace.lgTrace || trace.nTrConvg>=3 )
			fprintf(ioQQQ,"    lgCoolNetConverge %li calls to ConvEdenIoniz, converged? %c\n",
			loop , TorF( lgConvergeCoolHeat ) );
		++loop;
	}

	if( thermal.ConstTemp > 0 )
	{
		/* constant temperature model - this is trick so that temperature
		 * updater uses derivative to go to the set constant temperature */
		*CoolNet = phycon.te - thermal.ConstTemp;
		*dCoolNetDT = 1.f;
	}
	else if( phycon.te < 4000.f )
	{
		/* at low temperatures the analytical cooling-heating derivative
		 * is often of no value - the species densities change when the
		 * temperature changes.  Use simple dCnet/dT=(C-H)/T - this is 
		 * usually too large and results is smaller than optical dT, but
		 * we do eventually converge */
		*dCoolNetDT = thermal.ctot / phycon.te;
	}
	else if( thermal.htot / thermal.ctot > 2.f )
	{
		/* we are far from the solution and the temperature is much lower than
		 * equilibrium - analytical derivative from cooling evaluation is probably 
		 * wrong since coolants are not correct */
		*dCoolNetDT = thermal.ctot / phycon.te;
	}
	else 
	{
		*dCoolNetDT = thermal.dCooldT - thermal.dHeatdT;
		if( dCoolNetDTOld * *dCoolNetDT < 0. )
		{
			/* derivative is oscillating */
			fprintf(ioQQQ,"PROBLEM CoolNet derivative oscillating in initial solution "
					"Te=%.3e dCoolNetDT=%.3e CoolNet=%.3e dCooldT=%.3e dHeatdT=%.3e\n",
					phycon.te , *dCoolNetDT , *CoolNet,
					thermal.dCooldT , thermal.dHeatdT);
			*dCoolNetDT = *CoolNet / phycon.te;
		}
	}
	return lgConvergeCoolHeat;
}

/*ChemImportance find fraction of heavy elements in molecules, O in ices */
STATIC void ChemImportance( void )
{
	DEBUG_ENTRY( "ChemImportance()" );

	FracMoleMax = 0.;
	long int ipMax = -1;
	for( long i=0; i<mole_global.num_calc; ++i )
	{
		if (mole.species[i].location == NULL && !mole_global.list[i]->isMonatomic())
		{
			for( molecule::nNucsMap::iterator atom=mole_global.list[i]->nNuclide.begin();
				  atom != mole_global.list[i]->nNuclide.end(); ++atom)
			{
				long nelem = atom->first->el()->Z-1;
				if( dense.lgElmtOn[nelem] )
				{
					realnum dens_elemsp = (realnum) mole.species[i].den*atom->second;
					if ( FracMoleMax*dense.gas_phase[nelem] < dens_elemsp ) 
					{
						FracMoleMax = dens_elemsp/dense.gas_phase[nelem];
						ipMax = i;
					}
				}
			}
		}
	}

	/* fraction of all oxygen in ices on grains */
	OxyInGrains = 0.0;
	for(long i=0;i<mole_global.num_calc;++i)
	{
		if (! mole_global.list[i]->lgGas_Phase && mole_global.list[i]->isIsotopicTotalSpecies() )
		{
			OxyInGrains += mole.species[i].den*mole_global.list[i]->nElement(ipOXYGEN); 
		}
	}
	/* this is now fraction of O in ices */
	OxyInGrains /= SDIV(dense.gas_phase[ipOXYGEN]);

	{
		/* option to print out entire matrix */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"DEBUG fraction molecular %.2e species %li O ices %.2e\n" , 
				FracMoleMax , ipMax ,OxyInGrains );
		}
	}

	return;
}

STATIC double FindTempChangeFactor(void)
{
	double TeChangeFactor;

	DEBUG_ENTRY( "FindTempChangeFactor()" );

	/* find fraction of atoms in molecules - need small changes
	 * in temperature if fully molecular since chemistry
	 * network is complex and large changes in solution would
	 * cause linearization to fail */
	/* this evaluates the global variables FracMoleMax and
	 * OxyInGrains */
	ChemImportance();

	/* Following are series of chemical 
	 * tests that determine the factor by which
	 * which can change the temperature.  must be VERY small when 
	 * ice formation on grains is important due to exponential sublimation
	 * rate.  Also small when molecules are important, to keep chemistry
	 * within limits of linearized solver 
	 * the complete linearization that is implicit in all these solvers
	 * will not work when large Te jumps occur when molecules/ices
	 * are important - solvers can't return to solution if too far away
	 * keep temp steps small when mole/ice is important */
	if( OxyInGrains > 0.99  )
	{
		TeChangeFactor = 0.999;
	}
	else if( OxyInGrains > 0.9  )
	{
		TeChangeFactor = 0.99;
	}
	/* >>chng 97 mar 3, make TeChangeFactor small when close to lower bound */
	else if( phycon.te < 5. ||
		/*>>chng 06 jul 30, or if molecules/ices are important */
		OxyInGrains > 0.1  )
	{
		TeChangeFactor = 0.97;
	}
	/*>>chng 07 feb 23, add test of chemistry being very important */
	else if( phycon.te < 20. || FracMoleMax > 0.1 )
	{
		/* >>chng 07 feb 24, changed from 0,9 to 0.95 to converge 
		 * pdr_coolbb.in test */
		TeChangeFactor = 0.95;
	}
	else
	{
		/* this is the default largest step */
		TeChangeFactor = 0.8;
	}
	return TeChangeFactor;
}

/*ConvInitSolution drive search for initial temperature, for illuminated face */
void ConvInitSolution()
{
	long int i, 
	  ionlup, 
	  nelem ,
	  nflux_old,
	  nelem_reflux ,
	  ion_reflux;
	/* current value of Cooling - Heating */
	bool lgConvergeCoolHeat;

	double 
	  TeChangeFactor, 
	  TeBoundLow, 
	  TeBoundHigh;

	DEBUG_ENTRY( "ConvInitSolution()" );

	/* set NaN for safety */
	set_NaN( TeBoundLow );
	set_NaN( TeBoundHigh );

	/* this counts number of times ConvBase is called by PressureChange, in current zone */
	conv.nPres2Ioniz = 0;
	/* this counts how many times ConvBase is called in this iteration
	 * and is flag used by various routines to understand they are 
	 * being called the first time*/
	conv.nTotalIoniz = 0;
	conv.hist_pres_nzone = -1;
	conv.hist_temp_nzone = -1;
	/* ots rates not oscillating */
	conv.lgOscilOTS = false;

	dense.lgEdenBad = false;
	dense.nzEdenBad = 0;
	/* these are variables to remember the biggest error in the
	 * electron density, and the zone in which it occurred */
	conv.BigEdenError = 0.;
	conv.AverEdenError = 0.;
	conv.BigHeatCoolError = 0.;
	conv.AverHeatCoolError = 0.;
	conv.BigPressError = 0.;
	conv.AverPressError = 0.;

	/* these are equal if set dr was used, and we know the first dr */
	if( fp_equal( radius.sdrmin, radius.sdrmax ) )
	{
		// lgSdrmaxRel true if sdrmax is relative to current radius, false if limit in cm
		double rfac = radius.lgSdrmaxRel ? radius.Radius : 1.;
		radius.drad = MIN2( rfac*radius.sdrmax, radius.drad );
		radius.drad_mid_zone = radius.drad/2.;
		radius.drad_x_fillfac = radius.drad * geometry.FillFac;
	}

	/* the H+ density is zero in sims with no H-ionizing radiation and no cosmic rays,
	 * the code does work in this limit.  The H0 density can be zero in limit
         * of very high ionization where H0 underflows to zero. */
	ASSERT( dense.xIonDense[ipHYDROGEN][0] >=0 && dense.xIonDense[ipHYDROGEN][1]>= 0.);

	if( trace.nTrConvg )
		fprintf( ioQQQ, "\nConvInitSolution entered \n" );
	
	// Set to false to switch to using only standard temperature convergence method
	const bool lgDoInitConv = false;
	// Set false to disable search phase logic and test whether general solvers are sufficient
	const bool lgSearchPhaseLogicEnabled = true;

	/********************************************************************
	 *
	 * this is second or higher iteration, reestablish original temperature
	 *
	 *********************************************************************/
	if( iteration != 1 )
	{
		/* this is second or higher iteration on multi-iteration model */
		if( trace.lgTrace || trace.nTrConvg )
		{
			fprintf( ioQQQ, " ConvInitSolution called, ITER=%2ld resetting Te to %10.2e\n", 
				iteration, struc.testr[0] );
		}

		if( trace.lgTrace || trace.nTrConvg )
		{
			fprintf( ioQQQ, " search set true\n" );
		}

		/* search phase must be turned on so that variables such as the ots rates,
		 * secondary ionization, and auger yields, can converge more quickly to 
		 * proper values */
		conv.lgSearch = lgSearchPhaseLogicEnabled;
		conv.lgFirstSweepThisZone = true;
		conv.lgLastSweepThisZone = false;

		/* reset to the zone one temperature from the previous iteration */
		TempChange( struc.testr[0] , false);

		/* find current pressure - sets pressure.PresTotlCurr */
		PresTotCurrent();

		/* get new initial temperature and pressure */
		ConvPresTempEdenIoniz();

		if( trace.lgTrace || trace.nTrConvg )
		{
			fprintf( ioQQQ, " ========================================================================================================================\n");
			fprintf( ioQQQ, " ConvInitSolution: search set false 2 Te=%.3e========================================================================================\n" , phycon.te);
			fprintf( ioQQQ, " ========================================================================================================================\n");
		}
		conv.lgSearch = false;

		/* get temperature a second time */
		if( lgSearchPhaseLogicEnabled )
			ConvTempEdenIoniz();

		/* the initial pressure should now be valid 
		 * this sets values of pressure.PresTotlCurr */
		PresTotCurrent();
	}

	else
	{
		/********************************************************************
		 *
		 * do first te from scratch 
		 *
		 *********************************************************************/

		/* say that we are in search phase */
		conv.lgSearch = lgSearchPhaseLogicEnabled;
		conv.lgFirstSweepThisZone = true;
		conv.lgLastSweepThisZone = false;

		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " ConvInitSolution called, new temperature.\n" );
		}

		/* coming up to here Te is either 4,000 K (usually) or 10^6 */
		TeBoundLow = phycon.TEMP_LIMIT_LOW/1.001;

		/* set temperature floor option  - StopCalc.TeFloor is usually
		* zero but reset this this command - will go over to constant 
		* temperature calculation if temperature falls below floor */
		TeBoundLow = MAX2( TeBoundLow , StopCalc.TeFloor );

		/* highest possible temperature - must not get this high since
		 * parts of code will abort if value is reached. 
		 * divide by 1.2 to prevent checking on temp > 1e10 */
		TeBoundHigh = phycon.TEMP_LIMIT_HIGH/1.2;

		/* set initial temperature, options are constant temperature,
		 * or approach equilibrium from either high or low TE */
		double TeFirst;
		if( thermal.ConstTemp > 0 )
		{
			/* constant temperature , set 4000 K then relax to desired value
			 * for cold temperatures, to slowly approach fully molecular solution.
			 * if constant temperature is higher than 4000., go right to it */
			/* true allow ionization stage trimming, false blocks it */
			TeFirst = thermal.ConstTemp ;
			if (lgDoInitConv)
				TeFirst = MAX2( 4000. , TeFirst );

			/* override this if molecule deliberately disabled */
			if( mole_global.lgNoMole )
				TeFirst = thermal.ConstTemp;
		}
		else if( thermal.lgTeHigh )
		{
			/* approach from high TE */
			TeFirst = MIN2( 1e6 , TeBoundHigh );
		}
		else
		{
			/* approach from low TE */
			TeFirst = MAX2( 4000., TeBoundLow );
		}

		// initial kinetic temperature should be at or above the local
		// energy density temperature if the radiation field is well
		// coupled to the matter, but never let this overrule
		// constant temperature command
		if( !thermal.lgTemperatureConstant )
			TeFirst = max( TeFirst , phycon.TEnerDen );

		TempChange(TeFirst , false);
		if( trace.lgTrace || trace.nTrConvg>=2 )
			fprintf(ioQQQ,"ConvInitSolution doing initial solution with temp=%.2e\n", 
			phycon.te);

		if (lgDoInitConv)
		{
			/* this sets values of pressure.PresTotlCurr */
			PresTotCurrent();
			
			thermal.htot = 1.;
			thermal.ctot = 1.;
			
			/* call cooling, heating, opacity, loop to convergence
			 * this is very first call to it, by default is at 4000K */
			
			double CoolNet=0, dCoolNetDT=0;
			/* do ionization twice to get stable solution, evaluating initial dR each time */
			lgConvergeCoolHeat = false;
			for( ionlup=0; ionlup<2; ++ionlup )
			{
				if( trace.lgTrace || trace.nTrConvg>=2 )
					fprintf( ioQQQ, "ConvInitSolution calling lgCoolNetConverge "
								"initial setup %li with Te=%.3e C=%.2e H=%.2e d(C-H)/dT %.3e\n", 
								ionlup,phycon.te  , thermal.ctot , thermal.htot,dCoolNetDT );
				/* do not flag oscillating d(C-H)/dT until stable */
				dCoolNetDTOld = 0.;
				lgConvergeCoolHeat = lgCoolNetConverge( &CoolNet , &dCoolNetDT, true );
				/* set thickness of first zone, this affects the pumping rates due
				 * to correction for attenuation across zone, so must be known 
				 * for ConvEdenIoniz to get valid solution - this is why it
				 * is in this loop */
				radius_first();
			}
			if( !lgConvergeCoolHeat )
				fprintf(ioQQQ," PROBLEM ConvInitSolution did not converge the heating-cooling during initial search.\n");
			
			/* we now have error in C-H and its derivative - following is better
			 * derivative for global case where we may be quite far from C==H */
			if( thermal.ConstTemp<=0 )
				dCoolNetDT = thermal.ctot / phycon.te;
			bool lgConvergedLoop = false;
			long int LoopThermal = 0;
			/* initial temperature is usually 4000K, if the dTe scale factor is 0.95
			 * it will take 140 integrations to lower temperature to 3K,
			 * and many more if ices are important.  */
			const long int LIMIT_THERMAL_LOOP=300;
			double CoolMHeatSave=-1. , TempSave=-1. , TeNew=-1.,CoolSave=-1;
			while( !lgConvergedLoop && LoopThermal < LIMIT_THERMAL_LOOP )
			{
				/* change temperature until sign of C-H changes */
				CoolMHeatSave = CoolNet;
				dCoolNetDTOld = dCoolNetDT;
				CoolSave = thermal.ctot;
				TempSave = phycon.te;
				
				/* find temperature change factor, a number less than 1*/
				TeChangeFactor = FindTempChangeFactor();
				ASSERT( TeChangeFactor>0. && TeChangeFactor< 1. );
				
				TeNew = phycon.te - CoolNet / dCoolNetDT;
				
				TeNew = MAX2( phycon.te*TeChangeFactor , TeNew );
				TeNew = MIN2( phycon.te/TeChangeFactor , TeNew );
				/* change temperature */
				TempChange(TeNew , true);
				lgConvergeCoolHeat = lgCoolNetConverge( &CoolNet , & dCoolNetDT, false );
				
				if( trace.lgTrace || trace.nTrConvg>=2 )
					fprintf(ioQQQ,"ConvInitSolution %4li TeNnew=%.3e "
							  "TeChangeFactor=%.3e Cnet/H=%.2e dCnet/dT=%10.2e Ne=%.3e"
							  " Ograins %.2e FracMoleMax %.2e\n",
							  LoopThermal , TeNew , TeChangeFactor ,
							  CoolNet/SDIV(thermal.htot) ,  dCoolNetDT,
							  dense.eden , OxyInGrains , FracMoleMax );
				
				/* keep changing temperature until sign of CoolNet changes
				 * in constant temperature case CoolNet is delta Temperature
				 * so is zero for last pass in this loop 
				 * this is for constant temperature case */
				if( fabs(CoolNet)< SMALLDOUBLE )
					/* CoolNet is zero */
					lgConvergedLoop = true;
				else if( (CoolNet/fabs(CoolNet))*(CoolMHeatSave/fabs(CoolMHeatSave)) <= 0)
				/* change in sign of net cooling */
					lgConvergedLoop = true;
				else if( phycon.te <= TeBoundLow || phycon.te>=TeBoundHigh)
					lgConvergedLoop = true;
				++LoopThermal;
			}
			
			if( !lgConvergedLoop )
				fprintf(ioQQQ,"PROBLEM ConvInitSolution: temperature solution not "
						  "found in initial search, final Te=%.2e\n",
						  phycon.te );
			
			/* interpolate temperature where C=H if not constant temperature sim 
			 * will have set constant temperature mode above if we hit temperature
			 * floor */
			if( thermal.ConstTemp<=0 && 
				 ! fp_equal( TempSave , phycon.te ) )
			{
				if( trace.lgTrace || trace.nTrConvg>=2 )
					fprintf(ioQQQ," Change in sign of C-H found, Te brackets %.3e "
							  "%.3e Cool brackets %.3e %.3e ",
							  TempSave , phycon.te , CoolMHeatSave , CoolNet);
				/* interpolate new temperature assuming Cool = a T^alpha,
				 * first find alpha */
				double alpha = log(thermal.ctot/CoolSave) / log( phycon.te/TempSave);
				if( fabs(alpha) < SMALLFLOAT )
					/* alpha close to 0 means constant temperature */
					TeNew = phycon.te;
				else
				{
					/* next find log a - work in logs to avoid infinities */
					if( thermal.ctot<0. || thermal.htot<0. )
						TotalInsanity();
					double Alog = log( thermal.ctot ) - alpha * log( phycon.te );
					/* the interpolated temperature where heating equals cooling */
					double TeLn = (log( thermal.htot ) - Alog) / alpha ;
					
					/* a sanity check */
					if( TeLn < log( MIN2(phycon.te , TempSave )) )
						TeNew = MIN2(phycon.te , TempSave );
					else if( TeLn > log( MAX2(phycon.te , TempSave )) )
						TeNew = MAX2(phycon.te , TempSave );
					else
						TeNew = pow( EE , TeLn );
				}

				ASSERT( TeNew>=MIN2 ( TempSave , phycon.te ) ||
						  TeNew<=MAX2 ( TempSave , phycon.te ));
				
				if( trace.lgTrace || trace.nTrConvg>=2 )
					fprintf(ioQQQ," interpolating to Te %.3e \n\n",
							  TeNew);
				
				/* change temperature */
				TempChange(TeNew , true);
			}
		}

		if( lgSearchPhaseLogicEnabled )
			ConvTempEdenIoniz();

		conv.lgSearch = false;

		// Solve again without search phase lo-fi physics before starting on
		// anything which might have a long-term impact
		ConvTempEdenIoniz();

		/* this sets values of pressure.PresTotlCurr */
		PresTotCurrent();

		if( trace.lgTrace || trace.nTrConvg )
		{
			fprintf( ioQQQ, "\nConvInitSolution: end 1st iteration search phase "
				"finding Te=%.3e\n\n" , phycon.te);
		}

		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " ConvInitSolution return, TE:%10.2e==================\n", 
			  phycon.te );
		}
	}

	/* we now have a fairly good temperature and ionization 
	 * iteration is 1 on first iteration - remember current pressure
	 * on first iteration so we can hold this constant in constant 
	 * pressure simulation.
	 *
	 * flag dense.lgDenseInitConstant false if
	 * *constant pressure reset* is used - 
	 * default is true, after first iteration initial density is used for
	 * first zone no matter what pressure results.  Small changes occur due
	 * to radiative transfer convergence, 
	 * when set false with *reset* option the density on later iterations 
	 * can change to keep the pressure constant */
	if( iteration==1 || dense.lgDenseInitConstant )
	{
		double PresNew = pressure.PresTotlCurr;

		/* option to specify the pressure as option on constant pressure command */
		if( pressure.lgPressureInitialSpecified )
			/* this is log of nT product - if not present then set zero */
			PresNew = pressure.PressureInitialSpecified;
		if( trace.lgTrace )
		{
			fprintf( ioQQQ, 
				"     PresTotCurrent 1st zn old values of PresTotlInit:%.3e"
				" to %.3e \n",
				pressure.PresTotlInit,
				PresNew);
		}

		pressure.PresTotlInit = PresNew;
	}

	if( dynamics.lgTimeDependentStatic )
	{
		// some sort of time dependent sim
		static double PresTotalInitTime;
		if( dynamics.isInitialRelaxIteration( iteration ) )
		{
			PresTotalInitTime = pressure.PresTotlInit;
		}
		else if (dense.lgPressureVaryTime)
		{
	 		pressure.PresTotlInit = PresTotalInitTime *
				pow( 1.+(dynamics.time_elapsed/dense.PressureVaryTimeTimescale) ,
				dense.PressureVaryTimeIndex);
		}
//		fprintf(ioQQQ,"DEBUG conv_init_solution time dependent time=%.2e sets "
//		  "pressure to %.2e\n", dynamics.time_elapsed ,
//		  pressure.PresTotlInit);
	}

	/* now bring current pressure to the correct pressure - must do this
	 * before initial values are saved in iter start/end */
	ConvPresTempEdenIoniz();

	/* this counts number of times ConvBase is called by PressureChange, in
	 * current zone these are reset here, so that we count from first zone 
	 * not search */
	conv.nPres2Ioniz = 0;

	dense.lgEdenBad = false;
	dense.nzEdenBad = 0;

	/* save counter upon exit so we can see how many iterations were
	 * needed to do true zones */
	conv.nTotalIoniz_start = conv.nTotalIoniz;

	/* these are variables to remember the biggest error in the
	 * electron density, and the zone in which it occurred */
	conv.BigEdenError = 0.;
	conv.AverEdenError = 0.;
	conv.BigHeatCoolError = 0.;
	conv.AverHeatCoolError = 0.;
	conv.BigPressError = 0.;
	conv.AverPressError = 0.;

	/*>>chng 06 jun 09, only reset on first try - otherwise have stable solution? */
	if( iteration == 1 )
	{
		/* now remember some things we may need even in first zone, these
		* are normally set towards end of zone calculation in RT_tau_inc */
		struc.testr[0] = (realnum)phycon.te;
		/* >>chng 02 May 2001 rjrw: add hden for dilution */
		struc.hden[0] = dense.gas_phase[ipHYDROGEN];
		/* pden is the total number of particles per unit vol */
		struc.DenParticles[0] = dense.pden;
		struc.heatstr[0] = thermal.htot;
		struc.coolstr[0] = thermal.ctot;
		struc.volstr[0] = (realnum)radius.dVeffAper;
		struc.drad_x_fillfac[0] = (realnum)radius.drad_x_fillfac;
		struc.histr[0] = dense.xIonDense[ipHYDROGEN][0];
		struc.hiistr[0] = dense.xIonDense[ipHYDROGEN][1];
		struc.ednstr[0] = (realnum)dense.eden;
		struc.o3str[0] = dense.xIonDense[ipOXYGEN][2];
		struc.DenMass[0] = dense.xMassDensity;
		struc.drad[0] = (realnum)radius.drad;
	}

	/* check that nflux extends above IP of highest ionization species present.
	 * for collisional case it is possible for species to exist that are higher
	 * IP than the limit to the continuum.  Need continuum to encompass all
	 * possible emission - to account for diffuse emission
	 * NB 
	 * on second iteration of multi-iteration model this may result in rfield.nflux increasing
	 * which can change the solution */
	nflux_old = rfield.nflux;
	nelem_reflux = -1;
	ion_reflux = -1;
	for( nelem=2; nelem < LIMELM; nelem++ )
	{
		/* do not include hydrogenic species in following */
		for( i=0; i < nelem; i++ )
		{
			if( dense.xIonDense[nelem][i+1] > 0. )
			{
				if( Heavy.ipHeavy[nelem][i] > rfield.nflux )
				{
					rfield.nflux = Heavy.ipHeavy[nelem][i];
					nelem_reflux = nelem;
					ion_reflux = i;
				}
			}
		}
	}
	rfield.nPositive = rfield.nflux;

	/* was the upper limit to the continuum updated? if so need to define
	 * continuum variables to this new range */
	if( nflux_old != rfield.nflux )
	{
		/* zero out parts of rfield arrays that were previously undefined */
		rfield_opac_zero( nflux_old-1 , rfield.nflux );

		/* if continuum reset up, we need to define gaunt factors through high end */
		/*tfidle(false);*/
		/* this calls tfidle, among other things */
		/* this sets values of pressure.PresTotlCurr */
		PresTotCurrent();

		/* redo ionization and update opacities */
		ConvBase(1);

		/* need to update continuum opacities */
		if( trace.lgTrace )
		{
			fprintf(ioQQQ," nflux updated from %li to %li, anu from %g to %g \n",
				nflux_old , rfield.nflux ,
				rfield.anu(nflux_old) , rfield.anu(rfield.nflux) );
			fprintf(ioQQQ," caused by element %li ion %li \n",
			  nelem_reflux ,ion_reflux );
		}
	}

	/* dynamics solution - change density slightly
	 * and call pressure solver to see if it returns to original density */
	if( strcmp(dense.chDenseLaw,"DYNA") == 0 )
	{
		/* >>chng 04 apr 23, change pressure and let solver come back to correct
		 * pressure.  This trick makes sure we are correctly either super or sub sonic 
		 * since solver will jump from one branch to the other if necessary */
		static const double PCHNG = 0.98;
		/* this ignores return condition, assume must be ok if got this far */
		ConvPresTempEdenIoniz();

		pressure.PresTotlInit *= PCHNG;
		ConvPresTempEdenIoniz();

		pressure.PresTotlInit /= PCHNG;
		ConvPresTempEdenIoniz();
	}
}
