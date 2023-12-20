/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PressureChange called by ConvPresTempEdenIoniz
 * evaluate the current pressure, density change needed to get it to converge,
 * applies this correction factor to all gas constituents,
 * sets conv.lgConvPres true if good pressure, false if pressure change capped */
/*lgConvPres finds amount by which density should change to move towards pressure equilibrium
 * returns true if pressure is converged */
#include "cddefines.h"

#include "pressure_change.h"

#include "colden.h"
#include "conv.h"
#include "cosmology.h"
#include "dark_matter.h"
#include "dense.h"
#include "dynamics.h"
#include "geometry.h" 
#include "phycon.h"
#include "pressure.h"
#include "radius.h"
#include "struc.h"
#include "thermal.h"
#include "trace.h"
#include "wind.h"

/* zoneDensity finds density where prescribed */
double zoneDensity()
{
	double new_density;

	DEBUG_ENTRY( "zoneDensity()" );

	pressure.PresTotlError = 0.;

	/* evaluate a series of possible options, and set new_density */
	/* inside out globule */
	if( strcmp(dense.chDenseLaw,"GLOB") == 0 )
	{
		/* GLBDST is distance from globule, or glbrad-DEPTH */
		if( radius.glbdst < 0. )
		{
			fprintf( ioQQQ, " Globule distance is negative, internal overflow has occured,  sorry.\n" );
			fprintf( ioQQQ, " This is routine lgConvPres, GLBDST is%10.2e\n", 
			  radius.glbdst );
			cdEXIT(EXIT_FAILURE);
		}
		new_density = 
			(radius.glbden*pow(radius.glbrad/(radius.glbdst),radius.glbpow))*
			(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
	}
	else if( cosmology.lgDo )
	{
		/* cosmological - density varies because of expansion of universe */
		double dnew = GetDensity( cosmology.redshift_current );
		new_density = dnew*(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
	}
	else if( (strcmp(dense.chDenseLaw,"WIND") == 0) ) {

		if ( wind.lgStatic() )
		{
			fprintf( ioQQQ, " PROBLEM WIND called with zero velocity - this is impossible.\n Sorry.\n" );
			TotalInsanity();
		}
		
		/* this is positive wind velocity the outflowing wind beyond sonic point */
		else if( wind.lgBallistic() )
		{
			
			/* this is logic for supersonic wind solution,
			 * well above sonic point. */
			if( nzone > 1 )
			{
				/* Wind model */
				
				if( trace.lgTrace && trace.lgWind )
				{
					fprintf(ioQQQ," lgConvPres sets AccelGravity %.3e lgDisk?%c\n",
							  wind.AccelGravity , 
							  TorF(wind.lgDisk) );
				}
				
				/* following is form of constant acceleration equation, v^2 = u^2 + 2 a s
				 * struc.windv[nzone-2] is velocity of previous zone
				 * this increments that velocity to form square of new wind 
				 * velocity for outer edge of this zone */
				double term = POW2(struc.windv[nzone-2]) + 2.*(wind.AccelTotalOutward - wind.AccelGravity)* radius.drad;
				
				/* increment velocity if it is substantially positive */
				fixit("RHS of comparison should be ~ sound speed squared");
				if( term <= 1e3 )
				{
					/* wind velocity is well below sonic point, give up, 
					 * do not change velocity */
					wind.lgVelPos = false;
				}
				else
				{
					/* wind.windv is velocity at OUTER edge of this zone */
					term = sqrt(term);
					if (wind.windv > 0)
					{
						wind.windv = (realnum) term;
					}
					else
					{
						wind.windv = -(realnum) term;
					}
					wind.lgVelPos = true;
				}
				
				if( trace.lgTrace && trace.lgWind )
				{
					fprintf(ioQQQ," lgConvPres new wind V zn%li %.3e AccelTotalOutward %.3e AccelGravity %.3e\n",
							  nzone,wind.windv, wind.AccelTotalOutward, wind.AccelGravity );
				}
			}
			
			/* conservation of mass sets density here */
			new_density = wind.emdot/(wind.windv*radius.r1r0sq) * 
				(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
		}
		else
		{
			fprintf(ioQQQ,"chDenseLaw==\"WIND\" must now be ballistic or static\n");
			TotalInsanity();
		}
	}

	else if( strcmp(dense.chDenseLaw,"SINE") == 0 )
	{
		/* rapid density fluctuation */
		if( dense.lgDenFlucRadius )
		{
			new_density = 
				(dense.cfirst*cos(radius.depth*dense.flong+dense.flcPhase) + 
				 dense.csecnd)*(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
		}
		else
		{
			new_density = 
				(dense.cfirst*cos(colden.colden[ipCOL_HTOT]*dense.flong+dense.flcPhase) + 
				 dense.csecnd)*(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
		}
	}

	else if( strcmp(dense.chDenseLaw,"POWR") == 0 )
	{
		/* power law function of radius */
		double dnew = dense.den0*pow(radius.Radius/radius.rinner,(double)dense.DensityPower);
		new_density = dnew*(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
	}

	else if( strcmp(dense.chDenseLaw,"POWD") == 0 )
	{
		/* power law function of depth */
		double dnew = dense.den0*pow(1. + radius.depth/dense.rscale,(double)dense.DensityPower);
		new_density = dnew*(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
	}

	else if( strcmp(dense.chDenseLaw,"POWC") == 0 )
	{
		/* power law function of column density */
		double dnew = dense.den0*pow(1.f + colden.colden[ipCOL_HTOT]/
		  dense.rscale,dense.DensityPower);
		new_density = dnew*(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
	}


	else if( strncmp( dense.chDenseLaw ,"DLW" , 3) == 0 )
	{
		if( strcmp(dense.chDenseLaw,"DLW1") == 0 )
		{
			/* call ACF sub */
			new_density = 
				dense_fabden(radius.Radius,radius.depth)*
				(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
		}
		else if( strcmp(dense.chDenseLaw,"DLW2") == 0 )
		{
			/* call table interpolation subroutine
			 * >>chng 96 nov 29, added dense_tabden */
			new_density = dense.DLW.tabval(radius.Radius,radius.depth) *
				(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
		}
		else if( strcmp(dense.chDenseLaw,"DLW3") == 0 )
		{
			/* call parametrized wind subroutine */
			new_density = dense_parametric_wind(radius.Radius) *
				(scalingDensity()/dense.gas_phase[ipHYDROGEN]);
		}
		else
		{
			fprintf( ioQQQ, " Insanity, lgConvPres gets chCPres=%4.4s\n", 
			  dense.chDenseLaw );
			cdEXIT(EXIT_FAILURE);
		}

		if( new_density/scalingDensity() > 3. || new_density/scalingDensity()< 1./3 )
		{
			static bool lgWARN2BIG=false;
			if( !lgWARN2BIG )
			{
				lgWARN2BIG = true;
				fprintf(ioQQQ,"\n\n >========== Warning!  The tabulated or functional change in density as a function of depth was VERY large. This is zone %li.\n",nzone);
				fprintf(ioQQQ," >========== Warning!  This will cause convergence problems. \n");
				fprintf(ioQQQ," >========== Warning!  The current radius is %.3e. \n",radius.Radius);
				fprintf(ioQQQ," >========== Warning!  The current depth is %.3e. \n",radius.depth);
				fprintf(ioQQQ," >========== Warning!  Consider using a more modest change in density vs radius. \n\n\n");
			}
		}
	}

	else if( strcmp(dense.chDenseLaw,"CDEN") == 0 )
	{
		/* this is the default, constant density */
		new_density = scalingDensity();
	}

	else
	{
		fprintf( ioQQQ, " Unknown density law=%s= This is a major internal error.\n", 
		  dense.chDenseLaw );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	return new_density;
}

enum {
	CPRE,
	SUBSONIC,
	SUPERSONIC,
	STRONGD,
	ORIGINAL,
	SHOCK,
	ANTISHOCK,
	ANTISHOCK2
};

/** returns updated best guess for equilibrium scalingDensity */
STATIC double stepDensity(const PresMode & presmode, solverState &st);

STATIC void logPressureState()
{
	/* >>chng 04 feb 11, add option to remember current density and pressure */
	conv.hist_pres_density.push_back( dense.gas_phase[ipHYDROGEN] );
	conv.hist_pres_current.push_back( pressure.PresTotlCurr );
	conv.hist_pres_error.push_back( pressure.PresTotlError );
}

STATIC bool lgTestPressureConvergence( double new_density)
{
	double density_change_factor = new_density/scalingDensity();

	/* now see whether current pressure is within error bounds */
	if( density_change_factor > 1. + conv.PressureErrorAllowed || 
		 density_change_factor < 1. - conv.PressureErrorAllowed )
	{
		return false;
	}
	else
	{
		return true;
	}
}

STATIC double limitedDensityScaling( double new_density, double dP_chng_factor )
{
	double density_change_factor = new_density/scalingDensity()-1.;

	/* dP_chng_factor is initially 1, becomes smaller if sign of pressure change, changes */
	double pdelta = conv.MaxFractionalDensityStepPerIteration;

	/* make sure that change is not too extreme */
	density_change_factor *= dP_chng_factor;
	density_change_factor = MIN2(density_change_factor,+pdelta);
	density_change_factor = MAX2(density_change_factor,-pdelta);
	return 1.+density_change_factor;
}

/*PressureChange evaluate the current pressure, and change needed to
 * get it to PresTotlInit,
 * return value is true if density was changed, false if no changes were necessary */
void PressureChange( 
	/* this is change factor, 1 at first, becomes smaller as oscillations occur */
	double dP_chng_factor, const PresMode & presmode, solverState &st, bool &lgStable )
{
	DEBUG_ENTRY( "PressureChange()" );

	/* first evaluate total pressure for this location, and current conditions
	 * CurrentPressure is just sum of gas plus local line radiation pressure */
	/* this sets values of pressure.PresTotlCurr, also wind velocity for dynamics */
	PresTotCurrent();
	logPressureState();

	double old_density = scalingDensity();
	double old_pressure = pressure.PresTotlCurr;

	double new_density = stepDensity(presmode, st);

	conv.lgConvPres = lgTestPressureConvergence(new_density);

	double density_change_factor = limitedDensityScaling(new_density, dP_chng_factor);

	{
		/*@-redef@*/
		enum{DEBUG_LOC=false};
		static long int nsave=-1;
		/*@+redef@*/
		if( DEBUG_LOC /*&& nzone > 150 && iteration > 1*/ )
		{
			if( nsave-nzone ) fprintf(ioQQQ,"\n");
			nsave = nzone;
			fprintf(ioQQQ,"nnzzone\t%li\t%.2f%%\t%.3f\t%.2e\t%.2e\t%.2e\n", 
				nzone,
				density_change_factor,
				/* when this is negative we need to raise the density */
				pressure.PresTotlError*100.,
			    pressure.PresTotlCurr, 
				pressure.PresGasCurr,
				pressure.pres_radiation_lines_curr );
		}
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
		  " PressureChange called, changing HDEN from %10.3e to %10.3e Set fill fac to %10.3e\n", 
					dense.gas_phase[ipHYDROGEN], density_change_factor*dense.gas_phase[ipHYDROGEN], 
					geometry.FillFac );
	}

	bool lgChange = ( density_change_factor != 1. );

	lgStable = true;
	if( lgChange )
	{
		ScaleAllDensities((realnum) density_change_factor);

		/* must call TempChange since ionization has changed, there are some
		 * terms that affect collision rates (H0 term in electron collision) */
		TempChange(phycon.te , false);

		/* heating cooling balance while doing ionization,
		 * this is where the heavy lifting is done, this calls PresTotCurrent,
		 * which sets pressure.PresTotlCurr */
		ConvTempEdenIoniz();
		double new_pressure = pressure.PresTotlCurr;
		double dpdrho = (new_pressure-old_pressure)/(new_density-old_density);
		double dpdrhoScale = (new_pressure+old_pressure)/(new_density+old_density);
		if ( presmode.zone == SUBSONIC)
		{
			if ( dpdrho < -0.01*dpdrhoScale )
				lgStable = false;
		}
		else if ( presmode.zone == SUPERSONIC)
		{
			if ( dpdrho > 0.01*dpdrhoScale )
				lgStable = false;
		}
		if (0)
			fprintf(ioQQQ,"ConvPres STABLE??? nz %ld od %.2e nd %.2e eps %.2e op %.2e np %.2e dpdn %.2e cmp %.2e %c\n",
					  nzone,old_density,new_density,new_density/old_density-1.0,
					  old_pressure,new_pressure,dpdrho,dpdrhoScale,TorF(lgStable));
		if (0)
			if (!lgStable)
			{
				conv.lgConvPres = false; // If evidence of instability found, ensure at least another iteration
				fixit("Disabled test on stability"); 
				// log zone warning when pressure instability diagnosed
				// possibly reduce zone stepsize when close
			}
	}

	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && (nzone>215)/**/ )
		{
			fprintf( ioQQQ, 
				"%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%c\n", 
			  radius.depth, 
			  pressure.PresTotlCurr, 
			  pressure.PresTotlInit + pressure.PresInteg, 
			  pressure.PresTotlInit, 
			  pressure.PresGasCurr, 
			  pressure.PresRamCurr, 
			  pressure.pres_radiation_lines_curr, 
			  /* subtract continuum rad pres which has already been added on */
			  pressure.PresInteg - pressure.pinzon, 
			  wind.windv/1e5,
			  sqrt(5.*pressure.PresGasCurr/3./dense.xMassDensity)/1e5,
			  TorF(conv.lgConvPres) );
		}
	}
}

/* ============================================================================== */
/* DynaPresChngFactor, called from PressureChange to evaluate factor needed
 * to find new density needed for 
 * current conditions and wind solution, returns ratio of new to old density */

/* object is to get the local ram pressure
 * RamNeeded = pressure.PresTotlInit + pressure.PresGasCurr + pressure.PresInteg;
 * to equal the sum of the inital pressur at the illuminated face, the local gas pressure,
 * and the integrate radiative acceleration from the incident continuum 
 *
 * the local gas pressure is linear in density if the temperature is constant,
 *
 * the local ram pressure is inversely linear in density because of the relationship
 * between velocity and density introduced by the mass flux conservation 
 */

/* stepDensity finds a density which should move towards pressure equilibrium
 * sets pressure.PresTotlError */
STATIC double stepDensity(const PresMode &presmode, solverState &st)
{
	DEBUG_ENTRY( "stepDensity()" );

	double PresTotlCorrect = st.press;
	
	double densityCurrent = scalingDensity();

	double er = pressure.PresTotlCurr-PresTotlCorrect;
	pressure.PresTotlError = er/PresTotlCorrect;

	if( trace.lgTrace && presmode.global != CPRE )
	{
		fprintf( ioQQQ, 
					" DynaPresChngFactor set PresTotlCorrect=%.3e PresTotlInit=%.3e PresInteg=%.3e DivergePresInteg=%.3e\n", 
					PresTotlCorrect , pressure.PresTotlInit , pressure.PresInteg*pressure.lgContRadPresOn,
					dynamics.DivergePresInteg );
	}
	
	if( dynamics.lgTracePrint && presmode.global != CPRE)
		fprintf(ioQQQ,"Pressure: init %g +rad %g +diverge %g = %g cf %g\n",
				  pressure.PresTotlInit, pressure.PresInteg*pressure.lgContRadPresOn,
				  dynamics.DivergePresInteg, PresTotlCorrect, pressure.PresTotlCurr);
	/*fprintf(ioQQQ,"DEBUG\t%.2f\thden\t%.4e\tPtot\t%.4e\tPgas\t%.4e\n",
	  fnzone, scalingDensity(),pressure.PresTotlCurr,pressure.PresGasCurr );*/
	
	double rhohat, dpdrho;
	if( presmode.global == CPRE || presmode.global == ORIGINAL || 
		 st.lastzone != nzone || fabs(er-st.erp) < SMALLFLOAT
		 || fabs(densityCurrent-st.dp) < SMALLFLOAT) 
	{
		if ( presmode.zone == SUBSONIC ) 
		{
			dpdrho = pressure.PresTotlCurr/densityCurrent;
		}
		else
		{
			dpdrho = -pressure.PresTotlCurr/densityCurrent;
		}
		rhohat = densityCurrent - er/dpdrho;
	}
	else
	{
		/* second iteration on this zone, do linear fit to previous Pres - rho curve
		 * Linear approximation to guess root with two independent samples */
		dpdrho = (st.erp-er)/(st.dp-densityCurrent);
		rhohat = densityCurrent - er/dpdrho;

		/* Subsonic case: pressure ^ with density ^ => increase density further */
		/* Super "  case: pressure ^ with density v => decrease density further */
		
		/* Force the solution towards the required branch when it
		 * appears to have the "wrong" value of dP/drho */
		if(presmode.zone == SUBSONIC && dpdrho <= 0)
		{
			if (er < 0.0)
				rhohat = 1.03*densityCurrent;
			else
				rhohat = 0.97*densityCurrent;
		}
		else if(presmode.zone == SUPERSONIC && dpdrho >= 0)
		{
			if (er > 0.0)
				rhohat = 1.03*densityCurrent;
			else
				rhohat = 0.97*densityCurrent;
		}
	}
	
	st.dp = densityCurrent;
	st.erp = er;
	st.lastzone = nzone;
	
	if( presmode.global != CPRE && dynamics.lgTracePrint )
		fprintf(ioQQQ,"windv %li r %g v %g f %g\n",
				  nzone,radius.depth,wind.windv,DynaFlux(radius.depth));
	
	/* convergence trace at this level */
	if( trace.nTrConvg>=1  )
	{
		if ( presmode.global == CPRE )
		{
			fprintf( ioQQQ, 
						" PresChng %s density old %.3e new %.3e current %.3e correct %.3e dpdn %.3e\n",
						dynamics.chPresMode , densityCurrent, rhohat , pressure.PresTotlCurr,
						PresTotlCorrect, dpdrho
				);
		}
		else
		{
			fprintf( ioQQQ, 
						" DynaPresChngFactor mode %s new density %.3f vel %.3e\n",
						dynamics.chPresMode , rhohat , wind.windv );
		}
	}
	
	return rhohat;
}

void PresMode::set()
{
	DEBUG_ENTRY( "PresMode::set()" );
	/* dynamics.lgSetPresMode is flag to indicate sane value of dynamics.chPresMode.  
	 * If set true with SET DYNAMICS PRESSURE MODE
	 * then do not override that choice */
	if( !dynamics.lgSetPresMode )
	{
		/* above set true if pressure mode was set with
		 * SET DYNAMICS PRESSURE MODE - if we got here
		 * it was not set, and must make a guess */
		if(pressure.PresGasCurr < pressure.PresRamCurr)
			strcpy( dynamics.chPresMode , "supersonic" );
		else
			strcpy( dynamics.chPresMode , "subsonic" );
		/* clear the flag - pressure mode has been set */
		dynamics.lgSetPresMode = true;
	}
	
	/* if globally looking for transonic solution, then locally sometimes
	 * supersonic, sometimes subsonic - this branch sets global flag,
	 * which can also be set with SET DYNAMICS PRESSURE MODE.
	 * Under default conditions, ChPresMode was set in previous branch
	 * to sub or super sonic depending on current velocity on first time*/

	if (strcmp(dense.chDenseLaw,"CPRE") != 0)
	{
		ASSERT (strcmp(dense.chDenseLaw,"DYNA") == 0);
	}

	if (strcmp(dense.chDenseLaw,"CPRE") == 0)
	{
		global = CPRE;
		pressure.lgSonicPointAbortOK = false;
	}
	else if( strcmp( dynamics.chPresMode , "original" ) == 0 )
	{
		global = ORIGINAL;
		pressure.lgSonicPointAbortOK = true;
	}
	else if( strcmp( dynamics.chPresMode , "subsonic" ) == 0 )
	{
		global = SUBSONIC;
		pressure.lgSonicPointAbortOK = true;
	}
	/* supersonic */
	else if( strcmp( dynamics.chPresMode , "supersonic" ) == 0 )
	{
		global = SUPERSONIC;
		pressure.lgSonicPointAbortOK = true;
	}
	/* strong d */
	else if( strcmp( dynamics.chPresMode , "strongd" ) == 0 )
	{
		global = STRONGD;
		pressure.lgSonicPointAbortOK = false;
	}
	else if( strcmp( dynamics.chPresMode , "shock" ) == 0 )
	{
		global = SHOCK;
		pressure.lgSonicPointAbortOK = false;
	}
	else if( strcmp( dynamics.chPresMode , "antishock-by-mach" ) == 0 )
	{
		global = ANTISHOCK2;
		pressure.lgSonicPointAbortOK = false;
	}
	else if( strcmp( dynamics.chPresMode , "antishock" ) == 0 )
	{
		global = ANTISHOCK;
		pressure.lgSonicPointAbortOK = false;
	}
	
	/* this sets pressure mode for the current zone based on global mode
	 * and local conditions */
	if(global == CPRE)
	{
		zone = SUBSONIC;
	} 
	else if(global == ORIGINAL)
	{
		/* in this mode use comparison between ram and gas pressure to determine
		 * whether sub or super sonic */
		if(pressure.PresGasCurr > pressure.PresRamCurr)
		{
			zone = SUBSONIC;
		}
		else
		{
			zone = SUPERSONIC;
		}
	}
	else if(global == STRONGD)
	{
		//if(nzone <= 1)
		zone = SUPERSONIC;
	}
	else if(global == SUBSONIC)
	{
		zone = SUBSONIC;
	}
	else if(global == SUPERSONIC)
	{
		zone = SUPERSONIC;
	}
	else if(global == SHOCK)
	{
		if(radius.depth < dynamics.ShockDepth)
		{
			zone = SUBSONIC;
		}
		else
		{
			zone = SUPERSONIC;
		}
	}
	else if(global == ANTISHOCK)
	{
		if(radius.depth < dynamics.ShockDepth)
		{
			zone = SUPERSONIC;
		}
		else
		{
			zone = SUBSONIC;
		}
	}
	else if(global == ANTISHOCK2)
	{
		/* WJH 19 May 2004: This version of the antishock mode will
			insert the antishock at the point where the isothermal Mach
			number falls below a certain value, dynamics.ShockMach */
		if( fabs(wind.windv) > dynamics.ShockMach*sqrt(pressure.PresGasCurr/dense.xMassDensity))
		{
			zone = SUPERSONIC;
		}
		else
		{
			zone = SUBSONIC;
		}
	}
	else
	{
		printf("Need to set global pressure mode\n");
		cdEXIT( EXIT_FAILURE );
	}
}

double pressureZone(const PresMode &presmode)
{
	DEBUG_ENTRY( "pressureZone()" );
	double press;
	if( presmode.global == CPRE )
	{
		/* constant pressure */
		if( pressure.lgContRadPresOn )
		{
			/* >>chng 01 oct 31, replace pneed with CorretPressure */
			/* this has pressure due to incident continuum */
			press = pressure.PresTotlInit + pressure.PresInteg;
		}
		else
		{
			/* this does not have pressure due to incident continuum*/
			press = pressure.PresTotlInit*
				/* following term normally unity, power law set with option par on cmmnd*/
				pow(radius.Radius/radius.rinner,(double)pressure.PresPowerlaw);
		}

		if( dark.lgNFW_Set || pressure.gravity_symmetry >=0 )
		{
			fixit("check correct zone pressure for NFW");
			// doesn't this exclude current zone gravity pressure.RhoGravity?
			// and doesn't the above exclude current rad press, pressure.pinzon?
			// That would be a second order correction, which would need
			// a consistent second order treatment elsewhere to make it
			// worthwhile.
			press += pressure.IntegRhoGravity;
		}
	}
	else 
	{
		/* update the current desired pressure */
		press = pressure.PresTotlInit + pressure.PresInteg*pressure.lgContRadPresOn
			+ dynamics.DivergePresInteg;
	}
	return press;
}

