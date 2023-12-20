/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* DynaIterEnd called at end of iteration when advection is turned on */
/* DynaStartZone called at start of zone calculation when advection is turned on */
/* DynaEndZone called at end of zone calculation when advection is turned on */
/* DynaIonize, called from ionize to evaluate advective terms for current conditions */
/* DynaPresChngFactor, called from PressureChange to evaluate new density needed for
 * current conditions and wind solution, returns ratio of new to old density */
/* timestep_next - find next time step in time dependent case */
/* DynaZero zero some dynamics variables, called from zero.c */
/* DynaCreateArrays allocate some space needed to save the dynamics structure variables, 
 * called from DynaCreateArrays */
/* DynaPrtZone - called to print zone results */
/* DynaSave save info related to advection */
/* DynaSave, save output for dynamics solutions */
/* ParseDynaTime parse the time command, called from ParseCommands */
/* ParseDynaWind parse the wind command, called from ParseCommands */
/* CHANGES: (M. Salz 14.05.2013)
 *  - introduce the case table in ParseDynaWind(), which parse the
 *    tabulated values and sets the logic variables
 *  - changed p.getNumberCheck() to p.getNumberPlain() for wind
 *    velocity, because in the tabulated case there is no number in the
 *    first line
 *  - introduce the tabulated case in the flux calculation*/

#include "cddefines.h"
#include "cddrive.h"
#include "struc.h"
#include "colden.h"
#include "radius.h"
#include "stopcalc.h"
#include "hextra.h"
#include "iterations.h"
#include "conv.h"
#include "timesc.h"
#include "dense.h"
#include "mole.h"
#include "thermal.h"
#include "pressure.h"
#include "phycon.h"
#include "wind.h"
#include "iso.h"
#include "dynamics.h"
#include "cosmology.h"
#include "parser.h"
#include "rfield.h"
#include "container_classes.h"
#include "save.h"
#include "trace.h"

static const bool lgPrintDynamics = false;

t_dynamics dynamics;
static int ipUpstream=-1,iphUpstream=-1,ipyUpstream=-1;

bool t_dynamics::isInitialRelaxIteration( long int iteration )
{
	DEBUG_ENTRY( "t_dynamics::isInitialRelaxIteration()" );

	if( iteration <= n_initial_relax )
		return true;
	return false;
}

bool t_dynamics::doNonEquilibriumSolve( long int iteration )
{
	DEBUG_ENTRY( "t_dynamics::doNonEquilibriumSolve()" );

	if( iteration > n_initial_relax &&
		( lgAdvection || lgTimeDependentStatic ) &&
		Rate != 0.0 && !lgEquilibrium )
	{
		ASSERT(Rate > 0.0);
		return true;
	}

	return false;
}


/* 
 * >>chng 01 mar 16, incorporate advection within dynamical solutions
 * this file contains the routines that incorporate effects of dynamics and advection upon the
 * thermal and ionization solutions.  
 *
 * This code was originally developed in March 2001 during
 * a visit to UNAM Morelia, in collaboration with Robin Williams, Jane Arthur, and Will Henney.
 * Development was continued by email, and in a meeting in July/Aug 2001 at U Cardiff
 * Further development June 2002, UNAM Morelia (Cloudy and the Duendes Verdes)
 */

/* the interpolated upstream densities of all ionization stages,
 * [element][ion] */
static multi_arr<double,2> UpstreamIon;
static multi_arr<double,3> UpstreamStatesElem;
/* total abundance of each element per hydrogen */
static vector<double> UpstreamElem;

/* hydrogen molecules */
static vector<double> Upstream_molecules;

/* space to save continuum when time command is used 
static realnum *dyna_flux_save;*/

/* array of times and continuum values we will interpolate upon */
static vector<double> time_elapsed_time , 
	time_flux_ratio ,
	time_dt,
	time_dt_scale_factor;
static bool lgtime_dt_specified;
static vector<int> lgtime_Recom;
static const int NTIME = 200;

/* number of time steps actually read in */
static long int nTime_flux=0;

/* routine called at end of iteration to determine new step sizes */
STATIC void DynaNewStep(void);

/* routine called at end of iteration to save values in previous iteration */
STATIC void DynaSaveLast(void);

/* routine called to determine mass flux at given distance */
/* static realnum DynaFlux(double depth); */

/* lookahead distance, as derived in DynaIterEnd */
static double Dyn_dr;

/* advected work */
static double AdvecSpecificEnthalpy;

/* depth of position in previous iteration */
static vector<double> Old_depth/*[NZLIM]*/;

/* HI ionization structure from previous iteration */
static vector<realnum> Old_histr/*[NZLIM]*/ ,
	/* Lyman continuum optical depth from previous iteration */
	Old_xLyman_depth/*[NZLIM]*/,
	/* old n_p density from previous iteration */
	Old_hiistr/*[NZLIM]*/,
	/* old pressure from previous iteration */
	Old_pressure/*[NZLIM]*/,
	/* H density - particles per unit vol */
	Old_density/*[NZLIM]*/ ,
	/* density - total grams per unit vol */
	Old_DenMass/*[NZLIM]*/ ,
	/* sum of enthalpy and kinetic energy per gram */
	EnthalpyDensity/*[NZLIM]*/,
	/* old electron density from previous iteration */
	Old_ednstr/*[NZLIM]*/,
	/* sum of enthalpy and kinetic energy per gram */
	Old_EnthalpyDensity/*[NZLIM]*/;

static multi_arr<realnum,2> Old_molecules;

/* the ionization fractions from the previous iteration */
static multi_arr<realnum,3> Old_xIonDense;

/* the gas phase abundances from the previous iteration */
static multi_arr<realnum,2> Old_gas_phase;

/* the iso levels from the previous iteration */
static multi_arr<realnum,4> Old_StatesElem;

/* the number of zones that were saved in the previous iteration */
static long int nOld_zone;

STATIC bool lgNeedTimestep ();
STATIC void InitDynaTimestep();

STATIC void save_DT( double dtime, double this_fact, string &reason )
{
	bool lgDoPun = ( save.lgDTOn && !( save.lgDTPLst && !iterations.lgLastIt ) );
	if( (trace.lgTrace && trace.lgDrBug) || lgDoPun )
	{
		static bool lgPrtRelax = false;
		if( ! lgPrtRelax )
		{
			for( long iter = 1; iter <= dynamics.n_initial_relax; iter++ )
			{
				fprintf( save.ipDTout,
					"%.12e\t%.3e\t%.5e\t%.3e\tRELAX\n",
					0.0,
					0.0,
					0.0,
					0.0 );
			}
			lgPrtRelax = true;
		}

		double Temp;
		if( cdTemp(
			/* four char string, null terminzed, giving the element name */
			"hydr",
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number */
			2 ,
			/* will be temperature */
			&Temp,
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ) )
			TotalInsanity();

		fprintf( save.ipDTout,
			"%.12e\t%.3e\t%.5e\t%.3e\t%s\n",
			dynamics.time_elapsed,
			Temp,
			dtime,
			this_fact,
			reason.c_str() );

		/* create hash marks on second and later iterations */
		if( iteration > 1 && save.lgDTHash )
		{
			static bool lgPrtHash = false;
			if( lgPrtHash )
				fprintf( save.ipDTout, "%s\n",save.chHashString.c_str() );
			lgPrtHash = true;
		}
	}
}

/*timestep_next - find next time step in time dependent case */
STATIC double timestep_next( void )
{
	DEBUG_ENTRY( "timestep_next()" );

	static double te_old=-1;
	double timestep_return = dynamics.timestep;
	double factor = 0.0;
	string reason = "";

	if( dynamics.lgRecom )
	{
		double timestep_Hp_temp = -1.;
		double te_new;

		if( cdTemp(
			/* four char string, null terminzed, giving the element name */
			"hydr",
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number */
			2 ,
			/* will be temperature */
			&te_new,
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"VOLUME" ) )
			TotalInsanity();

		if( te_old>0 )
		{
			double dTdStep = fabs(te_new-te_old)/te_new;
			/* this is the change factor to get 0.1 frac change in mean temp */
			double dT_factor = 0.04/SDIV(dTdStep);
			if( dT_factor > 2.0 )
			{
				dT_factor = 2.0;
				reason = "DTMAX";
			}
			else if( dT_factor < 0.01 )
			{
				dT_factor = 0.01;
				reason = "DTMIN";
			}
			else
			{
				ostringstream oss;

				oss << "change in temperature weighted by H+,";
				oss << " old T = " << scientific << setprecision(3) << te_old;
				oss << " new T = " << scientific << setprecision(3) << te_new;

				reason = oss.str();
			}
			timestep_Hp_temp = dynamics.timestep * dT_factor;
			factor = dT_factor;
		}
		else
		{
			factor = 1.;
			reason = "initial";
		}
		te_old = te_new;
		if( timestep_Hp_temp > 0. )
			timestep_return = timestep_Hp_temp;
	}
	else
	{
		timestep_return = dynamics.timestep_init;
		factor = 1.;
		reason = "constant timestep";
	}

	if( lgPrintDynamics )
		fprintf(ioQQQ,"DEBUG timestep_next returns %.3e, old temp %.2e\n" , timestep_return, te_old );

	save_DT( timestep_return, factor, reason );

	return( timestep_return );
}

/* ============================================================================== */
/* DynaIonize, called from ConvBase to evaluate advective terms for current conditions,
 * calculates terms to add to ionization balance equation */
void DynaIonize(void)
{
	DEBUG_ENTRY( "DynaIonize()" );

	/* the time (s) needed for gas to move dynamics.Dyn_dr  */
	/* >>chng 02 dec 11 rjrw increase dynamics.dynamics.timestep when beyond end of previous zone, to allow -> eqm */
	if( !dynamics.lgTimeDependentStatic )
	{
		/* in time dependent model dynamics.timestep only changes once at end of iteration
		 * and is constant across a model */
		/* usual case - full dynamics */
		dynamics.timestep = -Dyn_dr/wind.windv;
	}
	/* printf("%d %g %g %g %g\n",nzone,radius.depth,Dyn_dr,radius.depth-Old_depth[nOld_zone-1],dynamics.timestep); */

	ASSERT(nzone<struc.nzlim );
	if( nzone > 0 )
		EnthalpyDensity[nzone-1] = (realnum)phycon.EnthalpyDensity;

	/* do nothing on first iteration or when looking beyond previous iteration */
	/* >>chng 05 jan 27, from hardwired "iteration < 2" to more general case,
	 * this is set with SET DYNAMICS RELAX command with the default of 2 */
	//For advection cases, switch to local equilibrium when depth is outside
	//the region of previous iteration.
	//Possibly should limit range of dynamical sources further by adding Dyn_dr???
	double depth = radius.depth; 
	if( dynamics.isInitialRelaxIteration( iteration ) ||
			( ! dynamics.lgTimeDependentStatic && 
				( depth < 0 || depth > dynamics.oldFullDepth ) ) )
	{
		/* first iteration, return zero */
		dynamics.Cool_r = 0.;
		dynamics.Heat_v = 0.;
		dynamics.dHeatdT = 0.;

		dynamics.Rate = 0.;

		for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			for( long ion=0; ion<nelem+2; ++ion )
			{
				dynamics.Source[nelem][ion] = 0.;
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( long nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso_sp[ipISO][nelem].numLevels_local; ++level )
					{
						dynamics.StatesElem[nelem][nelem-ipISO][level] = 0.;
					}
				}
			}
		}
		for( long mol=0;mol<mole_global.num_calc;mol++)
		{
			dynamics.molecules[mol] = 0.;
		}
		return;
	}

	if( dynamics.lgTracePrint )
	{
		fprintf(ioQQQ,"workwork\t%li\t%.3e\t%.3e\t%.3e\n",
			nzone,
			phycon.EnthalpyDensity,
			0.5*POW2(wind.windv)*dense.xMassDensity ,
			5./2.*pressure.PresGasCurr
			); /**/
	}

	/* net cooling due to advection */
	/* >>chng 01 aug 01, removed hden from dilution, put back into here */
	/* >>chng 01 sep 15, do not update this variable the very first time this routine
	 * is called at the new zone. */
	dynamics.Cool_r = 1./dynamics.timestep*dynamics.lgCoolHeat;
	dynamics.Heat_v = AdvecSpecificEnthalpy/dynamics.timestep*dynamics.lgCoolHeat;
	dynamics.dHeatdT = 0.*dynamics.lgCoolHeat;

#	if 0
	if( dynamics.lgTracePrint || nzone>17 && iteration == 10)
	{
		fprintf(ioQQQ,
			"dynamics cool-heat\t%li\t%.3e\t%i\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			nzone,
			phycon.te, 
			dynamics.lgCoolHeat,
			thermal.htot,
			phycon.EnthalpyDensity/dynamics.timestep,
			AdvecSpecificEnthalpy*scalingDensity()/dynamics.timestep,
			phycon.EnthalpyDensity,
			AdvecSpecificEnthalpy*scalingDensity(),
			scalingDensity(),
			dynamics.timestep);
	}
#	endif

	/* second or greater iteration, have advective terms */
	/* this will evaluate advective terms for current physical conditions */

	/* the rate (s^-1) atoms drift in from upstream, a source term for the ground */

	/* dynamics.Hatom/dynamics.timestep is the source (cm^-3 s^-1) of neutrals,
		 normalized to (s^-1) by the next higher ionization state as for all 
		 other recombination terms */

	/* dynamics.xIonDense[ipHYDROGEN][1]/dynamics.timestep is the sink (cm^-3 s^-1) of
		 ions, normalized to (s^-1) by the ionization state as for all other
		 ionization terms */

	dynamics.Rate = 1./dynamics.timestep;
	ASSERT(dynamics.timestep > 0);

	for( long mol=0;mol<mole_global.num_calc;mol++)
	{
		dynamics.molecules[mol] = Upstream_molecules[mol]*scalingDensity() / dynamics.timestep;
	}

	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* check that the total number of each element is conserved along the flow */
			if(fabs(UpstreamElem[nelem]*scalingDensity() -dense.gas_phase[nelem])/dense.gas_phase[nelem]>=1e-3) 
			{
				fprintf(ioQQQ,
					"PROBLEM conservation error: zn %li elem %li upstream %.8e abund %.8e (up-ab)/up %.2e\n",
					nzone ,
					nelem,
					UpstreamElem[nelem]*scalingDensity(),
					dense.gas_phase[nelem] ,
					(UpstreamElem[nelem]*scalingDensity()-dense.gas_phase[nelem]) /
						(UpstreamElem[nelem]*scalingDensity()) );
			}
			/* ASSERT( fabs(UpstreamElem[nelem]*scalingDensity() -dense.gas_phase[nelem])/dense.gas_phase[nelem]<1e-3); */
			for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
			{
				/* Normalize to next higher state in current zone, except at first iteration
				where upstream version may be a better estimate (required for
				convergence in the small dynamics.timestep limit) */

				dynamics.Source[nelem][ion] = 
					/* UpstreamIon is ion number per unit hydrogen because dilution is 1/hden 
					 * this forms the ratio of upstream atom over current ion, per dynamics.timestep,
					 * so Source has units cm-3 s-1 */
					UpstreamIon[nelem][ion] * scalingDensity() / dynamics.timestep;

			}
			for( long ion=0; ion<dense.IonLow[nelem]; ++ion )
			{
				dynamics.Source[nelem][ion] = 0.;
				dynamics.Source[nelem][dense.IonLow[nelem]] += 
					UpstreamIon[nelem][ion] * scalingDensity() / dynamics.timestep;				
			}
			for( long ion=dense.IonHigh[nelem]+1;ion<nelem+2; ++ion )
			{
				dynamics.Source[nelem][ion] = 0.;
				dynamics.Source[nelem][dense.IonHigh[nelem]] += 
					UpstreamIon[nelem][ion] * scalingDensity() / dynamics.timestep;
			}
		}
	}

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem<LIMELM; ++nelem)
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long level=0; level < iso_sp[ipISO][nelem].numLevels_local; ++level )
				{
					dynamics.StatesElem[nelem][nelem-ipISO][level] = 
						UpstreamStatesElem[nelem][nelem-ipISO][level] * scalingDensity() / dynamics.timestep;
				}
			}
		}
	}

#	if 0
	fprintf(ioQQQ,"dynamiccc\t%li\t%.2e\t%.2e\t%.2e\t%.2e\n",
		nzone,
		dynamics.Rate,
		dynamics.Source[ipHYDROGEN][0],
		dynamics.Rate,
		dynamics.Source[ipCARBON][3]);
#	endif
#if 0
	long nelem = ipCARBON;
	long ion = 3;
	/*if( nzone > 160 && iteration > 1 )*/
	fprintf(ioQQQ,"dynaionizeeee\t%li\t%i\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
		nzone, 
		ipUpstream,
		radius.depth ,
		Old_depth[ipUpstream],
		dense.xIonDense[nelem][ion], 
		UpstreamIon[nelem][ion]* scalingDensity(),
		Old_xIonDense[ipUpstream][nelem][ion] ,
		dense.xIonDense[nelem][ion+1], 
		UpstreamIon[nelem][ion+1]* scalingDensity(),
		Old_xIonDense[ipUpstream][nelem][ion+1] ,
		dynamics.timestep,
		dynamics.Source[nelem][ion]
		);
#endif
	if( dynamics.lgTracePrint )
	{
		fprintf(ioQQQ,"    DynaIonize, %4li photo=%.2e , H recom= %.2e \n",
			nzone,dynamics.Rate , dynamics.Source[0][0]  );
	}
	return;
}

/* ============================================================================== */
/* DynaStartZone called at start of zone calculation when advection is turned on */
void DynaStartZone(void)
{
	/* this routine is called at the start of a zone calculation, by ZoneStart:
	 *
	 * it sets deduced variables to zero if first iteration,
	 *
	 * if upstream depth is is outside the computed structure on previous iteration,
	 * return value at shielded face 
	 *
	 * Also calculates discretization_error, an estimate of the accuracy of the source terms.
	 *
	 * */

	/* this is index of upstream cell in structure stored from previous iteration */
	double upstream, dilution, dilutionleft, dilutionright, frac_next;

	/* Properties for cell half as far upstream, used to converge dynamics.timestep */
	double hupstream, hnextfrac=-BIGFLOAT, hion, hmol, hiso;

	/* Properties for cell at present depth, used to converge dynamics.timestep */
	double ynextfrac=-BIGFLOAT, yion, ymol, yiso;

	long int nelem , ion, mol;

	DEBUG_ENTRY( "DynaStartZone()" );

	/* do nothing on first iteration */
	if( iteration < 2 )
	{
		dynamics.Upstream_density = 0.;
		AdvecSpecificEnthalpy = 0.;
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			for( ion=0; ion<nelem+2; ++ion )
			{
				UpstreamIon[nelem][ion] = 0.;
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso_sp[ipISO][nelem].numLevels_max; ++level )
					{
						UpstreamStatesElem[nelem][nelem-ipISO][level] = 0.;
					}
				}
			}
		}
		/* hydrogen molecules */
		for(mol=0; mol<mole_global.num_calc;mol++) 
		{
			Upstream_molecules[mol] = 0;
		}

		ipUpstream = -1;
		iphUpstream = -1;
		ipyUpstream = -1;
		return;
	}

	/* radius.depth is distance from illuminated face of cloud to outer edge of
	 * current zone, which has thickness radius.drad */

	/* find where the down stream point is, in previous iteration: we
	 * are looking for point where gas in current cell was in previous
	 * iteration */

	/* true, both advection and wind solution */
	/* don't interpolate to the illuminated side of the first cell */
	upstream = MAX2(Old_depth[0] , radius.depth + Dyn_dr);
	hupstream = 0.5*(radius.depth + upstream);

	/* now locate upstream point in previous stored structure,
	 * will be at the same point or higher than we found previously */
	while( ipUpstream+1 < nOld_zone && Old_depth[ipUpstream+1] < upstream )
	{
		++ipUpstream;
	}
	ASSERT( ipUpstream+1 <= nOld_zone );

	/* above loop will give ipUpstream == nOld_zone-1 if computed structure has been overrun */

	if( (ipUpstream != -1) && (ipUpstream != nOld_zone-1)&& ((Old_depth[ipUpstream+1] - Old_depth[ipUpstream])> SMALLFLOAT) )
	{
		/* we have not overrun radius scale of previous iteration */
		ASSERT( Old_depth[ipUpstream] <= upstream && upstream <= Old_depth[ipUpstream+1] );

		frac_next = ( upstream - Old_depth[ipUpstream])/
			(Old_depth[ipUpstream+1] - Old_depth[ipUpstream]);

		if (! (frac_next >= 0.0 && frac_next <= 1.0))
		{
			fprintf(ioQQQ,"PROBLEM frac_next out of range %.16e %.16e %.16e %.16e\n",
					  Old_depth[ipUpstream],upstream,Old_depth[ipUpstream+1],frac_next);
		}

		ASSERT (frac_next >= 0.0 && frac_next <= 1.0);
		dynamics.Upstream_density = (realnum)(Old_density[ipUpstream] +
			(Old_density[ipUpstream+1] - Old_density[ipUpstream])*
			frac_next);
		dilutionleft = 1./Old_density[ipUpstream];
		dilutionright = 1./Old_density[ipUpstream+1];

		/* fractional changes in density from passive advection */
		/* >>chng 01 May 02 rjrw: use hden for dilution */
		/* >>chng 01 aug 01, remove hden here, put back into vars when used in DynaIonize */
		dilution = 1./dynamics.Upstream_density;

		// we have a mix of realnum and double here, so make sure double is used for the test...
		double lo = double(Old_EnthalpyDensity[ipUpstream])*dilutionleft,
			hi = double(Old_EnthalpyDensity[ipUpstream+1])*dilutionright;
		/* the advected enthalphy */
		AdvecSpecificEnthalpy = lo + (hi-lo)*frac_next;
		double x = AdvecSpecificEnthalpy;
		
		if( ! fp_bound(lo,x,hi) )
		{
			fprintf(ioQQQ,"PROBLEM interpolated enthalpy density is not within range  %.16e\t%.16e\t%.16e\t%e\t%e\n",
				lo, x, hi, (hi-x)/(hi-lo), (x-lo)/(hi-lo));
		}

		ASSERT( fp_bound(lo,x,hi) );

		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			UpstreamElem[nelem] = 0.;
			for( ion=0; ion<nelem+2; ++ion )
			{
				/* Easier to bring out the division by the old hydrogen
				 * density rather than putting in dilution and then
				 * looking for how dilution is defined.  The code is
				 * essentially equivalent, but slower */
				/* UpstreamIon is ion number per unit material (measured
				 * as defined in scalingdensity()), at the upstream
				 * position */
				UpstreamIon[nelem][ion] = 
					Old_xIonDense[ipUpstream][nelem][ion]/Old_density[ipUpstream] +
					(Old_xIonDense[ipUpstream+1][nelem][ion]/Old_density[ipUpstream+1] - 
					Old_xIonDense[ipUpstream][nelem][ion]/Old_density[ipUpstream])*
					frac_next;

				UpstreamElem[nelem] += UpstreamIon[nelem][ion];
			}
		}

		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso_sp[ipISO][nelem].numLevels_max; ++level )
					{
						UpstreamStatesElem[nelem][nelem-ipISO][level] = 
							Old_StatesElem[ipUpstream][nelem][nelem-ipISO][level]/Old_density[ipUpstream] +
							(Old_StatesElem[ipUpstream+1][nelem][nelem-ipISO][level]/Old_density[ipUpstream+1] - 
							 Old_StatesElem[ipUpstream][nelem][nelem-ipISO][level]/Old_density[ipUpstream])*
							frac_next;
						/* check for NaN */
						ASSERT( !isnan( UpstreamStatesElem[nelem][nelem-ipISO][level] ) );
					}
				}
			}
		}

		for(mol=0;mol<mole_global.num_calc;mol++)
		{
			Upstream_molecules[mol] = Old_molecules[ipUpstream][mol]/Old_density[ipUpstream] +
				(Old_molecules[ipUpstream+1][mol]/Old_density[ipUpstream+1] - 
				 Old_molecules[ipUpstream][mol]/Old_density[ipUpstream])*
				frac_next;
			/* Externally represented species are already counted in UpstreamElem */
			if(mole.species[mol].location == NULL && mole_global.list[mol]->isIsotopicTotalSpecies())
			{
				for(molecule::nNucsMap::iterator atom= mole_global.list[mol]->nNuclide.begin(); 
					 atom != mole_global.list[mol]->nNuclide.end(); ++atom)
				{
					UpstreamElem[atom->first->el()->Z-1] += 
						Upstream_molecules[mol]*atom->second;
				}
			}
		}
	}
	else
	{
		/* SPECIAL CASE - we have overrun the previous iteration's radius */
		long ipBound = ipUpstream;
		if (ipBound == -1)
			ipBound = 0;
		dynamics.Upstream_density = Old_density[ipBound];
		/* fractional changes in density from passive advection */
		dilution = 1./dynamics.Upstream_density;
		/* AdvecSpecificEnthalpy enters as heat term */
		AdvecSpecificEnthalpy = Old_EnthalpyDensity[ipBound]/Old_density[ipBound];
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			UpstreamElem[nelem] = 0.;
			for( ion=0; ion<nelem+2; ++ion )
			{
				/* UpstreamIon is ion number per unit hydrogen */
				UpstreamIon[nelem][ion] = 
					Old_xIonDense[ipBound][nelem][ion]/Old_density[ipBound];
				UpstreamElem[nelem] += UpstreamIon[nelem][ion];
			}
		}

		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso_sp[ipISO][nelem].numLevels_max; ++level )
					{
						UpstreamStatesElem[nelem][nelem-ipISO][level] = 
							Old_StatesElem[ipBound][nelem][nelem-ipISO][level]/Old_density[ipBound];
						/* check for NaN */
						ASSERT( !isnan( UpstreamStatesElem[nelem][nelem-ipISO][level] ) );
					}
				}
			}
		}

		for(mol=0;mol<mole_global.num_calc;mol++) 
		{
			Upstream_molecules[mol] = Old_molecules[ipBound][mol]/Old_density[ipBound];
			if(mole.species[mol].location == NULL && mole_global.list[mol]->isIsotopicTotalSpecies())
			{
				for(molecule::nNucsMap::iterator atom=mole_global.list[mol]->nNuclide.begin();
					 atom != mole_global.list[mol]->nNuclide.end(); ++atom)
				{
					UpstreamElem[atom->first->el()->Z-1] += 
						Upstream_molecules[mol]*atom->second;
				}
			}
		}
	}

	/* Repeat enough of the above for half-step and no-step to judge convergence:
	 * the result of this code is the increment of discretization_error */

	while( iphUpstream+1 < nOld_zone && Old_depth[iphUpstream+1] < hupstream )
	{
		++iphUpstream;
	}
	ASSERT( iphUpstream+1 <= nOld_zone );

	while( ipyUpstream+1 < nOld_zone && Old_depth[ipyUpstream+1] < radius.depth )
	{
		++ipyUpstream;
	}
	ASSERT( ipyUpstream+1 <= nOld_zone );

	dynamics.dRad = BIGFLOAT;

	if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && (Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
		hnextfrac = ( hupstream - Old_depth[iphUpstream])/
			(Old_depth[iphUpstream+1] - Old_depth[iphUpstream]);
	else
		hnextfrac = 0.;

	if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && (Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
		ynextfrac = ( radius.depth - Old_depth[ipyUpstream])/
			(Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]);
	else
		ynextfrac = 0.;

	// Shouldn't be jumping over large numbers of upstream cells
	if(ipUpstream != -1 && ipUpstream < nOld_zone-1)
		dynamics.dRad = MIN2(dynamics.dRad,
												 5*(Old_depth[ipUpstream+1] - Old_depth[ipUpstream]));
	

	// Value for scaling zonal changes to set zone width
	const double STEP_FACTOR=0.05;

	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		for( ion=0; ion<nelem+2; ++ion )
		{
			double f1;
			if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && (Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
			{
				yion = 
					Old_xIonDense[ipyUpstream][nelem][ion]/Old_density[ipyUpstream] +
					(Old_xIonDense[ipyUpstream+1][nelem][ion]/Old_density[ipyUpstream+1] - 
					 Old_xIonDense[ipyUpstream][nelem][ion]/Old_density[ipyUpstream])*
					ynextfrac;
			}
			else
			{		
				long ipBound = ipyUpstream;
				if (-1 == ipBound)
					ipBound = 0;
				yion = Old_xIonDense[ipBound][nelem][ion]/Old_density[ipBound];
			}
			if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && (Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
			{
				hion = 
					Old_xIonDense[iphUpstream][nelem][ion]/Old_density[iphUpstream] +
					(Old_xIonDense[iphUpstream+1][nelem][ion]/Old_density[iphUpstream+1] - 
					 Old_xIonDense[iphUpstream][nelem][ion]/Old_density[iphUpstream])*
					hnextfrac;
			}
			else
			{		
				long ipBound = iphUpstream;
				if (-1 == ipBound)
					ipBound = 0;
				hion = Old_xIonDense[ipBound][nelem][ion]/Old_density[ipBound];
			}

			/* the proposed thickness of the next zone */
			f1 = fabs(yion - UpstreamIon[nelem][ion] );
			if( f1 > SMALLFLOAT )
			{
				dynamics.dRad = MIN2(dynamics.dRad,STEP_FACTOR*fabs(Dyn_dr) * 
				/* don't pay attention to species with abundance relative to H less than 1e-6 */
					MAX2(yion + UpstreamIon[nelem][ion],1e-6 ) / f1);
			}

			/* Must be consistent with convergence_error below */
			/* this error is error due to the advection length not being zero - a finite
			 * advection length.  no need to bring convergence error to below
			 * discretization error.  when convergece error is lower than a fraction of this
			 * errror we reduce the advection length. */
			dynamics.discretization_error += POW2(yion-2.*hion+UpstreamIon[nelem][ion]); 
			dynamics.error_scale2 += POW2(UpstreamIon[nelem][ion]-yion);
		}
	}

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem<LIMELM; ++nelem)
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long level=0; level < iso_sp[ipISO][nelem].numLevels_max; ++level )
				{
					double f1;
					if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && 
						(Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
					{
						yiso = 
							Old_StatesElem[ipyUpstream][nelem][nelem-ipISO][level]/Old_density[ipyUpstream] +
							(Old_StatesElem[ipyUpstream+1][nelem][nelem-ipISO][level]/Old_density[ipyUpstream+1] - 
							 Old_StatesElem[ipyUpstream][nelem][nelem-ipISO][level]/Old_density[ipyUpstream])*
							ynextfrac;
					}
					else
					{		
						long ipBound = ipyUpstream;
						if (-1 == ipBound)
							ipBound = 0;
						yiso = Old_StatesElem[ipBound][nelem][nelem-ipISO][level]/Old_density[ipBound];
					}
					if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && 
						(Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
					{
						hiso = 
							Old_StatesElem[iphUpstream][nelem][nelem-ipISO][level]/Old_density[iphUpstream] +
							(Old_StatesElem[iphUpstream+1][nelem][nelem-ipISO][level]/Old_density[iphUpstream+1] - 
							 Old_StatesElem[iphUpstream][nelem][nelem-ipISO][level]/Old_density[iphUpstream])*
							hnextfrac;
					}
					else
					{		
						long ipBound = iphUpstream;
						if (-1 == ipBound)
							ipBound = 0;
						hiso = Old_StatesElem[ipBound][nelem][nelem-ipISO][level]/Old_density[ipBound];
					}

					/* the proposed thickness of the next zone */
					f1 = fabs(yiso - UpstreamStatesElem[nelem][nelem-ipISO][level] );
					if( f1 > SMALLFLOAT )
					{
						dynamics.dRad = MIN2(dynamics.dRad,fabs(Dyn_dr*STEP_FACTOR) * 
						/* don't pay attention to species with abundance relative to H less tghan 1e-6 */
						  MAX2(yiso + UpstreamStatesElem[nelem][nelem-ipISO][level],1e-6 ) / f1);
					}
					/* Must be consistent with convergence_error below */
					/* this error is error due to the advection length not being zero - a finite
					 * advection length.  no need to bring convergence error to below
					 * discretization error.  when convergece error is lower than a fraction of this
					 * error we reduce the advection length. */
					dynamics.discretization_error += POW2(yiso-2.*hiso+UpstreamStatesElem[nelem][nelem-ipISO][level]); 
					dynamics.error_scale2 += POW2(UpstreamStatesElem[nelem][nelem-ipISO][level]);
				}
			}
		}
	}

	for(mol=0; mol < mole_global.num_calc; mol++) 
	{
		double f1;
		if(ipyUpstream != -1 && ipyUpstream != nOld_zone-1 && (Old_depth[ipyUpstream+1] - Old_depth[ipyUpstream]>SMALLFLOAT))
		{
			ymol = 
				Old_molecules[ipyUpstream][mol]/Old_density[ipyUpstream] +
				(Old_molecules[ipyUpstream+1][mol]/Old_density[ipyUpstream+1] - 
				 Old_molecules[ipyUpstream][mol]/Old_density[ipyUpstream])*
				ynextfrac;
		}
		else
		{		
			long ipBound = ipyUpstream;
			if (-1 == ipBound)
				ipBound = 0;
			ymol = Old_molecules[ipBound][mol]/Old_density[ipBound];
		}
		if(iphUpstream != -1 && iphUpstream != nOld_zone-1 && (Old_depth[iphUpstream+1] - Old_depth[iphUpstream]>SMALLFLOAT))
		{
			hmol = 
				Old_molecules[iphUpstream][mol]/Old_density[iphUpstream] +
				(Old_molecules[iphUpstream+1][mol]/Old_density[iphUpstream+1] - 
				 Old_molecules[iphUpstream][mol]/Old_density[iphUpstream])*
				hnextfrac;
		}
		else
		{		
			long ipBound = iphUpstream;
			if (-1 == ipBound)
				ipBound = 0;
			hmol = Old_molecules[ipBound][mol]/Old_density[ipBound];
		}

		/* the proposed thickness of the next zone */
		f1 = fabs(ymol - Upstream_molecules[mol] );
		if( f1 > SMALLFLOAT )
		{
			dynamics.dRad = MIN2(dynamics.dRad,fabs(Dyn_dr*STEP_FACTOR) * 
				/* don't pay attention to species with abundance relative to H less than 1e-6 */
				MAX2(ymol + Upstream_molecules[mol],1e-6 ) / f1 );
		}

		/* Must be consistent with convergence_error below */
		/* >>chngf 01 aug 01, remove scalingDensity() from HAtom */
		/* >>chngf 02 aug 01, multiply by cell width */
		dynamics.discretization_error += POW2(ymol-2.*hmol+Upstream_molecules[mol]); 
		dynamics.error_scale2 += POW2(Upstream_molecules[mol]-ymol);
	}

	if( dynamics.lgTracePrint )
	{
		fprintf(ioQQQ," DynaStartZone, %4li photo=%.2e , H recom= %.2e dil %.2e \n",
			nzone,dynamics.Rate , dynamics.Source[0][0] , dilution*scalingDensity() );
	}
	return;
}

/* DynaEndZone called at end of zone calculation when advection is turned on */
void DynaEndZone(void)
{
	DEBUG_ENTRY( "DynaEndZone()" );

	/* this routine is called at the end of a zone calculation, by ZoneEnd */

	dynamics.DivergePresInteg += wind.windv*(DynaFlux(radius.depth)-DynaFlux(radius.depth-radius.drad)); 

	if(dynamics.lgTracePrint)
		fprintf(ioQQQ,"Check dp: %g %g mom %g %g mas %g\n",
			wind.windv*(DynaFlux(radius.depth)-DynaFlux(radius.depth-radius.drad)),
			2*wind.windv*DynaFlux(radius.depth)*radius.drad/(1e16-radius.depth),
			wind.windv*DynaFlux(radius.depth),
			wind.windv*DynaFlux(radius.depth)*(1e16-radius.depth)*(1e16-radius.depth),
			DynaFlux(radius.depth)*(1e16-radius.depth)*(1e16-radius.depth));
	return;
}


/* ============================================================================== */
/* DynaIterEnd called at end of iteration when advection is turned on */
void DynaIterEnd(void)
{
	/* this routine is called by IterRestart at the end of an iteration 
	 * when advection is turned on.  currently it only derives a 
	 * dynamics.timestep by looking at the spatial derivative of
	 * some stored quantities */
	long int i;
	static long int nTime_dt_array_element = 0;

	DEBUG_ENTRY( "DynaIterEnd()" );

	/* set stopping outer radius to current outer radius once we have let
	 * solution relax dynamics.n_initial_relax times
	 * note the off by one confusion
	 * want to do two static iterations then start dynamics 
	 * iteration was incremented before this call so iteration == 2 at
	 * end of first iteration */
	if( iteration == dynamics.n_initial_relax+1)
	{
		fprintf(ioQQQ,"DYNAMICS DynaIterEnd sets stop radius to %.2e after "
			"dynamics.n_initial_relax=%li iterations.\n",
			radius.depth,
			dynamics.n_initial_relax);
		for( i=iteration-1; i < iterations.iter_alloc; ++i )
			/* set stopping radius to current radius, this stops
			 * dynamical solutions from overrunning the upstream scale
			 * and extending to very large radius due to unphysical heat
			 * appearing to advect into region */
			iterations.StopThickness[i] = radius.depth;

		/* in the event time dependent cooling is required,
		 * but no timestep has been given,
		 * make a guess of the initial timestep based
		 * on the shortest cooling time in the model */
		if( lgNeedTimestep() )
		{
			InitDynaTimestep();
		}
	}

	dynamics.DivergePresInteg = 0.;

	/* This routine is only called if advection is turned on at the end of
	 * each iteration.  The test 
	 * is to check whether wind velocity is also set by dynamics code */

	/* !dynamics.lgStatic true - a true dynamical model */
	if( !dynamics.lgTimeDependentStatic )
	{
		if(iteration == dynamics.n_initial_relax+1 )
		{
			/* let model settle down for n_initial_relax iterations before we begin
			 * dynamics.n_initial_relax set with "set dynamics relax" 
			 * command - it gives the first iteration where we do dynamics 
			 * note the off by one confusion - relax is 2 by default,
			 * want to do two static iterations then start dynamics 
			 * iteration was incremented before this call so iteration == 2 
			 * at end of first iteration */
			if( dynamics.AdvecLengthInit> 0. )
			{
				Dyn_dr =  dynamics.AdvecLengthInit;
			}
			else
			{
				/* -ve dynamics.adveclengthlimit sets length as fraction of first iter total depth */
				Dyn_dr = -dynamics.AdvecLengthInit*radius.depth;
			}

			if (wind.windv0 > 0)
				Dyn_dr = -Dyn_dr;

			if( dynamics.lgTracePrint )
			{
				fprintf(ioQQQ," DynaIterEnd, dr=%.2e \n",
					Dyn_dr );
			}
		}
		else if( ! dynamics.isInitialRelaxIteration( iteration ) )
		{
			/* evaluate errors and update Dyn_dr */
			DynaNewStep();
		}
	}
	else
	{
		/* this is time-dependent static model */
		static double HeatInitial=-1. , HeatRadiated=-1. ,
			DensityInitial=-1;
		/* n_initial_relax is number of time-steady models before we start 
		 * to evolve, set with "set dynamics relax" command */
		Dyn_dr =  0.;
		if( lgPrintDynamics )
			fprintf(ioQQQ,
				"DEBUG times enter dynamics.timestep %.2e elapsed_time %.2e iteration %li relax %li \n",
				dynamics.timestep ,
				dynamics.time_elapsed,
				iteration , dynamics.n_initial_relax);
		if( ! dynamics.isInitialRelaxIteration( iteration ) )
		{
			/* evaluate errors */
			DynaNewStep();

			/* this is set true on CORONAL XXX TIME INIT command, says to use constant
			 * temperature for first n_initial_relax iterations, then let run free */
			if( dynamics.lg_coronal_time_init  )
			{
				thermal.lgTemperatureConstant = false;
				thermal.ConstTemp = 0.;
			}

			/* time variable branch, now without dynamics */
			/* elapsed time - don't increase dynamics.time_elapsed during first two
			 * two iterations since this sets static model */
			dynamics.time_elapsed += dynamics.timestep;
			/* dynamics.timestep_stop is -1 if not set with explicit stop time */
			if( dynamics.timestep_stop > 0 && dynamics.time_elapsed > dynamics.timestep_stop )
			{
				dynamics.lgStatic_completed = true;
			}

			/* stop lowest temperature time command */
			if( (phycon.te < StopCalc.TempLoStopIteration) ||
				(phycon.te > StopCalc.TempHiStopIteration ) )
				dynamics.lgStatic_completed = true;

			/* this is heat radiated, after correction for change of H density in constant
			 * pressure cloud */
			HeatRadiated += (thermal.ctot-dynamics.Cool()) * dynamics.timestep *
				(DensityInitial / scalingDensity());
		}
		else
		{
			/* this branch, during initial relaxation of solution */
			HeatInitial = 1.5 * pressure.PresGasCurr;
			HeatRadiated = 0.;
			DensityInitial = scalingDensity();
			if( lgPrintDynamics )
				fprintf(ioQQQ,"DEBUG relaxing times requested %li this is step %li\n", 
					dynamics.n_initial_relax , iteration);
		}
		if( lgPrintDynamics )
			fprintf(ioQQQ,"DEBUG heat conser HeatInitial=%.2e HeatRadiated=%.2e\n",
				HeatInitial , HeatRadiated );

		if( dynamics.lgStatic_completed )
		{
			// Do nothing but talk.
			if( lgPrintDynamics )
				fprintf(ioQQQ,"DEBUG lgtimes -- stop time reached.\n" );
		}
		/* at this point dynamics.time_elapsed is the time at the end of the 
		 * previous iteration.  We need dt for the next iteration */
		else if( time_elapsed_time.size() > 0 && dynamics.time_elapsed > time_elapsed_time[nTime_dt_array_element] )
		{
			/* time has advanced to next table point, 
			 * set dynamics.timestep to specified value */
			++nTime_dt_array_element;
			/* this is an assert since we already qualified the array
			 * element above */
			ASSERT( nTime_dt_array_element < nTime_flux );

			/* option to set flag for recombination logic */
			if( lgtime_Recom[nTime_dt_array_element] )
			{
				if( lgPrintDynamics )
					fprintf(ioQQQ,"DEBUG dynamics turn on recombination logic\n");
				dynamics.lgRecom = true;
				/* set largest possible zone thickness to value on previous
					* iteration when light was on - during recombination conditions
					* become much more homogeneous and dr can get too large,
					* crashing into H i-front */
				radius.sdrmax = radius.dr_max_last_iter/3.;
				radius.lgSdrmaxRel = false;
			}

			if( lgtime_dt_specified )
			{
				if( lgPrintDynamics )
					fprintf(ioQQQ,"DEBUG lgtimes increment Time to %li %.2e\n" ,nTime_dt_array_element,
						dynamics.timestep);
				/* this is the new dynamics.timestep */
				dynamics.timestep = time_dt[nTime_dt_array_element];
				/* option to change time step factor - default is 1.2 set in DynaZero */
				if( time_dt_scale_factor[nTime_dt_array_element] > 0. )
					dynamics.timestep_factor = time_dt_scale_factor[nTime_dt_array_element];
			}
		}
		else if( lgtime_dt_specified )
		{
			/* we are between two points in the table, increase dynamics.timestep */
			dynamics.timestep *= dynamics.timestep_factor;
			if( lgPrintDynamics )
				fprintf(ioQQQ,"DEBUG lgtimes increment Timeby dynamics.timestep_factor  to %li %.2e\n" ,
					nTime_dt_array_element,
					dynamics.timestep );
			if( time_elapsed_time.size() > 0 && dynamics.time_elapsed+dynamics.timestep > time_elapsed_time[nTime_dt_array_element] )
			{
				if( lgPrintDynamics )
					fprintf(ioQQQ,"DEBUG lgtimes but reset to %.2e\n" ,dynamics.timestep );
				dynamics.timestep = 1.0001*(time_elapsed_time[nTime_dt_array_element]-dynamics.time_elapsed);
			}
		}
		else
		{
			/* time not specified in third column, so use initial */
			dynamics.timestep = timestep_next();
		}

		if( cosmology.lgDo )
		{
			cosmology.redshift_step = -(1.f + cosmology.redshift_current ) *
				GetHubbleFactor( cosmology.redshift_current ) * (realnum)dynamics.timestep;
			cosmology.redshift_current += cosmology.redshift_step;
		}

		if( lgPrintDynamics )
			fprintf(ioQQQ,"DEBUG times exit dynamics.timestep %.2e elapsed_time %.2e scale %.2e ",
				dynamics.timestep ,
				dynamics.time_elapsed,
				rfield.time_continuum_scale);

		if( cosmology.lgDo )
		{
			fprintf(ioQQQ,"redshift %.3e ", cosmology.redshift_current );
		}

		fprintf(ioQQQ,"\n" );
	}

	/* Dyn_dr == 0 is for static time dependent cloud */
	ASSERT( dynamics.isInitialRelaxIteration( iteration ) ||
		Dyn_dr != 0. || (Dyn_dr == 0. && wind.lgStatic()) );

	/* reset the upstream counters */
	ipUpstream = iphUpstream = ipyUpstream = -1;
	dynamics.discretization_error = 0.;
	dynamics.error_scale2 = 0.;

	/* save results from previous iteration */
	DynaSaveLast();
	return;
}

/*DynaNewStep work out convergence errors */
STATIC void DynaNewStep(void)
{
	long int ilast = 0,
		i,
		nelem,
		ion,
		mol;

	double frac_next=-BIGFLOAT,
		Oldi_density,
		Oldi_ion,
		Oldi_iso,
		Oldi_mol;

	DEBUG_ENTRY( "DynaNewStep()" );

	/*n = MIN2(nzone, NZLIM-1);*/
	dynamics.convergence_error = 0;
	dynamics.error_scale1 = 0.;

	ASSERT( nzone < struc.nzlim);
	for(i=0;i<nzone;++i) 
	{
		/* Interpolate for present position in previous solution */
		while( ilast+1 < nOld_zone && Old_depth[ilast+1] < struc.depth[i] )
		{
			++ilast;
		}
		ASSERT( ilast+1 <= nOld_zone );

		if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
		{
			frac_next = ( struc.depth[i] - Old_depth[ilast])/
				(Old_depth[ilast+1] - Old_depth[ilast]);
			Oldi_density = Old_density[ilast] +
				(Old_density[ilast+1] - Old_density[ilast])*
				frac_next;
		}
		else
		{
			Oldi_density = Old_density[ilast];
		}
		/* Must be consistent with discretization_error above */
		/* >>chngf 02 aug 01, multiply by cell width */
		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			for( ion=0; ion<nelem+2; ++ion )
			{
				if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
				{
					Oldi_ion = (Old_xIonDense[ilast][nelem][ion] +
						(Old_xIonDense[ilast+1][nelem][ion]-Old_xIonDense[ilast][nelem][ion])*
						frac_next);
				}
				else
				{
					Oldi_ion = Old_xIonDense[ilast][nelem][ion];
				}
				dynamics.convergence_error += POW2((double)Oldi_ion/Oldi_density-struc.xIonDense[i][nelem][ion]/scalingZoneDensity(i)) /* *struc.dr[i] */;  

				/* >>chng 02 nov 11, add first error scale estimate from Robin */
				//fprintf(ioQQQ,"%g %g %g\n",dynamics.error_scale1,
				//				struc.xIonDense[nelem][ion][i],scalingZoneDensity(i));
				dynamics.error_scale1 += POW2((double)struc.xIonDense[i][nelem][ion]/scalingZoneDensity(i));			
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso_sp[ipISO][nelem].numLevels_local; ++level )
					{
						if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
						{
							Oldi_iso = (Old_StatesElem[ilast][nelem][nelem-ipISO][level] +
							  (Old_StatesElem[ilast+1][nelem][nelem-ipISO][level]-Old_StatesElem[ilast][nelem][nelem-ipISO][level])*
													frac_next);
						}
						else
						{
							Oldi_iso = Old_StatesElem[ilast][nelem][nelem-ipISO][level];
						}
						dynamics.convergence_error += POW2(Oldi_iso/Oldi_density-struc.StatesElem[i][nelem][nelem-ipISO][level]/scalingZoneDensity(i));	/* *struc.dr[i] */;  

						/* >>chng 02 nov 11, add first error scale estimate from Robin */
						dynamics.error_scale1 += POW2(struc.StatesElem[i][nelem][nelem-ipISO][level]/scalingZoneDensity(i));
					}
				}
			}
		}

		for( mol=0; mol < mole_global.num_calc; mol++)
		{
			if(ilast != nOld_zone-1 && ((Old_depth[ilast+1] - Old_depth[ilast])> SMALLFLOAT) )
			{
				Oldi_mol = (Old_molecules[ilast][mol] +
										(Old_molecules[ilast+1][mol]-Old_molecules[ilast][mol])*
										frac_next);
			}
			else
			{
				Oldi_mol = Old_molecules[ilast][mol];
			}
			dynamics.convergence_error += POW2((double)Oldi_mol/Oldi_density-struc.molecules[i][mol]/scalingZoneDensity(i)) /* *struc.dr[i] */;  

			/* >>chng 02 nov 11, add first error scale estimate from Robin
			 * used to normalize the above convergence_error */
			dynamics.error_scale1 += POW2((double)struc.molecules[i][mol]/scalingZoneDensity(i));			
		}
	}

	/* convergence_error is an estimate of the convergence of the solution from its change during the last iteration,
		 discretization_error is an estimate of the accuracy of the advective terms, calculated in DynaStartZone above:
		 if dominant error is from the advective terms, need to make them more accurate.
	*/

	/* report properties of previous iteration */
	fprintf(ioQQQ,"DYNAMICS DynaNewStep: Dyn_dr %.2e convergence_error %.2e discretization_error %.2e error_scale1 %.2e error_scale2 %.2e\n",
		Dyn_dr, dynamics.convergence_error , dynamics.discretization_error ,
		dynamics.error_scale1 , dynamics.error_scale2 
		);

	/* >>chng 02 nov 29, dynamics.convergence_tolerance is now set to 0.1 in init routine */
	if( dynamics.convergence_error < dynamics.convergence_tolerance*dynamics.discretization_error ) 
		Dyn_dr /= 1.5;
	return;
}

/*DynaSaveLast save results from previous iteration */
STATIC void DynaSaveLast(void)
{
	long int i,
		ion,
		nelem,
		mol;

	DEBUG_ENTRY( "DynaSaveLast()" );

	/* Save results from previous iteration */
	nOld_zone = nzone;
	dynamics.oldFullDepth = struc.depth[nzone-1];
	ASSERT( nzone < struc.nzlim );
	for( i=0; i<nzone; ++i )
	{
		Old_histr[i] = struc.histr[i];
		Old_depth[i] = struc.depth[i];
		Old_xLyman_depth[i] = struc.xLyman_depth[i];
		/* old n_p density from previous iteration */
		Old_hiistr[i] = struc.hiistr[i];
		/* old pressure from previous iteration */
		Old_pressure[i] = struc.pressure[i];
		/* old electron density from previous iteration */
		Old_ednstr[i] = struc.ednstr[i];
		/* energy term */
		Old_EnthalpyDensity[i] = EnthalpyDensity[i];
		Old_density[i] = scalingZoneDensity(i);
		Old_DenMass[i] = struc.DenMass[i];

		for(mol=0;mol<mole_global.num_calc;mol++)
		{
			Old_molecules[i][mol] = struc.molecules[i][mol];
		}

		for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
		{
			Old_gas_phase[i][nelem] = dense.gas_phase[nelem];
			for( ion=0; ion<nelem+2; ++ion )
			{
				Old_xIonDense[i][nelem][ion] = struc.xIonDense[i][nelem][ion];
			}
		}
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem<LIMELM; ++nelem)
			{
				if( dense.lgElmtOn[nelem] )
				{
					for( long level=0; level < iso_sp[ipISO][nelem].numLevels_max; ++level )
					{ 
						Old_StatesElem[i][nelem][nelem-ipISO][level] = struc.StatesElem[i][nelem][nelem-ipISO][level];
						ASSERT( !isnan( Old_StatesElem[i][nelem][nelem-ipISO][level] ) );
					}
				}
			}
		}		
	}
	return;
}

realnum DynaFlux(double depth)
{
	realnum flux;

	DEBUG_ENTRY( "DynaFlux()" );

	if(dynamics.FluxIndex == 0) 
	{
		flux = (realnum)dynamics.FluxScale;
	}	
	else 
	{
		flux = (realnum)(dynamics.FluxScale*pow(fabs(depth-dynamics.FluxCenter),dynamics.FluxIndex));
		if(depth < dynamics.FluxCenter)
			flux = -flux;
	}
	if(dynamics.lgFluxDScale) 
	{
		/*flux *= struc.DenMass[0]; */
		/* WJH 21 may 04, changed to use dense.xMassDensity0, which should be strictly constant */
		flux *= dense.xMassDensity0; 
	}

	if( wind.lgTabulated )
	{
	       flux = (realnum)wind.wlaw.tabval(radius.Radius, radius.depth);
               flux *= dense.xMassDensity;
	}
	return flux;
}

/* ============================================================================== */
/*DynaZero zero some dynamics variables, called from zero.c, 
 * before parsing commands */
void t_dynamics::zero( void )
{
	int ipISO;

	DEBUG_ENTRY( "t_dynamics::zero()" );

	/* the number of zones in the previous iteration */
	nOld_zone = 0;

	/* by default advection is turned off */
	lgAdvection = false;
	/*Velocity = 0.;*/
	AdvecSpecificEnthalpy = 0.;
	Cool_r = 0.;
	Heat_v = 0.;
	dHeatdT = 0.;
	HeatMax = 0.;
	CoolMax = 0.;
	Rate = 0.;

	/* sets recombination logic, keyword RECOMBINATION on a time step line */
	lgRecom = false;

	/* don't force populations to equilibrium levels */
	lgEquilibrium = false;

	/* set true if time dependent calculation is finished */
	lgStatic_completed = false;

	/* vars that determine whether time dependent soln only - set with time command */
	lgTimeDependentStatic = false;
	timestep_init = -1.;
	/* this factor multiplies the time step */
	timestep_factor = 1.2;
	time_elapsed = 0.;

	/* set the first iteration to include dynamics rather than constant pressure */
	/* iteration number, initial iteration is 1, default is 3 - changed with SET DYNAMICS FIRST command */
	n_initial_relax = 3;

	/* set initial value of the advection length,
	 * neg => fraction of depth of init model, + length cm */
	AdvecLengthInit = -0.1;

	/* this is a tolerance for determining whether dynamics has converged */
	convergence_tolerance = 0.1;

	/* this says that set dynamics pressure mode was set */
	lgSetPresMode = false;

	/* set default values for uniform mass flux */
	FluxScale = 0.;
	lgFluxDScale = false;
	FluxCenter = 0.;
	FluxIndex = 0.;
	dRad = BIGFLOAT;

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* factor to allow turning off advection for one of the iso seq, 
		 * this is done with command "no advection h-like" or "he-like" 
		 * only for testing */
		lgISO[ipISO] = true;
	}
	/* turn off advection for rest of ions, command "no advection metals" */
	lgMETALS = true;
	/* turn off thermal effects of advection, command "no advection cooling" */
	lgCoolHeat = true;
	DivergePresInteg = 0.;

	discretization_error = 0.;
	error_scale2 = 0.;
	Upstream_density = 0.;
	return;
}


/* ============================================================================== */
/* DynaCreateArrays allocate some space needed to save the dynamics structure variables, 
 * called from DynaCreateArrays */
void DynaCreateArrays()
{
	DEBUG_ENTRY( "DynaCreateArrays()" );

	Upstream_molecules.resize(mole_global.num_calc);
	dynamics.molecules.resize(mole_global.num_calc);

	UpstreamElem.resize(LIMELM);

	dynamics.Source.reserve(LIMELM);
	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
		dynamics.Source.reserve(nelem, nelem+2);
	dynamics.Source.alloc();
	UpstreamIon.alloc(dynamics.Source.clone());
	dynamics.Source = 0.;

	UpstreamStatesElem.reserve(LIMELM);
	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			UpstreamStatesElem.reserve(nelem, nelem+1);
			for( long ion=0; ion < nelem+1; ion++ )
			{
				long ipISO = nelem-ion;
				if( ipISO < NISO )
				{
					UpstreamStatesElem.reserve(nelem, nelem-ipISO, iso_sp[ipISO][nelem].numLevels_max);
				}
				else
				{
					fixit("for now, point non-iso ions to NULL");
				}
			}
		}
	}
	UpstreamStatesElem.alloc();
	dynamics.StatesElem.alloc(UpstreamStatesElem.clone());

	dynamics.Rate = 0.;

	Old_EnthalpyDensity.resize(struc.nzlim);
	Old_ednstr.resize(struc.nzlim);
	EnthalpyDensity.resize(struc.nzlim);
	Old_DenMass.resize(struc.nzlim);
	Old_density.resize(struc.nzlim);
	Old_pressure.resize(struc.nzlim);
	Old_histr.resize(struc.nzlim);
	Old_hiistr.resize(struc.nzlim);
	Old_depth.resize(struc.nzlim);
	Old_xLyman_depth.resize(struc.nzlim);

	Old_xIonDense.reserve(struc.nzlim);
	Old_StatesElem.reserve(struc.nzlim);
	Old_gas_phase.reserve(struc.nzlim);
	Old_molecules.reserve(struc.nzlim);
	/* now create diagonal of space for ionization arrays */
	for( long ns=0; ns < struc.nzlim; ++ns )
	{
		Old_xIonDense.reserve(ns, LIMELM);
		Old_StatesElem.reserve(ns, LIMELM);
		Old_gas_phase.reserve(ns, LIMELM);
		if (mole_global.num_calc != 0)
			Old_molecules.reserve(ns, mole_global.num_calc);
		for( long nelem=0; nelem < LIMELM; ++nelem )
		{
			Old_xIonDense.reserve(ns, nelem, LIMELM+1);
			if( dense.lgElmtOn[nelem] )
			{
				Old_StatesElem.reserve(ns, nelem, nelem+1);
				for( long ion=0; ion < nelem+1; ion++ )
				{
					long ipISO = nelem-ion;
					if( ipISO < NISO )
					{
						Old_StatesElem.reserve(ns, nelem, ion, iso_sp[ipISO][nelem].numLevels_max);
					}
					else
					{
						fixit("for now, point non-iso ions to NULL");
					}
				}
			}
		}
	}
	Old_xIonDense.alloc();
	Old_StatesElem.alloc();
	Old_gas_phase.alloc();
	Old_molecules.alloc();

	for( long i=0; i < struc.nzlim; i++ )
	{
		/* these are values if H0 and tau_912 from previous iteration */
		Old_histr[i] = 0.;
		Old_xLyman_depth[i] = 0.;
		Old_depth[i] = 0.;
		dynamics.oldFullDepth = 0.;
		/* old n_p density from previous iteration */
		Old_hiistr[i] = 0.;
		/* old pressure from previous iteration */
		Old_pressure[i] = 0.;
		/* old electron density from previous iteration */
		Old_ednstr[i] = 0.;
		Old_density[i] = 0.;
		Old_DenMass[i] = 0.;
		Old_EnthalpyDensity[i] = 0.;
	}
	Old_xIonDense = 0.;
	Old_StatesElem = 0.;
	Old_molecules = 0.;
	return;
}

/*advection_set_default - called to set default conditions
 * when time and wind commands are parsed,
 * lgWind is true if dynamics, false if time dependent */
STATIC void advection_set_default( bool lgWind )
{

	DEBUG_ENTRY( "advection_set_default()" );

	/* turn off prediction of next zone's temperature, as guessed in ZoneStart,
	 * also set with no tepredictor */
	thermal.lgPredNextTe = false;

	/* use the new temperature solver
	strcpy( conv.chSolverEden , "new" ); */

	/*  constant total pressure, gas+rad+incident continuum
	 *  turn on radiation pressure */
	pressure.lgPres_radiation_ON = true;
	pressure.lgPres_magnetic_ON = true;
	pressure.lgPres_ram_ON = true;

	/* we need to make the solvers much more exact when advection is in place */
	if( lgWind )
	{
		/* turn on advection */
		dynamics.lgAdvection = true;

		/* increase precision of solution */
		conv.EdenErrorAllowed = 1e-3;
		/* the actual relative error is relative to the total heating and cooling,
		 * which include the dynamics.heat() and .cool(), which are the advected heating/cooling.
		 * the two terms can be large and nearly cancel, what is written to the .heat() and cool files
		 * by save files has had the smaller of the two subtracted, leaving only the net advected 
		 * heating and cooling */
		conv.HeatCoolRelErrorAllowed = 3e-4f;
		conv.PressureErrorAllowed = 1e-3f;
	}

	return;
}

/* ============================================================================== */
/* lgNeedTimestep check if an initial timestep is required */
STATIC bool lgNeedTimestep ()
{
	DEBUG_ENTRY( "lgNeedTimestep()" );

	return  (dynamics.lgTimeDependentStatic && dynamics.lgRecom && dynamics.timestep_init <= 0.);
}

/* ============================================================================== */
/* InitDynaTimestep - initialize timestep based on shortest zone cooling time; used
 * with 'coronal init time' without 'time first timestep' */
STATIC void InitDynaTimestep()
{
	DEBUG_ENTRY( "InitDynaTimestep()" );

	/* set default flags - false says that time dependent, not dynamical solution */
	advection_set_default(false);

	wind.windv0 = 0.;
	wind.setStatic();
	wind.windv = wind.windv0;

	timesc.calc_therm_timesc( 1 );
	dynamics.timestep_init = dynamics.timestep = 1e-4 * timesc.time_therm_short;

	return;
}

/* ============================================================================== */
/* ParseDynaTime parse the time command, called from ParseCommands */
void ParseDynaTime( Parser &p )
{
	DEBUG_ENTRY( "ParseDynaTime()" );

	/*flag set true when time dependent only */
	dynamics.lgTimeDependentStatic = true;

	dynamics.timestep_init = p.getNumberCheckAlwaysLogLim("dynamics.timestep",30.);

	dynamics.timestep = dynamics.timestep_init;
	if( p.nMatch( "TRAC" ) )
		dynamics.lgTracePrint = true;

	/* this is the stop time and is optional */
	dynamics.timestep_stop = p.getNumberDefaultAlwaysLog("stop time", -1.);

	/* set default flags - false says that time dependent, not dynamical solution */
	advection_set_default(false);

	wind.windv0 = 0.;
	wind.setStatic();
	wind.windv = wind.windv0;

	/* create time step and flux arrays */
	time_elapsed_time.resize(NTIME);
	time_flux_ratio.resize(NTIME);
	time_dt.resize(NTIME);
	time_dt_scale_factor.resize(NTIME);
	lgtime_Recom.resize(NTIME);

	/* number of lines we will save */
	nTime_flux = 0;

	/* get the next line, and check for eof */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
			" Hit EOF while reading time-continuum list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* third column might set dt - if any third column is missing then
	 * this is set false and only time on command line is used */
	lgtime_dt_specified = true;

	while( ! p.hasCommand("END") )
	{
		if( nTime_flux >= NTIME )
		{
			fprintf( ioQQQ, 
				" Too many time points have been entered; the limit is %d.  Increase variable NTIME in dynamics.c.\n", 
			  NTIME );
			cdEXIT(EXIT_FAILURE);
		}

		if( p.nMatch("CYCLE") )
		{
			double period = p.getNumberCheckAlwaysLog("log time");
			ASSERT( period > time_elapsed_time[nTime_flux-1] );
			long pointsPerPeriod = nTime_flux;
			while( nTime_flux < NTIME - 1 )
			{
				time_elapsed_time[nTime_flux] = period + time_elapsed_time[nTime_flux-pointsPerPeriod];
				time_flux_ratio[nTime_flux] = time_flux_ratio[nTime_flux-pointsPerPeriod];
				time_dt[nTime_flux] = time_dt[nTime_flux-pointsPerPeriod];
				time_dt_scale_factor[nTime_flux] = time_dt_scale_factor[nTime_flux-pointsPerPeriod];
				nTime_flux++;
			}
			//Tell the code to continue cyclically by equating two named time points
			fprintf( ioQQQ, " Adding cycles with period = %e s.\n", period );

			/* get next line and check for eof */
			p.getline();
			if( p.m_lgEOF )
			{
				fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			continue;
		}

		time_elapsed_time[nTime_flux] = p.getNumberCheckAlwaysLog("log time");
		if( nTime_flux >= 1 )
			ASSERT( time_elapsed_time[nTime_flux] > time_elapsed_time[nTime_flux-1] );
		time_flux_ratio[nTime_flux] = p.getNumberCheckAlwaysLog("log flux ratio");

		/* this is optional dt to set time step - if not given then initial
		 * time step is always used */
		time_dt[nTime_flux] =  p.getNumberDefaultAlwaysLog("log time step",-1.);

		/* if any of these are not specified then do not use times array */
		if( time_dt[nTime_flux] < 0.0 )
			lgtime_dt_specified = false;

		/* this is optional scale factor to increase time */
		time_dt_scale_factor[nTime_flux] = p.getNumberDefaultAlwaysLog(
			"scale factor to increase time",-1.);

		/* turn on recombination front logic */
		if( p.nMatch("RECOMBIN") )
		{
			/* this sets flag dynamics.lgRecom true so that all of code knows recombination
			 * is taking place */
			lgtime_Recom[nTime_flux] = true;
		}
		else
		{
			lgtime_Recom[nTime_flux] = false;
		}

		/* this is total number stored so far */
		++nTime_flux;

		/* get next line and check for eof */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}

	}

	if( nTime_flux < 2 )
	{
		fprintf( ioQQQ, " At least two instances of time must be specified.  There is an implicit instance at t=0.\n"
			" The user must specify at least one additional time.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( lgPrintDynamics )
	{
		for( long i=0; i < nTime_flux; i++ )
		{
			fprintf( ioQQQ, "DEBUG time dep %.2e %.2e %.2e %.2e\n",
				time_elapsed_time[i],
				time_flux_ratio[i] ,
				time_dt[i],
				time_dt_scale_factor[i]);
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}
/* ============================================================================== */
/* ParseDynaWind parse the wind command, called from ParseCommands */
void ParseDynaWind( Parser &p)
{
	int iVelocity_Type;
	bool lgModeSet=false;
	/* compiler flagged possible paths where dfdr used but not set -
	 * this is for safety/keep it happy */
	double dfdr=-BIGDOUBLE;

	DEBUG_ENTRY( "ParseDynaWind()" );

	if( p.nMatch( "TRAC" ) )
		dynamics.lgTracePrint = true;

	/* Flag for type of velocity law: 
	 * 1 is original, give initial velocity at illuminated face
	 * 2 is face flux gradient (useful if face velocity is zero),
	 * set to zero, but will be reset if velocity specified */
	iVelocity_Type = 0;
	/* wind structure, parameters are initial velocity and optional mass
	 * v read in in km s-1 and convert to cm s-1, mass in solar masses */
	if( p.nMatch( "VELO" ) )
	{
		wind.windv0 = (realnum)(p.getNumberPlain("velocity")*1e5);
		wind.windv = wind.windv0;
		wind.setDefault();
		iVelocity_Type = 1;
	}

	if( p.nMatch( "BALL" ) )
	{
		wind.setBallistic();
		lgModeSet = true;
	}

	if( p.nMatch( "STAT" ) )
	{
		wind.windv0 = 0.;
		wind.setStatic();
		lgModeSet = true;
		iVelocity_Type = 1;
	}
	
 	if ( 1 == iVelocity_Type && !lgModeSet)
	{
		if (wind.windv0 > 0.)
		{
			fprintf(ioQQQ,"CAUTION, BALListic option needed to switch off pressure gradient terms\n");
		}
		else if (wind.windv0 == 0.)
		{
			fprintf(ioQQQ,"CAUTION, STATic option needed for zero speed solutions\n");
		}	
	}

	if( p.nMatch("DFDR") )
	{
		/* velocity not specified, rather mass flux gradient */
		dfdr = p.getNumberPlain("flux gradient");
		iVelocity_Type = 2;
	}

	/* center option, gives xxx */
	if( p.nMatch("CENT") )
	{
		/* physical length in cm, can be either sign */
		dynamics.FluxCenter = p.getNumberPlain(
			"centre of mass flux distribution");
	}

	/* flux index */
	if( p.nMatch("INDE") )
	{
		/* power law index */
		dynamics.FluxIndex = p.getNumberPlain(
			"power law index of mass flux distribution");
	}

	/* the case where velocity was set */
	if(iVelocity_Type == 1) 
	{
		/* was flux index also set? */
		if(dynamics.FluxIndex == 0)
		{
			dynamics.FluxScale = wind.windv0;
			dynamics.lgFluxDScale = true;
			/* Center doesn't mean much in this case -- make sure it's
			 * in front of grid so DynaFlux doesn't swap signs where
			 * it shouldn't */
			dynamics.FluxCenter = -1.;
		}
		else
		{
			/** \todo	2	Need to include divergence terms in pressure balance 
			 * if flux index is != 0 */
			/* velocity was set but flux index was not set - estimate it */
			dynamics.FluxScale = wind.windv0*
				pow(fabs(dynamics.FluxCenter),-dynamics.FluxIndex);

			dynamics.lgFluxDScale = true;
			if(dynamics.FluxCenter > 0)
			{
				dynamics.FluxScale = -dynamics.FluxScale;
			}
		}
	} 
	/* the case where flux gradient is set */
	else if(iVelocity_Type == 2)
	{
		if(dynamics.FluxIndex == 0)
		{
			fprintf(ioQQQ,"Can't specify gradient when flux is constant!\n");
			/* use this exit handler, which closes out MPI when multiprocessing */
			cdEXIT(EXIT_FAILURE);
		}
		/** \todo	2	Need to include divergence terms in pressure balance 
		 * if flux index is != 0 */
		/* Can't specify FluxScale from dvdr rather than dfdr, as
		 * d(rho)/dr != 0 */ 
		dynamics.FluxScale = dfdr/dynamics.FluxIndex*
			pow(fabs(dynamics.FluxCenter),1.-dynamics.FluxIndex);
		if(dynamics.FluxCenter > 0)
		{
			dynamics.FluxScale = -dynamics.FluxScale;
		}
		dynamics.lgFluxDScale = false;

		/* put in bogus value simply as flag -- assume that surface velocity
		 * is small or we wouldn't be using this to specify. */
		wind.windv0 = -0.01f;
		wind.setDefault();
	}
	else
	{
		/* assume old-style velocity-only specification */ 
		/* wind structure, parameters are initial velocity and optional mass
		 *  v in km/sec, mass in solar masses */
		wind.windv0 = (realnum)(p.getNumberPlain("wind velocity")*1e5);
		if (wind.windv0 < 0.) 
		{
			wind.setDefault();
		}
		else if (wind.windv0 > 0.)
		{
			wind.setBallistic();
		}
		else
		{
			wind.setStatic();
		}

		dynamics.FluxScale = wind.windv0;
		dynamics.FluxIndex = 0.;
		dynamics.lgFluxDScale = true;
		/* Center doesn't mean much in this case -- make sure it's
		 * in front of grid so DynaFlux doesn't swap signs where
		 * it shouldn't */
		dynamics.FluxCenter = -1.;
	}

	wind.windv = wind.windv0;

#	ifdef FOO
	fprintf(ioQQQ,"Scale %g (*%c) Index %g Center %g\n",
		dynamics.FluxScale,(dynamics.lgFluxDScale)?'D':'1',
		dynamics.FluxIndex,dynamics.FluxCenter);
#	endif

	/* option to include advection */
	if( p.nMatch( "ADVE" ) )
	{
		/* set default flags - true says dynamical solution */
		advection_set_default(true);
		strcpy( dense.chDenseLaw, "DYNA" );
	}

	else
	{
		/* this is usual hypersonic outflow */
		if( wind.windv0 <= 1.e6 )
		{
			/* speed of sound roughly 10 km/s */
			fprintf( ioQQQ, " >>>>Initial wind velocity should be greater than speed of sound; calculation only valid above sonic point.\n" );
			wind.lgWindOK = false;
		}

		/* set the central object mass, in solar masses */
		wind.comass = (realnum)p.getNumberDefault("central object mass",1.);
		/* default is one solar mass */

		/* option for rotating disk, keyword is disk */
		wind.lgDisk = false;
		if( p.nMatch( "DISK") )
			wind.lgDisk = true;

		strcpy( dense.chDenseLaw, "WIND" );
	}


	/*  option to turn off continuum radiative acceleration */
	if( p.nMatch("NO CO") )
	{
		pressure.lgContRadPresOn = false;
	}
	else
	{
		pressure.lgContRadPresOn = true;
	}

	/* read in velocity table at end, not to interfere with other options */
	if( p.nMatch("TABL") )
	{
               /* velocity is tabulated */
               wind.lgTabulated = true;
	       
               /* parse the tabulated values  */
	       p.readLaw(wind.wlaw);

               /* set initial velocity and static + ballistic = false */
               wind.windv0 = wind.wlaw.val[0];
               wind.windv = wind.windv0;
               wind.setDefault();

               dynamics.FluxScale = wind.windv0;
               dynamics.lgFluxDScale = false;
               /* Center doesn't mean much in this case -- make sure it's
                * in front of grid so DynaFlux doesn't swap signs where
                * it shouldn't */
               dynamics.FluxCenter = -1.;

               /* do not solve for pressure, so it is not important */
               pressure.lgPres_radiation_ON = false;
               pressure.lgContRadPresOn = false;
               pressure.lgPres_magnetic_ON = false;
               pressure.lgPres_ram_ON = true;
               
               /* do not include advective/expansion cooling or heating
                * --> must be solved outside of Cloudy */
               dynamics.lgCoolHeat = 0;

               /* tabulated velocity can only be used with tabulated density 
                * correct dynamics, except of advection must be solve outside
                * of Cloudy */
               strcpy( dense.chDenseLaw, "DLW2" );
       }
       else
       {
               /* not a tabulated velocity law */
               wind.lgTabulated = false;
       }

	return;
}

/*DynaPrtZone called to print zone results */
void DynaPrtZone( void )
{

	DEBUG_ENTRY( "DynaPrtZone()" );

	ASSERT( nzone>0 && nzone<struc.nzlim );

	if( nzone > 0 )
	{
		fprintf(ioQQQ," DYNAMICS Advection: Uad %.2f Uwd%.2e FRCcool: %4.2f Heat %4.2f\n",
			timesc.sound_speed_adiabatic/1e5 ,
			wind.windv/1e5 ,
			dynamics.Cool()/thermal.ctot,
			dynamics.Heat()/thermal.ctot);
	}

	ASSERT( EnthalpyDensity[nzone-1] > 0. );

	fprintf(ioQQQ," DYNAMICS Eexcit:%.4e Eion:%.4e Ebin:%.4e Ekin:%.4e ET+pdv %.4e EnthalpyDensity/rho%.4e AdvSpWork%.4e\n",
		phycon.EnergyExcitation,
		phycon.EnergyIonization,
		phycon.EnergyBinding,
  		0.5*POW2(wind.windv)*dense.xMassDensity,
		5./2.*pressure.PresGasCurr ,
		EnthalpyDensity[nzone-1]/scalingDensity() , AdvecSpecificEnthalpy
	);
	return;
}

/*DynaPunchTimeDep - save info about time dependent solution */
void DynaPunchTimeDep( FILE* ipPnunit , const char *chJob )
{

	DEBUG_ENTRY( "DynaPunchTimeDep()" );

	if( strcmp( chJob ,  "END" ) == 0 )
	{
		double te_mean,
			H2_mean,
			H0_mean,
			Hp_mean,
			Hep_mean;
		/* save info at end */
		if( cdTemp(
			/* four char string, null terminated, giving the element name */
			"HYDR", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 means that chLabel is a special case */
			2, 
			/* will be temperature */
			&te_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"HYDR", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			2, 
			/* will be fractional ionization */
			&Hp_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"HYDR", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			1, 
			/* will be fractional ionization */
			&H0_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"H2  ", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			0, 
			/* will be fractional ionization */
			&H2_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		if( cdIonFrac(
			/* four char string, null terminated, giving the element name */
			"HELI", 
			/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
			 * 0 says special case */
			2, 
			/* will be fractional ionization */
			&Hep_mean, 
			/* how to weight the average, must be "VOLUME" or "RADIUS" */
			"RADIUS" ,
			/* if true then weighting also has electron density, if false then only volume or radius */
			false ) )
		{
			TotalInsanity();
		}
		fprintf( ipPnunit , 
			"%.12e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n" , 
			dynamics.time_elapsed , 
			dynamics.timestep ,
			rfield.time_continuum_scale , 
			scalingDensity(),
			te_mean , 
			Hp_mean , 
			H0_mean , 
			H2_mean , 
			Hep_mean ,
			/* ratio of CO to total H column densities */
			findspecieslocal("CO")->column / SDIV( colden.colden[ipCOL_HTOT] ),
			cosmology.redshift_current,
			dense.eden/dense.gas_phase[ipHYDROGEN]
			);
	}
	else
		TotalInsanity();
	return;
}

/*DynaSave save dynamics - info related to advection */
void DynaSave(FILE* ipPnunit , char chJob )
{
	DEBUG_ENTRY( "DynaSave()" );

	if( chJob=='a' )
	{
		/* this is save dynamics advection, the only save dynamics */
		fprintf( ipPnunit , "%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			radius.depth_mid_zone,
			thermal.htot , 
			dynamics.Cool() , 
			dynamics.Heat() , 
			dynamics.dCooldT() ,
			dynamics.Source[ipHYDROGEN][ipHYDROGEN],
			dynamics.Rate,
			phycon.EnthalpyDensity/scalingDensity() ,
			AdvecSpecificEnthalpy
			);
	}
	else
		TotalInsanity();
	return;
}

#define MERGE 0
double t_dynamics::Heat()
{
	double heat = Heat_v*scalingDensity();
	if (MERGE)
	{
		double cool = Cool_r*phycon.EnthalpyDensity;
		if (heat > cool)
			return heat-cool;
		else
			return 0.;
	}
	return heat;
}

double t_dynamics::Cool()
{
	double cool = Cool_r*phycon.EnthalpyDensity;
	if (MERGE)
	{
		double heat = Heat_v*scalingDensity();
		if (heat < cool)
			return cool-heat;
		else
			return 0.;
	}
	return cool;
}
#undef MERGE

double t_dynamics::dCooldT()
{
	return Cool_r*5./2.*pressure.PresGasCurr/phycon.te;
}

void DynaIterStart(void)
{
	DEBUG_ENTRY( "DynaIterStart()" );

	if( 0 == nTime_flux || dynamics.isInitialRelaxIteration( iteration ) )
	{
		rfield.time_continuum_scale = 1.;
		return;
	}
	else if( dynamics.time_elapsed <= time_elapsed_time[0] )
	{
		/* if very early times not specified assume no flux variation yet */
		rfield.time_continuum_scale = (realnum)time_flux_ratio[0];
	}
	else if( dynamics.time_elapsed > time_elapsed_time[nTime_flux-1] )
	{
		if( dynamics.lgStatic_completed )
			rfield.time_continuum_scale = (realnum)time_flux_ratio[nTime_flux-1];
		else
		{
			fprintf( ioQQQ, " PROBLEM - DynaIterStart - I need the continuum at time %.2e but the table ends at %.2e.\n" ,
				dynamics.time_elapsed, time_elapsed_time[nTime_flux-1]);
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		rfield.time_continuum_scale = (realnum)linint(
			/* the times in seconds */
			&time_elapsed_time[0], 
			/* the rfield.time_continuum_scale factors */
			&time_flux_ratio[0], 
			/* the number of rfield.time_continuum_scale factors */
			nTime_flux,
			/* the desired time */
			dynamics.time_elapsed);
	}
	
	if( lgPrintDynamics )
		fprintf(ioQQQ,"DEBUG time dep reset continuum iter %ld dynamics.timestep %.2e elapsed time %.2e scale %.2e",
			  iteration,
			  dynamics.timestep ,
			  dynamics.time_elapsed, 
			  rfield.time_continuum_scale);
	if( dynamics.lgRecom )
	{
		fprintf(ioQQQ," recom");
	}
	fprintf(ioQQQ,"\n");
	
	/* make sure that at least one continuum source is variable */
	long int nTimeVary = 0;
	for( long int i=0; i < rfield.nShape; i++ )
	{
		/* this is set true if particular continuum source can vary with time 
		 * set true if TIME appears on intensity / luminosity command line */
		if( rfield.lgTimeVary[i] )
			++nTimeVary;
	}
	if( nTimeVary==0 )
		fprintf(ioQQQ," !No continuum sources vary with time. Was this intended?");
	
	
	if( hextra.lgTurbHeatVaryTime )
	{
		/* vary extra heating */
		hextra.TurbHeat = hextra.TurbHeatSave * rfield.time_continuum_scale;
		if( lgPrintDynamics )
			fprintf(ioQQQ,"DEBUG TurbHeat vary new heat %.2e\n",
					  hextra.TurbHeat);
	}
}
