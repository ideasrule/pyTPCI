/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include "module.h"
#include "container_classes.h"

/** routine called at start of iteration when advection is turned on */
void DynaIterStart(void);

/** routine called at end of iteration when advection is turned on */
void DynaIterEnd(void);

/** DynaStartZone called at start of iteration when advection is turned on */
void DynaStartZone(void);

/** DynaEndZone called at end of iteration when advection is turned on */
void DynaEndZone(void);

/** DynaIonize, called from ionize to evaluate advective terms for current conditions */
void DynaIonize(void);

/**DynaCreateArrays allocate some space needed to save the dynamics structure variables, called from atmdat_readin */
void DynaCreateArrays( void );

/** ParseDynaWind parse the wind command, called from ParseCommands 
\param *chCard
*/
class Parser;
void ParseDynaWind( Parser &p );

/** ParseDynaTime parse the time command, called from ParseCommands, in dynamics.c 
\param *chCard
*/
void ParseDynaTime( Parser &p );

/** DynaPrtZone - called to print zone results */
void DynaPrtZone( void );

/**DynaSave save info related to advection 
\param ipPnunit
\param chJob
*/
void DynaSave(FILE* ipPnunit , char chJob );

/**DynaPunchTimeDep - save info about time dependent solution 
\param ipPnunit
\param *chJob
*/
void DynaPunchTimeDep( FILE* ipPnunit , const char *chJob );

realnum DynaFlux(double depth);

/** all of these are initialized in zero */
struct t_dynamics : public module
{
	const char *chName() const
	{
		return "dynamics";
	}
	void zero();
	void comment(t_warnings&) {}

	/** @brief isInitialRelaxIteration - Check if initial iteration
	 * 	 
	 *  At the beginning of a dynamical run, a number of steady-state
	 *  iterations are performed to relax the solution before any
	 *  dynamical solves are done.  This method checks if the given
	 *  @p iteration is one of the initial relaxation iterations. 
	 *
	 * @param[in] iteration		iteration number
	 * 
	 * @return @c bool		true, if initial iteration; else false
	 */
	bool isInitialRelaxIteration( long int iteration );

	bool doNonEquilibriumSolve( long int iteration );


	/** is advection turned on ?, set to false in zero */
	bool lgAdvection;

	/** advective cooling minus heating (erg cm^-3 s^-1) */
	double Cool_r, Heat_v, dHeatdT;

	double Cool(), Heat(), dCooldT();

	/** largest fraction of cooling and heating */
	double CoolMax, HeatMax;

	/** the advection rate (s^-1) */
	double Rate;

	/** the advective recombination rate (cm^-3 s^-1)*/
	multi_arr<double,2> Source;

	/** the advective isolevel balance terms */
	multi_arr<double,3> StatesElem;

	/** save molecular network densities */
	vector<double> molecules;

	/** factor to turn off advection for H-like or He-like iso seq, no advection h-like, he-like*/
	bool lgISO[NISO];

	/** factor to turn off advection for rest of ions,  "no advection metals"*/
	bool lgMETALS;

	/** factor to turn off advective cooling */
	bool lgCoolHeat;

	/** var set with coronal time init - says to use constant temperature on first
	 * relax iterations then let temp run free */
	bool lg_coronal_time_init;

	/** set true if time dependent static simulation */
	bool lgTimeDependentStatic;

	/** elapsed time in time dependent static model */
	double time_elapsed;

	/** true if recombination logic in place */
	bool lgRecom;

	/** true if model ends since time dependent model is finished */
	bool lgStatic_completed;

	/** the initial value of the advection length, reset with set advection length */
	double AdvecLengthInit;

 	/** the center of the particle flux law */
	double FluxCenter;

	/** flag set by the "set dynamics pressure mode" command */
	char chPresMode[20];

	/** the shock depth in cm set with "set dynamics shock depth" command */
	double ShockDepth;

	/** the isothermal Mach number at which to insert an antishock
	   set with "set dynamics antishock Mach" command */
	double ShockMach;

	/** set how many iterations we will start with, before allowing
	 * changes.  This allows the solution to relax to an equilibrium
	 * set with "set dynamics relax" command */
	long int n_initial_relax;

 	/** the scale of the particle flux law */
	double FluxScale;

	/** whether we also need to scale by the face density */
	bool lgFluxDScale;

 	/** the power law index of the particle flux law */
	double FluxIndex;

	/** the proposed thickness for the next zone when advection is included */
	double dRad;

	/** the depth of the last iteration */
	double oldFullDepth;

	/** convergence_error and discretization_error give estimates of convergence: 
	 :: convergence_error -- change between the last iterations;
	 :: discretization_error -- error in the upstream interpolation.
	 When (and if) discretization_error >> convergence_error, the interpolation length should be decreased.

	 They should both be based on the same norm of the models, but what norm may be
	 experimented with -- at present, it's H+/Htot just weighted by cell number, 
	 which makes the estimates sensitive to the structure of the primary ionization front.
	 */
	/* the error from comparing this iteration with the previous one */
	double convergence_error;

	/** the allowed rel error, by default 0.1 */
	double convergence_tolerance;

	/** the error to be expected based on the coarseness of current advection length */
	double discretization_error;

	/** two ways of scaling the error estimate for convergence */
	double error_scale1 , error_scale2;

	/** flag set true if set dynamics flow type was set - this means
	 * to use the specified option, and not to derive one */
	bool lgSetPresMode;

	/** Upstream density */
	realnum Upstream_density;

	realnum DivergePresInteg;

	/** Enforce equilibrium populations */
	bool lgEquilibrium;

	/* set true with trace option on time command */
	bool lgTracePrint;

	 /* initial timestep (seconds) read off time command,
	  * each iteration accounts for this much time */
	double timestep_init,
		 timestep,
		 timestep_stop,
		 timestep_factor;
};
extern t_dynamics dynamics;

#endif /* DYNAMICS_H_ */
