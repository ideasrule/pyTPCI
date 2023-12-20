/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef STRUC_H_
#define STRUC_H_

#include "container_classes.h"

/**
 * these save EXTERN structure variables, like te vs depth
 * zeroed out in zero, incremented in RT_tau_inc
 */

struct t_struc {

	/** this is the new variable that replaces the old NZLIM,
	 * it gives the largest number of zones that will be needed
	 * in the current calculation, and is used to create
	 * space for the following variables, and also similar
	 * variables in dynamics.c */
	long int nzlim;

	/** nzone from previous iteration, -1 on first iteration */
	long int nzonePreviousIteration;

	/** these will all become vectors with length given by the limit to the number of zones */
	vector<realnum> testr, 
	  volstr, 
	  drad_x_fillfac, 
	  histr, 
	  hiistr, 
	  ednstr, 
	  o3str, 
	  /** the total pressure, including all terms (ram, radiation, integrated incident, gas)*/
	  pressure,
	  /** wind velocity, cm s-1 */
	  windv,
	  /** total outward acceleration cm s^-2 */
	  AccelTotalOutward,
	  /** inward acceleration due to gravity, a positive number */
	  AccelGravity,
	  /** just the gas pressure, nkT */
	  GasPressure,
	  /** line radiation pressure */
	  pres_radiation_lines_curr,
	  /** >>chng 02 May 2001 rjrw: add hden for dilution */
	  hden,
	  /** total particles per unit vol */
	  DenParticles,
	  /** density, gm/cm3 total grams per unit vol */
	  DenMass,
	  /** depth of this position */
	  depth,
	  /** the thickness of the current zone */
	  drad, 
	  /** Lyman continuum optical depth for current iteration */
	  xLyman_depth,
	  /** previous iteration's radius and dr scale */
	  depth_last,
	  drad_last;

	/** largest relative change in temperature in any zone between iterations */
	realnum TempChangeMax;

	/** save ionization balance array across model */
	multi_arr<realnum,3> xIonDense;

	/** save iso level array across model */
	multi_arr<realnum,4> StatesElem;

	/** the hydrogen molecules */
	multi_arr<realnum,2> molecules;
	vector<realnum>	H2_abund;

	/** total gas phase abundances */
	multi_arr<realnum,2> gas_phase;

	/** cooling and heating for each zone */
	vector<double> coolstr,
	  heatstr;

	/** this is the relative ionization that is the limit for choosing 
	 * zones using it, and for detecting it in prt_comment,
	 * default is 1e-3 */
	realnum dr_ionfrac_limit;

	t_struc()
	{
		/* limit on ionization we check for zoning and prtcomment */
		dr_ionfrac_limit = 1e-3f;
	}
};

extern t_struc struc;

#endif /* STRUC_H_ */
