/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PRESSURE_H_
#define PRESSURE_H_

#include "module.h"

class TransitionProxy;

/**PressureTotal determine the gas and line radiation pressures for current conditions,
 * this sets the variable pressure.PresTotlCurr */
void PresTotCurrent(void);

/**PressureRadiationLine calculate radiation pressure due to a particular line */
double PressureRadiationLine( const TransitionProxy &t, realnum DopplerWidth );

/** variables dealing with pressure across model */
struct t_pressure : public module {

	const char *chName() const
	{
		return "pressure";
	}

	void zero();
	void comment(t_warnings&) {}

	/** lowest and highest total current pressure to desired
	 * pressures in model, a measure of how well converged it is.  This
	 * includes total pressure, so in a constant gas pressure model will
	 * show a large excursion */
	realnum PresLow, 
	  PresHigh;

	realnum PresPowerlaw;

	/** the ram pressure */
	double PresRamCurr;

	/** current turbulent pressure */
	double PresTurbCurr;

	/** the current pressure, and the correct pressure, the ratio
	 * of the two is the error, PresTotlCurr is set in PressureTotal */
	double PresTotlCurr, 
	  PresTotlError,
	  /** PresGasCurr is gas pressure, nkT, set in PressureTotal */
	  PresGasCurr;

	/** PresTotlInit - total pressure at the illuminated face */
	double PresTotlInit;

	/** option to set pressure at illuminated face as number on 
	 * constant pressure command */
	bool lgPressureInitialSpecified;
	/** linear pressure in nkT units, zero if not set */
	double PressureInitialSpecified;

	/** pres_radiation_lines_curr is line radiation pressure for current zone */
	double pres_radiation_lines_curr;

	/** flag saying whether or not incident continuum should be included
	 * in total pressure, turned off with constant gas pressure */
	bool lgContRadPresOn;

	/** total pressure related variables
	 * PresInteg is integrated radiation pressure */
	realnum PresInteg, 
	  pinzon;

	/** integrated radiative acceleration, if no opacity present,
	 * only dilution due to 1/r^2
	 * this is the acceleration used in Eddington luminosity comparison */
	realnum PresIntegElecThin,
		pinzon_PresIntegElecThin;

	/** (Self-)gravity forces: Yago Ascasibar (UAM, Spring 2009) */
	double RhoGravity_dark;
	double RhoGravity_self;
	double RhoGravity_external;
	double RhoGravity;
	double IntegRhoGravity;
	int gravity_symmetry;
	double self_mass_factor;
	vector<double> external_mass[3];

	/** lgPres_radiation_ON says whether radiation pressure enabled, turned off with
	* constant density, constant gas pressure commands, on with constant pressure */
	bool lgPres_radiation_ON;
	bool lgPres_magnetic_ON;
	bool lgPres_ram_ON;

	realnum 
	  /** RadBetaMax is largest ratio of radiation to gas pressure */
	  RadBetaMax, 
	  /** pbeta is ratio of radition to gas pressure, evaluated in PressureTotal */
	  pbeta, 
	  /** PresMas is largest pressure that occurred in the calculation */
	  PresMax;

	/** pointer to line with greatest radiation pressure */
	long int ipPradMax_line;

	/** zone where greatest radiation pressure occurred */
	long int ipPradMax_nzone;

	/** string with label for line with greatest contributrion to pressure*/
	string chLineRadPres;

	/** lgPradCap true if radiation pressure capped on first iteration
	 * lgPradDen capped by thermalization length */
	bool lgPradCap, 
	  lgPradDen;

	/** flag true if radiation pressure is turned on, part of total pressure */
	bool lgLineRadPresOn;

	/** option to not abort on high radiation pressure - set true in initialization,
	 * set false with NO ABORT on constant radiation pressure command */
	bool lgRadPresAbortOK;

	/** option to abort when we reach the sonic point - set according to
	 * the dynamics pressure mode */ 
	bool lgSonicPointAbortOK; 

	/** we hit the sonic point */
	bool lgSonicPoint;

	/** True when we are in limbo while trying to find a strong-D
	 * solution. This is when there is no possible solution for the
	 * current target pressure. We just grit our teeth and plough
	 * onwards, hoping to come out the other side and to sort it all
	 * out on a later iteration of the dynamics */
	bool lgStrongDLimbo;

	};
extern t_pressure pressure;


#endif /* PRESSURE_H_ */
