/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DOPPVEL_H_
#define DOPPVEL_H_

#include "module.h"

/** doppvel.h */
struct t_DoppVel : public module {
	const char *chName() const
	{
		return "DoppVel";
	}

	void zero();
	void comment(t_warnings&) {}

	/** this is the turbulent velocity in cm/s, normally zero, set with
	 * turbulence command */
	realnum TurbVel;

	/** this is the turbulent velocity in cm/s at the beginning of the calculation,
	 * this is only important if turbulence law is set with vlaw command */
	realnum TurbVelZero;

	/** this is the turbulent velocity in cm/s, normally zero, set with
	 * turbulence command */
	realnum TurbVelLaw;

	/** this flag says that a turbulent velocity law (as a function of depth) 
	 * has been set with the vlaw command */
	bool lgTurbLawOn;

	/** this is the factor F in equation 34 of 
	 *>>refer	pressure	turb	Heiles, C. & Troland, T.H. 2005, 624, 773 */
	realnum Heiles_Troland_F;

	/** is TurbVel included in pressure? - can be done two ways, with the velocity
	 * being set of with equipartition - true when TurbVel set if not equipartition,
	 * false with NO PRESSURE option on turbulence command */
	bool lgTurb_pressure;

	/** equipartition option on turbulence command, to set turbulence from B 
	 * if true then TurbVel set from B */
	bool lgTurbEquiMag;

	/** The scale in cm over which the turbulence is dissipated.  Normally 0,
	 * only set if dissipate keyword appears on turbulence command, 
	 * entered as log of scale */
	realnum DispScale;

	};
extern t_DoppVel DoppVel;

/** GetDopplerWidth get the doppler width at current conditions for a given mass
 */
realnum GetDopplerWidth( realnum massAMU );
/** GetAveVelocity get the average particle velocity at current conditions for a given mass
 */
realnum GetAveVelocity( realnum massAMU );


#endif /* DOPPVEL_H_ */
