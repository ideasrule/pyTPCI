/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "doppvel.h"
t_DoppVel DoppVel;

void t_DoppVel::zero()
{
	DEBUG_ENTRY( "t_DoppVel::zero()" );

	/* velocity field information */
	/* the turbulent velocity at illuminated face, internally in cm/s,
	 * but entered with turbulence command in km/s */
	TurbVel = 0.;
	/* is a turbulent gradient imposed?  Default is no. */
	lgTurbLawOn = false;
	/* the log of the turbulence gradient power law. Default is zero. */
	TurbVelLaw = 0.;
	/* the parameter F in eq 34 of
	 *>>refer	pressure	turb	Heiles, C. & Troland, T.H. 2005, 624, 773 */
	Heiles_Troland_F = 0.;
	/* is TurbVel included in pressure? - can be done two ways, with the velocity
	 * being set of with equipartition - true when TurbVel set if not equipartition,
	 * false with NO PRESSURE option on turbulence command */
	lgTurb_pressure = true;
	/* The scale in cm over which the turbulence is dissipated.  Normally 0,
	 * only set if dissipate keyword appears on turbulence command */
	DispScale = 0.;
	/* equipartition option on turbulence command, to set turbulence from B */
	lgTurbEquiMag = false;
	
}
