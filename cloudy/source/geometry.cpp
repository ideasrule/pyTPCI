/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "geometry.h"
t_geometry geometry;

void t_geometry::zero()
{
	/* variable to do with geometry */
	nprint = 1000;
	lgZoneSet = false;
	lgZoneTrp = false;

	fiscal = 1.;
	FillFac = 1.;
	filpow = 0.;

	/* default is open geometry, not sphere */
	lgSphere = false;
	/* the radiative transport covering factor */
	covrt = 0.;
	/* the geometric covering factor */
	covgeo = 1.;
	/* default is expanding when geometry set */
	lgStatic = false;
	/* option to tell code not to complain when geometry static done without iterating,
	 * set with (OK) option on geometry command */
	lgStaticNoIt = false;
	/* this is exponent for emissivity contributing to observed luminosity, r^2.
	 * set to 1 with aperture slit, to 0 with aperture beam command */
	iEmissPower = 2;
}
