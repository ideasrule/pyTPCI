/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "hydrogenic.h"
t_hydro hydro;

void t_hydro::zero()
{
	DEBUG_ENTRY( "t_hydro::zero()" );
	cintot = 0.;

	/* option to print emissivity instead of intensity/luminosity */
	lgHiPop2 = false;
	pop2mx = 0.;

	/* flag for Lya masing */
	HCollIonMax = 0.;

	/* type of hydrogen atom top off, options are " add" and "scal" 
	 * in versions 90 this was " add", but was "scal" in 91
	 * >>chng 99 jan 16, changed back to " add"*/
	/*strcpy( chHTopType, "scal" );*/
	strcpy( chHTopType, " add" );

	/* Lya excitation temperature, counter for hotter than gas */
	TexcLya = 0.;
	TLyaMax = 0.;
	nLyaHot = 0;

	/* option to kill damping wings of Lya */
	DampOnFac = 1.;

	/* is FeII pumping by H Lyman lines included? */
	lgLyaFeIIPumpOn = true;

	/* is continuum pumping of H Lyman lines included?  yes, but turned off
	 * with atom h-like Lyman pumping off command */
	lgLymanPumping = true;

	/* multiplicative factor for all continuum pumping of H I Lyman lines,
	 * account for possible emission in the line */
	xLymanPumpingScaleFactor = 1.f;
}
