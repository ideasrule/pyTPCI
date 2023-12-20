/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "continuum.h"

t_continuum continuum;

void t_continuum::zero()
{
	lgCoStarInterpolationCaution = false;
	lgCon0 = false;

	/* upper limit to energies of inner shell opacities in Ryd
	 * this is 1 MeV by default */
	EnergyKshell = 7.35e4;

	lgPrtIsotropicCont = true;
}
