/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H_cross_section returns cross section (cm^-2), 
 * given EgammaRyd, the photon energy in Ryd,
 * ipLevel, the index of the level, 0 is ground,
 * nelem is charge, equal to 0 for Hydrogen */
#include "cddefines.h" 
#include "hydrogenic.h" 
#include "hydro_bauman.h"

/*H_cross_section returns cross section (cm^-2), 
 * given EgammaRyd, the photon energy in Ryd,
 * EthRyd, the threshold energy in Ryd,
 * quantum numbers n and l
 * nelem is charge, equal to 0 for Hydrogen */
double H_cross_section( double EgammaRyd , double EthRyd, long n, long l, long nelem )
{
	double cs;
	double rel_photon_energy;

	/* >>chng 02 apr 24, more protection against calling with too small an energy */
	/* evaluating H-like photo cs at He energies, may be below threshold -
	 * prevent this from happening */
	rel_photon_energy = EgammaRyd / EthRyd;
	rel_photon_energy = MAX2( rel_photon_energy , 1. + FLT_EPSILON*2. );

	cs = H_photo_cs(rel_photon_energy , n, l, nelem + 1 );

	ASSERT( cs > 0. && cs < 1.E-8 );

	return cs;
}
