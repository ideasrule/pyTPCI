/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCosmology parse cosmological parameters and options */
#include "cddefines.h"
#include "physconst.h"
#include "cosmology.h"

t_cosmology cosmology;

realnum GetHubbleFactor(realnum z)
{
	realnum H_z, H_z_squared;

	DEBUG_ENTRY( "GetHubbleFactor()" );

	/* NB - the factor of (1e5/MEGAPARSEC) converts from km/s/Mpc to km/s/km */
	H_z_squared = POW2(cosmology.H_0 * (realnum)(1e5/MEGAPARSEC)) * ( 
		cosmology.omega_lambda + 
		cosmology.omega_matter * POW3( 1.f + z ) +
		cosmology.omega_rad * POW4( 1.f + z ) + 
		cosmology.omega_k * POW2( 1.f + z ) );

	H_z = sqrt( H_z_squared );

	return H_z;
}

realnum GetDensity(realnum z)
{
	realnum density;

	DEBUG_ENTRY( "GetDensity()" );

	fixit("Cosmological He fraction should be abund.aprim[1] by default, but controlled by command line option");
	cosmology.f_He = 0.079f;

	/* from Switzer & Hirata 2007, equation 2 */
	density = 1.123e-5f * (cosmology.omega_baryon*cosmology.h*cosmology.h) / (1.f + 3.9715f*cosmology.f_He) * 
		pow3( 1.f + z );

	return density;
}
