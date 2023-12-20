/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "gravity.h"
#include "dense.h"
#include "pressure.h"
#include "radius.h"
#include "colden.h"
#include "dark_matter.h"
#include "cosmology.h"
#include "physconst.h"

void GravitationalPressure( void )
{
	DEBUG_ENTRY( "GravitationalPressure()" );

	double M_self, M_external, g_dark, g_self, g_external;
	double R = radius.Radius - radius.dRadSign * radius.drad / 2;

	M_self = M_external = g_dark = g_self = g_external = 0.;

	if( dark.lgNFW_Set )
	{
		double rho_crit = 3.*POW2(cosmology.H_0*1e5/MEGAPARSEC)/(8.*PI*GRAV_CONST);
		double c_200 = dark.r_200/dark.r_s;
		ASSERT( c_200 > 0. );
		double delta_c = (200./3.) * POW3( c_200 ) / ( log(1.+c_200) - c_200/(1.+c_200) );
		double M_dark = 4*PI*rho_crit*delta_c*POW3(dark.r_s);
		M_dark *= 1./(1.+R/dark.r_s) + log(1.+R/dark.r_s) - 1.;

		g_dark = -GRAV_CONST * M_dark / (R*R);
	}

	/* (self-)Gravity forces: Yago Ascasibar (UAM, Spring 2009) */

	// add all external mass components
	for( unsigned int i=0; i < pressure.external_mass[0].size(); i++ )
	{
		double M_i = pressure.external_mass[0][i];
		if( R < pressure.external_mass[1][i] )
			M_i *= pow( R / pressure.external_mass[1][i], pressure.external_mass[2][i] );
		M_external += M_i;
	}

	// evaluate gravitational acceleration
	switch( pressure.gravity_symmetry )
	{
		case -1: // no self-gravity
			break;
		
		case 0: // spherical symmetry
			M_self = dense.xMassTotal - dense.xMassDensity * radius.dVeffVol;
			// dense.xMassTotal has already been updated in radius_increment.cpp
			M_self *= 4*PI* radius.rinner*radius.rinner;
			M_self *= pressure.self_mass_factor;

			g_self     = - GRAV_CONST * M_self / (R*R);
			g_external = - GRAV_CONST * M_external * SOLAR_MASS / (R*R);
			break;
		
		case 1: // mid-plane symmetry
			M_self = colden.TotMassColl + dense.xMassDensity * radius.drad_x_fillfac / 2.;
			// colden.TotMassColl will be updated later on in radius_increment.cpp
			M_self *= pressure.self_mass_factor;
			M_self *= 2.; // mid-plane symmetry

			// Gravitational acceleration due to an infinite slab is independent of distance from
			// the slab and is given by g = - 2 PI G * sigma.
			g_self     = - 2*PI* GRAV_CONST * M_self;
			g_external = - 2*PI* GRAV_CONST * M_external * SOLAR_MASS/PARSEC/PARSEC;

			if( dark.lgNFW_Set )
		    		fprintf( ioQQQ, " WARNING: Setting both mid-plane baryonic gravity symmetry and an NFW dark matter halo is almost certainly unphysical!\n"  );
			break;
	    	
		default:
			fprintf( ioQQQ, " Unknown gravitational symmetry = %d !!!\n", pressure.gravity_symmetry );
			TotalInsanity();
	};

	pressure.RhoGravity_dark = dense.xMassDensity * g_dark * radius.drad_x_fillfac;
	pressure.RhoGravity_self = dense.xMassDensity * g_self * radius.drad_x_fillfac;
	pressure.RhoGravity_external = dense.xMassDensity * g_external * radius.drad_x_fillfac;
	pressure.RhoGravity = pressure.RhoGravity_dark + pressure.RhoGravity_self + pressure.RhoGravity_external;

	return;
}
