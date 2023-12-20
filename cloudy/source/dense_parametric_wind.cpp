/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "dense.h"
#include "radius.h"
#include "physconst.h"

/*dense_parametric_wind called by dlaw command, returns density for any density law */
double dense_parametric_wind(double rad)
{
	/* >> refer wind profile	Springmann, U., 1994, A&A, 289, 505 */

	// SOLAR_MASS / 3e7 converts from M_solar/year to g/s
	double Mdot = dense.DensityLaw[0] * SOLAR_MASS / 3e7;
	double v_inf = dense.DensityLaw[1] * 1e5;
	double Beta2 = dense.DensityLaw[2];
	double Beta1 = dense.DensityLaw[3];
	double v_0 = dense.DensityLaw[4] * 1e5;
	double v_star = dense.DensityLaw[5] * 1e5;

	double r_star = radius.rinner;
	double x = MIN2( 0.01, 1. - r_star/rad );
	double v_r = v_star + (v_inf - v_0) * sqrt( Beta1 * x + (1.-Beta1) * pow(x, Beta2) );
	double mu = 1.;
	if( dense.wmole > 0. )
		mu = dense.wmole;
	double density = Mdot / ( PI4 * ATOMIC_MASS_UNIT * mu * pow2(rad) * v_r);
	return density;
}
