/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COSMOLOGY_H_
#define COSMOLOGY_H_

/** this is the current temperature of the cosmic background radiation<BR>
 >>refer	CMB	temp	Mather, J.C., Fixsen, D.J., Shafer, R.A., Mosier, C., & <BR>
 * >>refercon	Wilkinson, D.T. 1999, ApJ, 512, 511 */
#define	CMB_TEMP	2.725

/**GetDensity find the baryonic density at the given redshift 
\param z
*/
realnum GetDensity(realnum z);

/**GetHubbleFactor find the Hubble factor at the given redshift 
\param z
*/
realnum GetHubbleFactor(realnum z);

/** cosmology.h saves options and parameters relating to cosmology */
struct t_cosmology {

	realnum redshift_current;
	realnum redshift_start;
	realnum redshift_step;

	realnum omega_baryon;
	realnum omega_rad;
	realnum omega_lambda;
	realnum omega_matter;
	realnum omega_k;

	realnum h;
	realnum H_0;

	realnum f_He;

	bool lgDo;

	t_cosmology()
	{
		redshift_current = 0.f;
		redshift_start = 0.f;
		redshift_step = 0.f;
		omega_baryon = 0.04592f;
		omega_rad = 8.23e-5f;
		omega_lambda = 0.7299177f;
		omega_matter = 0.27f;
		omega_k = 0.f;
		/* the Hubble parameter in 100 km/s/Mpc */
		h = 0.71f;
		/* the Hubble parameter in km/s/Mpc */
		H_0 = 100.f*h;
		lgDo = false;
	}
};
extern t_cosmology cosmology;


#endif /* COSMOLOGY_H_ */
