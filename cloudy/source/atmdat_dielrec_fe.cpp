/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atmdat_dielrec_fe Dielectronic recombination rates for Fe from Arnaud & Raymond 1992 */
#include "cddefines.h"
#include "atmdat.h"

/*atmdat_dielrec_fe Dielectronic recombination rates for Fe from Arnaud & Raymond 1992 */
double atmdat_dielrec_fe(
		// ion on physics scale
		long int ion,
		// electron temperture
		double t)
{
	static const double dfe[25][8] = {
		{5.120e00,1.29e01,0.00e00,0.00e00,2.20e-04,1.00e-04,0.00e00,0.00e00},
		{1.670e01,3.14e01,0.00e00,0.00e00,2.30e-03,2.70e-03,0.00e00,0.00e00},
		{2.860e01,5.21e01,0.00e00,0.00e00,1.50e-02,4.70e-03,0.00e00,0.00e00},
		{3.730e01,6.74e01,0.00e00,0.00e00,3.80e-02,1.60e-02,0.00e00,0.00e00},
		{5.420e01,1.00e02,0.00e00,0.00e00,8.00e-02,2.40e-02,0.00e00,0.00e00},
		{4.550e01,3.60e02,0.00e00,0.00e00,9.20e-02,4.10e-02,0.00e00,0.00e00},
		{6.670e01,1.23e02,0.00e00,0.00e00,1.60e-01,3.60e-02,0.00e00,0.00e00},
		{6.610e01,1.29e02,0.00e00,0.00e00,1.80e-01,7.00e-02,0.00e00,0.00e00},
		{2.160e01,1.36e02,0.00e00,0.00e00,1.40e-01,2.60e-01,0.00e00,0.00e00},
		{2.220e01,1.44e02,0.00e00,0.00e00,1.00e-01,2.80e-01,0.00e00,0.00e00},
		{5.960e01,3.62e02,0.00e00,0.00e00,2.25e-01,2.31e-01,0.00e00,0.00e00},
		{7.500e01,2.05e02,0.00e00,0.00e00,2.40e-01,1.70e-01,0.00e00,0.00e00},
		{3.630e01,1.93e02,0.00e00,0.00e00,2.60e-01,1.60e-01,0.00e00,0.00e00},
		{3.940e01,1.98e02,0.00e00,0.00e00,1.90e-01,9.00e-02,0.00e00,0.00e00},
		{2.460e01,2.48e02,5.60e02,0.00e00,1.20e-01,1.20e-01,6.00e-01,0.00e00},
		{5.600e02,0.00e00,0.00e00,0.00e00,1.23e00,0.00e00,0.00e00,0.00e00},
		{2.250e01,1.17e02,3.41e02,6.83e02,2.53e-03,3.36e-02,1.81e-01,1.92e00},
		{1.620e01,9.60e01,3.30e02,7.29e02,5.67e-03,7.82e-02,3.18e-02,1.26e00},
		{2.370e01,8.51e01,3.29e02,7.87e02,1.60e-02,7.17e-02,9.06e-02,7.39e-01},
		{1.320e01,6.66e01,2.97e02,7.14e02,1.85e-02,9.53e-02,7.90e-02,1.23e00},
		{3.910e01,8.03e01,3.92e02,9.19e02,9.20e-04,1.29e-01,1.92e-01,9.12e-01},
		{7.320e01,3.16e02,8.77e02,0.00e00,1.31e-01,8.49e-02,6.13e-01,0.00e00},
		{1.000e-01,3.62e01,3.06e02,9.28e02,1.10e-02,4.88e-02,8.01e-02,5.29e-01},
		{4.625e03,6.00e03,0.00e00,0.00e00,2.56e-01,4.52e-01,0.00e00,0.00e00},
		{5.300e03,0.00e00,0.00e00,0.00e00,4.30e-01,0.00e00,0.00e00,0.00e00}
	};

	double rate, te;

	DEBUG_ENTRY( "atmdat_dielrec_fe()" );
	/*Dielectronic recombination rates for Fe from 
	 * >>refer	Fe	rec	Arnaud, M. & Raymond, J. 1992, ApJ, 398, 394 */

	/* ion - spectroscopic symbol of final ion
	 * t - temperature, K
	 * d - rate coefficient, cm^3 s^-1 */

	if( ion > 26 )
	{
		fprintf( ioQQQ, " atmdat_dielrec_fe invalid ion%10ld\n", ion );
		cdEXIT(EXIT_FAILURE);
	}

	else if( ion == 26 )
	{
		/* d is the rate */
		rate = 0.0;
	}

	else
	{
		te = t*EVRYD/TE1RYD;
		rate = 0.0;
		for( int j=0; j < 4; j++ )
		{
			int k = j + 4;
			rate += (dfe[ion-1][k]*sexp(dfe[ion-1][j]/te));
		}
		rate = rate/powpq(t,3,2);
	}
	return rate;
}
