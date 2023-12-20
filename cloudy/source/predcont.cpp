/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "predcont.h"
#include "rfield.h"

t_PredCont::t_PredCont()
{
	/* energies where diffuse continuum is to be entered into line array 
	 * NB - if these numbers change, the wavelength in the printout will change 
	 * too, since the wavelength is derived form this */
	/* >>>chng 99 mar 23, adjusted energies so that wavelength line list is
	 * the same as it was in C90 - small changes were caused by going over
	 * to proper Rydberg constant */

	p_val.reserve(100);

	add(1275.,"MHz"); /* VLA */
	add(1465.,"MHz"); /* VLA */
	add(4535.,"MHz"); /* VLA */
	add(4885.,"MHz"); /* VLA */
	add(8435.,"MHz"); /* VLA */
	add(8735.,"MHz"); /* VLA */
	add(3.4,"cm");
	add(14965.,"MHz"); /* VLA */
	add(22460.,"MHz"); /* VLA */
	add(30.,"GHz"); /* OCRA */
	add(43340.,"MHz"); /* VLA */
	add(7.445e-04,"Ryd"); 
	add(1.498e-03,"Ryd"); 
	add(2.211e-03,"Ryd"); 
	add(2.952e-03,"Ryd"); 
	add(3.677e-03,"Ryd"); 
	add(3.7501e-03,"Ryd"); /* Ney-Allen */
	add(3.9915e-03,"Ryd"); /* Ney-Allen */
	add(4.2543e-03,"Ryd"); /* Ney-Allen */
	add(4.314e-03,"Ryd"); 
	add(4.6446e-03,"Ryd"); /* Ney-Allen */
	add(5.162e-03,"Ryd"); 
	add(5.2462e-03,"Ryd"); /* Ney-Allen */
	add(5.8079e-03,"Ryd"); /* Ney-Allen */
	add(6.240e-03,"Ryd"); 
	add(7.3312e-03,"Ryd"); /* Ney-Allen */
	add(7.9936e-03,"Ryd"); /* Ney-Allen */
	add(8.7119e-03,"Ryd"); /* Ney-Allen */
	add(9.6125e-03,"Ryd"); /* Ney-Allen */
	add(9.77243e-03,"Ryd");
	add(1.1099e-02,"Ryd"); /* Ney-Allen */
	add(1.2022e-02,"Ryd"); /* Ney-Allen */
	add(1.29253e-02,"Ryd"); 
	add(2.2152e-02,"Ryd"); 
	add(3.92044e-02,"Ryd"); 
	add(5.54593e-02,"Ryd"); 
	/* next two either side of n=4 edge of hydrogen, set to 1.5% off either direction*/
	/* >>chng 00 sep 18, had been too close in energy */
	add(6.1563e-02,"Ryd"); 
	add(6.3437e-02,"Ryd"); 
	add(8.1460e-02,"Ryd"); 
	/* >>chng 00 sep 14, changed energies of paschen jump to be farther away as
	 * per note on BJ */
	add(0.1094,"Ryd"); 
	add(0.1128,"Ryd"); 
	add(0.14675,"Ryd"); 
	add(0.18653,"Ryd"); 
	/* >>chng 00 sep 14, next two energies changed since they were too close to BJ
	 * and so both ended up shortward of limit*/
	/* these two are the Balmer jump, below and above. */
	/* continuum binning not much better than 1% so offset energies by more */
	add(0.246,"Ryd");
	add(0.254,"Ryd");
	add(0.375,"Ryd");  /* peak on two photon continuum */
	add(0.38096,"Ryd"); 
	add(0.43994,"Ryd"); 
	add(0.44394,"Ryd"); 
	add(0.50811,"Ryd"); 
	add(0.57489,"Ryd");
	add(0.62487,"Ryd"); 
	add(0.67155,"Ryd"); 
	add(0.70244,"Ryd"); 
	add(0.72163,"Ryd"); 
	add(0.74812,"Ryd"); 
	add(0.76172,"Ryd"); 
	add(0.77551,"Ryd"); 
	add(0.79681,"Ryd"); 
	add(0.81859,"Ryd"); 
	add(0.8260,"Ryd"); 
	add(0.84859,"Ryd"); 
	add(0.85618,"Ryd"); 
	add(0.87967,"Ryd"); 
	add(1000.,"A");
	/* points on either side of Lyman jump,
	 * energies changed to be robust when energy grid changes,
	 * grid resolution is about 1%, so change from 0.99783 and 1.000
	 * to 1 +/- 1.5%
	 * >>chng 00 sep 23 change wavelength points for next two */
	add(0.985,"Ryd"); 
	add(1.015,"Ryd"); 
	add(1.199,"Ryd"); 
	add(1.299,"Ryd"); 
	add(1.4984,"Ryd"); 
	add(1.58441,"Ryd"); 
	/* points on either side of Lyman jump,
	 * energies changed to be robust when energy grid changes,
	 * grid resolution is about 1%, so change from 1.80433 and 1.809
	 * to 1.807 +/- 1.5%
	 * >>chng 00 sep 23 change wavelength points for next two */
	add(1.780,"Ryd"); 
	add(1.834,"Ryd"); 
	add(2.283,"Ryd");
}

long t_PredCont::find(double energy, const char* unit) const
{
	DEBUG_ENTRY( "t_PredCont::find()" );

	for( size_t i=0; i < p_val.size(); ++i )
		if( fp_equal( p_val[i].get(unit), energy ) )
			return i;
	return -1;
}

long t_PredCont::add(double energy, const char* unit)
{
	DEBUG_ENTRY( "t_PredCont::add()" );

	long ind = find(energy, unit);
	if( ind < 0 )
	{
		p_val.push_back( EnergyEntry(energy, unit) );
		ind = p_val.size()-1;
	}
	double eRyd = p_val[ind].Ryd();
	if( eRyd < rfield.emm() || eRyd > rfield.egamry() )
	{
		fprintf( ioQQQ, " The energy %g Ryd (%g %s) is not within the default Cloudy range\n",
			 eRyd, energy, unit );
		fprintf( ioQQQ, " The energy must be between %g and %g Ryd\n",
			 rfield.emm(), rfield.egamry() );
		cdEXIT(EXIT_FAILURE);
	}
	return ind;
}
