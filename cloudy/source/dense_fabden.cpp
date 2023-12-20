/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*dense_fabden called by dlaw command, returns density for any density law */
#include "cddefines.h"
#include "rfield.h"
#include "dense.h"

/*dense_fabden implements the dlaw command, returns density using
 * current position and up to ten parameters on dlaw command line */
double dense_fabden(double radius, 
  double depth)
{
	double fabden_v = exp10(dense.DensityLaw[0]);
	if( rfield.lgUSphON )
		fabden_v *= pow(radius/rfield.rstrom,dense.DensityLaw[1]);
	else
		fabden_v *= pow(depth/rfield.rstrom,dense.DensityLaw[1]);

	return fabden_v; 
}
