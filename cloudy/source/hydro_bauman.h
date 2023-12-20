/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYDRO_BAUMAN_H_
#define HYDRO_BAUMAN_H_

/**\file hydro_bauman.h these are the hydrogenic routines written by Robert Bauman 
     For references, see h_bauman.c                                  
     \author Robert Bauman
     */

/************************************************************************/
/*  IN THE FOLLOWING WE HAVE  n > n'                                    */
/************************************************************************/

/** returns hydrogenic photoionization cross section in cm^2             
\param photon_energy photon energy relative to threshold
\param n principal quantum number, 1 for ground
\param l angular momentum, 0 for s
\param iz charge, 1 for H+, 2 for He++, etc
*/
double H_photo_cs( 
	double photon_energy, 
	long int n, 
	long int l, 
	long int iz
);

/**\verbatim Calculates the Einstein A's for hydrogen                           
   for the transition n,l --> n',l'                                   
   units of sec^(-1)
                                                                      
   In the following, we have n > n' \endverbatim                                   
   \param n principal quantum number, 1 for ground, upper level
   \param l angular momentum, 0 for s
   \param np principal quantum number, 1 for ground, lower level
   \param lp angular momentum, 0 for s
   \param iz Nuclear charge, 1 for H+, 2 for He++, etc
*/
realnum H_Einstein_A(
	long int n,
	long int l,
	long int np,                
	long int lp,               
	long int iz
);

#endif /* HYDRO_BAUMAN_H_ */
