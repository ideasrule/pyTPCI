/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_RECOM_H_
#define HELIKE_RECOM_H_

double DielecRecombRate( long ipISO, long nelem, long ipLo );

/** Recomb_Seaton59 - computes total recombination into levels greater than nmax = n.
 \param nelem
 \param temp
 \param n
*/
double Recomb_Seaton59( long nelem, double temp, long n );

/**He_cross_section returns cross section (cm^-2), is a wrapper for cross_section 
 \param EgammaRyd, the photon energy in Ryd,
 \param EthRyd, the threshold energy in Ryd,
 \param ipLevel, the index of the level, 0 is ground, 3 within 2 3P,
 \param n
 \param l
 \param S
 \param ipZ is charge, equal to 1 for Helium,
*/
double He_cross_section( double EgammaRyd , double EthRyd, long n, long l, long S, long nelem );



#endif /* HELIKE_RECOM_H_ */
