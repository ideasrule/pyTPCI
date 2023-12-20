/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_EINSTA_H_
#define HELIKE_EINSTA_H_

const int N_HE1_TRANS_PROB = 651;

void HelikeTransProbSetup();

/** compute energy diffference in wn and Aul for given line
 * return is 0 for success, 1 for failure 
\param nelem charge on the C scale, 1 is helium
\param Eff_nupper upper quantum numbers
\param lHi
\param sHi
\param jHi
\param Eff_nlower lower quantum numbers
\param lLo
\param sLo
\param jLo
*/
double he_1trans(  
		long nelem, 
		double Eff_nupper, long lHi, long sHi, long jHi,
		double Eff_nlower, long lLo, long sLo, long jLo);

#endif /* HELIKE_EINSTA_H_ */
