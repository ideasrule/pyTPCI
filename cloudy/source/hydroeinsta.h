/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYDROEINSTA_H_
#define HYDROEINSTA_H_

/**HydroEinstA calculates Einstein A's from  osillator strengths
\param lower
\param iupper
*/
double HydroEinstA(long int lower, 
	  long int iupper);

realnum hydro_transprob( long nelem, long ipHi, long ipLo );

#endif /* HYDROEINSTA_H_ */
