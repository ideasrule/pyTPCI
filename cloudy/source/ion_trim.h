/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ION_TRIM_H_
#define ION_TRIM_H_

void ion_trim_init();

void ion_trim_untrim(long nelem);
void ion_trim_invalidate( long nelem );

void ion_trim_from_set();

void ion_trim_small (long nelem, double abund_total);

void ion_widen(long nelem);
void ion_trim_validate(long nelem, bool lgIonizTrimCalled);

void mole_ion_trim(void);

/**ion_trim raise or lower most extreme stages of ionization considered 
\param nelem element number on the C scale, 5 for C
*/
void ion_trim(
	long int nelem );

void ion_trim2(
	long int nelem,
	long int oldIonRange[2] );

#endif /* ION_TRIM_H_ */

