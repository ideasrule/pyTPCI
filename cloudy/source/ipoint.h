/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef IPOINT_H_
#define IPOINT_H_

/** ipoint - the basic routine that generates an index for the continuum array 
\return array index on the Fortran or Physics scale, [1] would be first one 
\param energy - photon energy in Ryd
*/
long ipoint(double energy);

/**ipLineEnergy generates a safe pointer to energy in line array with energy
 * given by the first argument (Rydbergs).  The second is a label that
 * will be placed at that energy if none have yet been set.  The last
 * is the energy of the next higher continuum.  the code will make sure
 * that the pointer is not within this continuum - that would create
 * serious energy problems if, for instance, a Lyman line could ionize H
 \return array index on the Fortran or Physics scale, [1] would be first one 
 \param energy - photon energy in Ryd
 \param chString - 4 char + null string giving label for the line
 \param ipIonEnergy - ipIonEnergy - if <0 ignored, if >0 will make sure that
return array index is less than this - used to make sure we don't overflow into
a higher continuum */
long ipLineEnergy(double energy , const char* chString , long ipIonEnergy );

/**ipContEnergy generates a safe pointer to energy in continuum array with energy
 * given by the first argument (Rydbergs).  The second is a label that
 * will be placed at that energy if none have yet been set.  The last
 * is the energy of the next higher continuum.  the code will make sure
 * that the pointer is not within this continuum - that would create
 * serious energy problems 
 \return array index on the Fortran or Physics scale, [1] would be first one 
 \param energy - photon energy in Ryd
 \param chString - 4 char + null string giving label for the continuum edge */
long ipContEnergy(double energy , const char* chString );

/**ipFineCont returns array index within fine energy mesh 
\return array index on the Fortran or Physics scale, [1] would be first one 
\param energy in Ryd */
long ipFineCont(
	double energy );

#endif /* IPOINT_H_ */
