/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_H_
#define HELIKE_H_

/** \verbatim
  Set this flag to one of the following values
 	0	don't print
 	1	print As
 	2	print only forbidden As
 	3	print Es
 	4	print threshold photoionization cross-sections
 	5	print radiative recombination coefficients.
        6	print photoionization cross-section grids. \endverbatim <BR>
  Out files are "Helike___.txt" with the blank being, respectively,
 *	"As","fAs","Es","ThPCS","RR", and "PCSgrid." */
const int printflag = 0;

/** the magic number for the table of He-like transition probabilities, YYMMDD */
const int TRANSPROBMAGIC = 190102;
/** the magic number for the table of He1 effective collision strengths, YYMMDD */
const int COLLISMAGIC = 190101;
/** the magic number for the table of He1 case B emissivities, YYMMDD */
const int CASEBMAGIC = 190107;

/**HeCollidSetup read in helium collision data files */
void HeCollidSetup( void );

/**helike_energy get helike level energy in cm-1 
\param nelem element number
\param n quantum number n
\param l quantum number l
\param s quantum number s
\param j quantum number j
*/
double helike_energy(long nelem, long n, long l, long s, long j);

/**helike_quantum_defect get quantum defect for helium-like level
\param nelem
\param n
\param lqn
\param S
\param j
*/
double helike_quantum_defect( long nelem, long n, long lqn, long S, long j );

/**helike_transprob get transition probability for helium-like transition [s-1]
\param nelem
\param ipHi
\param ipLo
*/
realnum helike_transprob( long nelem, long ipHi, long ipLo );

/**AGN_He1_CS routine to save table needed for AGN3 - collision strengths of HeI 
\param *ioPun
*/
void AGN_He1_CS( FILE *ioPun );

#endif /* HELIKE_H_ */
