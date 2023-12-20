/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HEAVY_H_
#define HEAVY_H_

/* heavy.h */
struct t_Heavy {

	/** array index for valence shell photoionization */
	long int ipHeavy[LIMELM][LIMELM], 
		/** array index for average lya of this ion */
	  ipLyHeavy[LIMELM][LIMELM-1],
	  /** array index for "Balmer" recom of this ion */
	  ipBalHeavy[LIMELM][LIMELM-1];

	/** radiative recombination continuum to ground state, erg cm-3 s-1 */
	double RadRecCon[LIMELM][LIMELM];
	
	/**xLyaHeavy is Lya for nelem,ion */
	realnum xLyaHeavy[LIMELM][LIMELM];

	/** ionization pot of valence shell, rydbergs, [nelem][ion] */
	double Valence_IP_Ryd[LIMELM][LIMELM];

	/** save number of shells for each element and ionization stage
	 * first dimension is atomic weight, second is ionization stage */
	long int nsShells[LIMELM][LIMELM];

	/** character representation of shells 1 to 7 of atoms and ions, '1s', '2s', etc */
	char chShell[7][3];

	t_Heavy()
	{
		/* list of shells, 1 to 7 */
		strcpy( chShell[0], "1s" );
		strcpy( chShell[1], "2s" );
		strcpy( chShell[2], "2p" );
		strcpy( chShell[3], "3s" );
		strcpy( chShell[4], "3p" );
		strcpy( chShell[5], "3d" );
		strcpy( chShell[6], "4s" );
	}
};
extern t_Heavy	Heavy;


#endif /* HEAVY_H_ */
