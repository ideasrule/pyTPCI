/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef SECONDARIES_H_
#define SECONDARIES_H_

#include "module.h"
#include "container_classes.h"

/** these are the global variables dealing with the effects of hot secondary electrons */
struct t_secondaries : public module {
	
	const char *chName() const
	{
		return "secondaries";
	}

	void zero();
	void comment(t_warnings&) {}

	void alloc();

	/** heating efficiency, unity for ionized gas, 0 < HeatEfficPrimary < 1 for
	 * partially neutral medium */
	realnum HeatEfficPrimary;
	/** number of secondary ionizations per primary erg of 
	 * photo-electron energy - It multiplies the high-energy heating rate.  
	 * zero for a fully ionized gas and -> 0.3908 / erg-Ryd for very neutral gas */
	realnum SecIon2PrimaryErg;
	/** secondary excitations of Lya per primary erg of energy */
	realnum SecExcitLya2PrimaryErg; 

	/** secondary ionization [nelem][ion] rate [s-1]*/
	multi_arr<realnum,2> csupra;
	multi_arr<realnum,2> csupra_effic;

	/** sec2total is ratio of secondary to total ionizations of ground state H,
	 * evaluated in hydrogen, used in ionization to determine how fine a 
	 * convergence on the secondary rates is needed */
	realnum sec2total;

	/** max ratio of sec2total, recall to identify secondary ionization models */
	realnum SecHIonMax;

	/**SetCsupra is H^0 secondary ionization rate set with set csupra command */
	realnum SetCsupra;

	/**lgCSetOn set with set csupra command */
	bool lgCSetOn;

	/**flag saying that secondary electron have been turned off, 
	 * set true with no secondaries command */
	bool lgSecOFF;

	/** pointer to continuum energies where secondary ionization can occur, 100eV */
	long int ipSecIon;

	/** this block of variables save the state of the code in startr,
	 * in reset, resets this state for the start of the next zone */
	realnum /*supsav, */
	  hetsav, 
	  savefi, 
	  x12sav;

	/** x12tot rate for excitation of Lya, units excitations s-1 per atom */
	realnum x12tot; 

	};

extern t_secondaries secondaries;


#endif /* SECONDARIES_H_ */
