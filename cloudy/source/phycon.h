/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PHYCON_H_
#define PHYCON_H_

/** phycon.h */

#include "module.h"

struct t_phycon : public module {
	const char *chName() const
	{
		return "phycon";
	}

	void zero();
	void comment(t_warnings&) {}

	/**te, electron temperature K */
	double te;

	/** the current electron temperature in eV */
	double te_eV;

	/** the current electron temperature in Ryd */
	double te_ryd;

	/** the current electron temperature in wavenumbers */
	double te_wn;

	/** 1/TE */
	double teinv;

	/** T^2 */
	double tesqrd;

	/** total ionization energy of gas, erg/cm^3, evaluated in PressureTotal,
	 * this is the amount of energy needed to go from pure atoms to the current
	 * ionization of the gas. */
	double EnergyIonization;

	/** the total internal energy of atoms and molecules within the gas, erg/cm^3.  This is
	 * the amount of energy to excite the current level populations from the ground.
	 * this includes electronic excitations of atoms and rotation excitations of molecules 
	 */
	double EnergyExcitation;

	/** the enthalpy per unit vol, updated with pressure updated */
	double EnthalpyDensity;

	/** this is the total binding energy of the molecules, and is negative, the energy
	 * need to get back to free atoms */
	double EnergyBinding;

	/** these are simple powers of the electron temperature, are
	 * evaluated in tfidle, and can be used to avoid exponentials*/
	double sqrte, 
	  te32, 
	  te90, 
	  te70, 
	  te40, 
	  te30, 
	  te20, 
	  te10, 
	  te07,
	  te05, 
	  te04, 
	  te03, 
	  te02, 
	  te01, 
	  te007, 
	  te005, 
	  te004, 
	  te003, 
	  te002, 
	  te001,
	/*>>chng 06 June 30 :Added te0001,0002,0003,0004,0005,0007- Humeshkar Nemala*/
		te0001,
		te0002,
		te0003,
		te0004,
		te0005,
		te0007;
	/**1 is log Te, 2 is (log T)^2, etc*/
	double telogn[7] ,
		/** this is used to rescale telogn array for other charge,
		 * ith element is log(i+1)^(i+1)*/
		sqlogz[LIMELM];

	/** alogte is base 10 log of temperature */
	double alogte;

	/** alnte is natural log of temperature */
	double alnte;

	double 
	  /** initial temperature */
	  TeInit, 
	  /** proposed temperature */
	  TeProp,
	  /** initial electron density */
	  EdenInit ,
	  /** proposed electron density */
	  EdenProp;

	/** energy density temperature */
	double TEnerDen;

	/** lag set if a physical condition has been disabled */
	bool lgPhysOK;

	/** largest relative changes in Te, ne, H+, H2, and CO in structure
	 * this is computed as part of prtcomment so does not exist when code not talking,
	 * set to zero in zero and still zero if prtcomment not called */
	realnum BigJumpTe , BigJumpne , BigJumpH2 , BigJumpCO;

	/** The default value of the stopping temperature */
	const double TEMP_STOP_DEFAULT;
	/** lowest temperature to ever allow */
	const double TEMP_LIMIT_LOW;
	/** highest temperature to ever allow */
	const double TEMP_LIMIT_HIGH, TEMP_LIMIT_HIGH_LOG;

	/** define constructor to set initial values for these constant variables 
	 * within class.  Order matters here - from PvH:
	 * This is what Stroustrup says:
	 * "The constructors [i.e. the initializers for TEMP_LIMIT_LOW, etc.] are 
	 * called in the order in which the members are declared in the class 
	 * rather than the order in which the members appear in the initializer 
	 * list. To avoid confusion, it is best to specify the initializers in 
	 * the member declaration order." */
	t_phycon() : TEMP_STOP_DEFAULT(4000.) , TEMP_LIMIT_LOW(2.8), 
		TEMP_LIMIT_HIGH(1.001e10) , TEMP_LIMIT_HIGH_LOG(10.0004340775) {}

	};

extern t_phycon phycon;


#endif /* PHYCON_H_ */
