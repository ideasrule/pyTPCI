/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HEXTRA_H_
#define HEXTRA_H_

/* hextra.h */

#include "module.h"

/** galactic background energy density of 1.8 eV cm-3 is from
 *>>refer	cr	background	Webber, W.R. 1998, ApJ, 506, 329 */ 
#define	CR_EDEN_GAL_BACK_EV_CMM3	1.8

struct t_hextra : public module {
	const char *chName() const
	{
		return "hextra";
	}
	void zero();
	void comment(t_warnings&) {}

	/** heat due to cosmic rays*/
	realnum cryden, 
	  crpowr, 
	  crtemp;

	/** true for cosmic rays in equipartition with magnetic field */
	bool lg_CR_B_equipartition;

	/** cosmic ray energy density erg cm-3 */
	double cr_energydensity;

	/** current cosmic ray density divided by default galactic background */
	realnum cryden_ov_background;

	/** default cosmic ray background density and rate */
	realnum background_density;
	realnum background_rate;

	/** extra heating set with hextra command, first the heating rate */
	realnum TurbHeat, 

	  /** save the initial value in case TurbHeat varies with time */
	  TurbHeatSave;

	/** options for heating to depth on depth, true is depth occurs on
	 * hextra command */
	bool lgHextraDepth;
	/** the scale radius for the heating */
	realnum turrad,
	  /** the scale radius from the back of the cloud */
	  turback;

	/** options for extra heating the depends on density, as set with 
	 * hextra command */
	bool lgHextraDensity;

	/** the scale density */
	realnum HextraScaleDensity;

  	/** options for extra heating is from SS model, as set with 
	 * hextra command */
	bool lgHextraSS;

	/** the parameter alpha of alpha model, dimensionless */
	realnum HextraSSalpha;

	/** mass of the black hole in grams */
	double HextraSS_M;

	/** radius from center in cm */
	realnum HextraSSradius;
  
	/** set true if extra heat varies with time in time dependent sims */
	bool lgTurbHeatVaryTime;

	/** totneu is neutron energy flux, erg cm-2 s-1	*/
	realnum totneu;
	/** flag lgNeutrnHeatOn says heating due to neutrons is enabled */
	bool lgNeutrnHeatOn;
	/** frcneu is fraction of total luminosity in neutrons, dimensionless */
	realnum frcneu;
	/** effneu is efficiency */
	realnum effneu;
	/** cross section for stopping relativistic neutrons */
	double CrsSecNeutron;

	t_hextra()
	{
		/* default cosmic ray background */
		/* >>chng 99 jun 24, slight change to value
		 * quoted by 
		 * >>refer	cosmic ray	ionization rate	McKee, C.M., 1999, astro-ph 9901370
		 * this will produce a total
		 * secondary ionization rate of 2.5e-17 s^-1, as tested in 
		 * test suite cosmicray.in.  If each ionization produces 2.4 eV of heat,
		 * the background heating rate should be 9.6e-29 * n*/
		/* >>chng 04 jan 26, update cosmic ray ionization rate for H0 to
		 * >>refer	cosmic ray	ionization	Williams, J.P., Bergin, E.A., Caselli, P., 
		 * >>refercon	Myers, P.C., & Plume, R. 1998, ApJ, 503, 689,
		 * H0 ionization rate of 2.5e-17 s-1 and a H2 ionization rate twice this
		 * >>chng 04 mar 15, comment said 2.5e-17 which is correct, but code produce 8e-17,
		 * fix back to correct value 
		 */
		/* NB - the rate is derived from the density.  these two are related by the secondary
		 * ionization efficiency problem.  background_rate is only here to provide the relationship
		 * for predominantly neutral gas.  the background_density is the real rate. 
		 hextra.background_density = 1.99e-9f;*/
		/* >>chng 05 apr 16, to get proper ionization rate in ism_set_cr_rate, where
		 * H is forced to be fully atomic, no molecules, density from 1.99 to 2.15 */
		/* >>chng 12 apr 05, update to
		 * >>refer	cosmic ray 	ionization	Indriolo, N., Geballe, T., Oka, T., & McCall, B.J. 2007, ApJ, 671, 1736
		 */
		background_density = 2.15e-9f*7.9f;
		background_rate = 2.5e-17f*7.9f;
	}
};
extern t_hextra hextra;

#endif /* HEXTRA_H_ */
