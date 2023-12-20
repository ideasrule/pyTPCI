/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef FREEBOUND_H_
#define FREEBOUND_H_

#define NUM_DR_TEMPS	19

class freeBound
{
public:
	/** a set of array indices for all atoms on the iso sequences */
	long int ipIsoLevNIonCon;

	/** ionization potential of level N in Ryd */
	double xIsoLevNIonRyd;

	/** radiative recombination rate coefficient, RadRecomb[fcn]
	 * RadRecomb[ipRecEsc] escape prob
	 * RadRecomb[ipRecNetEsc] net escape prob, accounting for absorption
	 * RadRecomb[ipRecRad] rate coef, cm^3 s^-1
	 * */
	double RadRecomb[3];

	/** total radiative recombination continuum, units erg cm-3 s-1 */
	double RadRecCon;
	/** radiative recombination cooling coefficient, units erg cm3 s-1 */
	double RadRecCoolCoef;

	/** state specific dielectronic recombination rate coefs, units cm^3 s^-1 */
	double DielecRecomb;

	/** state specific dielectronic recombination rate coefs as a function of temperature
	 * DielecRecombVsTemp[Temp] rate coef, units cm^3 s^-1 */
	double DielecRecombVsTemp[NUM_DR_TEMPS];

	/** all processes from level n to the continuum, units s-1 */
	double RateLevel2Cont;

	/** all processes from the continuum to level n, units s-1 */
	double RateCont2Level;

	/** ipOpac pointers for photoionization cross sections */
	long int ipOpac;

	/** continuum to total opacity factors for each level */
	double ConOpacRatio;

	/** lte population of each level, cm^3 */
	double PopLTE;

	/** collisional ionization rate coefficient from each level (cm3 s-1) */
	double ColIoniz;

	/** photoionization rate, s-1 */
	double gamnc;

	/** RecomInducRate will become induced recombination rate coefficient
	 * when multiplied by lte population.
	 * integral of photorate times exp(-hu/kt)
	 * for ind rec, produced by gamma routine needs to be multiplied
	 * by lte pop to become real rate  */
	double RecomInducRate;

	/** RecomInducCool_Coef becomes rate coef for induced recombination cooling,
	 * when multiplied by lte population.
	 * this times hnu-hnuo0 to get cooling,
	 * evaluated in gamma routine and saved */
	double RecomInducCool_Coef;

	/** photoelectric heating rate */
	double PhotoHeat;

	/** error in sum of As. */
	double SigmaAtot;

	/** effective recombination and standard deviation in it */
	double RadEffec, SigmaRadEffec;
	
	void Reset();
};


#endif /* FREEBOUND_H_ */

