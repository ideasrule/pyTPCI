/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef THERMAL_H_
#define THERMAL_H_

/* thermal.h */

#include "module.h"
#include "depth_table.h"

#define	NCOLNT	10000L

/**TempChange change kinetic temperature, calls tfidle
*/
void TempChange(
			 double TempNew ,
			 /* option to force update of all variables */
			 bool lgForceUpdate);

/**TempChange change kinetic temperature, calls tfidle but does not
 * check on temperature floor or update all variables
*/
void TempChange(
			 double TempNew );

class t_thermal : public module
{
public:
	const char *chName() const
	{
		return "thermal";
	}

	void zero();
	void comment(t_warnings&) {}

	/** flag saying whether to use temperature predictor for next zone, 
	 * in constant density models */
	bool lgPredNextTe;

	/** normally false, is set true if constant temperature model is 
	 * assumed, this can be because the temperature floor was hit */
	bool lgTemperatureConstant;

	/** set true when command parser sees constant temperature assumption,
	 * not set true when temperature floor is hit.  distinguishes
	 * between true constant temperature sim, and one where floor
	 * was hit */
	bool lgTemperatureConstantCommandParsed;

	/** three uses, 
	 * temperature set with constant temperature command,
	 * initial forced temperature with force temperature command
	 * also set if temperature floor is hit */
	realnum ConstTemp;

	/** constant grain temperature */
	realnum ConstGrainTemp;

	/** keep track of possibly thermally unstable models
	 * nUnstable is number of zones that were possibly thermally unstable
	 * lgUnstable says that current conditions possibly unstable */
	long int nUnstable;
	bool lgUnstable;

	/** remember the highest and lowest temperature that occurs in the model */
	realnum thist, 
	  tlowst;

	/** flag set if temperature map is from hot to cool, set with high command */
	bool lgTeHigh;

	/** flag set if energy density of rad field greater tha
	* compton temp - this is unphysical but could be set
	* by users */
	bool lgEdnGTcm;

	/** flag saying that temperature law has been specified with tlaw command */
	bool lgTLaw;

	/** flag to do Berltodi & Draine simple temperature law,
	 * set with tlaw bd96 */
	bool lgTeBD96;
	/** the initial temperature in their equation */
	realnum T0BD96,
		/** the coefficient on column density for temp drop off */
		SigmaBD96;

	/** these incorporate the Sternberg & Neufeld density/temperature relationship */
	realnum T0SN99;
	bool lgTeSN99;

	/** Depth table for temperature */
	DepthTable tlaw;
	bool lgTeTLaw;

	/** remember strongest coolants
	 * these save arrays of line heating and cooling
	 * CoolAdd is called by n level cooling routines, and CoolAdd
	 * fills in cooling (cooling) or heating (heatnt)
	 * heatnt is always positive, although it is negative cooling */
	realnum collam[NCOLNT];
	double cooling[NCOLNT], 
	  heatnt[NCOLNT];
	long int ncltot;
#	define	NCOLNT_LAB_LEN	15
	char chClntLab[NCOLNT][NCOLNT_LAB_LEN+1];

	/* element coolants, the last one is molecular */
	double elementcool[LIMELM + 1];
	/* heavy element collisional ionization cooling
	 * zero with CoolHeavy.colmet, add to element cooling
	 * at cool_save.cpp */ 
	double heavycollcool[LIMELM];

	/* cooling due to level 2 lines */
	double dima;

	/** flag set true during cooling map, saying to keep cool even
	 * if cooling is negative */
	bool lgCNegChk;

	/** max of negative coolants, and a pointer to it */
	realnum CoolHeatMax;
	realnum wlCoolHeatMax;
	char chCoolHeatMax[NCOLNT_LAB_LEN+1];

	/** integrated cooling over model */
	double totcol, 
	  /** cooling in this zone */
	  ctot, 
	  /** heatl is total line heating, t(ipLnHeat) */
	  heatl,
	  /**coolheat is other coolants that were heat sources */
	  coolheat;

	/** derivative of cooling wrt temperature */
	double dCooldT;

	/** derivative of cooling, 1/te^2, .5/T */
	double tsq1, 
	  halfte;

	/** this is set to phycon.te in tfidle, is used to insure that all temp
	 * vars are properly updated when conv_ionizeopacitydo is called 
	 * NB must be same type as phycon.te */
	double te_update;

	/** info about 'extra' cooling, lgCextOn says it is on */
	bool lgCExtraOn;
	realnum CoolExtra, 
	  cextpw;

	/** this flag indicates (true) that we are between when cooling was set to
	 * zero with call to CoolZero, and when final sum was used.  Any call
	 * after final summation in CoolSum, where set (false), would be 
	 * ignored and so is fatal error */
	bool lgCoolEvalOK;

	/** value of, and pointer to, strongest g-bar cooling line */
	realnum GBarMax;
	long int ipMaxExtra;

	/** heating - cooling due to charge transfer ionization / recombination */
	double char_tran_heat , char_tran_cool;
	
	/** total heat input to this zone */
	double	htot,

		/** total energy input over calculated structure, updated in lines */
		power, 

		/** derivative of total heating in this zone, evaluated in SumHeat*/
		dHeatdT;

	/** total free free heating integrated over model */
	double FreeFreeTotHeat;

	/** HeatLineMax is largest fractional heating due to lines */
	realnum HeatLineMax;
private:
	/** heating per unit vol, erg cm^-3 s^-1, heating[nelem][ion]  */
	double m_heating[LIMELM][LIMELM];
public:
	double heating( long nelem, long ion )
	{
		return m_heating[nelem][ion];
	}
	void setHeating( long nelem, long ion, double heating )
	{
		m_heating[nelem][ion] = heating;
	}
	void AddHeating( long nelem, long ion, double heating )
	{
		m_heating[nelem][ion] += heating;
	}
};
extern t_thermal thermal;

struct t_phoHeat
{
	/**HeatNet is heating due to individual species */
	double HeatNet;
	/** this is the part of the heating that cannot do secondary ionizations */
	double HeatLowEnr;
	/** this is the part of the heating that does secondaries, but without efficiency */
	double HeatHiEnr;
};

#endif /* THERMAL_H_ */
