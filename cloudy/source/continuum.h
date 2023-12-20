/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CONTINUUM_H_
#define CONTINUUM_H_

#include "module.h"

/**ContCreatePointers create pointers for lines and continua, one time per coreload */
void ContCreatePointers();

/**ContSetIntensity derive intensity of incident continuum */
void ContSetIntensity();

/**IncidentContinuumHere derive intensity of incident continuum*/
void IncidentContinuumHere();

/** set up continuum energy mesh if first call, otherwise reset to original mesh */
void ContCreateMesh();

/**ContNegative sanity check for negative continuum intensities */
void ContNegative();

/**ffun evaluate total flux for sum of all continuum sources 
 \param anu photon energy (Rydberg) where continuum is evaluated 
 \param frac_beam_time fraction of beamed continuum that is varies with time
 \param frac_beam_const fraction of beamed continuum that is constant
 \param frac_isotropic fraction of continuum that is isotropic
*/ 
double ffun(
			/* the energy in Rydbergs where the continuum will be evaluated */
			double anu , 
			/* fraction of beamed continuum that is varies with time */
			double *frac_beam_time,
			/* fraction of beamed continuum that is constant */
			double *frac_beam_const,
			/* fraction of continuum that is isotropic */
			double *frac_isotropic );

/**ffun version without fractions */
double ffun(double anu);

/**ffun1 derive flux at a specific energy, for one continuum 
\param anu photon energy (Rydberg) where continuum is evaluated 
 */
double ffun1(double xnu);

/*outsum sum outward continuum beams */
void outsum(double *outtot, double *outin, double *outout);

/**cont_gaunt_calc do table look up of gaunt factor 
\param temp
\param z
\param photon
*/
double cont_gaunt_calc(double, double, double);

struct t_continuum : public module {
	const char *chName () const
	{
		return "continuum";
	}
	void zero();
	void comment(t_warnings&) {}

	/** flag saying that parts of continuum are zero */
	bool lgCon0,
	  lgCoStarInterpolationCaution;

	/** TotalLumin is total intensity in incident continuum erg cm-2 s-1 */
	double TotalLumin, 
	  totlsv;

	/** the incident continuum at Hb, FUV out of damped Lya, and La */
	realnum cn4861, 
	  cn1216, 
	  cn1367,
	  cn2066,
	  sv4861, 
	  sv1216,
	  sv2066,
	  sv1367;

	realnum fluxv;
	realnum fbeta;

	/** these are number, labels, and bounds of continuum bands
	 * they are specified in continuum_bands.ini in the data dir */
	long int nContBand;
	vector<string> chContBandLabels;
	vector<realnum> ContBandWavelength;
	vector<long> ipContBandLow, ipContBandHi;
	/** these are fractions of first and last bin to include in the 
	 * band */
	vector<realnum> BandEdgeCorrLow, BandEdgeCorrHi;

	/** option to include isotropic continua in output */
	bool lgPrtIsotropicCont;

	/** this is highest energy where k-shell opacities are counted
	 * can be adjusted with "set kshell" command */
	long int KshellLimit;
	realnum EnergyKshell;

	/* set check energy every zone to check energy balance, slow */
	bool lgCheckEnergyEveryZone;
};

extern t_continuum continuum;


/** SpeciesPseudoContCreate - initialize requested pseudo-continua */
void SpeciesPseudoContCreate();

/** SpeciesPseudoContAccum - accumulate pseudo-continua */
void SpeciesPseudoContAccum();

/** addUniqueSpeciesBand - add unique combination of species and bands file to list
 *
 *  Note that this command requires that the atomic data be read in.
 *
 *  \param [in] filename     name of bands file
 *  \param [in[ speciesLabel species whose emission to accumulate
 */
void addUniqueSpeciesBand( const string &filename, const string &speciesLabel );

/** SpeciesBandsCreate - initialize requested species bands files */
void SpeciesBandsCreate();

/** SpeciesBandsAccum - accumulate emission in species bands */
void SpeciesBandsAccum();

#endif /* CONTINUUM_H_ */
