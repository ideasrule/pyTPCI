/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYDROGENIC_H_
#define HYDROGENIC_H_

/**
 * \file hydrogenic.h                                                       
 * this file contains the variables for the model hydrogen ions,            
 * and prototypes for the series of routines that drive the atom            
 * the EXTERN structure hydro is defined here -                             
 * all H variables should migrate here  <BR>                                    
 *                                                                          
 *                                                                          */

#include "module.h"

/**HydroCSInterp calculate collision strengths for all transitions of h-like iso sequence, all colliders 
\param nelem
\param ipHi
\param ipLo
\param Collider
*/
realnum HydroCSInterp( long int nelem, long int ipHi, long int ipLo, long int Collider );

/**GetHlikeCollisionStrength calculate collision strengths for any transition of H-like iso sequence
\param nelem
\param ipCollider
\param nHi
\param lHi
\param sHi
\param gHi
\param IP_Ryd_Hi
\param nLo
\param lLo
\param sLo
\param IP_Ryd_Lo
\param Aul
\param tauLo
\param EnerWN
\param EnerErg
\param where
*/
realnum GetHlikeCollisionStrength( long nelem, long ipCollider,
					long nHi, long lHi, long sHi, long gHi, double IP_Ryd_Hi,
					long nLo, long lLo, long sLo,/*long gLo,*/ double IP_Ryd_Lo,
					double Aul, double tauLo, double EnerWN, double EnerErg, const char **where );

/**HydroLevel calls iso_level to solve for ionization balance level populations of model hydrogen atom 
\param ipZ
*/
void HydroLevel(long ipISO, long int ipZ);

/**HydroRecCool hydrogen recombination cooling 
\param n
\param ipZ
*/
double HydroRecCool(long int n, long int ipZ);

/** returns the ratio of recombination cooling to recombination coefficient 
 \param t the scaled temperature, T * n^2 / Z^2, n is prin quant number, Z is charge, 1 for H
*/
double HCoolRatio( 
	double t );

/**H_cross_section - get Hydrogenic cross section 
 \param EgammaRyd
 \param EthRyd
 \param n
 \param l
 \param nelem
*/
double H_cross_section( double EgammaRyd , double EthRyd, long n, long l, long nelem );

/** all of these are initialized in zero */
struct t_hydro : public module {

	const char *chName() const
	{
		return "hydro";
	}

	void zero();
	void comment(t_warnings&) {}

	/** lgHiPop2 flag set if H n=2 population gets large relative to ground
	 * pop2mx is maximum population of n=2 relative to ground */
	bool lgHiPop2;
	realnum pop2mx;

	/** says whether (true) or not (false) pumping of the FeII model atom by HI Lya is included
	 *  normally true, set false with the 'NO FEII' command */
	bool lgLyaFeIIPumpOn;

	/** dstfe2lya is destruction probability for Lya onto FeII,
	 * net deexcitation of Lya but not ots destruction */
	realnum dstfe2lya;

	/** width of Lya */
	realnum HLineWidth;

	/** TexcLya is the excitation temperature of Lya */
	realnum TexcLya;

	/** nLyaHot is counts how ofter Lya hotter than gas */
	long int nLyaHot;

	/** TLyaMax is hottest */
	realnum TLyaMax, 
	/** TeLyaMax is electron temp at point where Lya max  */
	  TeLyaMax;

	/** nZTLaMax is the zone where this happened */
	long int nZTLaMax;

	/** chHTopType is the method.used to top off the H atom */
	char chHTopType[5];

	/** relative importance of photo ioniz from n=2 of H */
	realnum H_ion_frac_photo;

	/** largest fraction of ground state H destruction due to collisional ionization */
	realnum HCollIonMax;

	/** fraction of H ionizations due to ground collisions */
	realnum H_ion_frac_collis;

	/** cintot is total induced cooling over model */
	double cintot;

	/** lgHInducImp says whether or not induced recombination is important*/
	bool lgHInducImp;

	/** usually 1, set to 0 with hydrogen damping off command, scales rayleigh scat */
	realnum DampOnFac;

	/** remember induced fractions for hydrogen  */
	realnum FracInd;
	long int ndclev;
	realnum fbul;
	long int nbul;

	/** is continuum pumping of H lyman lines included?  yes, but turned off
	 * with atom h-like lyman pumping off command */
	bool lgLymanPumping;

	/** multiplicative scale factor for HI lyman line pump rate, takes into account
	 * possible emission lines - NB test against equal to 1.f in rt_lines_all to
	 * see if it has been set */
	realnum xLymanPumpingScaleFactor;

	};
extern t_hydro hydro;

#endif /* HYDROGENIC_H_ */
