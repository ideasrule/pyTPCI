/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef OPACITY_H_
#define OPACITY_H_

/**\file opacity.h
 routines dealing with creation and evaluration of opacities */

#include "module.h"
#include "vectorize.h"

/** set true when allocated, init to false */
extern bool lgOpacAllocated;

/**OpacityCreateAll compute initial set of opacities for all species */
void OpacityCreateAll(void);

/**OpacityAdd1Subshell add opacity due to single species to main opacity array
\param ipOpac ipOpac is opacity index within opac opacity offset for this species 
\param ipLowLim lower freq limit to opacity range on energy mesh 
\param ipUpLim upper limit to opacity range on energy mesh
\param abundance abundance, we bail if zero
\param chStat either static 's' or volitile 'v'
*/
void OpacityAdd1Subshell(
	long int ipOpac, 
	long int ipLowLim, 
	long int ipUpLim, 
	realnum abundance,
	char chStat );

/**OpacityAddTotal derive total opacity for this position */
void OpacityAddTotal(void);

/**OpacityAdd1Element enter total photo cross section for all subshells 
 * of a single element into opacity array 
 \param ipZ is 0 for H, 1 for He, etc 
 */
void OpacityAdd1Element(
		long int ipZ);

/** OpacityZero - called by OpacityAddTotal to zero opacity array after saving
 * it into oldopac */
void OpacityZero(void);

/** OpacityZeroOld - only set old opac to current value during search phase */
void OpacityZeroOld(void);

/**OpacityAdd1SubshellInduc add opacity of individual species, including stimulated emission 
\param ipOpac pointer to opacity offset with stack
\param low low energy limit to opacity bound 
\param ihi high-energy limit to opacity bound
\param a the abundance of the species in this level
\param b the departure coefficient
\param chStat either 's' for static opacities, or 'v' for volitile 
*/
void OpacityAdd1SubshellInduc(
  long int ipOpac, 
  long int low, 
  long int ihi, 
  double a, 
  double b,
  char chStat );

struct t_opac : public module {

	const char *chName() const
	{
		return "opac";
	}
	void zero();
	void comment(t_warnings&) {}
	/** [0] is optical depth for current position,
	 * [1] is total optical depth from previous iteration */

	/**TauAbsGeo, TauScatGeo, absorption, scattering, optical depths,
	 * at current position, includes through far side of slab*/
	vector<realnum> TauAbsGeo[2];
	vector<realnum> TauScatGeo[2];

	/** TauTotalGeo is total optical depth at each energy, at current position, 
	 * includes through far side of slab */
	vector<realnum> TauTotalGeo[2];

	/** these are the integrated absorption and scattering optical depths
	 * to the illuminated face of the cloud */
	vector<realnum> TauAbsFace, TauScatFace;

	/**local absorption opacity, cm-1*/
	vector<double> opacity_abs;
	 
	/** local scattering opacity, cm-1 */
	vector<double> opacity_sct;

	/**save previous opacity */
	vector<double> OldOpacSave;

	/** albedo is local gas albedo*/
	vector<double> albedo;

	/** initial opacities from zone 1, used to reset opacity at restart */
	/** and the saved value of local opacity*/
	vector<double> opacity_abs_savzon1;
	/** opacity_sct_savzon1 save local opacity at start of calculations*/
	vector<double> opacity_sct_savzon1;

	/** these static opacities are only evaluated 
	 * one time per zone*/
	vector<double> OpacStatic;

	/** density/temp factors needed for free-free opacity */
	vector<double> FreeFreeOpacity;
	vector<double> eeFreeFreeOpacity;

	/** temperature used for the last ee brems evaluation */
	double eeFreeFreeTemp;

	/** exp(-dTau) for the current zone, evaluated in radinc */
	vector_avx<double> ExpZone;

	/** saves value of E2(tau), exponential integral,
	 * where tau is optical depth to illuminated face */
	vector<realnum> E2TauAbsFace;
	/** same things but for outward direction.  Only defined on second and later iteratios */
	vector<realnum> E2TauAbsTotal;
	vector<realnum> E2TauAbsOut;
	/** total absorption optical depth across computed structure */
	vector<realnum> TauAbsTotal;

	/** exp(-tau) to illuminated face */
	vector<realnum> ExpmTau;

	/** factors that account for attenuation of light across this zone
	 * should be nearly unity */
	vector<realnum> tmn;

	/** true by default, and set false with the *no static opacities* command.  
	 * When false always update all opacities */
	bool lgOpacStatic;

	/** this flag is set true in OpacityZero 
	 * when the OpacStatic array is zeroed, 
	 * and is false if the array has been left alone. 
	 * all later opacities must be reevaluated when this
	 * flag is true */
	bool lgRedoStatic;

	/**< this is the stack used to hold opacities - entered one time when code
	 * is initialized, in routine OpacityCreateAll */
	vector<double> OpacStack;

	/** taumin is the smallest optical depths allowed, */
	realnum taumin;

	/** tlamin is smallest Lya optical depth, 
	 * modified to large value if case b used */
	realnum tlamin;

	/** this flag says we are doing case b */
	bool lgCaseB;

	/** this flag turns off n=2 collisional excitations for comparison 
	 * with their paper, keyword HUMMER on caseb command */
	bool lgCaseB_HummerStorey;

	/** this flag turns off excited state photoionization, 
	 * keyword NO PHOTO on caseb command */
	bool lgCaseB_no_photo;

	/** another case b option, turn off background opacities, no Pdest */
	bool lgCaseB_no_pdest;

	/** optical depths to Compton and H- */
	realnum telec, 
	  thmin;

	/** flag set if negative opacities occured */
	bool lgOpacNeg;

	/** flag saying whether (true) or not (false) scattering opacity
	 * is enabled */
	bool lgScatON;

	/** IO unit to save negative opacities */
	bool lgNegOpacIO;

	/** variable dealing with the option to generate and use an
	 * ancillary file containing a stored opacity set.  These
	 * are all initialized in zero1  */
	/**lgCompileOpac flag saying to compile opacity only */
	bool lgCompileOpac;

	/** set false if no file opacity command entered, option to
	 * ignore opacity file */
	bool lgUseFileOpac;

	/** total number of opacity cells used in opacity stack
	 * in OpacityCreateAll used as a counter to remember where
	 * next opacity goes */
	long int nOpacTot;

	/**
	 *
	 * NBNBNBthis must exactly parallel the read/write statements
	 * in OpacityCreateAll
	 */

	/** ipRayScat opacity pointer for Rayleigh scattering*/
	long int ipRayScat, 

	/** iopcom compton scatterin, total recoil*/
	iopcom, 

	/** ippr is threshold for pair production, ioppr is opacity offset*/
	ippr, 
	ioppr, 

	/** ipBrems opacity offset pointer for brems (free-free)*/
	ipBrems, 

	/** iphmra ratio of h minus to neut h ff brems opacity*/
	iphmra, 

	/** iphmop H- bound free opacity*/
	iphmop, 

	/** ih2pnt lower, upper limits to bound, ih2pof, opacity offset*/
	ih2pnt[2], 
	ih2pof, 
	/* excited state versions of above H2+ variables */
	ih2pnt_ex[2], 
	ih2pof_ex, 

	/** iophe1 points for photo from singlet levels*/
	iophe1, 

	/** lowest levels of he triplets*/
	ioptri, 

	/** \verbatim 
	 * ipElement[nelem][ioniz][shell][purpose],
	 * array index on fortran scale, to energy range are set in ipShells, called by ContCreatePointers
	 * pointers to opacity stack offset defined in OpacityCreate1Element where
	 * opacities are set
	 *
	 * first dim 
	 * [0] is pointer to low-energy threshold in energy array
	 * [1] is highest energy for shell, set by LimitSh
	 * [2] is opacity offset within large opacity stack
	 *
	 * if shell does not exist, set upper limit
	 * to less than lower limit so this never looped upon
	 * these are used as flags by LimitSh to check whether
	 * this is a real shell - 
	 * following code is in ipShells for this non-existant case 
	 * OpacPoint.ipElement[nelem][ion][nshell][0] = 2;
	 * OpacPoint.ipElement[nelem][ion][nshell][1] = 1;
	 * all routines must protect against this since opacities for these
	 * undefined shells is also undefined
	 *
	 * second dim is shell, 0 for k shell, up to 6,
	 * (Shell.chShell[ns] give the label for each of these shells)
	 * These quantities are only defined for the number of shells that exist. 
	 * The number of shells is given by Heavy.nsShells[nelem][ion]
	 *
	 * third dimension is ion stage, 0 for atom
	 * 
	 * last dim is atomic number of element, 0 for H
	 *\endverbatim
	 */
	ipElement[LIMELM][LIMELM][7][3], 

	/** in1 is [NI] excited state*/
	in1[3], 

	/** pointers to oxy excited states*/
	ipo3exc[3], 
	ipo3exc3[3], 
	ipo1exc[3];

	/** threshold energies corresponding to ipo3exc[0] and ipo3exc3[0] */
	double o3exc;
	double o3exc3;

	/** photo to excited O+ levels*/
	long iopo2d, 

	/** photoionization from upper lev of Mg II 2798*/
	ipmgex, 
	ipOpMgEx, 

	/** calcium excited states*/
	ica2ex[2], 
	ica2op;

	/** index of carbon k-shell threshold in energy array*/
	long int ipCKshell;

	/** smallest ots rate, set to 0 in scalar */
	realnum otsmin;

	/** max correction for stim emission in continuum opacities at Lyman and Balmer edges */
	realnum stimax[2];

	};
extern t_opac opac;



#endif /* OPACITY_H_ */
