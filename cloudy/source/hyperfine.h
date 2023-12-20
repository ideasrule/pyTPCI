/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYPERFINE_H_
#define HYPERFINE_H_

#include "module.h"

/** HyperfineCreate - read in data files and create space for hyperfine lines,
 * called by atmdat_readin at start of calculation */
void HyperfineCreate(void);

/** HyperfineCS - returns interpolated collision strength for element nelem and ion ion 
\param i
*/
/*double HyperfineCS( long nelem , long ion );*/
double HyperfineCS(  size_t i  );

class TransitionProxy;

/** HyperfineTspin
 *
 * Compute spin (excitation) temperature for HFS line.
 * The spin temperature is always finite (> 0, or < 0 if masing).
 * The kinetic temperature is used if the spin temperature evaluates to 0,
 * or if the upper limit has a population that is too small (~< 2e-36 cm^-3)
 * to reliably compute the spin temperature.
 *
 * \param t [in]   Transition
 * \returns Tspin
 */
double HyperfineTspin( const TransitionProxy &t );

/** H21_cm_pops - fine level populations for 21 cm with Lya pumping included*/
void H21_cm_pops( void );

/**H21cm_H_atom computes rate for H 21 cm from upper to lower excitation by atomic hydrogen 
 * from 
 * >>refer	H1	cs	Allison, A.C. & Dalgarno, A., 1969, ApJ 158, 423 
 \param temp temperature
 */
double H21cm_H_atom( double temp );

/**H21cm_proton - evaluate proton spin changing H atom collision rate, from
*>>refer	21cm	p coll	Furlanetto, S.R. & Furlanetto, M.R. 2007,
*>>refcon	MNRAS, doi:10.1111/j.1365-2966.2007.11921.x 
\param temp temperature
*/
double H21cm_proton( double temp );

/**H21cm_electron computes rate for H 21 cm from upper to lower excitation by electrons 
 * >>refer	H1	cs	Smith, F.J., 1966, Planet. Space Sci 14, 929 
 \param temp
 */
double H21cm_electron( double temp );

struct t_hyperfine : public module {

	const char *chName() const
	{
		return "hyperfine";
	}

	void zero();
	void comment(t_warnings&) {}

	/** the isotope abundances relative to main species abundances, as read in from the hyperfine.dat file */
	vector<realnum> HFLabundance;

	/** this is the computed 21cm spin temperature */
	double Tspin21cm;

	/** option to turn off Lya pumping of 21 cm */
	bool lgLya_pump_21cm;

	/** total cooling due to all hyperfine lines */
	double cooling_total;

	/** largest relative cooling due to hyperfine structure lines */
	realnum cooling_max;

	/** type and enum for determining shape of Lya source function at line center, for 21 cm pumping.
	* EXCITATION: 2p / 1s excitation temperature
	* KINETIC: gas kinetic temperature
	* CONSTANT: S_nu = constant */
	typedef enum { EXCITATION, KINETIC, CONSTANT } LyaSourceFunctionShape;
	LyaSourceFunctionShape LyaSourceFunctionShape_assumed;

	};
extern t_hyperfine hyperfine;

#endif /* HYPERFINE_H_ */
