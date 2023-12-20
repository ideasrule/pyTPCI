/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HELIKE_CS_H_
#define HELIKE_CS_H_

#include "iso.h"


/**reduced_amu - Calculate reduced mass in AMU
 *
 * \param nelen   [in]   element index (0 for H)
 * \param Collider[in]   collision partner
 *
 * \returns  reduced mass
 */
double reduced_amu( long nelem, long Collider );

/**HeCollid evaluate collisional rates 
\param nelem
*/
void HeCollid( long int nelem);

/**HeCSInterp interpolate on He1 collision strengths 
\param nelem
\param ipHi
\param ipLo
\param Collider
*/
realnum HeCSInterp( long int nelem,
				 long int ipHi,
				 long int ipLo,
				 long int Collider );

/**GetHelikeCollisionStrength calculate collision strengths for any transition of He-like iso sequence
\param nelem
\param ipCollider
\param nHi
\param lHi
\param sHi
\param jHi
\param gHi
\param IP_Ryd_Hi
\param nLo
\param lLo
\param sLo
\param jLo
\param gLo
\param IP_Ryd_Lo
\param Aul
\param tauLo
\param EnerWN
\param EnerErg
\param where
*/
realnum GetHelikeCollisionStrength( long nelem, long Collider,
					long nHi, long lHi, long sHi, long jHi, long gHi, double IP_Ryd_Hi,
					long nLo, long lLo, long sLo, long jLo, long gLo, double IP_Ryd_Lo,
					double Aul, double tauLo, double EnerWN, double EnerErg, const char **where );

/* Three different collision treatments, based on
 * Seaton 1962;
 * Pengelly and Seaton 1964; and
 * Vrinceanu and Flannery 2001.
 */

/**CS_l_mixing_PS64 Collision treatment based on Pengelly and Seaton 1964
\param nelem, the chemical element, 1 for He
\param ipISO,
\param tau,
\param target_charge,
\param n,
\param l,
\param gLo,
\param deltaE_eV
\param lp,

\param Collider
*/
double CS_l_mixing_PS64(
	long nelem,
	long ipISO,
	double tau,
	double target_charge,
	long int n,
	long int l,
	double gLo,
	long int lp,
	double deltaE_eV,
	long int Collider);

/**CS_l_mixing_PS64_expI Collision treatment based on Pengelly and Seaton 1964,
 * using the exponential integral
\param nelem, the chemical element, 1 for He
\param ipISO,
\param tau,
\param target_charge,
\param n,
\param l,
\param gLo,
\param lp,
\param s,
\param deltaE_eV
\param lp,
\param Collider
\param lgPSM20
*/
double CS_l_mixing_PS64_expI(
	long nelem,
	long ipISO,
	double tau,
	double target_charge,
	long int n,
	long int l,
	double gLo,
	long int lp,
	//double s,
	double deltaE_eV,
	long int Collider,
	bool lgPSM20);

/**CS_l_mixing_VF01 Collision treatment based on Vrinceanu and Flannery 2001
\param ipISO
\param nelem
\param n
\param l
\param lp
\param s
\param gLo
\param tauLo
\param IP_Ryd_Hi
\param IP_Ryd_Lo
\param temp
\param Collider
*/
double CS_l_mixing_VF01(
	long ipISO,
	long nelem,
	long n,
	long l,
	long lp,
	long s,
	long gLo,
	double tauLo,
	double IP_Ryd_Hi,
	double IP_Ryd_Lo,
	double temp,
	long Collider );

/**CS_l_mixing_VOS12 Collision treatment based on Vrinceanu Onofrio & Sadeghpour 2012
 * corresponding to semiclassical treatment: equation (7)
 http://adsabs.harvard.edu/abs/2012ApJ...747...56V
 */
double CS_l_mixing_VOS12(long n, long l, long lp,
		long nelem, double gLo, long Ztarget, long Collider, double sqrte);

/**CS_l_mixing_VOS12 Collision treatment based on Vrinceanu Onofrio & Sadeghpour 2012
 * corresponding to quantal treatment: equation (2)
 http://adsabs.harvard.edu/abs/2012ApJ...747...56V
 */
double CS_l_mixing_VOS12QM(
	long ipISO,
	long nelem,
	long n,
	long l,
	long lp,
	long s,
	long gLo,
	double tauLo,
	double IP_Ryd_Hi,
	double IP_Ryd_Lo,
	double temp,
	long Collider );

/** CSresolver - this function averages collapsed-to-collapsed collision strengths into
                 collapsed-to-resolved and resolved-to-resolved
\param ipISO
\param nHi
\param lHi
\param sHi
\param nLo
\param lLo
\param sLo
*/
double CSresolver(long ipISO, long nHi,long lHi,long sHi,long nLo,
		long lLo, long sLo, long nHighestResolved);


#endif /* HELIKE_CS_H_ */
