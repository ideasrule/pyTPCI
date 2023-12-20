/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYDRO_VS_RATES_H_
#define HYDRO_VS_RATES_H_


/** VS80 stands for Vriens and Smeets 1980<BR>
    This routine calculates thermally-averaged collision strengths.<BR>
    \param ipISO
    \param nHi
    \param IP_Ryd_Hi
    \param nLo
    \param IP_Ryd_Lo
    \param Aul
    \param nelem
    \param Collider
    \param temp
*/
double CS_VS80( long ipISO, long nHi, double IP_Ryd_Hi, long nLo, double IP_Ryd_Lo,
		double Aul, long nelem, long Collider, double temp );

/**hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients 
 \param ionization_energy_Ryd
 \param Te
 \param stat_level
 \param stat_ion
 */
double hydro_vs_coll_recomb( double ionization_energy_Ryd, double Te, double stat_level, double stat_ion );

/**hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients 
 \param ionization_energy_Ryd
 \param Te
 */
double hydro_vs_ioniz( double ionization_energy_Ryd, double Te );


/**Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for all n, and Z
 \param ipISO the isoelectronic sequence 
 \param nelem element, >=1 since only used for ions<BR>
              nelem = 1 is helium the least possible charge
 \param n 	 principal quantum number, > 1<BR>
		 since only used for excited states<BR> 
 \param ionization_energy_Ryd
 \param temperature
 */
double Hion_coll_ioniz_ratecoef(
		long int ipISO ,
		long int nelem,
		long int n,
		double ionization_energy_Ryd,
		double temperature );

/**hydro_vs_deexcit generate hydrogenic collisional ionization rate coefficients 
 * for quantum number n 
 \param nHi
 \param gHi
 \param IP_Ryd_Hi 
 \param nLo
 \param gLo
 \param IP_Ryd_Lo 
 \param Aul
 */
double hydro_vs_deexcit( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul );

/**hydro_vs_deexcit generate hydrogenic collisional ionization rate coefficients
 * for quantum number n
 \param nelem
 \param ipISO
 \param nHi
 \param nLo
 \param IP_Ryd_Lo
 */
double hydro_Lebedev_deexcit(long ipISO, long nelem, long nHi, long nLo, double IP_Ryd_Lo);

/**hydro_Fujimoto_deexcit - Collisional de-excitation rates from original equation (6-12)
 * from Fujimoto (1978) IPPJ-AM-8, Institute of Plasma Physics, Nagoya University, Nagoya.
 *
 * \param ipISO    [in]  iso-sequence (0 for H-like, 1 for He-like)
 * \param nHi      [in]  upper level principal quantum number
 * \param lHi      [in]  upper level orbital angular momentum
 * \param Aul      [in]  Einstein A for transition between levels
 * \param ip_Ryd_Hi[in]  upper level ionization potential, in Rydberg
 * \param ip_Ryd_Lo[in]  lower level ionization potential, in Rydberg
 *
 * \returns n-n' collision strength
 */
double hydro_Fujimoto_deexcit(long ipISO, long nHi, long nLo, double Aul, double ip_Ryd_Hi, double ip_Ryd_Lo);

/** hydro_vanRegemorter_deexcit - Collisional deexcitation rates according to
 * Van Regemorter, ApJ 136 (1962) 906
 *
 * \param ipISO     [in]  iso-sequence (0 for H-like, 1 for He-like)
 * \param nelem     [in]  element (0 for H)
 * \param nHi       [in]  upper level principal quantum number
 * \param lHi       [in]  upper level orbital angular momentum
 * \param lLo       [in]  lower level orbital angular momentum
 * \param EnerWN    [in]  energy difference between levels, in wavenumbers
 * \param deltaE_eV [in]  energy difference between levels, in eV
 * \param Aul       [in]  Einstein A for transition between levels
 * \param ipCollider[in]  collider type
 *
 * \returns n-n' collision strength
 */
double hydro_vanRegemorter_deexcit(long ipISO,  const long nelem,
					const long nHi, const long lHi, const long lLo,
					const double EnerWN, const double deltaE_eV,
					const double Aul,
					const long ipCollider );

/**CS_ThermAve_PR78 - Calculate thermal averaged Percival and Richard excitation
 * rates
 *
 * \param ipISO [in]  iso-sequence (0 for H-like, 1 for He-like)
 * \param nelem [in]  element (0 for H)
 * \param nHi   [in]  upper level principal quantum number
 * \param deltaE[in]  energy difference between levels, in Rydberg
 * \param temp  [in]  gas temparature
 *
 * \returns n-n' collision strength
*/
double CS_ThermAve_PR78(long ipISO, long nelem, long nHi, long nLo, double deltaE, double temp );

/**Therm_ave_coll_str_int_PR78 The integrand for calculating the thermal average
 * of collision strengths
\param EOverKT
*/
double Therm_ave_coll_str_int_PR78( double EOverKT );

/**C2_PR78
 \param x
 \param y
 */
inline double C2_PR78(double x, double y);

/** CS_PercivalRichards78 calculates collision strength from Percival & Richards
 * (1978) MNRAS 183, 329
 \param Ebar
 */
double CS_PercivalRichards78( double Ebar );


#endif /* HYDRO_VS_RATES_H_ */
