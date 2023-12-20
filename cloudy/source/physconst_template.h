/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

// >>>>>>>> WARNING: do not #define PHYSCONST_TEMPLATE_H_

// >>>>>>>> WARNING: This is a header template file, designed to be included
// on more than one occasion, with the macro NEW_CONSTANT defined to achieve
// different actions on each include, and no include guard
//
// Do not include anything here which isn't required for this single purpose

/**\file physconst_template.h \verbatim
 * physical constants used by Cloudy, mostly taken from
 * >>refer	phys	const	CODATA 2018, https://physics.nist.gov/cuu/Constants/
 * <BR><BR>
 * NB - these are all printed with the "print constants" command, 
 * which is in parse_print.cpp. Using the NEW_CONSTANT macro guarantees
 * that any constant added here is printed there as well \endverbatim
 */

/*********************************************************************
 * first come math constants                                         *
 *********************************************************************/

// Define macros for mathematical operations, undefined at foot of file
// Can be converted back when constexpr() can be relied on (C++11)

#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))

/** the number e */
NEW_CONSTANT( EE, 2.718281828459045235360287 );

/** the Euler constant (aka Euler-Mascheroni constant or gamma) */
NEW_CONSTANT( EULER, 0.577215664901532860606512090082 );

/** the Euler constant parameter exp(-gamma/2) */
NEW_CONSTANT( EXPEULER2, 0.74930600128844902360587 );

/** pi */
NEW_CONSTANT( PI, 3.141592653589793238462643 );

/** 2*pi */
NEW_CONSTANT( PI2, 6.283185307179586476925287 );

/** 4*pi */
NEW_CONSTANT( PI4, 12.56637061435917295385057 );

/** 8*pi */
NEW_CONSTANT( PI8, 25.13274122871834590770115 );

/** sqrt(2) */
NEW_CONSTANT( SQRT2, 1.414213562373095048801689 );

/** sqrt(pi) */
NEW_CONSTANT( SQRTPI, 1.772453850905516027298167 );

/** sqrt(pi/2) */
NEW_CONSTANT( SQRTPIBY2, 1.253314137315500251207883 );

/** ln(2) */
NEW_CONSTANT( LN_TWO, 0.6931471805599453094172321 );

/** ln(10) */
NEW_CONSTANT( LN_TEN, 2.302585092994045684017991 );

/** log(e) */
NEW_CONSTANT( LOG10_E, 0.4342944819032518276511289 );

/** factor that converts optical depth into extinction in mags,
 * 2.5 log e */
NEW_CONSTANT( OPTDEP2EXTIN, 1.085736204758129569127822 );

/** 180/pi */
NEW_CONSTANT( RADIAN, 57.29577951308232087679815 );

/*********************************************************************
 * astronomical constants go here                                    *
 *********************************************************************/

/** solar mass in gram
 * >>refer	phys	const	http://pdg.lbl.gov/2019/reviews/rpp2018-rev-astrophysical-constants.pdf */
NEW_CONSTANT( SOLAR_MASS, 1.98848e33 );

/** solar luminosity erg s-1
 * >>refer	phys	const	http://pdg.lbl.gov/2019/reviews/rpp2018-rev-astrophysical-constants.pdf */
NEW_CONSTANT( SOLAR_LUMINOSITY, 3.828e33 );

/** luminosity of a star with absolute bolometric magnitude 0, in erg s-1
 * >>refer	phys	const	http://pdg.lbl.gov/2019/reviews/rpp2018-rev-astrophysical-constants.pdf */
NEW_CONSTANT( MBOL_ZERO_LUMINOSITY, 3.0128e35 );

/** astronomical unit, cm, nearly the length of the semimajor
 * axis of the Earth's elliptical orbit around the sun */
/* >>refer	phys	const	http://pdg.lbl.gov/2019/reviews/rpp2018-rev-astrophysical-constants.pdf */
NEW_CONSTANT( AU, 1.49597870700e13 );

/*********************************************************************
 * fundamental constants go next, eventually rest should be defined  *
 * in terms of these, these are CODATA 2018 values.                  *
 *********************************************************************/

/** electron mass, gram */
NEW_CONSTANT( ELECTRON_MASS, 9.1093837015e-28 );

/** electron mass, in u */
NEW_CONSTANT( ELECTRON_MASS_U, 5.48579909065e-4 );

/** proton mass, gram */
NEW_CONSTANT( PROTON_MASS, 1.67262192369e-24 );

/** proton mass, in u */
NEW_CONSTANT( PROTON_MASS_U, 1.007276466621 );

/** mass of the alpha particle, in u */
NEW_CONSTANT( ALPHA_MASS_U, 4.001506179127 );

/** this is the Boltzmann factor, erg/K, exact */
NEW_CONSTANT( BOLTZMANN, 1.380649e-16 );

/** speed of light, cm/s, exact */
NEW_CONSTANT( SPEEDLIGHT, 2.99792458e10 );

/** Planck's constant, exact */
NEW_CONSTANT( HPLANCK, 6.62607015e-27 );

/** Avogadro constant, exact */
NEW_CONSTANT( AVOGADRO, 6.02214076e23 );

/** Gravitational constant, cm^3/g/s^2 */
NEW_CONSTANT( GRAV_CONST, 6.67430e-8 );

/** elementary charge, in C (SI units), exact, to use this must convert to cgs */
NEW_CONSTANT( ELEM_CHARGE, 1.602176634e-19 );

/** infinite mass rydberg constant, in cm^-1 */
NEW_CONSTANT( RYD_INF, 1.0973731568160e5 );

/** ionization potential of real hydrogen atom, in inf mass ryd, based on CODATA 2006,
 * uncertainty 10e-12, calculated by Peter van Hoof */
NEW_CONSTANT( HIONPOT, 0.999466508345 );

/** ionization potential of real He^+ ion, in inf mass ryd, based on CODATA 2006,
 * uncertainty 13e-11, calculated by Peter van Hoof */
NEW_CONSTANT( HE2IONPOT, 3.99963199547 );

/** ionization potential of H- in inf mass ryd,
 * >>refer	phys	const	T. Anderson, 2004, Phys. Reports, 394, 157 */
NEW_CONSTANT( HMINUSIONPOT, 0.055432956 );

/** atomic mass unit, gram */
NEW_CONSTANT( ATOMIC_MASS_UNIT, 1.66053906660e-24 );

/*********************************************************************
 * below here should be derived constants                            *
 *                                                                   *
 * NB - explicit values in comments are approximate                  *
 *      and are not maintained !                                     *
 *********************************************************************/

/** molar mass constant, g/mol */
NEW_CONSTANT( MOL_MASS_CONST, AVOGADRO*ATOMIC_MASS_UNIT );

/** number of arcsec in 1 radian, 206264.806 */
NEW_CONSTANT( AS1RAD, RADIAN*3600. );

/** number of square arcsec in 1 steradian, 4.254517e10 */
NEW_CONSTANT( SQAS1SR, pow2(AS1RAD) );

/** number of square arcsec in the whole sky, 5.3463838e11 */
NEW_CONSTANT( SQAS_SKY, PI4*SQAS1SR );

/** parsec in cm, 3.085678e18 */
NEW_CONSTANT( PARSEC, AU*AS1RAD );

/** megaparsec in cm, 3.085678e24 */
NEW_CONSTANT( MEGAPARSEC, 1.e6*PARSEC );

/** h/2pi, 1.05457e-27 */
NEW_CONSTANT( H_BAR, HPLANCK/(2.*PI) );

/** elementary charge, in ESU, 4.8032e-10 */
NEW_CONSTANT( ELEM_CHARGE_ESU, ELEM_CHARGE*SPEEDLIGHT/10. );

/** electric constant, in F/m, 8.854e-12 */
NEW_CONSTANT( ELECTRIC_CONST, 1.e11/(PI4*pow2(SPEEDLIGHT)) );

/** this is the factor that appears in front of Boltzmann factor to get
 * LTE level populations for hydrogenic ions. It is given in the
 * first parts of section 5 of part 2 of hazy, and is actually
 * ( planck^2 / (2 pi m_e k ) )^3/2, but cannot evaluate powers here,
 * must raise this to 3/2 when used, HION_LTE_POP, 5.556e-11 cm^2 K */
NEW_CONSTANT( HION_LTE_POP, pow2(HPLANCK)/(PI2*BOLTZMANN*ELECTRON_MASS) );

/** SAHA is ( h^2/2/pi/m/k )^3/2, is correct constant for free electron
	 SAHA, 4.14132e-16 cm^3 K^(3/2) */
#ifdef HAVE_CONSTEXPR
NEW_CONSTANT( SAHA, sqrt(pow3(HION_LTE_POP)) );
#else
// Need to use explicit constant rather than formula as sqrt() isn't
// guaranteed to be evaluated at compile-time.  Checked by an ASSERT
// in t_physconst::t_physconst() in physconst.cpp
NEW_CONSTANT( SAHA, 4.1413302848114741e-16 );
#endif

/** number of ergs per wavenumber, 1.9864e-16 */
NEW_CONSTANT( ERG1CM, HPLANCK*SPEEDLIGHT );

/** degrees kelvin per unit wavenumber, 1.4388 */
NEW_CONSTANT( T1CM, HPLANCK*SPEEDLIGHT/BOLTZMANN );

/** kJ/mol per unit wavenumber, 0.01196266 */
NEW_CONSTANT( KJMOL1CM, ERG1CM*AVOGADRO/1e10 );

/** ratio of Rydberg constants R_H/R_inf, 0.999456 */
NEW_CONSTANT( H_RYD_FACTOR, 1./(1.+ELECTRON_MASS_U/PROTON_MASS_U) );

/** ratio of Rydberg constants R_He/R_inf, 0.999863 */
NEW_CONSTANT( HE_RYD_FACTOR, 1./(1.+ELECTRON_MASS_U/ALPHA_MASS_U) );

/** number of Ryd per wavenumber, 9.11267e-6 */
NEW_CONSTANT( WAVNRYD, 1./RYD_INF );

/** Angstrom per infinite mass Ryd, 911.2671 */
NEW_CONSTANT( RYDLAM, 1.e8/RYD_INF );

/** ergs per inf mass Ryd, 2.180e-11 */
NEW_CONSTANT( EN1RYD, HPLANCK*SPEEDLIGHT*RYD_INF );

/** the temperature of 1 Rydberg
 te1ryd is h/k is temp of 1 Rydberg, 1.579e5 */
NEW_CONSTANT( TE1RYD, HPLANCK*SPEEDLIGHT*RYD_INF/BOLTZMANN );

/** Kelvins per eV, 1.1604e4 */
NEW_CONSTANT( EVDEGK, ELEM_CHARGE*1.e7/BOLTZMANN );

/** eV per inf mass Ryd, 13.606 */
NEW_CONSTANT( EVRYD, HPLANCK*SPEEDLIGHT*RYD_INF/ELEM_CHARGE*1.e-7 );

/** ergs per eV, 1.602176e-012 */
NEW_CONSTANT( EN1EV, EN1RYD/EVRYD );

/** frequency of one Ryd for infinite mass nuclei, 3.289842e15 */
NEW_CONSTANT( FR1RYD, SPEEDLIGHT*RYD_INF );

/**2 h FR1RYD^3 / c^2 for infinite mass nucleus, 0.5250 */
NEW_CONSTANT( HNU3C2, 2.*HPLANCK*SPEEDLIGHT*pow3(RYD_INF) );

/** frequency of ionization potential of H (not inf mass), 3.288087e15 - never used */
NEW_CONSTANT( FR1RYDHYD, SPEEDLIGHT*RYD_INF*HIONPOT );

/** H_BAR in eV sec, 6.582e-16 */
NEW_CONSTANT( HBAReV, H_BAR/EN1EV );  

/** wavelength (A) of ionization potential of Hydrogen, 911.7535 - never used */
NEW_CONSTANT( RYDLAMHYD, RYDLAM/HIONPOT );

/** Stefan-Boltzmann constant, 5.6704e-5 */
NEW_CONSTANT( STEFAN_BOLTZ, pow2(PI*pow2(BOLTZMANN))/(60.*pow3(H_BAR)*pow2(SPEEDLIGHT)) );

/** the frequency of one eV, 2.418e14 */
NEW_CONSTANT( FREQ_1EV, SPEEDLIGHT*RYD_INF/EVRYD );

/** the fine-structure constant a= 2pi e^2/hc 7.297 352 533 x 10-3 */
NEW_CONSTANT( FINE_STRUCTURE, pow2(ELEM_CHARGE_ESU)/SPEEDLIGHT/H_BAR );

/** the square of the fine-structure constant */
NEW_CONSTANT( FINE_STRUCTURE2, pow2(FINE_STRUCTURE) );

/** Bohr radius in cm, 5.29177249e-9 */
NEW_CONSTANT( BOHR_RADIUS_CM, FINE_STRUCTURE/(PI4*RYD_INF) );

/** the two photon constant as defined by Breit & Teller, as in equation 4 of Spitzer & Greenstein 51, 2.18313 */
NEW_CONSTANT( TWO_PHOT_CONST, 9.*pow3(FINE_STRUCTURE2)*FR1RYD/2048. );

/** this is the square of the value roughly equal to 8.629e-6 that appears in converting 
 * collision strengths to rates. The constant is h^2/((2PI*me)^3/2 * k^1/2). */
NEW_CONSTANT( COLL_CONST, SAHA*BOLTZMANN/HPLANCK );

/** this is the square of the value roughly equal to 4.123e11 that appears in the integration
 * of photoionization cross-sections to obtain recombination coefficients. */
#ifdef HAVE_CONSTEXPR
NEW_CONSTANT( MILNE_CONST, SPEEDLIGHT*sqrt(pow3(FINE_STRUCTURE2)*pow3(TE1RYD)/PI) );
#else
// Need to use explicit constant rather than formula as sqrt() isn't
// guaranteed to be evaluated at compile-time.  Checked by an ASSERT
// in t_physconst::t_physconst() in physconst.cpp
NEW_CONSTANT( MILNE_CONST, 4.1234755895831189e+11 );
#endif

/** This is the constant used in converting oscillator strengths to As. The formula is
 * Aul, TRANS_PROB_CONST * f(u,l) * wavenumber^2. TRANS_PROB_CONST is 0.667025 */ 
NEW_CONSTANT( TRANS_PROB_CONST, PI4*HPLANCK*FINE_STRUCTURE/ELECTRON_MASS );

/** This is the atomic physics constant that enters the calculation of the absorption coefficient
 * given the energy separation between two levels, and the oscillator strength.
 * Numerical value: 0.01497361235435992
 * Units: esu^2 / (g cm/s) = cm^2 / s
 * The absorption coefficient is defined as in Eq. (9-33) of Mihalas, Stellar Atmospheres, 2nd Edition */
NEW_CONSTANT( ABSOR_COEFF_CONST, SQRTPI*pow2(ELEM_CHARGE_ESU)/(ELECTRON_MASS*SPEEDLIGHT) );

/** Thomson cross-section, cm^2 */
NEW_CONSTANT( SIGMA_THOMSON, PI8/3.*pow2(FINE_STRUCTURE*H_BAR/(ELECTRON_MASS*SPEEDLIGHT)) );

/** Express hc in erg Angstroms, approximately equal to 1.9864e-8
 * Used to get h*nu from wavelength in Angstroms*/
NEW_CONSTANT( HC_ERG_ANG, HPLANCK*SPEEDLIGHT*1e8 );

/** Jeans constant as in Hazy 3 -- other constant factors are quoted */
NEW_CONSTANT( JEANS, PI*BOLTZMANN/(GRAV_CONST*ATOMIC_MASS_UNIT) );

/* Free-free emission constant.
 * This is defined as the constant in the FF cooling equation (5.14a) in Rybicki & Lightman
 * divided by the Planck constant.  This is obtained by integrating equation (5.14a) over
 * frequency and changing variables to energy in units of R_inf.  A factor R_inf*c comes out
 * of the integral, which is accounted for by multiplying the constant by EN1RYD = R_inf * h * c,
 * but we need to correct for the extra h.
 * The exact equation is:
 *     32*pi * e^6 / (3 * h * m_e * c^3) * sqrt(2*pi / (3 * k_Boltzmann * m_e))
 */
#ifdef HAVE_CONSTEXPR
NEW_CONSTANT( FREE_FREE_EMIS,
	32. * PI * pow(ELEM_CHARGE_ESU, 6.) /
	(3. * pow3(SPEEDLIGHT) * ELECTRON_MASS * HPLANCK ) *
	sqrt(2. * PI / (3. * BOLTZMANN * ELECTRON_MASS ) ) );
#else
// Need to use explicit constant rather than formula as sqrt() isn't
// guaranteed to be evaluated at compile-time.  Checked by an ASSERT
// in t_physconst::t_physconst() in physconst.cpp
NEW_CONSTANT( FREE_FREE_EMIS, 1.0325265080202923e-11 );
#endif

/* Free-free absorption constant.
 * The constant is obtained by considering the integral of absorption over frequency
 * to compute the total heating rate.  The logic is similar to that for the emission
 * constant above.  NB  The radiation field is already expressed as dnu*Inu, so that
 * leaves 3 powers of (R_inf * c) in the denominator, from the nu^-3 dependence --
 * see Rybicki & Lightman eqn (5.18a).
 *
 * The exact constant is:
 *     4 * e^6 / (3 * h * m * c) * sqrt(2*pi / (3 * k_Boltzmann * m_e)) / (R_inf * c)^3
 */
#ifdef HAVE_CONSTEXPR
NEW_CONSTANT( FREE_FREE_ABS,
	 4. * pow(ELEM_CHARGE_ESU, 6.) /
	(3. * HPLANCK * ELECTRON_MASS * SPEEDLIGHT ) *
	sqrt(2. * PI / (3. * BOLTZMANN * ELECTRON_MASS ) ) /
	pow3(SPEEDLIGHT * RYD_INF) );
#else
// Need to use explicit constant rather than formula as sqrt() isn't
// guaranteed to be evaluated at compile-time.  Checked by an ASSERT
// in t_physconst::t_physconst() in physconst.cpp
NEW_CONSTANT( FREE_FREE_ABS, 1.0369973575937702e-38 );
#endif
