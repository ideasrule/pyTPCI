/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* CS_VS80 compute thermally averaged collision strength for collisional deexcitation 
 * of hydrogenic atoms, from
 * >>refer	H1	collision	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
 *hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients */
/*Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for ions 
 * with all n, and Z*/
#include "cddefines.h"
#include "dense.h"
#include "phycon.h"
#include "iso.h"
#include "helike_cs.h"
#include "hydro_vs_rates.h"
#include "thirdparty.h"
#include "lines_service.h"
#include "integrate.h"

inline double get_iso_statw( const long ipISO, const long n )
{
	if( ipISO == ipH_LIKE )
	{
		return ipDOUBLET * double( pow2(n) );
	}
	else if( ipISO == ipHE_LIKE )
	{
		if( n == 1 )
		{
			return 1.;  // ipSINGLET * pow2( 1 );
		}
		else
		{
			return (ipSINGLET+ipTRIPLET) * double( pow2(n) );
		}
	}
	else
		TotalInsanity();
}

long global_ipISO, global_nelem, global_nHi, global_nLo;
double kTRyd, global_deltaE;


STATIC double hydro_vs_coll_str( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul, long nelem, long Collider, double energy );

STATIC double expe1( double x);

namespace {
	class my_Integrand : public std::unary_function<double, double> 
	{
	public:
		long nelem, Collider;
		double Aul;
		double temp;
		long nHi, gHi, nLo, gLo;
		double IP_Ryd_Hi, IP_Ryd_Lo;
		
		double operator() (double EOverKT) const
			{
				double col_str = hydro_vs_coll_str( nHi, gHi, IP_Ryd_Hi, nLo, gLo, IP_Ryd_Lo, 
					Aul, nelem, Collider, EOverKT * EVRYD * temp/TE1RYD );
				return exp( -1.*EOverKT ) * col_str;
			}
	};
}

/*
 Neither of these rates can be modified to account for impact by non-electrons because they 
 are fits to thermally-averaged rates for electron impact...it can not be unravelled to 
 expose a projectile energy that can then be scaled to account for other projectiles.
 Instead, we have included the original cross section formula (eq 14) from 
 Vriens & Smeets (1980) below.
*/

/* VS80 stands for Vriens and Smeets 1980 */
/* This routine calculates thermally-averaged collision strengths. */
double CS_VS80( long ipISO, long nHi, double IP_Ryd_Hi, long nLo, double IP_Ryd_Lo, double Aul, long nelem, long Collider, double temp )
{
	double coll_str;

	double gLo = get_iso_statw( ipISO, nLo );
	double gHi = get_iso_statw( ipISO, nHi );

	if( Collider == ipELECTRON )
	{
		coll_str = hydro_vs_deexcit( nHi, gHi, IP_Ryd_Hi, nLo, gLo, IP_Ryd_Lo, Aul );
	}
	else
	{
		/* evaluate collision strength, two options,
		 * do thermal average (very slow) if set with
		 * SET COLLISION STRENGTH AVERAGE command,
		 * default just do point at kT
		 * tests show no impact on test suite, much faster */
		if( iso_ctrl.lgCollStrenThermAver )
		{
			my_Integrand func;

			func.nHi = nHi;
			func.gHi = gHi;
			func.IP_Ryd_Hi = IP_Ryd_Hi;
			func.nLo = nLo;
			func.gLo = gLo;
			func.IP_Ryd_Lo = IP_Ryd_Lo;
			func.nelem = nelem;
			func.temp = temp;
			func.Collider = Collider;
			func.Aul = Aul;

			Integrator<my_Integrand,Gaussian32> VS80( func );
			/* do average over Maxwellian */
			coll_str = VS80.sum( 0., 1. );
			coll_str += VS80.sum( 1., 10. );
		}
		else
		{
			/* the argument to the function is equivalent to evaluating the function at
			 * hnu = kT */
			coll_str = hydro_vs_coll_str( nHi, gHi, IP_Ryd_Hi, nLo, gLo, IP_Ryd_Lo, 
					Aul, nelem, Collider, EVRYD*temp/TE1RYD );
		}
	}

	ASSERT( coll_str >= 0. );
	return coll_str;
}

/*hydro_vs_coll_str compute collision strength for collisional deexcitation for hydrogen atom, 
 * from Vriens and Smeets */
STATIC double hydro_vs_coll_str( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul, long nelem, long Collider, double energy )
{
	DEBUG_ENTRY( "hydro_vs_coll_str()" );

	// number of colliders is 4 in def, should be macro
	ASSERT( Collider >= 0.&& Collider <4 );
	double reduced_mass = dense.AtomicWeight[nelem]*colliders.list[Collider].mass_amu/
		(dense.AtomicWeight[nelem]+colliders.list[Collider].mass_amu)*ATOMIC_MASS_UNIT;

	/* This comes from equations 14,15, and 16 of Vriens and Smeets. */
	/* >>refer he-like cs Vriens, L. \& Smeets, A. H. M. Phys. Rev. A, 22, 940 */ 
	/* Computes the Vriens and Smeets
	 * rate coeff. for collisional dexcitation between any two levels of H.
	 * valid for all nhigh, nlow
	 * at end converts this excitation rate to collision strength	*/

	double n = (double)nHi;
	double p = (double)nLo;

	double g_n = (double)gHi;
	double g_p = (double)gLo;

	double ryd = EVRYD;
	double s = fabs(n-p);
	ASSERT( s > 0. );

	double Epn = EVRYD * (IP_Ryd_Lo - IP_Ryd_Hi);
	double Epi = EVRYD * IP_Ryd_Lo;

	/* This is an absorption oscillator strength. */
	double abs_osc_str = GetGF( Aul, Epn*RYD_INF/EVRYD, g_n)/g_p;
	double rEpn = 1./Epn;
	double Apn = 2.*ryd*rEpn*abs_osc_str;

	double rp = 1./p;
	double bp = rp*(1.4*log(p) - .7 + rp*(- .51 + rp*(1.16 - 0.55*rp)));

	double Bpn = 4.*ryd*ryd/pow3(n)*rEpn*rEpn*(1. + Epi*rEpn*(4./3. + bp*Epi*rEpn));

	double delta = exp(-safe_div(Bpn,Apn)) - 0.4*Epn/ryd;


	/* Scale the energy to get an equivalent electron energy. */
	energy *= colliders.list[ipELECTRON].mass_amu/colliders.list[Collider].mass_amu;

	double cross_section;
	/* cross section in units of PI*a_o^2 */
	if( energy/2./ryd+delta <= 0 /*|| energy < Epn*/ )
	{
		cross_section = 0.;
	}
	else
	{
		double gamma = ryd*(8. + 23.*POW2(s*rp))*s*s/
			( s*s*(8. + 1.1*n*s) + 0.8 + 0.4*powpq(n*s,3,2)*fabs(s-1.0) );
		cross_section = 2.*ryd/(energy + gamma)*(Apn*log(energy/2./ryd+delta) + Bpn);
		cross_section = MAX2( cross_section, 0. );
	}

	/* convert to collision strength */
	double coll_str = ConvCrossSect2CollStr( cross_section * PI*BOHR_RADIUS_CM*BOHR_RADIUS_CM, g_n, energy/EVRYD, reduced_mass );

	ASSERT( coll_str >= 0. );

	return( coll_str );
}

/*hydro_vs_coll_recomb generate hydrogenic collisional recombination rate coefficients */
double hydro_vs_coll_recomb( double ionization_energy_Ryd, double Te, double stat_level, double stat_ion )
{
	double coef, 
	  denom, 
	  epi, 
	  t_eV;

	DEBUG_ENTRY( "hydro_vs_coll_recomb()" );

	/* This is equation 9 of
	 * >>refer	H1	collision recomb	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940 */

	/* this is kT in eV */
	t_eV = Te/EVDEGK;

	/* this is the ionization energy relative to kT, dimensionless */
	epi = ionization_energy_Ryd * EVRYD / t_eV;

	/* this is the denominator of equation 8 of VS80. */
	denom = pow(epi,2.33) + 4.38*pow(epi,1.72) + 1.32*epi;

	/* this is equation 9 of VS80 */
	coef = 3.17e-27 / pow3(t_eV) * stat_level / stat_ion / denom;

	ASSERT( coef >= 0. );
	return( coef );
}


/*hydro_vs_ioniz generate hydrogenic collisional ionization rate coefficients */
double hydro_vs_ioniz( double ionization_energy_Ryd, double Te )
{
	double coef, 
	  denom, 
	  epi, 
	  t_eV;

	DEBUG_ENTRY( "hydro_vs_ioniz()" );

	/* a function written to calculate the rate coefficients 
	 * for hydrogen collisional ionizations from
	 * Jason Ferguson, summer 94
	 * valid for all n
	 * >>refer	H1	collision	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
	 * */

	/* this is kT in eV */
	t_eV = Te/EVDEGK;

	/* this is the ionization energy relative to kT, dimensionless */
	epi = ionization_energy_Ryd * EVRYD / t_eV;

	/* this is the denominator of equation 8 of VS80. */
	denom = pow(epi,2.33) + 4.38*pow(epi,1.72) + 1.32*epi;

	/* this is equation 8 of VS80 */
	coef = 9.56e-6 / powpq(t_eV,3,2) * dsexp( epi ) / denom;

	ASSERT( coef >= 0. );
	return( coef );
}

/*Hion_coll_ioniz_ratecoef calculate hydrogenic ionization rates for all n, and Z*/
double Hion_coll_ioniz_ratecoef(
		/* the isoelectronic sequence */
		long int ipISO ,
		/* element, >=1 since only used for ions 
		 * nelem = 1 is helium the least possible charge */
		long int nelem,
		/* principal quantum number, > 1
		 * since only used for excited states */
		long int n,
		double ionization_energy_Ryd,
		double Te )
{
	double H, 
	  HydColIon_v, 
	  Rnp, 
	  charge,
	  chim, 
	  eone, 
	  etwo, 
	  efour,
	  g, 
	  rate, 
	  rate2, 
	  boltz,
	  t1, 
	  t2, 
	  t3, 
	  t4, 
	  tev, 
	  xn, 
	  y;
	static const double arrH[4] = {1.48,3.64,5.93,8.32};
	static const double arrRnp[8] = {2.20,1.90,1.73,1.65,1.60,1.56,1.54,1.52};
	static const double arrg[10] = {0.8675,0.932,0.952,0.960,0.965,0.969,0.972,0.975,0.978,0.981};

	static double small = 0.;

	DEBUG_ENTRY( "Hion_coll_ioniz_ratecoef()" );
	/*calculate hydrogenic ionization rates for all n, and Z
	 * >>refer	HI	cs	Allen 1973, Astro. Quan. for low Te.
	 * >>refer	HI	cs	Sampson and Zhang 1988, ApJ, 335, 516 for High Te.
	 * */

	charge = nelem - ipISO;
	/* this routine only for ions, nelem=0 is H, nelem=1 he, etc */
	ASSERT( charge > 0);
	ASSERT( n>1 );

	if( n > 4 )
	{
		H = 2.15*n;
	}
	else
	{
		H = arrH[n-1];
	}

	if( n > 8 )
	{
		Rnp = 1.52;
	}
	else
	{
		Rnp = arrRnp[n-1];
	}

	if( n > 10 )
	{
		g = arrg[9];
	}
	else
	{
		g = arrg[n-1];
	}

	tev = EVRYD/TE1RYD*Te;
	xn = (double)n;
	chim = EVRYD * ionization_energy_Ryd;
	y = chim/tev;
	boltz = dsexp( chim/tev );

	eone = e1(y);
	etwo = e2(y);
	efour = expn(4,y);

	t1 = 1./xn*eone;
	t2 = 1./xn*efour;
	t3 = 3.*H/xn/(3. - Rnp)*(y*etwo - 2.*y*eone + boltz);
	t4 = 3.36*y*(eone - etwo);
	rate = 7.69415e-9*sqrt(Te)*9.28278e-3*powi(xn/(charge+1),4)*g*y;
	rate *= t1 - t2 + t3 + t4;
	rate2 = 2.1e-8*sqrt(Te)/chim/chim*dsexp(2.302585*5040.*chim/Te);

	/* don't let the rates go negative */
	rate = MAX2(rate,small);
	rate2 = MAX2(rate2,small);

	/* Take the lowest of the two, they fit nicely together... */
	/* >>chng 10 Sept 02, sometimes one of these is zero and the other is positive.
	 * in that case take the bigger one.	*/
	if( rate==0. || rate2==0. )
		HydColIon_v = MAX2(rate,rate2);
	else
		HydColIon_v = MIN2(rate,rate2);

	ASSERT( HydColIon_v >= 0. );
	return( HydColIon_v );
}

/*hydro_vs_deexcit compute collisional deexcitation rates for hydrogen atom, 
 * from Vriens and Smeets 1980 */
double hydro_vs_deexcit( long nHi, long gHi, double IP_Ryd_Hi, long nLo, long gLo, double IP_Ryd_Lo, double Aul )
{
	DEBUG_ENTRY( "hydro_vs_deexcit()" );

	double kT_eV = EVRYD * phycon.te / TE1RYD;

	/* This comes from equations 24 of Vriens and Smeets. */
	/* >>refer he-like cs Vriens, L. \& Smeets, A. H. M. Phys. Rev. A, 22, 940 */ 
	/* Computes the Vriens and Smeets
	 * rate coeff. for collisional dexcitation between any two levels of H.
	 * valid for all nhigh, nlow
	 * at end converts this excitation rate to collision strength	*/

	double n = (double)nLo;
	double p = (double)nHi;

	ASSERT( n!=p );

	double g_n = (double)gLo;
	double g_p = (double)gHi;

	double ryd = EVRYD;
	double s = fabs(n-p);

	double Enp = EVRYD * (IP_Ryd_Lo - IP_Ryd_Hi);
	double rEnp = 1./Enp;
	double Eni = EVRYD * IP_Ryd_Hi;

	ASSERT( Enp > 0. );

	/* This is an absorption oscillator strength. */
	double abs_osc_str = GetGF( Aul, Enp*RYD_INF/EVRYD, g_p)/g_n;
	double Anp = 2.*ryd*rEnp*abs_osc_str;

	double rn = 1./n;
	double bn = rn*( 1.4*log(n) - .7 +rn*(- .51 + rn*(1.16 - 0.55*rn)));
	
	double Bnp = 4.*ryd*ryd/pow3(p)*rEnp*rEnp*(1. + Eni*rEnp*(4./3. + bn*Eni*rEnp));

	double delta_np = exp(-safe_div(Bnp,Anp)) + 0.06 * s*s / (p*n*n);

	double rate;

	if( 0.3*kT_eV/ryd+delta_np <= 0 )
	{
		rate = 0.;
	}
	else
	{
		double Gamma_np = ryd*log(1. + n*n*n*kT_eV/ryd) * (3. + 11.*s*s*rn*rn) * s * s /
			( s*s* (6. + 1.6*p*s) + 0.3 + 0.8*powpq(p*s,3,2)*fabs(s-0.6) );

		rate = 1.6E-7 * sqrt(kT_eV) * g_n / (g_p * ( kT_eV + Gamma_np ) ) *
			( Anp * log(0.3*kT_eV/ryd + delta_np) + Bnp );
	}

	double col_str = rate / COLL_CONST * phycon.sqrte * g_p;

	return col_str ;
}

/* used original equation (6-12) from Fujimoto (1978) IPPJ-AM-8, Institute of Plasma Physics,
* Nagoya University, Nagoya.*/
double hydro_Fujimoto_deexcit(long ipISO, long nHi, long nLo, double Aul,
				double ip_Ryd_Hi, double ip_Ryd_Lo)
{
	DEBUG_ENTRY( "hydro_Fujimoto_deexcit()" );

	double ecollryd = phycon.te/TE1RYD;
	double deltaE_Ryd = abs(ip_Ryd_Hi - ip_Ryd_Lo);
	double osc_strength = Aul / (TRANS_PROB_CONST*POW2(deltaE_Ryd/WAVNRYD));
	double u = deltaE_Ryd/ecollryd;
	double beta = 0.5;
	int gamma = 7;
	double delta = 0.2;

	double onepdelta = 1.+delta;
	double oscgamma = powpq(osc_strength/deltaE_Ryd,-gamma,10);

	double k = 2.19e-10*pow2(1./ecollryd)*phycon.sqrte*osc_strength;

	double factor1 = (log(onepdelta)+expe1(onepdelta*u))/u;
	double factor2 = exp(-2.*beta*oscgamma)/(beta*oscgamma+u);
	double factor3 = log(onepdelta)+expe1(onepdelta*(beta*oscgamma +u));

	double gLo = get_iso_statw( ipISO, nLo );
	double gHi = get_iso_statw( ipISO, nHi );
	double rate = k * gLo/gHi * (factor1 - factor2*factor3 );
	double col_str = rate/COLL_CONST*phycon.sqrte*gHi;

	ASSERT (col_str >= 0);

	return col_str;

}


/* Beigman and Lebedev cross sections from Lebedev and Beigman (1988)
 * Physics of Highly-Excited Atoms and Ions, Springer Series on Atoms and Plasmas, vol. 22, Springer chapter 8 eq. 8.30
 * Rates are excitation, but easily converted to collision strengths
 */
double hydro_Lebedev_deexcit(long ipISO, long nelem, long nHi, long nLo, double IP_Ryd_Lo)
{
	DEBUG_ENTRY( "hydro_Fujimoto_deexcit()" );

	/*core charge*/
	double Z = (double)nelem-(double)ipISO + 1.;
	double deltan = (double)nHi-(double)nLo;
	double npn = (double)nHi + (double)nLo;
	double npn2 = npn*npn;
	double xn = IP_Ryd_Lo*TE1RYD/phycon.te;
	double e1expxn = expe1(xn);

	//
	// these are the two terms in the curly brackets in Eq. (8.30)
	//
	double phi = 2.*pow2(nHi*nLo/npn2)/deltan*(4*deltan - 1.)*e1expxn/deltan;
	phi += pow3(2.*nLo/npn)*npn*(deltan - 0.6)*(1.33+pow2(nLo)*deltan)*(1.-xn*e1expxn)/deltan/pow2(nLo)/pow2(nHi);

	double theta = phycon.te / TE1RYD/pow2(Z);
	double ftheta = log(1.+nLo*theta/(Z*deltan*sqrt(theta)+2.5));
	ftheta /= log(1.+nLo*sqrt(theta)/Z/deltan);

	//
	// NB NB
	//
	// Eq. (8.30) gives the excitation rate coefficient, and involves the
	// term exp(-DE/kT).  Omit that term here to obtain its product with the
	// excitation rate coeff, q(exc)*exp(+DE/KT), stored in rate_x_exp.
	//
	// Omission gives the impression we compute the CS incorrectly.
	// Tripped over this twice already.  Added comment for clarity.
	//
	double rate_x_exp = 2.* pow2(BOHR_RADIUS_CM) * FINE_STRUCTURE * SPEEDLIGHT * SQRTPI;
	rate_x_exp *= double(nLo)*pow3((double)nHi/(deltan*Z))/sqrt(theta)*phi*ftheta;

	double gLo = get_iso_statw( ipISO, nLo );

	double col_str = rate_x_exp / COLL_CONST * phycon.sqrte * gLo;

	ASSERT (col_str >= 0);

	return col_str;
}

/* Van regemorter formula for allowed transitions. Van Regemorter, 1962, ApJ, 136, 906
 * The interval 0.005 < y < infty is interpolated from the results of Table 2
 * from Van Regemorter paper and adjusted at high energies to avoid discontinuities. */
double hydro_vanRegemorter_deexcit( const long ipISO, const long nelem, const long nHi, const long lHi, const long lLo,
					const double EnerWN, const double deltaE_eV, const double Aul,
					const long ipCollider )
{
	DEBUG_ENTRY( "hydro_vanRegemorter_deexcit()" );

	/* ensure that the transition is allowed */
	if ( lHi > 0 && lLo >0 && abs(lHi - lLo) !=1 )
		return 0.;

	double Py = 1.;
	double y = deltaE_eV*EVDEGK/phycon.te;
	static const double valy[11] =
		{
			log10(0.005), -2., log10(0.02), log10(0.04), -1., log10(0.2),
			log10(0.4), 0., log10(2.), log10(4.), 1.
		};

	double a1 = sqrt(3.)/(2.*PI) * e1(0.005);

	if ( nelem == ipHYDROGEN || nelem == ipHELIUM )
	{
		if (y <= 0.005)
			Py = sqrt(3.)/(2.*PI) * e1(y);
		else if (y <= 10.)
		{
			//Py = 0.128384/sqrt(y);

			static const double val[11] =
				{
					log10(a1), log10(1.16), log10(0.956), log10(0.758),
					log10(0.493), log10(0.331), log10(0.209), -1.,
					log10(0.063), log10(0.040), log10(0.023)
				};
			Py = linint(valy,val,11,log10(y));
			Py = exp10(Py);
		}
		else
			Py = 0.066/sqrt(y);

		/*if(nHi==nLo+1)
		 	fprintf(ioQQQ,"vrgm nhi %li, nlo %li, y %g\n",nHi,nLo,y);*/
	}
	else
	{
		if (y <= 0.005)
			Py = sqrt(3.)/(2.*PI) * e1(y);
		else if (y <= 10.)
		{
			static const double val[11] =
				{
					log10(a1), log10(1.16), log10(0.977), log10(0.788),
					log10(0.554), log10(0.403), log10(0.290), log10(0.214),
					log10(0.201), log10(0.2), log10(0.2)
				};
			Py = linint(valy,val,11,log10(y));
			Py = exp10(Py);
			//Py = 0.154023 + 0.1099165/sqrt(y);
		}
		else
			Py = 0.200;
	}

	double massratio = reduced_amu( nelem, ipCollider ) / ELECTRON_MASS;

	double CStemp = 20.6*Aul/pow3(EnerWN)/phycon.sqrte*Py;

	double gHi = get_iso_statw( ipISO, nHi );
	double factor = ( COLL_CONST * powpq(massratio, -3, 2) ) / phycon.sqrte / gHi;

	/* convert to collision strength */
	CStemp /= factor;

	return CStemp;
}

/* exp(xn)*e1(xn) is a function which slowly goes to 0 for big xn,
 * to avoid rounding errors and crashes, it is approximated to an asymptotic
 * form for xn >300.
 *
 * 	xn  e1(xn)		 exp(xn)	 exp(xn)*e1(xn)
 * 	--	------		---------	--------------
 * 10  4.15697e-06  22026.5      9.156e-2
 * 50  3.78326e-24  5.18471e+21  1.962e-2
 * 100 3.6836e-46   2.68812e+43  9.901e-3
 * 200 6.88523e-90  7.22597e+86  4.975e-3
 * 300 1.71038e-133 1.94243e+130 3.322e-3
 * 400 4.77601e-177 5.22147e+173 2.494e-3
 * 500 1.42208e-220 1.40359e+217 1.996e-3
 */
STATIC double expe1( double x)
{

	double e1exp=0.;
	if (x < 300.)
		e1exp = exp(x)*e1(x);
	else
		e1exp = (x*x+4.0634*x+1.15198)/(x*x+5.03637*x+4.19160)/x;

	return e1exp;
}

double CS_ThermAve_PR78(long ipISO, long nelem, long nHi, long nLo, double deltaE, double temp )
{

	double coll_str;

	DEBUG_ENTRY( "CS_ThermAve_PR78()" );

	global_ipISO = ipISO;
	global_nelem = nelem;
	global_nHi = nHi;
	global_nLo = nLo;
	global_deltaE = deltaE;

	kTRyd = temp / TE1RYD;

	if( !iso_ctrl.lgCS_therm_ave[ipISO] )
	{
		/* Do NOT average over Maxwellian */
		coll_str = CS_PercivalRichards78( kTRyd );
	}
	else
	{
		/* DO average over Maxwellian */
		coll_str =  qg32( 0.0, 1.0, Therm_ave_coll_str_int_PR78);
		coll_str += qg32( 1.0, 10.0, Therm_ave_coll_str_int_PR78);
	}

	return coll_str;
}

/* The integrand for calculating the thermal average of collision strengths */
double Therm_ave_coll_str_int_PR78( double EOverKT )
{
	double integrand;

	DEBUG_ENTRY( "Therm_ave_coll_str_int_PR78()" );

	integrand = exp( -1.*EOverKT ) * CS_PercivalRichards78( EOverKT * kTRyd );

	return integrand;
}

inline double C2_PR78(double x, double y)
{
	return pow2(x)*log1p(2./3.*x)/(2.*y + 1.5*x);
}

double CS_PercivalRichards78( double Ebar )
{
	DEBUG_ENTRY( "CS_PercivalRichards78()" );

	if( Ebar < global_deltaE )
		return 0.;

	long ipISO = global_ipISO;
	long nelem = global_nelem;
	double np = (double)global_nHi;
	double n  = (double)global_nLo;
	double n2 = pow2(n);
	double nnp = n*np;
	double Ebar2 = pow2(Ebar);
	double s = np - n;
	long is = global_nHi - global_nLo ;

	ASSERT( s > 0. );

	double A = 8./(3.*s) * pow3(np/(s*n)) * (0.184 - 0.04/pow2(cbrt(s))) * powi(1.-0.2*s/nnp, 1+2*is);

	double Z3 = (double)(nelem + 1 - ipISO);
	double Z32 = pow2(Z3);

	double D = exp(-Z32/(nnp*Ebar2));

	double L = log((1.+0.53*Ebar2*nnp/Z32) / (1.+0.4*Ebar));

	double F = powi(1.-0.3*s*D/nnp, 1+2*is);

	double G = 0.5*pow3(Ebar*n2/(Z3*np));

	double h1 = sqrt(2. - pow2(n/np));
	double h2 = 2.*Z3/(n2*Ebar);
	double xPlus = h2/(h1+1.);
	double xMinus = h2/(h1-1.);

	double y = 1./(1.-D*log(18.*s)/(4.*s));

	double H = C2_PR78(xMinus,y) - C2_PR78(xPlus,y);

	/* this is the LHS of equation 1 of PR78 */
	double cross_section  = A*D*L + F*G*H;
	/* this is the result after solving equation 1 for the cross section */
	cross_section *= PI*pow2(n2*BOHR_RADIUS_CM/Z3) / Ebar;


	double stat_weight = get_iso_statw( ipISO, n );

	/* convert to collision strength */
	return cross_section*stat_weight*Ebar / (PI*pow2(BOHR_RADIUS_CM));
}

#if	0
STATIC void TestPercivalRichards( void )
{
	double CStemp;

	/* this reproduces Table 1 of PR78 */
	for( long i=0; i<5; i++ )
	{
		double Ebar[5] = {0.1, 0.4, 0.8, 1.0, 10.};

		CStemp = CS_PercivalRichards78( 0, 2, 12, 10, Ebar[i] );
	}

	/* this reproduces Table 2 of PR78 */
	for( long i=0; i<5; i++ )
	{
		double Ebar[5] = {0.1, 0.4, 0.8, 1.0, 10.};

		CStemp = CS_ThermAve_PR78( ipISO, 0, N_(ipHi), N_(ipLo), phycon.te );
	}

	return;
}
#endif

