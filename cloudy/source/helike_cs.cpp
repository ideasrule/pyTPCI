/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HeCollid evaluate collisional rates */
/*HeCSInterp interpolate on He1 collision strengths */
/*AtomCSInterp do the atom	*/
/*IonCSInterp do the ions	*/
/*CS_l_mixing_PS64 - find rate for l-mixing collisions by protons, for neutrals */
#include "cddefines.h"
#include "dense.h"
#include "helike.h"
#include "helike_cs.h"
#include "hydro_vs_rates.h"
#include "iso.h"
#include "phycon.h"
#include "thirdparty.h"
#include "thirdparty_quadpack.h"
#include "trace.h"
#include "freebound.h"
#include "lines_service.h"
#include "integrate.h"
#include "vectorize.h"
#include "hydroeinsta.h"
#include "parser.h"

/** vector of temperatures corresponding to collision strengths stuffed into HeCS.	*/
const unsigned int NCSTEMP = 9U;
realnum CSTemp[NCSTEMP];
/** array of collision strengths read from data file...this is interpolated upon.	*/
map<QNPair, array<realnum,NCSTEMP>> HeCS;

STATIC realnum HeCSTableInterp( long nelem, long Collider, 
			long nHi, long lHi, long sHi, long jHi, 
			long nLo, long lLo, long sLo, long jLo );

STATIC double CS_l_mixing_S62( double deltaE_eV, double IP_Ryd_ground, long gLo, double Aul, long nelem, long Collider, double temp );

/* all of these are used in the calculation of Stark collision strengths
 * following the algorithms in Vrinceanu & Flannery (2001). */
template<class P>
class my_Integrand_VF01_E;

template<class P>
STATIC double collision_strength_VF01( long ipISO, double velOrEner,
												   const my_Integrand_VF01_E<P>& vf );

/* These are masses relative to the proton mass of the electron, proton, he+, and alpha particle. */
static const double ColliderCharge[4] = {1.0, 1.0, 1.0, 2.0};

double reduced_amu( long nelem, long Collider )
{
	return dense.AtomicWeight[nelem]*colliders.list[Collider].mass_amu/
		(dense.AtomicWeight[nelem]+colliders.list[Collider].mass_amu)*ATOMIC_MASS_UNIT;	
}

namespace
{
	void compare_VOS12();
	class StarkCollTransProb_VF01;
	class StarkCollTransProb_VOS12QM;
}

template<class P>
class my_Integrand_VF01_E
{
public:
	long ipISO, nelem, n, l, lp, gLo;
	double tauLo, IP_Ryd_Hi, IP_Ryd_Lo;
	long Collider;
	double ColliderCharge, temp;
	double reduced_mass, aveRadius, RMSv, reduced_b_max, ColliderMass;

	my_Integrand_VF01_E( long ipISO_1, long nelem_1, long n_1, long l_1, long lp_1, long s, long gLo_1,
								double tauLo_1, double IP_Ryd_Hi_1, double IP_Ryd_Lo_1, long Collider_1, double temp_1) :
		ipISO(ipISO_1), nelem(nelem_1), n(n_1), l(l_1), lp(lp_1), gLo(gLo_1),
		tauLo(tauLo_1), IP_Ryd_Hi(IP_Ryd_Hi_1), IP_Ryd_Lo(IP_Ryd_Lo_1), Collider(Collider_1), temp(temp_1)
		{
			ColliderCharge = ::ColliderCharge[Collider];
			ASSERT( ColliderCharge > 0. );
			reduced_mass = reduced_amu(nelem, Collider);
			ASSERT( reduced_mass > 0. );
			// Don't police max value of n here: this isn't a validity
			// constraint *for this function*, and it may be useful to
			// use higher values for testing
			if( ! ( s > 0 && n > 0 ) ) // && n <= iso_sp[ipISO][nelem].n_HighestResolved_max ) )
			{
				fprintf( ioQQQ, "invalid parameter for my_Integrand_VF01_E\n" );
				cdEXIT(EXIT_FAILURE);
			}
			/* this is root mean squared velocity */
			/* use this as projectile velocity for thermally-averaged cross-section */
			double vred = sqrt(2.*BOLTZMANN*temp/reduced_mass);
			// This is the n-scaled Bohr radius
			aveRadius = (BOHR_RADIUS_CM/((double)nelem+1.-(double)ipISO))*POW2((double)n);
			ASSERT( aveRadius < 1.e-4 );
			/* >>chng 05 jul 14, as per exchange with Ryan Porter & Peter van Hoof, avoid
			 * roundoff error and give ability to go beyond zinc */
			/*ASSERT( aveRadius >=  BOHR_RADIUS_CM );*/
			ASSERT( aveRadius > 3.9/LIMELM * BOHR_RADIUS_CM );	

			/* vn = n * H_BAR/ m / r = Z * e^2 / n / H_BAR 
			 * where Z is the effective charge. */
			RMSv = ((double)nelem+1.-(double)ipISO)*POW2(ELEM_CHARGE_ESU)/(double)n/H_BAR;
			ASSERT( RMSv > 0. );

			/* From here Pengelly and Seaton cut-offs MNRAS (1964),127-165 are applied
			 * to H and He
			 */

			/* Those are variables used for cut-off calculation */
			double meanb;
			double deltaE_eV= abs(IP_Ryd_Hi - IP_Ryd_Lo)*EVRYD;
			double reduced_b_maxa, reduced_b_maxb;
			if( ipISO == ipH_LIKE )
			{
				/* Debye radius: appears to be too large, results in 1/v^2 variation. */
				/* Reduced here means in units of aveRadius: */

				reduced_b_max = sqrt( BOLTZMANN*temp/ColliderCharge/dense.eden/PI)
					/ (2.*ELEM_CHARGE_ESU)/aveRadius;

				/* lifetime Radius in units of aveRadius */
				reduced_b_maxb = 0.72*vred*tauLo/aveRadius;

				reduced_b_max = MIN2( reduced_b_max, reduced_b_maxb );

			}
			else if( ipISO == ipHE_LIKE )
			{
				double omega_qd;

				/* Keep the precession criterion cut-off for Vrinceanu calculation
				 * as it seems to be large compared with PS cut-offs for all
				 * calculations. If other cut-offs applied it could not be enough for
				 * low temperatures and high densities
				 */
				if(iso_ctrl.lgCS_Vrinceanu[ipISO] && (l != 0 && lp != 0 ))
				{
					double quantum_defect1  = (double)n- (double)nelem/sqrt(IP_Ryd_Lo);
					double quantum_defect2  = (double)n- (double)nelem/sqrt(IP_Ryd_Hi);
				
					/* The magnitude of each quantum defect must be between zero and one. */
					ASSERT( fabs(quantum_defect1)  < 1.0 );
					ASSERT( fabs(quantum_defect1)  > 0.0 );
					ASSERT( fabs(quantum_defect2)  < 1.0 );
					ASSERT( fabs(quantum_defect2)  > 0.0 );
				
					/* The quantum defect precession frequencies */
					double omega_qd1 = fabs( 5.* quantum_defect1 * (1.-0.6*POW2((double)l/(double)n)) / POW3( (double)n ) / (double)l );
					/* >>chng 06 may 30, this needs lp not l. */
					double omega_qd2 = fabs( 5.* quantum_defect2 * (1.-0.6*POW2((double)lp/(double)n)) / POW3( (double)n ) / (double)lp );
					/* Take the average for the two levels, for reciprocity. */
					omega_qd = 0.5*( omega_qd1 + omega_qd2 );
				
					ASSERT( omega_qd > 0. );

					reduced_b_max = sqrt( 1.5 * ColliderCharge * n / omega_qd )/aveRadius;
				}
				else
				{
				// Debye Radius, reduced in units of aveRadius
					reduced_b_maxa = sqrt( BOLTZMANN*temp/ColliderCharge/dense.eden/PI )
											/ (2.*ELEM_CHARGE_ESU)/aveRadius;
				// Lifetime cut-off
					reduced_b_maxb = 0.72*vred*tauLo/aveRadius;

				// Criteria for degeneration of levels
					if (iso_ctrl.lgCS_PSdeg[ipISO] && iso_ctrl.lgCS_VOS12QM[ipISO])
					{
					/*PS64 criterion for degeneration
					 * if PS_crit = 1.55 Rc = 0.72*v*tau
					 */
						double PS_crit = deltaE_eV*tauLo/HBAReV;
						if (PS_crit > 1.55)
							reduced_b_maxb = 1.12*HBAReV*vred/deltaE_eV/aveRadius;
					}
					else if (iso_ctrl.lgCS_B72[ipISO] && iso_ctrl.lgCS_VOS12QM[ipISO])
					{
					/*B72 criterion cut off, Brocklehurst M., 1972, MNRAS, 157, 211
					 * beta = <L>DE/hbar/2/<W> < 4
					 */
						if (dense.xIonDense[0][0] > 0 )
						{
							// mean impact parameter is the minimum of 1/(ion density) or Debye radius
							meanb= powpq(1/dense.xIonDense[0][0],1,3);
							meanb = MIN2(meanb,reduced_b_maxa);
						}
						else
							meanb = reduced_b_maxa;
						double beta_brock= meanb*abs(deltaE_eV)/(HBAReV*4.*PI*vred);
						if (beta_brock >= 0.4)
							reduced_b_maxb = 1.12*HBAReV*vred/deltaE_eV/aveRadius;
					}

					reduced_b_max = min( reduced_b_maxb ,reduced_b_maxa);

					//we need to ensure a large enough integral limit for Vrinceanu
					if (iso_ctrl.lgCS_Vrinceanu[ipISO])
						reduced_b_max = max( reduced_b_maxb ,reduced_b_maxa);

				}
			}
			else
				/* rethink this before using on other iso sequences. */
				TotalInsanity();

			ColliderMass = colliders.list[Collider].mass_amu;
			ASSERT( ColliderMass > 0. );
		}

	double operator() (double EOverKT) const
	{
		/* This reduced mass is in grams.	*/
		double col_str = collision_strength_VF01<P>( ipISO, EOverKT * temp / TE1RYD, *this );
		return exp(-EOverKT) * col_str;
	}
};

// Version of my_Integrand_VF01_E designed to be integrated over y = exp(-EOverkT) in the range (0,1]
// Changing the integrand allows the E integral to be formally extended to E=infty, and means the
// integrand should be more uniform
template<class P>
class my_Integrand_VF01_E_log
{
public:
	my_Integrand_VF01_E<P> data;

	my_Integrand_VF01_E_log( long ipISO_1, long nelem_1, long n_1, long l_1, long lp_1, long s, long gLo_1,
								double tauLo_1, double IP_Ryd_Hi_1, double IP_Ryd_Lo_1, long Collider_1, double temp_1) :
		data(ipISO_1, nelem_1, n_1, l_1, lp_1, s, gLo_1,
			  tauLo_1, IP_Ryd_Hi_1, IP_Ryd_Lo_1, Collider_1, temp_1) {}

	double operator() (double y) const
	{
		double EOverKT = -log(y);
		/* This reduced mass is in grams.	*/
		double col_str = collision_strength_VF01<P>( data.ipISO, EOverKT * data.temp / TE1RYD,
																data );
		return col_str;
	}
};

void HeCollidSetup()
{
	DEBUG_ENTRY( "HeCollidSetup()" );

	if (0)
	{
		compare_VOS12();
		cdEXIT(EXIT_FAILURE);
	}

	/* get the collision strength data for the He 1 lines */
	DataParser d("he1_cs.dat", ES_NONE);

	/* check that magic number is ok */
	d.getline();
	d.checkMagic(COLLISMAGIC);

	/* get the array of temperatures */
	d.getline();
	d.getToken(CSTemp, NCSTEMP);
	d.checkEOL();

	/* now read in the CS data */
	long nelem = ipHELIUM;
	while( d.getline() )
	{
		/* get lower and upper level index */
		long ipLo, ipHi;
		d.getToken(ipLo);
		d.getToken(ipHi);
		if( ipLo >= ipHi )
			d.errorAbort("invalid level indices");

		if( ipHi >= iso_sp[ipHE_LIKE][nelem].numLevels_max - iso_sp[ipHE_LIKE][nelem].nCollapsed_max )
			continue;
		else
		{
			int nLo, jLo, nHi, jHi;
			string confLo, confHi;
			d.getToken(nLo);
			if( nLo <= 0 )
				d.errorAbort("invalid principal quantum number");
			d.getToken(confLo);
			if( confLo.length() != 2 )
				d.errorAbort("invalid configuration string");
			int sLo = confLo[0]-'0';
			if( sLo != 1 && sLo != 3 )
				d.errorAbort("invalid value for S in this configuration");
			int lLo = getL(confLo[1]);
			if( lLo < 0 )
				d.errorAbort("invalid value for L in this configuration");
			d.getToken(jLo);
			d.getToken(nHi);
			if( nHi <= 0 )
				d.errorAbort("invalid principal quantum number");
			d.getToken(confHi);
			if( confHi.length() != 2 )
				d.errorAbort("invalid configuration string");
			int sHi = confHi[0]-'0';
			if( sHi != 1 && sHi != 3 )
				d.errorAbort("invalid value for S in this configuration");
			int lHi = getL(confHi[1]);
			if( lHi < 0 )
				d.errorAbort("invalid value for L in this configuration");
			d.getToken(jHi);
			QNPair inPair(nHi, lHi, sHi, 2*jHi+1, nLo, lLo, sLo, 2*jLo+1);
			realnum* p = HeCS[inPair].data();
			d.getToken(p, NCSTEMP);
			d.checkEOL();
		}
	}
}

/* Choose either AtomCSInterp or IonCSInterp */
realnum HeCSInterp(long nelem,
			long ipHi,
			long ipLo,
			long Collider )
{
	DEBUG_ENTRY( "HeCSInterp()" );

	ASSERT( nelem >= ipHELIUM );
	ASSERT( nelem < LIMELM );

	long nHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].n();
	long lHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].l();
	long sHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].S();
	long jHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].j();
	long gHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].g();
	double IP_Ryd_Hi = iso_sp[ipHE_LIKE][nelem].fb[ipHi].xIsoLevNIonRyd;
	long nLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].n();
	long lLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].l();
	long sLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].S();
	long jLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].j();
	long gLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].g();
	double IP_Ryd_Lo = iso_sp[ipHE_LIKE][nelem].fb[ipLo].xIsoLevNIonRyd;
	double Aul = iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).Emis().Aul();
	// collisions are from high to low level, then initial level lifetime is from higher level
	double tauLo = iso_sp[ipHE_LIKE][nelem].st[ipHi].lifetime();
	double EnerWN = iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).EnergyWN();
	double EnerErg = iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).EnergyErg();

	if( !iso_ctrl.lgColl_excite[ipHE_LIKE] ||
		( nHi==nLo && !iso_ctrl.lgColl_l_mixing[ipHE_LIKE] ) )
	{
		return 0.f;
	}
	const char *where="      ";
	realnum cs = GetHelikeCollisionStrength( nelem, Collider, 
				nHi, lHi, sHi, jHi, gHi, IP_Ryd_Hi,
				nLo, lLo, sLo, jLo, gLo, IP_Ryd_Lo,
				Aul, tauLo, EnerWN, EnerErg, &where );

	ASSERT( cs >= 0.f );

	return cs;
}

realnum GetHelikeCollisionStrength( long nelem, long Collider,
					long nHi, long lHi, long sHi, long jHi, long gHi, double IP_Ryd_Hi,
					long nLo, long lLo, long sLo, long jLo, long gLo, double IP_Ryd_Lo,
					double Aul, double tauLo, double EnerWN, double EnerErg, const char **where )
{
	DEBUG_ENTRY( "GetHelikeCollisionStrength()" );



	bool lgResolvedData = true;

	/* a string used in the output to designate where each cs comes from.	*/
	*where = "      ";

	/* init values, better be positive when we exit */
	realnum cs = -1.f; 

	/* this may be used for splitting up the collision strength within 2^3P */
	realnum j_resolve = 1.f;

	/* Energy difference in eV */
	double deltaE_eV = EnerErg/EN1EV;

	/* lower level is within 2^3P */
	if( nLo==2 && lLo==1 && sLo==3 )
	{
		j_resolve *= (2.f*jLo+1.f) / 9.f;
	}

	/* upper level is within 2^3P */
	if( nHi==2 && lHi==1 && sHi==3 )
	{
		j_resolve *= (2.f*jHi+1.f) / 9.f;
	}

	/* for most of the helium iso sequence, the order of the J levels within 2 3P 
	 * in increasing energy, is 0 1 2 - the exception is atomic helium itself,
	 * which is swapped, 2 1 0 */

	/* this branch is for tabulated data, mostly from Bray et al 2000. */
	if( (cs = HeCSTableInterp( nelem, Collider, nHi, lHi, sHi, jHi, nLo, lLo, sLo, jLo )) >= 0.f ) 
	{
		// These are already j-resolved, so return this to unity.
		if( nLo==2 && lLo==1 && sLo==3 && nHi==2 && lHi==1 && sHi==3 ) 
		{
			j_resolve = 1.f;
		}

		*where = "table";

		ASSERT( cs >= 0.f );
		/* statistical weights included */
	}
	/* this branch is ground to n=2 or from n=2 to n=2, for ions only	*/
	/*>>refer Helike	CS	Zhang, Honglin, & Sampson, Douglas H. 1987, ApJS 63, 487	*/
	/*>>refer Helike	CS Fe	Si et al., 2017, A&A, 600, A85 */
	else if( nelem!=ipHELIUM && nHi==2 && nLo<=2 && Collider==ipELECTRON)
	{
		*where = "Zhang ";
		if( nelem == ipIRON )
			*where = "Si+2017";
		j_resolve = 1.;

		/* Collisions from ground	*/
		if( nLo == 1 )
		{
			if( lHi==0 && sHi==3 ) // to 2tripS
			{
				if( nelem == ipIRON )
					cs = 3.16e-4;
				else
					cs = 0.25f/POW2(nelem+1.f);
			}
			else if( lHi==0 && sHi==1 ) // to 2singS
			{
				if( nelem == ipIRON )
					cs = 1.05e-3;
				else
					cs = 0.4f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==0 ) // to 2tripP0
			{
				if( nelem == ipIRON )
					cs =1.51e-4;
				else
					cs = 0.15f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==1 ) // to 2tripP1
			{
				if( nelem == ipIRON )
					cs = 7.90e-4;
				else
					cs = 0.45f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==2 ) // to 2tripP2
			{
				if( nelem == ipIRON )
					cs = 7.13e-4;
				else
					cs = 0.75f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==1 ) // to 2singP
			{
				if( nelem == ipIRON )
					cs = 4.18e-3;
				else
					cs = 1.3f/POW2(nelem+1.f);
			}
			else
				TotalInsanity();
		}
		/* collisions from 2tripS to n=2	*/
		else if( nLo==2 && lLo==0 && sLo==3 )
		{
			if( lHi==0 && sHi==1 ) // to 2singS
			{
				if( nelem == ipIRON )
					cs = 2.16e-3;
				else
					cs = 2.75f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==0 ) // to 2tripP0
			{
				if( nelem == ipIRON )
					cs = 1.43e-1;
				else
					cs = 60.f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==1 ) // to 2tripP1
			{
				if( nelem == ipIRON )
					cs = 3.84e-1;
				else
					cs = 180.f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==2 ) // to 2tripP2
			{
				if( nelem == ipIRON )
					cs = 6.28e-1;
				else
					cs = 300.f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==1 ) // to 2singP
			{
				if( nelem == ipIRON )
					cs = 2.95e-4;
				else
					cs = 5.8f/POW2(nelem+1.f);
			}
			else
				TotalInsanity();
		}
		/* collisions from 2singS to n=2	*/
		else if( nLo==2 && lLo==0 && sLo==1 )
		{
			if( lHi==1 && sHi==3 && jHi==0 ) // to 2tripP0
			{
				if( nelem == ipIRON )
					cs = 3.99e-4;
				else
					cs = 0.56f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==1 ) // to 2tripP1
			{
				if( nelem == ipIRON )
					cs = 7.23e-2;
				else
					cs = 1.74f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==2 ) // to 2tripP2
			{
				if( nelem == ipIRON )
					cs = 2.05e-3;
				else
					cs = 2.81f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==1 ) // to 2singP
			{
				if( nelem == ipIRON )
					cs = 3.68e-1;
				else
					cs = 190.f/POW2(nelem+1.f);
			}
			else
				TotalInsanity();
		}
		/* collisions from 2tripP0 to n=2	*/
		else if( nLo==2 && lLo==1 && sLo==3 && jLo==0 )
		{
			if( lHi==1 && sHi==3 && jHi==1 ) // to 2tripP1
			{
				if( nelem == ipIRON )
					cs = 5.61e-3;
				else
					cs = 8.1f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==3 && jHi==2 ) // to 2tripP2
			{
				if( nelem == ipIRON )
					cs = 8.2e-3;
				else
					cs = 8.2f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==1 ) // to 2singP
			{
				if( nelem == ipIRON )
					cs = 2.73e-3;
				else
					cs = 3.9f/POW2(nelem+1.f);
			}
			else
				TotalInsanity();
		}
		/* collisions from 2tripP1 to n=2	*/
		else if( nLo==2 && lLo==1 && sLo==3 && jLo==1 )
		{
			if( lHi==1 && sHi==3 && jHi==2 ) // to 2tripP2
			{
				if( nelem == ipIRON )
					cs = 2.46e-2;
				else
					cs = 30.f/POW2(nelem+1.f);
			}
			else if( lHi==1 && sHi==1 ) // to 2singP
			{
				if( nelem == ipIRON )
					cs = 9.39e-3;
				else
					cs = 11.7f/POW2(nelem+1.f);
			}
			else
				TotalInsanity();
		}
		/* collisions from 2tripP2 to n=2	*/
		else if( nLo==2 && lLo==1 && sLo==3 && jLo==2 )
		{
			ASSERT( lHi==1 && sHi==1 );
			/* to 2singP	*/
			if( nelem == ipIRON )
				cs =1.22e-3;
			else
				cs = 19.5f/POW2(nelem+1.f) * (realnum)iso_ctrl.lgColl_l_mixing[ipHE_LIKE];
		}
		else
			TotalInsanity();

		/* statistical weights included */
	}

	/* this branch, n-same, l-changing collisions, with no data */
	else if( (nHi==nLo) && (sHi==sLo) )
	{
		/*This is always resolved, hence the asserts */
		ASSERT( nHi <= iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max );
		ASSERT( lLo >= 0);
		ASSERT( lHi >= 0);


		if( iso_ctrl.lgCS_Seaton[ipHE_LIKE] && Collider == ipELECTRON )
		{
			/* Use the method given in 
			 * >>refer He	CS	Seaton, M. J. 1962, Proc. Phys. Soc. 79, 1105 
			 * statistical weights included. This must be default for non-degenerate case and l<3 */
			/* FG 01/10/2020: Seaton does not deal well with protons, and so it is degraded to deal
			 * with only electrons. The default for l<3 is now PSM.
			 */
			double IP_Ryd_ground = iso_sp[ipHE_LIKE][nelem].fb[0].xIsoLevNIonRyd;
			cs = (realnum)CS_l_mixing_S62( deltaE_eV, IP_Ryd_ground, gHi, Aul, nelem, Collider, (double)phycon.te );
			*where = "S62   ";
		}
		else if( iso_ctrl.lgCS_Vrinceanu[ipHE_LIKE] )
		{

			if( (lLo>=3 && lHi>=3))
			{
				if (lHi != lLo )
				{
					/* Use the method given in
					 * >>refer He	CS	Vrinceanu, D. \& Flannery, M. R. 2001, PhysRevA 63, 032701
					 * statistical weights included.
					 * All multipoles included */
					cs = (realnum)CS_l_mixing_VF01( ipHE_LIKE,
							nelem,
							nLo,
							lHi,
							lLo,
							sLo,
							gHi,
							tauLo,
							IP_Ryd_Hi,
							IP_Ryd_Lo,
							(double)phycon.te,
							Collider );
					*where = "VF01  ";

				}
				else
				{
					cs = 0.f;
					*where = "none  ";
				}

			}
			else
			{
				cs = 0.f;
				*where = "none  ";
			}
		}
		else if( iso_ctrl.lgCS_VOS12[ipHE_LIKE] )
		{
			if(lLo>=3 && lHi>=3)
			{
				/* That avoids problems when dl=0 in case of j changing
				 * collisions as 2^3P_0 -> 2^3P_1 transition
				 * (non degenerate case)
				 */
				if ( lLo != lHi )
				{
					/* Use the method given in
					 * >>refer He	CS	 Vrinceau etal ApJ 747, 56 (2012) equation (7) semiclassical treatment
					 * statistical weights included */
					cs = (realnum)CS_l_mixing_VOS12(
							nHi,
							lHi,
							lLo,
							nelem,
							gHi,
							nelem+1.
							-ipHE_LIKE,
							Collider,
							phycon.sqrte);
					*where = "VOS12 ";
				}
				else
				{
					cs = 0.f;
					*where = "none  ";
				}
			}
			else
			{
				cs = 0.f;
				*where = "none  ";
			}
		}
		else if( iso_ctrl.lgCS_VOS12QM[ipHE_LIKE] )
		{

			if ( lLo != lHi )
			{
				/* Use the method given in
				 * >>refer He	CS	Vrinceau et al ApJ 747, 56 (2012) eq. (2) quantal treatment
				 * statistical weights included */
				cs = (realnum)CS_l_mixing_VOS12QM( ipHE_LIKE,
						nelem,
						nLo,
						lHi,
						lLo,
						sLo,
						gHi,
						tauLo,
						IP_Ryd_Hi,
						IP_Ryd_Lo,
						(double)phycon.te,
						Collider );
				*where = "VOS12Q";
			}
			else
			{
				cs = 0.f;
				*where = "none  ";
			}
			if (0)
			{
				/* This	is for diagnostic, in order to compare different methods */
				long nLo1 = nLo, nHi1 = nHi;
				long lLo1 = lLo;
				long lHi1 = lHi;
				long ipISO1 = ipHE_LIKE;
				double cs1 = (realnum)CS_l_mixing_VF01( ipISO1,
						nelem,
						nLo1,
						lLo1,
						lHi1,
						sLo,
						gLo,
						tauLo,
						IP_Ryd_Hi,
						IP_Ryd_Lo,
						(double)phycon.te,
						Collider );
				double cs2 = (realnum)CS_l_mixing_VF01( ipISO1,
						nelem,
						nHi1,
						lHi1,
						lLo1,
						sHi,
						gHi,
						tauLo,
						IP_Ryd_Lo,
						IP_Ryd_Hi,
						(double)phycon.te,
						Collider );
				double cs3 = (realnum)CS_l_mixing_VOS12QM( ipISO1,
						nelem,
						nHi1,
						lHi1,
						lLo1,
						sHi,
						gHi,
						tauLo,
						IP_Ryd_Lo,
						IP_Ryd_Hi,
						(double)phycon.te,
						Collider );
				fprintf(ioQQQ,"Check VF01 %ld %ld %ld %ld %ld %ld %g: %g (%g) VOS12 %g (%g) VOS12QM %g PS64 %g\n",
						Collider,nLo1,lLo1,lHi1,sLo,nelem,phycon.te,
						cs1,cs2,
						CS_l_mixing_VOS12(nLo1,lLo1,lHi1,nelem,gLo,nelem+1.-ipHE_LIKE,Collider,phycon.sqrte),
						CS_l_mixing_VOS12(nHi1,lHi1,lLo1,nelem,gHi,nelem+1.-ipHE_LIKE,Collider,phycon.sqrte),
						cs3,
						CS_l_mixing_PS64(nelem,ipHE_LIKE,tauLo,nelem+1.-ipISO1,nLo1,lLo1,gHi,lHi1,deltaE_eV,Collider));
				double rate_t1, rate_t2, rate_t;
				double oHi=1./(double)gHi;
				double ogLo=1./(double)gLo;
				double reduced_mass_collider_system = dense.AtomicWeight[nelem]*colliders.list[Collider].mass_amu/
						(dense.AtomicWeight[nelem]+colliders.list[Collider].mass_amu)*ATOMIC_MASS_UNIT;
				double ratef = powpq(ELECTRON_MASS/reduced_mass_collider_system,3,2) * COLL_CONST/phycon.sqrte;

				rate_t1 = cs1*ratef*oHi;
				rate_t2 = cs2*ratef*ogLo;
				rate_t = rate_t1+rate_t2;
				fprintf(ioQQQ,"Rates for H %ld %ld %ld %ld %ld %ld %ld %ld %g: %g %g \n",
						nLo1,lLo1,lHi1,sLo,sHi,jLo,jHi,nelem,phycon.te,rate_t, dense.eden);
				cdEXIT(EXIT_FAILURE);
			}
		}
		/* this branch, l changing by one */
		else if( abs(lHi-lLo)==1 && iso_ctrl.lgCS_PS64[ipHE_LIKE])
		{
			/* >>refer	He	cs	Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165 */
			/* P & S 1964 only consider dipole l-changing */
			/* statistical weights included */
			/*FG 01/10/2020: PS64 and PSM deal with l<3 by default */
			if (iso_ctrl.lgCS_PSClassic[ipHE_LIKE])
			{
				/* Classical Pengelly and Seaton */
				cs = (realnum)CS_l_mixing_PS64(
						nelem,
						ipHE_LIKE,
						tauLo,
						nelem+1.-ipHE_LIKE,
						nLo,
						lHi,
						gHi,
						lLo,
						deltaE_eV,
						Collider);
				*where = "PS64  ";
			}
			else
			{
				/* PS-M: Modified PS method
				 * Refer to F. Guzman et al. MNRAS (2016) 464, 312
				 *
				 * 11/01/22 changed default to PSM20 formulas from
				 * N. Badnell et al. MNRAS (2021) 507, 2922
				 */

				cs = CS_l_mixing_PS64_expI(
						nelem,
						ipHE_LIKE,
						tauLo,
						nelem+1.-ipHE_LIKE,
						nHi,
						lHi,
						gHi,
						lLo,
						deltaE_eV,
						Collider,
						iso_ctrl.lgCS_PSM20[ipHE_LIKE]);

				if(iso_ctrl.lgCS_PSM20[ipHE_LIKE])
					*where = "PSM   ";
				else
					*where = "PSM20 ";

			}
		}
		else
		{
				/* l changes by more than 1, but same-n collision */
				cs = 0.f;
				*where = "none  ";
		}
	}
	/* n-changing collision */
	else if( nHi != nLo )
	{
		if( iso_ctrl.lgCS_None[ipHE_LIKE] )
		{
			cs = 0.f;
			*where = "no gb ";
		}
		/* collisions using Vriens & Smeets only for heavy atoms if not specified */
		/* V&S do not have Z dependence, so is only used for neutrals */
		else if( nelem == ipHELIUM && (iso_ctrl.lgCS_Vriens[ipHE_LIKE] || Collider != ipELECTRON))
		{
			/* >>refer He CS	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940
			 * statistical weight IS included in the routine */
			cs = (realnum)CS_VS80( ipHE_LIKE, nHi, IP_Ryd_Hi, nLo, IP_Ryd_Lo, Aul, nelem, Collider, phycon.te );
			*where = "Vriens";
			lgResolvedData = false;
		}
		/* only electron impact collisions */
		else if (Collider == ipELECTRON)
		{
			if(iso_ctrl.lgCS_Lebedev[ipHE_LIKE])
			{
				/* Lebedev and Beigman (1998) Phys. Highly excited atoms and ions p. 225 eq. 8.30
				 */
				cs = hydro_Lebedev_deexcit(ipHE_LIKE, nelem, nHi, nLo, IP_Ryd_Lo);
				*where = "lebed";
				lgResolvedData = false;
			}
			else if(iso_ctrl.lgCS_Fujim[ipHE_LIKE])
			{

				cs = hydro_Fujimoto_deexcit(ipHE_LIKE, nHi, nLo, Aul, IP_Ryd_Hi, IP_Ryd_Lo);
				*where = "Fuji ";
				lgResolvedData = false;
			}
			else if( iso_ctrl.lgCS_vrgm[ipHE_LIKE])
			{
				cs = hydro_vanRegemorter_deexcit( ipHE_LIKE, nelem, nHi, lHi, lLo,
									EnerWN, deltaE_eV, Aul,
									Collider );
				*where = "vrgm ";
				lgResolvedData = false;
			}
			else if( iso_ctrl.nCS_new[ipHE_LIKE] && nelem==ipHELIUM )
			{
				/* Don't know if stat weights are included in this, but they're probably
				 * wrong anyway since they are based in part on the former (incorrect)
				 * implementation of Vriens and Smeets rates */

				/* two different fits, allowed and forbidden */
				if( Aul > 1. )
				{
					/* permitted lines - large A */
					double x = log10(MAX2(34.7,EnerWN));

					if( iso_ctrl.nCS_new[ipHE_LIKE] == 1 )
					{
					/* this is the broken power law fit, passing through both quantal
					 * calcs at high energy and asymptotically goes to VS at low energies
					 * THIS IS DISCONTINUOUS*/
						if( x < 4.5 )
						{
							/* low energy fit for permitted transitions */
							cs = (realnum)exp10(  -1.45*x + 6.75);
						}
						else
						{
							/* higher energy fit for permitted transitions */
							cs = (realnum)exp10(  -3.33*x+15.15);
						}
					}
					else if( iso_ctrl.nCS_new[ipHE_LIKE] == 2 )
					{
						/* single parallel fit for permitted transitions, runs parallel to VS */
						cs = (realnum)exp10(  -2.3*x+10.3);
					}
					else
						TotalInsanity();
				}
				else
				{
					/* fit for forbidden transitions */
					if( EnerWN < 25119.f )
					{
						cs = 0.631f;
					}
					else
					{
						cs = (realnum)exp10( -3.*log10(EnerWN)+12.8);
					}
				}
				lgResolvedData = false;

				*where = "newgb ";
			}
			else
			{
				/* F. Guzman 07/18/19 Heavy ions impact are 2 to three orders of magnitude lower than electron impacts at nebular typical
				 * temperatures. The heavy ion impact approach using Vriens & Smeets 1980 used here is valid for high energies
				 * and is based in a mass scalling of the collision strength (Burgess and Tully (2005).
				 * Cloudy implementation of Percival & Richards 1978 does not include this treatment. Furthermore the Coulomb
				 * repulsion between ions will make these cross sections smaller only to be relevant at even higher energies.  */

				/* Percival and Richards (1978) have got a Z dependence so their rates are preferred */
				cs = CS_ThermAve_PR78( ipHE_LIKE, nelem, nHi, nLo,
							EnerErg / EN1RYD, phycon.te );
				*where = "PR78  ";

				lgResolvedData = false;
			}
		}
		else
			cs = 0.;
	}
	else if (sHi != sLo)
	{
		/* what's left are deltaN=0, spin changing collisions.
		 * These have not been accounted for.	*/
		/* Make sure what comes here is what we think it is.	*/
		ASSERT( nHi == nLo );
		ASSERT( sHi != sLo );
		cs = 0.f;
		*where = "spin  ";
	}
	else
		TotalInsanity();


	/*
	 * Resolve collision strengths for n-changing collisions
	 * from collapsed to resolved levels
	 */
	if (!lgResolvedData && nLo <= iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max)
	{
		realnum l_resolve = 
				realnum( CSresolver(ipHE_LIKE, nHi, lHi, sHi, nLo, lLo, sLo,
						iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max) );

		enum {DEBUG_LOC = false};
		if( DEBUG_LOC )
		{ 
			if( nelem == ipIRON &&
				//(nHi == iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max + 15) &&
				// nLo == iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
				nLo == 1 && nHi == 3 ) 
				// nHi == 3 && sHi == 1 && (nLo == 1 || nLo == 2) )
				//	sLo == 1 )
			{
				fprintf( ioQQQ, "nelem: %ld"
						"  (nHi, lHi, sHi): (%ld, %ld, %ld) ->"
						"  (nLo, lLo, sLo): (%ld, %ld, %ld)"
						"\t nn' CS: %.4e"
						"\t l-reslv: %.4e"
						"\t j-reslv: %.4e"
						"\t lj-reslv: %.4e"
						"\t final CS: %.4e"
						"\t '%s'\n",
						nelem,
						nHi, lHi, sHi,
						nLo, lLo, sLo,
						cs, l_resolve, j_resolve,
						l_resolve * j_resolve,
						cs * l_resolve * j_resolve,
						*where );
			}
		}

		cs *= l_resolve;
	}

	/* take factor into account, usually 1, ratio of stat weights if within 2 3P 
	 * and with collisions from collapsed to resolved levels */
	cs *= j_resolve;


	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/

		if( DEBUG_LOC && /* ( nelem==ipOXYGEN ) &&*/  (cs > 0.f) ) //&& (iso_sp[ipHE_LIKE][nelem].st[ipHi].n() == 2) 
			//&& ( iso_sp[ipHE_LIKE][nelem].st[ipLo].n() <= 2 ) )
			fprintf(ioQQQ,"DEBUGGG HeCSInterp %li\t%li\t%li\t%li\t-\t%li\t%li\t%li\t%e\t%s\n",
				nelem,
				nLo, sLo, lLo,
				nHi, sHi, lHi,
				cs, *where );
	}

	ASSERT( cs >= 0.f );

	return cs;
}

STATIC realnum HeCSTableInterp( long nelem, long Collider, 
			long nHi, long lHi, long sHi, long jHi, 
			long nLo, long lLo, long sLo, long jLo )
{
	DEBUG_ENTRY( "HeCSTableInterp()" );

	if( nelem!=ipHELIUM )
		return -1_r;
	else if( Collider!= ipELECTRON )
		return -1_r;
	else if( nHi > 5 )
		return -1_r;
	/*else if( nHi > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
	{
		fixit( "HeCS allocation must be changed in order to remove this and do collapsed levels here." );
		return -1.f;
	}
	// Cannot do collapsed levels with unspecified l here.
	else if( lLo < 0 || lHi < 0 )
		return -1_r;*/
	
	realnum cs = 0.;
	long topLo =1;
	long topHi =1;
	long l = lHi;
	long lp = lLo;
	long s = sHi;
	long sp = sLo;
	long j = jHi;
	long jp = jLo;
	long topsHi = 1;
	long topsLo = 1;

	if( nLo > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max)
	{
		topLo = nLo;
		topsLo = 2;
	}
	
	if (nHi > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max)
	{
		topHi = nHi;
		topsHi = 2;
	}

	long gHi = 2*j+1;
	long gLo = 2*jp+1;



	/* The loops ensure the addition of coll. strengths when the routine is called for collapsed levels */
	for (long Hi=0; Hi< topHi; Hi ++)
	{
		realnum csHi=0.;
		for (long su = 0; su<topsHi ; su++)
		{
			realnum csLo = 0.;
			if (nHi > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max)
			{
				l = Hi;
				s = 2*su+1;
			}
			/*The pairs have j=-1 and so g=-1, however this might not be true in future versions of the table */
			/*if (jHi<0)
				gHi = (2*l+1)*s;*/

			for (long Lo=0; Lo < topLo; Lo++)
			{
				realnum csLo_s =0;
				for (long sl = 0 ; sl < topsLo ; sl++)
				{
					realnum cs0 = -1.f;

					if (nLo > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max)
					{
						lp = Lo;
						sp = 2*sl+1;
					}
					/*The pairs have j=-1 and so g=-1, however this might not be true in future versions of the table */
					/*if(jLo < 0)
						gLo = (2*lp+1)*sp;*/

					// if this fails, the upper and lower level were reversed in the call
					ASSERT( HeCS.find(QNPair(nLo, lp, sp, gLo, nHi, l, s, gHi)) == HeCS.end() );

					auto p = HeCS.find(QNPair(nHi, l, s, gHi, nLo, lp, sp, gLo));
					if( p == HeCS.end() || p->second[0] < 0_r )
						return -1_r;

					/* this is the case where we have quantal calculations */
					/* >>refer	He1	cs	Bray, I., Burgess, A., Fursa, D.V., & Tully, J.A., 2000, A&AS, 146, 481-498 */
					/* check whether we are outside temperature array bounds,
					 * and use extreme value if we are */
					cs0 = linint(CSTemp, p->second.data(), NCSTEMP, realnum(phycon.alogte));

					if (nLo > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max && cs0 >= 0.)
						csLo_s += cs0;
					else
					{
						csLo_s = cs0;
						break;
					}

				}

				if (nLo > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max && csLo_s >= 0.)
					csLo += csLo_s;
				else
				{
					csLo = csLo_s;
					break;
				}
			}

			if (nHi > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max && csLo >= 0.)
				csHi += csLo;
			else
			{
				csHi = csLo;
				break;
			}
		}


		if (nHi > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max && csHi >= 0.)
			cs += csHi;
		else
		{
			cs = csHi;
			break;
		}
	}

	ASSERT( cs >= 0_r );
	return cs;
}

class my_Integrand_S62
{
	double deltaE, osc_strength, temp, I_energy_eV;
public:
	my_Integrand_S62(double deltaE, double osc_strength, double temp, double I_energy_eV) :
		deltaE(deltaE), osc_strength(osc_strength), temp(temp), I_energy_eV(I_energy_eV)
	{}
	void operator() (const double proj_energy_overKT[], double res[], long n) const;
};

/*CS_l_mixing_S62 - find rate for l-mixing collisions by protons, for neutrals */
/* The S62 stands for Seaton 1962 */
STATIC double CS_l_mixing_S62( double deltaE_eV, double IP_Ryd_ground, long gLo, double Aul, long nelem, long Collider, double temp )
{
	/* >>refer	He	l-mixing	Seaton, M.J., 1962, Proc. Phys. Soc. */
	DEBUG_ENTRY( "CS_l_mixing_S62()" );

	if( Aul <= iso_ctrl.SmallA )
	{
		return 0.;
	}

	ASSERT( TRANS_PROB_CONST*POW2(deltaE_eV/WAVNRYD/EVRYD) > 0. );

	double osc_strength = Aul / (TRANS_PROB_CONST*POW2(deltaE_eV/WAVNRYD/EVRYD));
	double I_energy_eV = EVRYD * IP_Ryd_ground;
	double reduced_mass(reduced_amu(nelem,Collider));

	// Integral is over [0,infty) in principle
	// Dimensionless parameters for integral can then be taken as
	// 1. kT/I_energy_eV; 2. osc_strength; 3. deltaE/I_energy_eV
	// Tabulation may be possible?  NB 1. is independent of atomic physics 
	// parameters specific to the system.  At present, the integral is independent of
	// the collider type, but seems likely that this will need to be corrected.

	// Major cost is very accurate Bessel function evaluations, which
	// doesn't seem consistent with accuracy of interpolation in
	// deriving zeta
	my_Integrand_S62 func(deltaE_eV, osc_strength, temp, I_energy_eV);

	double energy_factor = phycon.te/EVDEGK * colliders.list[ipELECTRON].mass_amu/colliders.list[Collider].mass_amu;
	double elow = 0., emid = 1.*energy_factor, ehigh = 10.*energy_factor;
	VecIntegrator<my_Integrand_S62,VecGaussian32> S62(func);
	/* This returns a thermally averaged collision strength */
	double coll_str = S62.sum( elow, emid );
	coll_str += S62.sum( emid, ehigh );
	ASSERT( coll_str > 0. );
	coll_str /= energy_factor;

	// Moved constant factors in conversion from cross-section to collision strength out of integral,
	// hence 1.0 for energy argument here
	coll_str = ConvCrossSect2CollStr(coll_str*PI*BOHR_RADIUS_CM*BOHR_RADIUS_CM, (double)gLo, 1.0, reduced_mass);
	return coll_str;
}

//* Calculate beta(F), where F(beta) = K0(beta)^2 + K1(beta)^2 
STATIC double S62BesselInvert(double zOverB2)
{
	DEBUG_ENTRY( "S62BesselInvert()" );
	double betaone;
	if( zOverB2 > 100. )
	{
		betaone = sqrt( 1./zOverB2 );
	}
	else if( zOverB2 < 0.54 )
	{
		/* Low betaone approximation */
		double logz = log(zOverB2/PI);
		betaone = (1./3.)*(1.28-logz);
		if( betaone > 2.38 )
		{
			/* average with this over approximation */
			betaone += -0.5*logz;
			betaone *= 0.5;
		}
	}
	else
	{
		long ip_zOverB2 = 0;
		static const double zetaOVERbeta2[101] = {
			1.030E+02,9.840E+01,9.402E+01,8.983E+01,8.583E+01,8.200E+01,7.835E+01,7.485E+01,
			7.151E+01,6.832E+01,6.527E+01,6.236E+01,5.957E+01,5.691E+01,5.436E+01,5.193E+01,
			4.961E+01,4.738E+01,4.526E+01,4.323E+01,4.129E+01,3.943E+01,3.766E+01,3.596E+01,
			3.434E+01,3.279E+01,3.131E+01,2.989E+01,2.854E+01,2.724E+01,2.601E+01,2.482E+01,
			2.369E+01,2.261E+01,2.158E+01,2.059E+01,1.964E+01,1.874E+01,1.787E+01,1.705E+01,
			1.626E+01,1.550E+01,1.478E+01,1.409E+01,1.343E+01,1.280E+01,1.219E+01,1.162E+01,
			1.107E+01,1.054E+01,1.004E+01,9.554E+00,9.094E+00,8.655E+00,8.234E+00,7.833E+00,
			7.449E+00,7.082E+00,6.732E+00,6.397E+00,6.078E+00,5.772E+00,5.481E+00,5.202E+00,
			4.937E+00,4.683E+00,4.441E+00,4.210E+00,3.989E+00,3.779E+00,3.578E+00,3.387E+00,
			3.204E+00,3.031E+00,2.865E+00,2.707E+00,2.557E+00,2.414E+00,2.277E+00,2.148E+00,
			2.024E+00,1.907E+00,1.795E+00,1.689E+00,1.589E+00,1.493E+00,1.402E+00,1.316E+00,
			1.235E+00,1.157E+00,1.084E+00,1.015E+00,9.491E-01,8.870E-01,8.283E-01,7.729E-01,
			7.206E-01,6.712E-01,6.247E-01,5.808E-01,5.396E-01};

		ASSERT( zOverB2 >= zetaOVERbeta2[100] );

		/* find beta in the table */
		long ilo=0, ihi = 100;
		for(;;)
		{
			long imid = (ilo+ihi)/2.;
			if ( zetaOVERbeta2[imid] > zOverB2 )
				ilo = imid;
			else
				ihi = imid;
			if (ihi == ilo+1)
				break;
		}
		ASSERT( ( zOverB2 < zetaOVERbeta2[ilo] ) && ( zOverB2 >= zetaOVERbeta2[ilo+1] ) );

		ip_zOverB2 = ilo;

		ASSERT( (ip_zOverB2 >=0) && (ip_zOverB2 < 100) );

		const double fp = 1.023292992280754131; // exp10(1./100.);
		betaone = exp10( (double)ip_zOverB2/100. - 1.)*
			( (zOverB2 - zetaOVERbeta2[ip_zOverB2]) * fp 
			 -zOverB2 + zetaOVERbeta2[ip_zOverB2+1] )/ 
			 (zetaOVERbeta2[ip_zOverB2+1] - zetaOVERbeta2[ip_zOverB2]);
	}
	//fprintf(ioQQQ,"%g %g %g\n",betaone,zOverB2,POW2(bessel_k0(betaone))+POW2(bessel_k1(betaone)));
	return betaone;
}

/* The integrand for calculating the thermal average of collision strengths */
void my_Integrand_S62::operator() (const double proj_energy[], double res[], long n) const
{
	DEBUG_ENTRY( "S62_Therm_ave_coll_str()" );

	ASSERT( deltaE > 0. );

	avx_ptr<double> x(n), val(n);
	for( long i=0; i < n; ++i )
		x[i] = -proj_energy[i]*EVDEGK/temp;
	vexp(x.ptr0(), val.ptr0(), 0, n);

	for( long i=0; i < n; ++i )
	{
		/* projectile energy in eV */
		/* Rnot = 1.1229*H_BAR/sqrt(ELECTRON_MASS*deltaE*EN1EV)/Bohr_radius; in units of Bohr_radius */
		/* The deltaE here is to make sure that the collider has no less
		 * than the energy difference between the initial and final levels. */
		double Dubya = proj_energy[i] + 0.5*deltaE;
		ASSERT( Dubya > 0. );

		/* betanot = sqrt((proj_energy[i]+deltaE)/I_energy_eV)*(deltaE/2./Dubya)*Rnot; */

		ASSERT( I_energy_eV > 0. );
		ASSERT( osc_strength > 0. );

		// z/b1**2 = 0.5*W**2/(deltaE**2*Phi_ij)  -- (33)
		// Phi_ij = I_energy_eV*osc_strength/deltaE -- (19)
		/* from equation 33 */
		double zOverB2 = 0.5*Dubya*Dubya/(deltaE*I_energy_eV*osc_strength);

		// zOverB2 = K_0^2+K_1^2 -- (11)
		ASSERT( zOverB2 > 0. );

		double betaone = S62BesselInvert(zOverB2);

		double zeta_of_betaone = zOverB2 * POW2(betaone);

		// Equation (39)
		/* cs1 = betanot * bessel_k0(betanot) * bessel_k1(betanot); */
		double k0val, k1val;
		bessel_k0_k1(betaone, &k0val, &k1val);
		double cs2 = 0.5*zeta_of_betaone + betaone * k0val * k1val;

		/* cross_section = MIN2(cs1, cs2); */
		double cross_section = cs2;

		/* cross section in units of PI * a_o^2 */
		cross_section *= 8. * (I_energy_eV/deltaE) * osc_strength * (I_energy_eV/(proj_energy[i]+deltaE));

		// Maxwell weighting cf Burgess & Tully 1992 A&A 254, 436 eq (21), Seaton (1953)
		// Moved constant factors in conversion from cross-section to collision strength out of integral
		res[i] = val[i] * (proj_energy[i]+deltaE)/EVRYD * cross_section;
	}
}

/*CS_l_mixing_PS64 - find rate for l-mixing collisions by protons, for neutrals */
/* This version does not assume that E_min/kt is small and works with the exponential integral */
/* Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165
 * assume Emin/KT << 1 which limites the range of application of the treatment.
 * This routine is based on N. Badnell developments on PS formulation.
 */
double CS_l_mixing_PS64_expI(
		long nelem,
		long ipISO,
		double tau,
		double target_charge,
		long n,
		long l,
		double g,
		long lp,
		double deltaE_eV,
		long Collider,
		bool lgPSM20)
{
	double cs;
	double RD,R12,Sij,Plowb,RC,RC1,/*R1,*/EC,ED,
	reduced_mass, reduced_mass_2_emass,bracket,contr, rate;
	double fb1, fb2;
	double eEm,eED,eEC,eEmt1Em;

	DEBUG_ENTRY( "CS_l_mixing_PS64_expI()" );

	const double ChargIncoming = ColliderCharge[Collider];
	// Bohr time
	const double tau_zero = 2.41889e-17;
	long lmax = MAX2(l,lp);
	double n2 = (double)n*(double)n;
	double lmax2 = (double)lmax*(double)lmax;
	reduced_mass = reduced_amu(nelem, Collider);
	reduced_mass_2_emass = reduced_mass / ELECTRON_MASS;

	//Projectile energy
	double tempryd=phycon.te/TE1RYD;
	// density in au
	double dens_au = dense.eden*pow(BOHR_RADIUS_CM,3);
	// projectile velocity in au
	double vred = sqrt(tempryd/reduced_mass_2_emass);
	// lifetime in au
	double tau_ua = tau/tau_zero;

	Plowb = 0.5;
	fb1 = 1.;
	fb2 = 1.;
	/*line strength multiplied by 2 (this factor comes from the spin of the electron) */
	Sij = 9.*n2*lmax*(n2-lmax2)/(2.*target_charge*target_charge);

	/* set uf cut offs* and Emin*/
	/*First, R1, low impact parameter cut-off. R12 = R1^2*/
	R12 = 2. * pow2(ChargIncoming)/(3.*Plowb*vred*vred*(2.*l+1));

	R12 *= Sij;
	//R1 = sqrt(R12);

	/* Debye cut-off in au*/
	RD = sqrt(tempryd/(8.*PI*dens_au));

	/*Lifetime cut off in au*/
	RC1 = 0.72*tau_ua*vred;

	//degeneration energy dependent cuts off
	if (ipISO > 0)
	{
		double meanb;
		double PS_crit = deltaE_eV*tau/HBAReV;

		//double aveRadius = (BOHR_RADIUS_CM/((double)nelem+1.-(double)ipHE_LIKE))*POW2((double)n);
		double bmax = sqrt( BOLTZMANN*phycon.te/dense.eden )
								/ (PI2*ELEM_CHARGE_ESU);
		if (dense.xIonDense[0][0] > 0 )
		{
			// mean impact parameter is the minimum of 1/(ion density) or Debye radius
			meanb= powpq(1/dense.xIonDense[0][0],1,3);
			meanb = MIN2(meanb,bmax);
		}
		else
			meanb = bmax;
		double v = sqrt(2.*BOLTZMANN*phycon.te/reduced_mass);
		double beta_brock= meanb*abs(deltaE_eV)/(HBAReV*4.*PI*v);

		/*PS64 criterion for degeneration
		 * if PS_crit = 1.55 Rc = 0.72*v*tau.
		 * PS64 eq. 29 */
		if (PS_crit > 1.55 && iso_ctrl.lgCS_PSdeg[ipHE_LIKE])
			RC1 = 1.12*HBAReV/*EXPEULER2*/*vred/deltaE_eV/tau_zero;

		// Brocklehurst criterion on eq. 3.12 MNRAS (1972) 157, 211
		if (beta_brock >= 0.4 && iso_ctrl.lgCS_B72[ipHE_LIKE])
			RC1 = 1.12*HBAReV/*EXPEULER2*/*vred/deltaE_eV/tau_zero;

			/* Emin for energy dependent cut off */
		//			Emin = R1/RC1;
	}
	RC = min(RC1,RD); /* use the minimum cut-off*/

	/* E/Emin=RC1**2/R12 or
	 * E/Emin = Rc1/R1
	 */
	//if (RC < sqrt(R12) )
	//{
		fb1 = 2./3.;
		R12 = R12*2.;
		//R1 = sqrt(R12);
		//RC = R1;
	//}
	double Emin = R12/(RC*RC); //Um

	if (RC == RD)
			fb2=1.;
	else
		{
			fb2=2.;
			Emin = sqrt(Emin);
		}

	//Emin *= Eproj;

	/* Resolved Dnl depending on Sij */
	double Dnl = 2. * pow2(ChargIncoming)*Sij/(3.*(2.*l+1.));

	ASSERT( phycon.te  / Dnl / reduced_mass_2_emass > 0. );
	ASSERT( Dnl > 0. );

	EC = RD*RD/(RC1*RC1); //Uc
	ED = R12/(RD*RD); //Um
	eEm = exp(-1.*Emin);
	eED = exp(-1.*ED);
	eEC = exp(-1.*EC);
	eEmt1Em = eEm*(1.+Emin);

	if(lgPSM20) 
	/* N. Badnell et al. MNRAS (2021) 507, 2922 */
	{
		double Um = ED;

		double Umc = sqrt(Um*EC);
		//Dnl = 4. * pow2(ChargIncoming)*Sij/(3.*(2.*l+1.));


		if(fb2 == 1)
		{
			/* this the asymptote of sqrt(pi/Um^3)*erf(Um)/2 - exp(-Um)/Um when Um->0*/
			double nf = 2./3.;
			/* we have looked for a value of Um low enough where the error is low
			 * when x=1e-6 rel error = 5.99e-7 */
			bracket = nf+e1(Um);		
			if (Um > 1e-6)
				bracket = sqrt(PI/pow3(Um))*erf(sqrt(Um))/2.-exp(-1.*Um)/Um + e1(Um);

		}
		else
		{
			/* Use of the incomplete gamma function P(n,x)= 1 - e_(n-1)(x)e^-x for (n E Z) to
			 * avoid cancellation errors
			 *
			 * where
			 *
			 *
			 *            n        j
			 *            -       x
			 *  e_n       >      ----
			 *            -       j!
			 *    	      j=0
			 *
			 * Functions 6.5.3, 6.5.11, and 6.5.13 in
			 * Abramowitz & Stegun Tenth Ed. 1972
			 */
			if (EC<Umc)
			/* this is unlikely and correspond to contribution of eq (122) in Nigel's talk
			 * it kicks on extreme high density cases
			 */
			{
				/* this the asymptote of sqrt(pi/Um^3)*erf(Um)/2 - exp(-Um)/Um */
				double nf = 2./3.;
				/* we have looked for a value of Um low enough where the error is low
                         	* when x=1e-6 rel error = 5.99e-7 */
				if (Um > 1e-6)
					nf = sqrt(PI/pow3(Um))*(erf(sqrt(Um))/2.) - exp(-1.*Um)/Um;

				bracket = 4.*(1.-igamc(3,EC)+0.25*exp(-1.*EC)*EC*EC)/pow3(Umc) -
						sqrt(PI/pow3(Um))*erf(sqrt(EC))/2. + nf;
				bracket += e1(Um);

			}
			else //eq 12 Badnell et al. 2021 MNRAS 507, 2922
			{
				bracket = 4.*(1.-igamc(3,Umc))/pow3(Umc);
				bracket += 2.*e1(Umc)-e1(EC);
			}


		}
	}
	else /* Guzman et al. MNRAS (2016) 464, 312 */
	{
		/* First exponential integral is used as the analitical solution of Maxwell averaged PS64 cross sections
		 * Different cases are used depending on the cut-off
		 */
		if ( fb1 == 1 && fb2 == 1)
			bracket = eEm + e1(Emin);
		else
		{
			if (fb2 ==1 )
				bracket = fb1*eEm+ fb2*e1(Emin);
			else
			{
				if (EC > Emin)
					bracket = fb1*eEm + 2.*e1(Emin) - e1(EC);
				else
					bracket = fb1*eED + e1(ED);
			}
		}

		//contribution 0<E<Emin, important for Emin/kt >>
		contr = 0.;
		if(fb1 != 1 )
		{
			if ( fb2 == 1 )
				contr = (1.-eEmt1Em)/Emin;
			else if (fb2 !=1 )
			{
				if (EC >= Emin)
				{
					contr = 2.*( 1. -eEmt1Em)/pow2(Emin);
					contr -= eEm;
				}
				else
				{
					contr = ( 2. -eEC*(2.+EC))/pow2(Emin);
					contr -= eEm*(1.+1./ED);
				}
			}

			contr *= 2./3.;

			bracket += contr;

		}
	}
	ASSERT( bracket >= 0.);


	if (bracket == 0. )
		return SMALLFLOAT;

	/* This is the rate coefficient.   Units: cm^3 s-1	*/

	double units = 2.*pow(BOHR_RADIUS_CM,3)*sqrt(PI)/vred/tau_zero;


	rate = units * Dnl* bracket;

	/* convert rate to collision strength */
	/* NB - the term in parentheses corrects for the fact that COLL_CONST is only appropriate
	 * for electron colliders and is off by reduced_mass_2_emass^-1.5 */
	cs = rate / ( COLL_CONST * powpq(reduced_mass_2_emass, -3, 2) ) * phycon.sqrte * g;


	ASSERT( cs > 0. );

	return cs;
}

/* Classical PS64, refer to eq. 43 Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165
	CS_l_mixing_PS64 - find rate for l-mixing collisions by protons, for neutrals */
double CS_l_mixing_PS64(
	long nelem, /* the chemical element, 1 for He */
	long ipISO, /* iso-sequence index */
	double tau, /* the radiative lifetime of the initial level.	*/
	double target_charge,
	long n,
	long l,
	double gLo,
	long lp,
	double deltaE_eV,
	long Collider)
{
	/* >>refer	H-like	l-mixing	Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165 */
	/* >>refer	He-like	l-mixing	Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165 */
	double cs, 
		TwoLogDebye, TwoLogRc1, 
		factor1, factor2, 
		bestfactor,	factorpart,
		reduced_mass, reduced_mass_2_emass,
		rate;
	const double ChargIncoming = ColliderCharge[Collider];

	DEBUG_ENTRY( "CS_l_mixing_PS64()" );

	/* In this routine, two different cutoff radii are calculated, and as per PS64,
	 * the least of these is selected.  We take the least positive result.	*/

	/* Only dipole l-changing is considered in this approach */
	ASSERT( abs(l-lp) == 1);
	/* This reduced mass is in grams.	*/
	reduced_mass = reduced_amu(nelem, Collider);
	/* this mass always appears relative to the electron mass, so define it that way */
	reduced_mass_2_emass = reduced_mass / ELECTRON_MASS;

	/* equation 46 of PS64 */
	/* min on density added to prevent this from becoming large and negative
	 * at very high n_e - Pengelly & Seaton did not apply this above 1e12 cm^-3 */
	/* This is 2 times the log of the Debye radius.	*/
	TwoLogDebye = 1.68 + log10( phycon.te / MIN2(1e11 , dense.eden ) );
	/* Brocklehurst (1971, equation 3.40) has 1.181 instead of 1.68.  This appears to be due to 
	 * Pengelly and Seaton neglecting one of the two factors of PI in their Equation 33 */
	//TwoLogDebye = 1.181 + log10( phycon.te / MIN2(1e10 , dense.eden ) );

	/* This is 2 times the log of cutoff = 0.72v(tau), where tau is the lifetime of the level nl.	
	 * This is PS64 equation 45 (same as Brocklehurst 1971 equation 3.41)  */
	TwoLogRc1 = 10.95 + log10( phycon.te * tau * tau / reduced_mass_2_emass );

	/* non-degenerate case */
	if (ipISO > 0)
	{
		double meanb;
		double PS_crit = deltaE_eV*tau/HBAReV;

		// This is Debye Radius
		double bmax = sqrt( BOLTZMANN*phycon.te/dense.eden/PI )
								/ (2*ELEM_CHARGE_ESU);
		double vred = sqrt(2.*BOLTZMANN*phycon.te/reduced_mass);

		if (dense.xIonDense[0][0] > 0 )
		{
		// mean impact parameter is the minimum of 1/(ion density) or Debye radius
		meanb= powpq(1/dense.xIonDense[0][0],1,3);
		meanb = MIN2(meanb,bmax);
		}
		else
			meanb = bmax;

		double beta_brock= meanb*abs(deltaE_eV)/(HBAReV*4.*PI*vred);
		double deltaE_cm = deltaE_eV*EN1EV/ERG1CM;

		/*PS64 criterion for degeneration
		 * if PS_crit = 1.55 Rc = 0.72*v*tau.
		 * PS64 eq. 29
		 */
		if (PS_crit > 1.55 && iso_ctrl.lgCS_PSdeg[ipHE_LIKE])
			/* Direct calculation */
			//TwoLogRc1 = 2.*log10(1.12*HBAReV*vred/**EXPEULER2*//deltaE_eV);
			/* Eq. 3.13 Brocklehurst 1971 for degenerate case */
			TwoLogRc1 = log10(phycon.te*ELECTRON_MASS/pow2(deltaE_cm)/reduced_mass)-11.22;

		// Brocklehurst criterion on eq. 3.12 MNRAS (1972) 157, 211
		if (beta_brock >= 0.4 && iso_ctrl.lgCS_B72[ipHE_LIKE])
			/* Direct calculation */
			//TwoLogRc1 = 2.*log10(1.12*HBAReV*vred*EXPEULER2/deltaE_eV);
			/* Eq. 3.13 Brocklehurst 1971 for degenerate case */
			TwoLogRc1 = log10(phycon.te*ELECTRON_MASS/pow2(deltaE_cm)/reduced_mass)-11.22;

	}

	/* this is equation 44 of PS64
	 * unresolved Dnl
	 */
	double Dnl = POW2( ChargIncoming / target_charge) * 6. * POW2( (double)n) *
		( POW2((double)n) - POW2((double)l) - l - 1);


	/***** NB NB NB NB
	 * Brocklehurst (1971) has a factor of electron density in the rate coefficient (equation 3.38).
	 * This is NOT a proper rate, as can be seen in his equations 2.2 and 2.4.  This differs
	 * from the formulism given by PS64. */
	//rate *= dense.eden;


	// Separate up and down rates in PS64 expression using
	// (2l+1) Q_{l,l-1} = (2l-1) Q_{l-1,l} = n^2 l - l^3
	// So (2l+1) Dnl' = (n^2 l - l^3) + (n^2 [l+1] - [l+1]^3) = (2l+1) (n^2 - l^2 - l - 1)

	// l, gLo are value for *initial* state.  It is likely (but not
	// certain) that l < lp.

	long lmax = MAX2(l,lp);
	/* This is resolved l->lp Dnl */
	//	double Dnlup = POW2( ChargIncoming / target_charge) * 6. * POW2( (double)n) *
	//			lmax * (n*n-lmax*lmax)  / double( 2*l+1 );
	double Dnlup = POW2( ChargIncoming / target_charge) * 6. * POW2( (double)n) * lmax * (n*n-lmax*lmax)  / double( 2*l+1 );


	ASSERT( Dnl > 0. );
	ASSERT( phycon.te  / Dnl / reduced_mass_2_emass > 0. );

	/* In the current configuration is the resolved/unresolved Dnl what get's inside the logarithm.
	 * That depends on the implementation. Resolved Dnl has been used in F. Guzman MNRAS (2016) 459, 3498
	 * Summers MNRAS (1977) uses unresolved Dnl. That allows to use resolved results as a branching ratio
	 * of total unresolved PS64 crates.
	 */
	factorpart = (11.54 + log10( phycon.te  / Dnl / reduced_mass_2_emass ) );

	if( (factor1 = factorpart + TwoLogDebye) <= 0.)
		factor1 = BIGDOUBLE;

	if( (factor2 = factorpart + TwoLogRc1) <= 0.)
		factor2 = BIGDOUBLE;

	/* Now we find the least positive result.	*/
	bestfactor = MIN2(factor1,factor2);

	ASSERT( bestfactor > 0. );

	/* If both factors are bigger than 100, toss out the result, and return SMALLFLOAT. */
	if( bestfactor > 100. )
	{
		return SMALLFLOAT;
	}

	/* This is the rate coefficient.   Units: cm^3 s-1	*/
	rate = 9.93e-6 * sqrt( reduced_mass_2_emass  ) * Dnlup / phycon.sqrte * bestfactor;

	/* convert rate to collision strength */
	/* NB - the term in parentheses corrects for the fact that COLL_CONST is only appropriate 
	 * for electron colliders and is off by reduced_mass_2_emass^-1.5 */
	cs = rate / ( COLL_CONST * powpq(reduced_mass_2_emass, -3, 2) ) *
		phycon.sqrte * gLo;

	ASSERT( cs > 0. );

	return cs;
}

/*CS_l_mixing - find collision strength for l-mixing collisions for neutrals */
/* The VF stands for Vrinceanu & Flannery 2001 */
/* >>refer	He	l-mixing	Vrinceanu, D. & Flannery, M. R. 2001, PhysRevA 63, 032701	*/
/* >>refer	He	l-mixing	Hezel, T. P., Burkhardt, C. E., Ciocca, M., He, L-W., */
/* >>refercon	Leventhal, J. J. 1992, Am. J. Phys. 60, 329 */

template<class P>
inline double CS_l_mixing(long ipISO,
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
				long Collider )
{
	double coll_str;

	DEBUG_ENTRY( "CS_l_mixing()" );

	typedef my_Integrand_VF01_E<P> integType;
	integType func(ipISO,nelem,n,l,lp,s,gLo,tauLo,IP_Ryd_Hi,IP_Ryd_Lo,Collider,temp);

	/* no need to do this for h-like */
	if( iso_ctrl.lgCS_Seaton[ipISO] && ipISO > ipH_LIKE )
	{
		ASSERT( l != 0 );
		ASSERT( lp != 0 );
	}

	if( !iso_ctrl.lgCS_VOS_thermal[ipISO] )
	{
		/* Do NOT average over Maxwellian */
		/* Compensate for the exponential weighting factor applied in
		 * the integrand */
		coll_str = func(1.0)*exp(1.0);
	}
	else
	{
		/* DO average over Maxwellian */
		if (1)
		{
			Integrator<integType, Gaussian32> VF01_E(func);
			
			coll_str = VF01_E.sum( 0.0, 1.0 );
			coll_str += VF01_E.sum( 1.0, 10.0 );
		}
		else
		{		  
			double xn=0., xl = 10.;
			if (func(xl) != 0.)
			{
				while (xl < 100. && func(xl/2.) == 0.)
				{
					xn = xl/2.;
					xl *= 2.;
				}
			}
			
			class integrate::Romberg<integrate::Midpoint<integType> >
				VF01_ER(integrate::Midpoint<integType>(func,xn,xl));
			VF01_ER.update(1e-3);
			coll_str = VF01_ER.sum();

			if (0)
			{
				if (0)
				{
					typedef my_Integrand_VF01_E_log<P> integLogType;
					integLogType funcl(ipISO,nelem,n,l,lp,s,gLo,tauLo,IP_Ryd_Hi,IP_Ryd_Lo,Collider,temp);
					double x1n = 0.0, x1x = 1.0;
					do {
						double x1m = 0.5*(x1n+x1x);
						if (funcl(x1m) > 0.)
							x1n = x1m;
						else
							x1x = x1m;
					} while ((x1x-x1n) > 0.001);
					class integrate::Simple<integrate::Midpoint<integLogType> >
						VF01_ESL(integrate::Midpoint<integLogType>(funcl,0.0,x1x));
					VF01_ESL.update(1e-3);

					Integrator<integType, Gaussian32> VF01_E(func);
					double coll_strg = VF01_E.sum( 0.0, 1.0 );
					coll_strg += VF01_E.sum( 1.0, 10.0 );
					
					class integrate::Romberg<integrate::Midpoint<integLogType> >
						VF01_ERL(integrate::Midpoint<integLogType>(funcl,0.0,x1x));
					VF01_ERL.update(1e-3);
					fprintf(ioQQQ,"Int %ld %ld->%ld ER %g(%g) %ld ESL %g(%g) %ld ERL %g(%g) %ld G %g S %g\n",
							  n,l,lp,
							  VF01_ER.sum(),VF01_ER.error(),VF01_ER.evals(),
							  VF01_ESL.sum(),VF01_ESL.error(),VF01_ESL.evals(),
							  VF01_ERL.sum(),VF01_ERL.error(),VF01_ERL.evals(),
							  coll_strg,func(1.0)
						);
					fflush(ioQQQ);
				}
				else
				{
					fprintf(ioQQQ,"Int %ld %ld->%ld ER %g(%g) %ld\n",
							  n,l,lp,
							  coll_str,VF01_ER.error(),VF01_ER.evals()
						);
					fflush(ioQQQ);
				}
			}
		}
	}

	return coll_str;
}

double CS_l_mixing_VF01(long ipISO,
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
				long Collider )
{

		return CS_l_mixing<StarkCollTransProb_VF01>(ipISO,nelem,n,l,lp,s,gLo,tauLo,
				IP_Ryd_Hi,IP_Ryd_Lo,temp,Collider);

}

double CS_l_mixing_VOS12QM(long ipISO,
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
				long Collider )
{
		return CS_l_mixing<StarkCollTransProb_VOS12QM>(ipISO,nelem,n,l,lp,s,gLo,
				tauLo,IP_Ryd_Hi,IP_Ryd_Lo,temp,Collider);
}

namespace
{
	class Chi_VF01 
	{
		double m_alpha, m_alphamin, m_sinfac;
	public:
		Chi_VF01(double alpha, double alphamin) : m_alpha(alpha), m_alphamin(alphamin) 
			{
				ASSERT( m_alpha >= 0. );
				/* deltaPhi is the effective angle swept out by the projector as viewed by the target.  */
				/* deltaPhi is PI for a complete collision in a straight trajectory
				 * Kazansky et al. JPhysB: At. Mol. Opt. Phys. 29, 3651 proposed a core effect correction
				 * This has an unwanted effect in the probability as it decreases at high impact parameters
				 * That shouldn't happen as deflection should be stronger at low impact parameters*/
				
				//double b_over_bmax = m_alphamin / m_alpha;		
				double deltaPhi = -PI;//2.*asin(b_over_bmax)-PI;
				m_sinfac = pow2(sin(0.5*deltaPhi*sqrt(1+m_alpha*m_alpha)));
			}
		double sinChiOver2sq() const
			{	
				/* >>refer He l-mixing Kazansky, A. K. & Ostrovsky, V. N. 1996, JPhysB: At. Mol. Opt. Phys. 29, 3651 */				
				if( m_alpha <= m_alphamin)
				{
					return 0.;
				}
				double denom = (1.+m_alpha*m_alpha), ratio = m_alpha/denom;
				return 4.*ratio*ratio*m_sinfac*(denom-m_alpha*m_alpha*m_sinfac);
			}
		double cosChi() const
			{
				double alphasq = pow2(m_alpha);
				return 1.-2.*alphasq*m_sinfac/(1.+alphasq);
			}
	};

	class StarkCollTransProb_VF01
	{
		long int n, l, lp;
		double cosU1, cosU2, sinU1sinU2;		
	public:
		StarkCollTransProb_VF01( long int n1, long int l1, long int lp1) : n(n1), l(l1), lp(lp1)
			{
				/* These are defined on page 11 of VF01 */ 
				cosU1 = 2.*POW2((double)l/(double)n) - 1.;
				cosU2 = 2.*POW2((double)lp/(double)n) - 1.;
				sinU1sinU2 = sqrt( ( 1. - cosU1*cosU1 )*( 1. - cosU2*cosU2 ) );
			}
		double operator()(double alpha, double alphamin) const;
		double bfun(double alpha, double alphamin) const
			{
				double sinChiOver2sq = Chi_VF01(alpha, alphamin).sinChiOver2sq();
				double cosChi = 1. - 2*sinChiOver2sq;
				return (cosU1*cosU2 + sinU1sinU2 - cosChi);
			}
	};

	double StarkCollTransProb_VF01::operator() (double alpha, double alphamin) const
	{
		DEBUG_ENTRY( "StarkCollTransProb_VF01::operator()()" );
		double probability;
		ASSERT( alpha >= 1e-30 );

		double sinChiOver2sq = Chi_VF01(alpha, alphamin).sinChiOver2sq();

		/*if ( lp == 0 )
			 // Refer to Vrinceanu and Flanery (2001)PRA 63, 032701
			 // before eq. (35) "for l->lp=0 transitions, the probability is zero."
			probability = 0.;
		else */
		if( l == 0 || lp == 0 )
		{
			long lf = max(l,lp);
			if( n*n*sinChiOver2sq  <= lf*lf ) //lf*lf ) //lp*lp )
			{
				probability = 0.;
			}
			else
			{
				/* Here is the initial state S case */
				/*ASSERT( sinChiOver2sq > 0. );
				ASSERT( n*n*sinChiOver2sq > lp*lp ); //lf*lf ); //lp*lp );*/
				/* This is equation 35 of VF01.  There is a factor of hbar missing in the denominator of the
				 * paper, but it's okay if you use atomic units (hbar = 1). */
				//probability = lp+0.5 / sqrt( n*n*sinChiOver2sq * (n*n*sinChiOver2sq - lp*lp) );
				probability = (lf+0.5) / sqrt( n*n*sinChiOver2sq * (n*n*sinChiOver2sq - lf*lf) );
			}

			if (lp == 0)
			{
				probability /= (2.*l+1);
			}
		}
		else
		{
#if 0
			// Rearrangement of cosChi above means that this is no longer necessary.
			if( OneMinusCosChi == 0. )
			{
				double hold = sin( deltaPhi / 2. );
				/* From approximation at bottom of page 10 of VF01. */
				OneMinusCosChi = 8. * alpha * alpha * POW2( hold );
			}
#endif	
		
			if( sinChiOver2sq == 0. )
			{
				probability = 0.;
			}
			else
			{
				double cosChi = 1. - 2*sinChiOver2sq;
				/* Here is the general case, with factor of OneMinusCosChi cancelled for efficiency */
				double A = (cosU1*cosU2 - sinU1sinU2 - cosChi);
				double B = (cosU1*cosU2 + sinU1sinU2 - cosChi);
				
				/* These are the three cases of equation 34. */
				if( B <= 0. )
				{
					probability = 0.;
				}
				else
				{
					ASSERT( sinChiOver2sq > 0. );

					probability = 2.*lp/(PI* /*H_BAR* */ n*n* sinChiOver2sq );

					/*lp factor in eq (34) of VF 2001 has been changed to lp+1/2
					 * in order to cope with statistical factors
					 */
					//probability = (2.*lp+1)/(PI* /*H_BAR* */ n*n* sinChiOver2sq );
					
					if( A < 0. )
					{
						ASSERT( B > A );
						probability *= ellpk( -A/(B-A) ) * sqrt( 2*sinChiOver2sq / (B - A) );
					}
					else
					{
						probability *= ellpk( A/B ) * sqrt( 2*sinChiOver2sq / B );
					}
				}
			}
		}

		return probability;	
	}

	template<class P>
	class my_Integrand_VF01_alpha : public D_fp
	{
		P s;
		double alphamin;
	public:
		my_Integrand_VF01_alpha(long int n, long int l, long int lp, double alphamin) : 
			s( n, l, lp), alphamin(alphamin) {}
		
		double operator() (double alpha) const
			{
				ASSERT( alpha >= 1.e-30 );
				
				return s( alpha, alphamin )/(alpha*alpha*alpha);
			}
		double bfun (double alpha) const
			{
				return s.bfun(alpha,alphamin);
			}
	};

	template<class P>
	class my_Integrand_VF01_beta : public D_fp
	{
		P s;
		double alphamin;
	public:
		// Change variables to beta = log(alpha) to improve accuracy & convergence
		// rate of Romberg integration
		my_Integrand_VF01_beta(long int n, long int l, long int lp, double alphamin) : 
			s( n, l, lp), alphamin(alphamin) {}
		
		double operator() (double beta) const
			{
				double alpha = exp(beta);
				ASSERT( alpha >= 1.e-30 );
				
				return s( alpha, alphamin )/(alpha*alpha);
			}
		double bfun (double alpha) const
			{
				return s.bfun(alpha,alphamin);
			}
	};

	class StarkCollTransProb_VOS12QM
	{
		long int n, l, lp;
		double cosU1, cosU2, sinU1sinU2;
		vector<double> common;
	public:
		StarkCollTransProb_VOS12QM( long int n1, long int l1, long int lp1) : n(n1), l(l1), lp(lp1), common(n)
		{
			/* These are defined on page 11 of VF01 */ 
			cosU1 = 2.*POW2((double)l/(double)n) - 1.;
			cosU2 = 2.*POW2((double)lp/(double)n) - 1.;
				
			sinU1sinU2 = sqrt( ( 1. - cosU1*cosU1 )*( 1. - cosU2*cosU2 ) );
			long Lmax = MIN2(n-1,lp+l), Lmin = abs(lp-l);
			double ffac = 1./sqrt(double(n));
			long nsq = n*n;
			for (long L=1; L < Lmin; ++L)
			{
				ffac *= L/sqrt(double(nsq-L*L));
			}
			if (ffac <= 0.)
			{
				fprintf(ioQQQ," PROBLEM: Catastrophic underflow in VOS12 QM transition probability\n"
						  " Try reducing resolved levels\n");
				cdEXIT(EXIT_FAILURE);
			}

			const int sjtyp = 1;
			if (sjtyp == 3)
			{
				double l1min,l1max,lmatch;
				long ier;
				rec6j(get_ptr(common)+Lmin,lp,l,0.5*(n-1),0.5*(n-1),0.5*(n-1),
						&l1min,&l1max,&lmatch,n-Lmin,&ier);
				ASSERT(Lmin == int(l1min) && ier >= 0);
			}
			for (long L=Lmin; L <= Lmax; ++L)
			{
				if (L > 0)
					ffac *= L/sqrt(double(nsq-L*L));
				if (0)
				{
					double factL = factorial(L);
					fprintf(ioQQQ,"TEST %ld %g %g\n",
							  L,ffac,factL*factL*factorial(n-L-1)/factorial(n+L));
				}
				double sixj1;
				// pass 2*j arguments to allow for half-integer values
				if (sjtyp == 1)
					sixj1 = sjs(2*lp,2*l,2*L,n-1,n-1,n-1);
				else if (sjtyp == 2)
					sixj1 = SixJFull(2*lp,2*l,2*L,n-1,n-1,n-1);
				else
					sixj1 = common[L];
				common[L] = sixj1*ffac*sqrt(2.*double(L)+1.);
			}
		}
		double operator()(double alpha, double alphamin) const;
		double bfun (double, double ) const
		{
			return 0.;
		}
	};

	double StarkCollTransProb_VOS12QM::operator() (double alpha, double alphamin) const
	{
		DEBUG_ENTRY( "StarkCollTransProb_VOS12QM::operator()()" );
		ASSERT( alpha >= 1e-30 );
		
		Chi_VF01 chi(alpha, alphamin);
		double sinChiOver2sq = chi.sinChiOver2sq();
		
		// Implement equation (2) of VOS12 -- QM probability
		// Definition of chi appears inconsistent between VF01 &
		// VOS12.  Compare the last display expression in the LH
		// column of VF01 032701-10 with VOS12 eq (3)
		//ASSERT (sinChiOver2sq <= 1.);

		double cosChiVOS12QM = chi.cosChi();
		ASSERT( cosChiVOS12QM <= 1);
		//double sinChiOver2sq = 1 - pow2(cosChiVOS12QM);
		long Lmax = MIN2(n-1,lp+l);
		UltraGen u(n,cosChiVOS12QM);
		for (long L=n-1; L > Lmax; --L)
		{
			u.step();
		}
		double frac = 4*sinChiOver2sq;
		double pqm = 0.;
		for (long L=Lmax; L >= abs(lp-l); --L)
		{
			double gegen = u.val();
			u.step();
			if (0)
			{
				// Test relative accuracy of different
				// Gegenbauer/ultraspherical implementations
				double gegen1 = gegenbauer(n-L-1,L+1,cosChiVOS12QM);
				if (! fp_equal_tol(gegen1,gegen,
										 1000*DBL_EPSILON+1e-10*fabs(gegen1)))
					fprintf(ioQQQ,"DEBUG %ld %ld %g: %g %g %g\n",n,L,cosChiVOS12QM,
							  gegen1,gegen,gegen/gegen1-1.);
			}
			pqm *= frac;
			double fac = common[L]*gegen;
			pqm += fac*fac;
		}
		pqm *= powi(frac,abs(lp-l))*(2*lp+1);
		ASSERT (pqm >= 0);

		return pqm;	
	}

	void compare_VOS12()
	{
		//This function allows to compare probabilities of VOS12 approaches
		DEBUG_ENTRY("compare_VOS12()");
		//VOS12, Fig 1.
		{
			long n=30, l = 4, lp = 3, npt = 1000;
			double alphamin=0.;
			StarkCollTransProb_VOS12QM qm( n, l, lp);
			StarkCollTransProb_VF01 cl( n, l, lp);
			double v = 0.25;
			for (long i=0; i<npt; ++i)
			{
				double ban = 900*(i+0.5)/double(npt);
				double alpha = 1.5/(ban*v);
				fprintf(ioQQQ,"%g %g %g\n",ban,
						  cl(alpha,alphamin),qm(alpha,alphamin));
			}
		}
		fprintf(ioQQQ,"\n\n");
		//VF01, Fig 16. -- one panel
		{
			long n=28, l = 5, lp = 13, npt = 1000;
			double alphamin=0.;
			StarkCollTransProb_VOS12QM qm( n, l, lp);
			StarkCollTransProb_VF01 cl( n, l, lp);
			for (long i=0; i<npt; ++i)
			{
				double alpha = (i+0.5)/double(npt);
				fprintf(ioQQQ,"%g %g %g\n",alpha,
						  cl(alpha,alphamin),qm(alpha,alphamin));
			}
		}
		fprintf(ioQQQ,"\n\n");

		{
			long n=8, l = 5, lp = 6, npt = 1000;
			double alphamin=0.;
			vector<StarkCollTransProb_VOS12QM> qm;
			for (long i=0; i<4; ++i)
			{
				long scale = 1<<i;
				qm.push_back(StarkCollTransProb_VOS12QM( 
									 scale*n, scale*l, scale*lp));
			}
			StarkCollTransProb_VF01 cl( n, l, lp);
			for (long i=0; i<npt; ++i)
			{
				double alpha = (i+0.5)/double(npt);
				fprintf(ioQQQ,"%g %g",alpha,cl(alpha,alphamin));
				long scale = 1;
				for (vector<StarkCollTransProb_VOS12QM>::iterator it=qm.begin();
					  it != qm.end(); ++it)
				{
					fprintf(ioQQQ," %g",scale*(*it)(alpha,alphamin));
					scale <<= 1;
				}
				fprintf(ioQQQ,"\n");
			}
		}
		fprintf(ioQQQ,"\n\n");
		//Compare all rates out of level
		{
			long n=30, l = 29, npt = 10000;
			double alphamin=0.;
			vector<StarkCollTransProb_VOS12QM> qm;
			for (long lp=0; lp<n; ++lp)
				qm.push_back(StarkCollTransProb_VOS12QM( n, l, lp));
			double v = 0.25;
			for (long i=0; i<npt; ++i)
			{
				double ban = 400*(i+0.5)/double(npt);
				double alpha = 1.5/(ban*v);
				fprintf(ioQQQ,"PC %g",ban);
				for (long lp=0; lp<n; ++lp)
					fprintf(ioQQQ," %g",
							  qm[lp](alpha,alphamin));
				fprintf(ioQQQ,"\n");
			}
		}
	}
}

/**CS_l_mixing_VOS12 Collision treatment based on Vrinceanu Onofrio & Sadeghpour 2012
 http://adsabs.harvard.edu/abs/2012ApJ...747...56V
 */
double CS_l_mixing_VOS12(long n, long l, long lp,
								 long nelem, double gLo, long Ztarget, long Collider, double sqrte)
{
	// Equation (9) of Vrinceau etal ApJ 747, 56 (2012)
	long lmin = (lp < l) ? lp : l, dl = abs(l-lp);
	double nlfac = double(n*n*(n*n*(l+lp)-lmin*lmin*(l+lp+2*dl)))/
		((l+0.5)*dl*dl*dl);
	long Z = ColliderCharge[Collider];
	double massratio = reduced_amu(nelem,Collider)/ELECTRON_MASS;
	double rate = 1.294e-5*sqrt(massratio)*Z*Z/(Ztarget*Ztarget*sqrte) * nlfac;
	double cs = rate / ( COLL_CONST * powpq(massratio, -3, 2) ) * sqrte * gLo;
	return cs;
}

template<class P>
STATIC double CSIntegral_QG32(const my_Integrand_VF01_alpha<P>& func, 
												 double alphamin_int, double alphamax_int)
{
	if (alphamin_int >= alphamax_int)
		return 0.;
	double CSIntegral = 0.;
	Integrator<my_Integrand_VF01_alpha<P>, Gaussian32> VF01_alpha(func);
	double step = (alphamax_int - alphamin_int)/5.;
	double alpha1 = alphamin_int;
	CSIntegral += VF01_alpha.sum( alpha1, alpha1+step );
	CSIntegral += VF01_alpha.sum( alpha1+step, alpha1+4.*step );
	return CSIntegral;
}

template<class P>
STATIC double CSIntegral_Romberg(long ipISO, const my_Integrand_VF01_beta<P>& func,
													 double alphamin_int, double alphamax_int, double eps)
{
	if (alphamin_int >= alphamax_int)
		return 0.;

	// Bisection search for lower limit to integral, should be
	// improved to e.g. secant method
	double alphalo = alphamin_int;
	if (! iso_ctrl.lgCS_VOS12QM[ipISO] && func.bfun(alphalo) <= 0.)
	{
		double alphahi = alphamax_int;
		while (alphahi-alphalo > 0.5*eps*(alphahi+alphalo))
		{
			double alphamid = 0.5*(alphahi+alphalo);
			if (func.bfun(alphamid) <= 0.)
				alphalo = alphamid;
			else
				alphahi = alphamid;
		}
		// Move edge *just* inside integration range
		alphalo = alphahi;
	}
	class integrate::Romberg<integrate::Midpoint<my_Integrand_VF01_beta<P> > >
		VF01_alphaR(integrate::Midpoint<my_Integrand_VF01_beta<P> >(func,log(alphalo),log(alphamax_int)));
	VF01_alphaR.update(eps);
	return VF01_alphaR.sum();
}

template<class P>
STATIC double CSIntegral_Romberg_alpha(long ipISO, const my_Integrand_VF01_alpha<P>& func,
													 double alphamin_int, double alphamax_int, double eps)
{
	if (alphamin_int >= alphamax_int)
		return 0.;
	// Bisection search for lower limit to integral, should be
	// improved to e.g. secant method
	double alphalo = alphamin_int;
	if (! iso_ctrl.lgCS_VOS12QM[ipISO] && func.bfun(alphalo) <= 0.)
	{
		double alphahi = alphamax_int;
		while (alphahi-alphalo > 0.5*eps*(alphahi+alphalo))
		{
			double alphamid = 0.5*(alphahi+alphalo);
			if (func.bfun(alphamid) <= 0.)
				alphalo = alphamid;
			else
				alphahi = alphamid;
		}
		// Move edge *just* inside integration range
		alphalo = alphahi;
	}
	class integrate::Romberg<integrate::Midpoint<my_Integrand_VF01_alpha<P> > >
		VF01_alphaR(integrate::Midpoint<my_Integrand_VF01_alpha<P> >(func,alphalo,alphamax_int));
	VF01_alphaR.update(eps);
	return VF01_alphaR.sum();
}

template<class P>
STATIC double collision_strength_VF01( long ipISO, double E_Proj_Ryd,
													const my_Integrand_VF01_E<P>& vf)
{

	DEBUG_ENTRY( "collision_strength_VF01()" );

	double reduced_b_max1;
	/* The proton mass is needed here because the ColliderMass array is
	 * expressed in units of the proton mass, but here we need absolute mass. */
	double reduced_vel = sqrt( 2.*E_Proj_Ryd*EN1RYD/vf.reduced_mass )/vf.RMSv;

	/* put limits on the reduced velocity.   These limits should be more than fair. */
	ASSERT( reduced_vel > 1.e-10 );
	ASSERT( reduced_vel < 1.e10 );

 	/* Reduced here means in units of vf.aveRadius: */
	//double reduced_b_min = 1.5 * vf.ColliderCharge / reduced_vel;
	/* alphamax has been increased in order to account for low cut-offs
	 * needed at low temperatures T<100 and high densities
	 */
	//if (iso_ctrl.lgCS_VOS12QM[ipISO])
	double alphamax = 15000./pow2(vf.n);
	/* alphamax can be fixed in different ways. Next is done assuming that cross section is zero at scaled
	 * temperatures T=1K are zero so the lower limit= higher limit*/
	/*double alphamax = 1.5*vf.ColliderCharge/(reduced_vel * vf.reduced_b_max)*
			sqrt(E_Proj_Ryd*phycon.te/BOLTZMANN);*/
	// before was 1.5*ColliderCharge/(reduced_vel * reduced_b_min); but above
	// definition of reduced_b_min meant this is always 1. to machine
	// accuracy
	// this alphamax was not enough to cover low temperature cross sections
	if (iso_ctrl.lgCS_Vrinceanu[ipISO])
		alphamax = 1.;
	//


	/*reduced_b_min is used only to be compared to reduced_b_max*/
	double reduced_b_min = 1.5 * vf.ColliderCharge /reduced_vel/alphamax;
	ASSERT( reduced_b_min > 1.e-10 );
	ASSERT( reduced_b_min < 1.e10 );

	/* Here cut-off, set for temperature, is scaled to real energy for Maxwell calculation */
	// Debye Radius, reduced in units of aveRadius
	double reduced_debye = sqrt( BOLTZMANN*vf.temp/vf.ColliderCharge/dense.eden/PI )
											/ (2.*ELEM_CHARGE_ESU)/vf.aveRadius;
	if (vf.reduced_b_max < reduced_debye && iso_ctrl.lgCS_VOS12QM[ipISO])
		reduced_b_max1 = vf.reduced_b_max*sqrt(E_Proj_Ryd*EN1RYD/BOLTZMANN/vf.temp);
	else
		reduced_b_max1 = vf.reduced_b_max;
	
	// Return strictly zero when there is no allowed range: falling
	// through with MAX2 generates junk values due to rounding error
	if (reduced_b_max1 <= reduced_b_min)
		return 0.;

	reduced_b_max1 = MAX2( reduced_b_max1, reduced_b_min );
	ASSERT( reduced_b_max1 > 0. );

	double alphamin = 1.5*vf.ColliderCharge/(reduced_vel * reduced_b_max1);
	ASSERT( alphamin > 0. );
	ASSERT( alphamax > 0. );

	// set up the integrator.

	double alphamin_int = MAX2( alphamin, 1.e-30 );
	double alphamax_int = MAX2( alphamax, 1.e-20 );

	double CSIntegral = 0.;

	if (0)
	{
		my_Integrand_VF01_alpha<P> func(vf.n, vf.l, vf.lp, alphamin);
		CSIntegral = CSIntegral_QG32(func,alphamin_int,alphamax_int);
	}
	else
	{
		my_Integrand_VF01_beta<P> funcb(vf.n, vf.l, vf.lp, alphamin);
		if (1)
		{
			CSIntegral = CSIntegral_Romberg(ipISO,funcb,alphamin_int,alphamax_int,1e-4);
		}
		else
		{
			double epsabs=0., epsrel=1e-5, abserr, qqq;
			const long limit=25, lenw=4*limit;
			long neval, ier,last, iwork[limit];
			double work[lenw];
			double lalphamin_int = log(alphamin_int),
				lalphamax_int = log(alphamax_int);
			dqags_(funcb,&lalphamin_int,&lalphamax_int,&epsabs,&epsrel,&qqq,
					 &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
			CSIntegral = qqq;
		}

		if (0)
		{
			double qqq;
			my_Integrand_VF01_alpha<P> func(vf.n, vf.l, vf.lp, alphamin);
			const char *ccc;
			if (0)
			{
				ccc = "G32";
			   qqq = CSIntegral_QG32(func,alphamin_int,alphamax_int);
			}
			else if (0)
			{
				ccc = "Unirom";
				qqq = CSIntegral_Romberg_alpha(ipISO,func,alphamin_int,alphamax_int,1e-5);
			}
			else
			{
				ccc = "QAGS";
				double epsabs=0., epsrel=1e-5, abserr;
				const long limit=25, lenw=4*limit;
				long neval, ier,last, iwork[limit];
				double work[lenw];
				double lalphamin_int = log(alphamin_int),
					lalphamax_int = log(alphamax_int);
				dqags_(funcb,&lalphamin_int,&lalphamax_int,&epsabs,&epsrel,&qqq,
						&abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);
				fprintf(ioQQQ,"Check VF QAGS1 err=%g neval=%ld ier=%ld\n",abserr,neval,ier);
			}
			if (CSIntegral > 0. || qqq > 0.)
			{
				fprintf(ioQQQ,"Check VF[%ld %ld->%ld %g %g]: %s %g var %g err=%g\n",
						  vf.n,vf.l,vf.lp,E_Proj_Ryd,alphamin,
						  ccc,qqq,CSIntegral,(qqq-CSIntegral)/SDIV(CSIntegral));
			}
		}
	}


	/* Factors outside integral	*/
	double ConstantFactors = 4.5*PI*POW2(vf.ColliderCharge*vf.aveRadius/reduced_vel);

	/* Calculate cross section */
	double cross_section = ConstantFactors * CSIntegral;

	/* convert to collision strength, cross section is already in cm^2 */
	double coll_str = ConvCrossSect2CollStr( cross_section, vf.gLo, E_Proj_Ryd, vf.reduced_mass );

	coll_str = MAX2( SMALLFLOAT, coll_str);

	return coll_str;
}

/* resolved to resolved and collapsed to resolved modification 
 *
 * convert resolved to collapsed:
 * rates: q(n->n') = \sum_l' \sum_l (2l+1)(2s+1) q(nl->nl')/ 2n^2
 * upsilons: Y(nn') = \sum_l' \sum_l Y(nln'l')
 * example: q(n=3->n=2) = 1/18(6q(3p->2s) + 2q(3s->2p) + 10q(3d->2p) + ...
 *
 * convert collapsed to collapsed to resolved:
 * rates: q(nl-n'l') = (2s'+1)(2l'+1)q(n->n')*n2/[\sum\sum(2s+1)(2l+1)(2s'+1)(2l'+1)]
 * with the constrains [{|s-s'|=0}]
 * upsilons: Y(nln'l') =(2s+1)(2l+1)(2s'+1)(2l'+1)Y(nn')/S
 * with the constrains {|s-s'|=0}
 *
 *where: S=[\sum\sum(2s+1)(2l+1)(2s'+1)(2l'+1)] = (2s+1)^2 n^2n'^2
 *
 * If i want a collapsed to resolved coefficient:
 *
 * 1) If the routine returns a collapsed value:
 *
 * rates: q(n -> n'l') = (2s+1)^2\sum(2l+1)(2l'+1)q(n->n')*n^2/S
 * upsilon: Y(nn'l') = (2s+1)^2\sum(2l+1)(2l'+1)Y(nn')/S
 *
 *
 * 2) If it is a resolved value
 *
 * rates: q(n->n'l') = \sum_l (2l+1)(2s+1)q(nl->n'l')/2n^2
 * upsilons: Y(nn'l') = \sum_l Y(nln'l')
 *
 */
double CSresolver(long ipISO, long nHi,long lHi,long sHi,long nLo,
		long lLo, long sLo, long n_HighestResolved)
{
	DEBUG_ENTRY( "CSresolver()" );

	if (sHi != sLo && nHi <= n_HighestResolved )
		return 0.;

	double factor=1.;

	/* S= \sum\sum (2s+1)(2l+1)(2s'+1)(2l'+1)
	 * with |s-s'|=0
	 */


	double S = pow2(nHi*nLo);

	if (ipISO==ipH_LIKE)
	{
		ASSERT( sLo == ipDOUBLET );
		S *= pow2( ipDOUBLET );
	}
	else if(ipISO==ipHE_LIKE)
	{
		if( nLo == 1 )
		{
			ASSERT( sLo == ipSINGLET );
			S *= pow2( ipSINGLET );
		}
		else
			S *= pow2( ipSINGLET ) + pow2( ipTRIPLET );
	}


	/* resolved to resolved */
	if ( nHi <= n_HighestResolved )
		factor *= sHi*sLo*(2.*lHi+1)*(2.*lLo+1)/S;
				//sHi*(2.*lHi+1)/16/pow4(nHi);
	/* collapsed to resolved */
	else if (nLo <= n_HighestResolved )
	{
		double factor1 = pow2(nHi)*(2.*lLo+1);
		factor *= sLo*sLo*factor1/S;
	}
	else
		TotalInsanity();

	return factor;
}
