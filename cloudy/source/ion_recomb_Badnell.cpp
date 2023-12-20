/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_recom_calculate calculate radiative and dielectronic recombination rate coefficients */
/*Badnell_rec_init This code is written by Terry Yun, 2005 *
 * It reads rate coefficient fits into 3D arrays and output array.out for testing *
 * The testing can be commented out */
/*Badnell_DR_rate_eval This code is written by Terry Yun, 2005 *
 * It interpolates the rate coefficients in a given temperature.*
   It receives ATOMIC_NUM_BIG, NELECTRONS values, temperature and returns the rate coefficient*
   It returns
        '-2': initial <= final
              init < 0 or init >302 or final < 0 or final > 302
        '-1': the transition is not defined
        '99': unknown invalid entries */ 
#include "cddefines.h"
#include "phycon.h"
#include "elementnames.h"
#include "atmdat.h"
#include "atmdat_adfa.h"
#include "iso.h"
#include "ionbal.h"
#include "save.h"
#include "freebound.h"
#include "dense.h"
#include "ran.h"
#include "thirdparty.h"

static const int MAX_FIT_PAR_DR = 9;
static multi_arr<double,3> DRFitParPart1; 
static multi_arr<double,3> DRFitParPart2;
static multi_arr<int,2> nDRFitPar;

static const int MAX_FIT_PAR_RR = 6;
static multi_arr<double,3> RRFitPar; 

/* flags to recall that we have read the fits from the main data files */
static multi_arr<bool,2> lgDRBadnellDefined,
	lgDRBadnellDefinedPart2,
	lgRRBadnellDefined,
	lgDR_BadWeb_exist;
static bool lgMustAllocateRec=true;
static double RecNoise[LIMELM],
	DR_Badnell_rate_coef_mean_ion[LIMELM];

static char chDRDataSource[LIMELM][LIMELM][10];
static char chRRDataSource[LIMELM][LIMELM][10];

/* these enable certain debugging print statements */
/* #define PRINT_DR */
/* #define PRINT_RR */

#if defined(PRINT_DR) || defined(PRINT_RR)
static const char FILE_NAME_OUT[] = "array.out";
#endif

/* This function computes the standard electron density-dependent
 * suppression factor of the collisional DR rate coefficient of H-like
 * to Cl-like ions, based on Hugh P. Summers' 1979 report AL-R-5 report.
 * It is then scalable for other choices of ionic charge and temperature.
 */
STATIC double CollisSuppres(
	/* atomic_number on physics scale = nuclear charge - 6 for C */
	long int atomic_number,
	/* ionic_charge = charge before recombination, physics scale, 3 for C+3 -> C+2 */
	long int ionic_charge,
	/* eden = electron density */
	double eden,
	/*T = temperature (K)*/
	double T )
{
	DEBUG_ENTRY( "CollisSuppres()" );

	// >>refer	DR	Nikolic D., Gorczyca T.W., Korista K.T., Ferland G.J., Badnell N.R., 2013, ApJ 768, 82
	// >>refer	DR	Nikolic D., Gorczyca T.W., Korista K.T., Ferland G.J., van Hoof, P.A.M., Williams, R.J.R., Badnell N.R., 2015, ApJ in preparation

	ASSERT( ionic_charge > 0 );

	/* fitting constants to compute nominal suppression factor function */
	const double mu = 0.000;	/* pseudo Voigt Lorentzian mixture */
	const double w = 5.64586;	/* suppression decay rate */
	const double x_a0 = 10.1821; 	/* log10 of the electron density fitting parameter for H-like ions */

	/* a fitting constant to compute the suppression factor corrected for an
	 * estimate of surviving DR based on the lowest dipole allowed core
	 * excitation energy
	 */
	const double c = 10.0;	/* smaller c means larger fraction will survive, and vice versa */

	double s, snew, x_a, E_c;
	double T_0, q_0, A_N;   /* seed temperature, charge, and sequence selector*/
	long int iso_sequence, N_1, N_2;

	eden = log10(eden);
	double Tlog = log10(T);

	/* the isoelectronic sequence number, iso_sequence=3 for Li-like, etc */
	iso_sequence = atomic_number - ionic_charge;
	ASSERT( iso_sequence > 0 );

	/* Temporarily save ionic_charge/10 needed later for excitation energy fits*/
	double c10 = double(ionic_charge) / 10.;

	N_1 = -1;
	N_2 = -1;
	/* initiate sequence-wise charge-dependent seed charge */
	if( (iso_sequence >= 1) && (iso_sequence <= 2) )       /* 1st row sequences */
	{
		N_1 = 1;
		N_2 = 2;
	}
	else if( (iso_sequence >= 3) && (iso_sequence <= 10) ) /* 2nd row sequences */
	{
		N_1 = 3;
		N_2 = 10;
	}
	else if( (iso_sequence >= 11) && (iso_sequence <= 18) ) /* 3rd row sequences */
	{
		N_1 = 11;
		N_2 = 18;
	}
	else if( (iso_sequence >= 19) && (iso_sequence <= 36) ) /* 4th row sequences */
	{
		N_1 = 19;
		N_2 = 36;
	}
	else if( (iso_sequence >= 37) && (iso_sequence <= 54) ) /* 5th row sequences */
	{
		N_1 = 37;
		N_2 = 54;
	}
	else if( (iso_sequence >= 55) && (iso_sequence <= 86) ) /* 6th row sequences */
	{
		N_1 = 55;
		N_2 = 86;
	}
	else if( iso_sequence >= 87 )                           /* 7th row sequences */
	{
		N_1 = 87;
		N_2 = 118;
	}
	//fprintf(ioQQQ, "DEBUGGG %li %li %.2e %.2e %li %li %li\n",
	//		atomic_number, ionic_charge, eden, T ,
	//		N_1 , N_2 , iso_sequence);
	ASSERT( N_1>0 && N_2>0 );

	/* initiate zig-zag approximation for A(N), which must be enveloped by two asymptotes:
	 * Amax = 12 + 10*N   and  Amin = 12 + 2*N
	 */
	A_N = 12.0 + 10.0 * N_1 + (10.0 * N_1 - 2.0 * N_2) * (iso_sequence - N_1) / (N_1 - N_2);
	ASSERT( A_N >= 16.0 );

	/* Now loop through specific sequences to assign computational estimates
	 * of lowest dipole allowed core excitation energies (these are fits to
	 * NIST statistical weighted energies) and readjust A(N) values for N<=5 sequences.
	 */

	if( iso_sequence == 1 )      /* H-like ions */
	{
		E_c = 0.0;
		A_N = 16.0;
	}
	else if( iso_sequence == 2 ) /* He-like ions */
	{
		E_c = 0.0;
		A_N = 18.0;
	}
	else if( iso_sequence == 3 ) /* Li-like ions */
	{
		E_c = 1.96274 + c10*(20.30014 + c10*(-0.97103 + c10*( 0.85453 + c10*( 0.13547 + 0.02401*c10))));
		A_N = 66.0;
	}
	else if( iso_sequence == 4 ) /* Be-like ions */
	{
		E_c = 5.78908 + c10*(34.08270 + c10*( 1.51729 + c10*(-1.21227 + c10*( 0.77559 - 0.00410*c10))));
		A_N = 66.0;
	}
	else if( iso_sequence == 5 ) /* B-like ions */
	{
		E_c = 0.0;
		A_N = 52.0;
	}
	else if( iso_sequence == 7 ) /* N-like ions */
	{
		E_c = 11.37092 + c10*(36.22053 + c10*( 7.08448 + c10*(-5.16840 + c10*( 2.45056 - 0.16961*c10))));
	}
	else if( iso_sequence == 11 ) /* Na-like ions */
	{
		E_c =  2.24809 + c10*(22.27768 + c10*(-1.12285 + c10*( 0.90267 + c10*(-0.03860 + 0.01468*c10))));
	}
	else if( iso_sequence == 12 ) /* Mg-like ions */         
	{
		E_c =  2.74508 + c10*(19.18623 + c10*(-0.54317 + c10*( 0.78685 + c10*(-0.04249 + 0.01357*c10))));
	}
	else if( iso_sequence == 15 ) /* P-like ions */
	{
		E_c =  1.42762 + c10*( 3.90778 + c10*( 0.73119 + c10*(-1.91404 + c10*( 1.05059 - 0.08992*c10))));
	}
	else
	{
		E_c = 0.0;	/* forces suppression factor to s for all T */
	}

	/* special case S2+ -> S+ recombination
	 * >>refer	DR	Badnell N.R., Ferland G.J., Gorczyca T.W., Nikolic D., Wagle G.A., 2015, ApJ 804, 100 */
	if( iso_sequence == 14 && ionic_charge == 2 )
		E_c = 17.6874;

	double xic = double(ionic_charge);

	/* check for low temperatures in sequences below Carbon-like*/
	if( iso_sequence <= 5 ) /* H-like to B-like ions */
	{
		double pp1,pp2,pp3,pp4,pp5,pp6,AN0; /* pi adjustment coefficients */

		double Tscal = T/pow2(xic);
		A_N *= 2.0 / ( 1.0 + exp(-sqrt(25000.0 / Tscal))); /* remove temperature cusp for Tscal < 25000 */ 
  
		if( iso_sequence == 1 ) /* H-like */
		{
			pp1 =  4.7902  + 0.32456 * pow(xic, 0.97838) * exp( -xic /  24.78084 ); // xc1
			pp2 = -0.0327  + 0.13265 * pow(xic, 0.29226);                           // w1
			pp3 = -0.66855 + 0.28711 * pow(xic, 0.29083) * exp( -xic /   6.65275 ); // A1
			pp4 =  6.23776 + 0.11389 * pow(xic, 1.24036) * exp( -xic /  25.79559 ); // xc2
			pp5 =  0.33302 + 0.00654 * pow(xic, 5.67945) * exp( -xic /   0.92602 ); // w2
			pp6 = -0.75788 + 1.75669 * pow(xic,-0.63105) * exp( -xic / 184.82361 ); // A2
		}                                                                                           
		else if( iso_sequence == 2 ) /* He-like */
		{
			pp1 =  4.82857 + 0.3     * pow(xic, 1.04558) * exp( -xic / 19.6508  ); // xc1
			pp2 = -0.50889 + 0.6     * pow(xic, 0.17187) * exp( -xic / 47.19496 ); // w1
			pp3 = -1.03044 + 0.35    * pow(xic, 0.3586 ) * exp( -xic / 39.4083  ); // A1
			pp4 =  6.14046 + 0.15    * pow(xic, 1.46561) * exp( -xic / 10.17565 ); // xc2
			pp5 =  0.08316 + 0.08    * pow(xic, 1.37478) * exp( -xic /  8.54111 ); // w2
			pp6 = -0.19804 + 0.4     * pow(xic, 0.74012) * exp( -xic /  2.54024 ); // A2
		}
		else if( iso_sequence == 3 ) /* Li-like */
		{
			pp1 =  4.55441 + 0.08    * pow(xic, 1.11864);                          // xc1
			pp2 =  0.3     + 2.0     * powi(xic, -2)     * exp( -xic / 67.36368 ); // w1
			pp3 = -0.4     + 0.38    * pow(xic, 1.62248) * exp( -xic /  2.78841 ); // A1
			pp4 =  4.00192 + 0.58    * pow(xic, 0.93519) * exp( -xic / 21.28094 ); // xc2
			pp5 =  0.00198 + 0.32    * pow(xic, 0.84436) * exp( -xic /  9.73494 ); // w2
			pp6 =  0.55031 - 0.32251 * pow(xic, 0.75493) * exp( -xic / 19.89169 ); // A2
		}
		else if( iso_sequence == 4 ) /* Be-like */
		{
			pp1 =  2.79861 + 1.0     * pow(xic, 0.82983) * exp( -xic / 18.05422 ); // xc1
			pp2 = -0.01897 + 0.05    * pow(xic, 1.34569) * exp( -xic / 10.82096 ); // w1
			pp3 = -0.56934 + 0.68    * pow(xic, 0.78839) * exp( -xic /  2.77582 ); // A1
			pp4 =  4.07101 + 1.0     * pow(xic, 0.7175 ) * exp( -xic / 25.89966 ); // xc2
			pp5 =  0.44352 + 0.05    * pow(xic, 3.54877) * exp( -xic /  0.94416 ); // w2
			pp6 = -0.57838 + 0.68    * pow(xic, 0.08484) * exp( -xic /  6.70076 ); // A2
		}
		else if( iso_sequence == 5 ) /* B-like */
		{
			pp1 =  6.75706 - 3.77435 *                     exp( -xic /  4.59785 ); // xc1
			pp2 =  0.0     + 0.08    * pow(xic, 1.34923) * exp( -xic /  7.36394 ); // w1
			pp3 = -0.63    + 0.06    * pow(xic, 2.65736) * exp( -xic /  2.11946 ); // A1
			pp4 =  7.74115 - 4.82142 *                     exp( -xic /  4.04344 ); // xc2
			pp5 =  0.26595 + 0.09    * pow(xic, 1.29301) * exp( -xic /  6.81342 ); // w2
			pp6 = -0.39209 + 0.07    * pow(xic, 2.27233) * exp( -xic /  1.9958  ); // A2
		}
		else
		{
			TotalInsanity();
		}
		/* define additional charge-temperature dependent portion of AN */
		AN0 = pp3 * exp( -0.5*pow2((Tlog - pp1)/pp2) ) + pp6 * exp( -0.5*pow2((Tlog - pp4)/pp5) );
		A_N *= (1. + AN0);

		if (A_N <= 0)
		{
			//fprintf(ioQQQ, "DEBUGDN an=%li q=%li logNe=%.2e Te=%.2e N_1=%li N_2=%li Niso=%li AN0=%.2e\n",
			//	atomic_number, ionic_charge, eden, T ,
			//	N_1 , N_2 , iso_sequence, AN0);
		}
		ASSERT( A_N > 0.0 );
	}

	double gg1,gg2,gg3,gg4,gg5,gg6,BN0=0.; /* gamma adjustment coefficients */

	if( iso_sequence == 5 ) /* B-like */
	{
		A_N = 52.0; /* IMPORTANT : rollback the changes introduced by the Psi - "detailed" factor above */
		gg1 =  6.91078 - 1.6385  * pow(xic, 2.18197) * exp( -xic / 1.45091 ); // xc1
		gg2 =  0.4959  - 0.08348 * pow(xic, 1.24745) * exp( -xic / 8.55397 ); // w1
		gg3 = -0.27525 + 0.132   * pow(xic, 1.15443) * exp( -xic / 3.79949 ); // A1
		gg4 =  7.45975 - 2.6722  * pow(xic, 1.7423 ) * exp( -xic / 1.19649 ); // xc2
		gg5 =  0.51285 - 0.60987 * pow(xic, 5.15431) * exp( -xic / 0.49095 ); // w2
		gg6 = -0.24818 + 0.125   * pow(xic, 0.59971) * exp( -xic / 8.34052 ); // A2        
	}
	else if( iso_sequence == 6 ) /* C-like */
	{
		gg1 =  5.90184 - 1.2997  * pow(xic, 1.32018) * exp( -xic /  2.10442 ); // xc1
		gg2 =  0.12606 + 0.009   * pow(xic, 8.33887) * exp( -xic /  0.44742 ); // w1
		gg3 = -0.28222 + 0.018   * pow(xic, 2.50307) * exp( -xic /  3.83303 ); // A1
		gg4 =  6.96615 - 0.41775 * pow(xic, 2.75045) * exp( -xic /  1.32394 ); // xc2
		gg5 =  0.55843 + 0.45    *                     exp( -xic /  2.06664 ); // w2
		gg6 = -0.17208 - 0.17353 *                     exp( -xic /  2.57406 ); // A2
	}
	else if( iso_sequence == 13 ) /* Al-like */
	{
		gg1 =  6.59628 - 3.03115 *                     exp( -xic / 10.519821 ); // xc1        
		gg2 =  1.20824 - 0.85509 * pow(xic, 0.21258) * exp( -xic / 25.56     ); // w1        
		gg3 = -0.34292 - 0.06013 * pow(xic, 4.09344) * exp( -xic /  0.90604  ); // A1        
		gg4 =  7.92025 - 3.38912 *                     exp( -xic / 10.02741  ); // xc2        
		gg5 =  0.06976 + 0.6453  * pow(xic, 0.24827) * exp( -xic / 20.94907  ); // w2        
		gg6 = -0.34108 - 0.17353 *                     exp( -xic /  6.0384   ); // A2
	}
	else if( iso_sequence == 14 ) /* Si-like */
	{
		gg1 =  5.54172 - 1.54639 * pow(xic, 0.01056) * exp( -xic /  3.24604  ); // xc1
		gg2 =  0.39649 + 0.8     * pow(xic, 3.19571) * exp( -xic /  0.642068 ); // w1 
		gg3 = -0.35475 - 0.08912 * pow(xic, 3.55401) * exp( -xic /  0.73491  ); // A1 
		gg4 =  6.88765 - 1.93088 * pow(xic, 0.23469) * exp( -xic /  3.23495  ); // xc2
		gg5 =  0.58577 - 0.31007 * pow(xic, 3.30137) * exp( -xic /  0.83096  ); // w2 
		gg6 = -0.14762 - 0.16941 *                     exp( -xic / 18.53007  ); // A2        
	}
	else /* leave it is as it is, set gamma values so that BN0 = 0 */
	{
		gg1 =  0.; // xc1
		gg2 =  1.; // w1 
		gg3 =  0.; // A1 
		gg4 =  0.; // xc2
		gg5 =  1.; // w2 
		gg6 =  0.; // A2  
	}
	if( gg3 != 0. || gg6 != 0. )
	{
		/* define additional charge-temperature dependent portion of AN relevant for secondary autoionization */
		BN0 = gg3 * exp( -0.5*pow2((Tlog - gg1)/gg2) ) + gg6 * exp( -0.5*pow2((Tlog - gg4)/gg5) );
		A_N *= (1. + BN0);
	}
	
	if (A_N <= 0)
	{
		fprintf(ioQQQ, "DEBUGDN an=%li q=%li logNe=%.2e Te=%.2e N_1=%li N_2=%li Niso=%li BN0=%.2e\n",	
			atomic_number, ionic_charge, eden, T ,
			N_1 , N_2 , iso_sequence, BN0);
	}
	ASSERT( A_N > 0.0 );

	/* initiate charge- and sequence-dependent seed charge qo
	 * qo = (1-sqrt(2/3q))*A(N)/sqrt(q)
	 */
	q_0 = 1.0 / sqrt((double)ionic_charge);
	q_0 = A_N * q_0 * (1.0 - 0.816497 * q_0);
	ASSERT( q_0 > 0.0 );

	/* initiate charge-dependent seed temperature in K */
	T_0 = 50000.0 * pow2( q_0 );

	/* scale log activation density to current charge and temperature */
	x_a = x_a0 + log10( powi( ((double)ionic_charge/q_0), 7 ) * sqrt( T/T_0 ) );

	/* Now we're going to modify this standard suppression factor curve to
	 * allow for the survival of some fraction of the total DR rate at
	 * generally lower temperatures T, when appropriate.
	 */

	/* here we compute the standard suppression factor function, s( n_e, T, ionic_charge ) */
	if( eden >= x_a )
	{
		s = ( mu/( 1. + pow2((eden-x_a)/w) ) +
		    (1. - mu) * exp( -LN_TWO * pow2((eden-x_a)/w) ) );
	}
	else
	{
		s = 1.;
	}
	/* converting the standard curve to the revised one allowing for
	 * survival at lower energies, eqn 14 of Nikolic+13
	 */
	if( iso_sequence == 1 || iso_sequence == 2 || iso_sequence == 10 )
		snew = 1. + (s-1.)*exp(-(20.0*erfc( 2.*(eden - x_a0) ) * EVDEGK)/(c*T));
	else
		snew = 1. + (s-1.)*exp(-(E_c*EVDEGK)/(c*T));

	ASSERT( snew >= 0. && snew <= 1. );
	return snew;
}

/**\verbatim Badnell_DR_rate_eval This code is written by Terry Yun, 2005 
It interpolates the rate coefficients in a given temperature.
It receives atomic number on Physics scale, with H = 1, 
and the number of core electrons before recombination, and returns the rate coefficient*
It returns
'-2': initial <= final
init < 0 or init >302 or final < 0 or final > 302
'-1': the transition is not defined
'99': unknown invalid entries                         
\endverbatim
\param z_val atomic number on C scale - He is 1
\param n_val number of core electrons before capture of free electron 
*/ 
STATIC double Badnell_DR_rate_eval(
	/* atomic number on C scale - He - 1 */
	int nAtomicNumberCScale, 
	/* number of core electrons before capture of free electron */
	int n_core_e_before_recomb )
{

	double RateCoefficient, sum;
	int i;

	DEBUG_ENTRY( "Badnell_DR_rate_eval()" );
	ASSERT( nAtomicNumberCScale>=0 && nAtomicNumberCScale<LIMELM );

	if( nAtomicNumberCScale==ipIRON && n_core_e_before_recomb>=12 && 
		n_core_e_before_recomb<=18 )
	{
		/* these data are from table 1 of
		*>>refer	Fe	DR	Badnell, N., 2006, ApJ, 651, L73
		* Fe 8+ to Fe 12+, but also include Fe13+ and Fe 14+,
		* so these are 26-8=18 to 26-14=12 
		* increasing number of bound electrons, 0 is 14 elec, 1 is 15 elec 
		* Fe 3p^q, q=2-6
		* these are not in badnell large dat file as of 2011 apr 24 */
		static const double cFe_q[7][8] =
		{
			{5.636e-4, 7.390e-3, 3.635e-2, 1.693e-1, 3.315e-2, 2.288e-1, 7.316e-2, 0.},
			{1.090e-3, 7.801e-3, 1.132e-2, 4.740e-2, 1.990e-1, 3.379e-2, 1.140e-1, 1.250e-1},
			{3.266e-3, 7.637e-3, 1.005e-2, 2.527e-2, 6.389e-2, 1.564e-1, 0.,       0.},
			{1.074e-3, 6.080e-3, 1.887e-2, 2.540e-2, 7.580e-2, 2.773e-1, 0.,       0.},
			{9.073e-4, 3.777e-3, 1.027e-2, 3.321e-2, 8.529e-2, 2.778e-1, 0.,       0.},
			{5.335e-4, 1.827e-3, 4.851e-3, 2.710e-2, 8.226e-2, 3.147e-1, 0.,       0.},
			{7.421e-4, 2.526e-3, 4.605e-3, 1.489e-2, 5.891e-2, 2.318e-1, 0.,       0.}
		};

		/* Table 2 of Badnell 06 */
		static const double EFe_q[7][8] =
		{
			{3.628e3, 2.432e4, 1.226e5, 4.351e5, 1.411e6, 6.589e6, 1.030e7, 0},
			{1.246e3, 1.063e4, 4.719e4, 1.952e5, 5.637e5, 2.248e6, 7.202e6, 3.999e9},
			{1.242e3, 1.001e4, 4.466e4, 1.497e5, 3.919e5, 6.853e5, 0.     , 0.},
			{1.387e3, 1.048e4, 3.955e4, 1.461e5, 4.010e5, 7.208e5, 0.     , 0.},
			{1.525e3, 1.071e4, 4.033e4, 1.564e5, 4.196e5, 7.580e5, 0.     , 0.},
			{2.032e3, 1.018e4, 4.638e4, 1.698e5, 4.499e5, 7.880e5, 0.     , 0.},
			{3.468e3, 1.353e4, 3.690e4, 1.957e5, 4.630e5, 8.202e5, 0.     , 0.}
		};
		/* nion is for the above block of numbers */
		long int nion = n_core_e_before_recomb - 12;
		ASSERT( nion>=0 && nion <=6 );

		sum = 0;
		/* loop over all non-zero terms */
		for( i=0; i < 8; ++i )
		{
			sum += (cFe_q[nion][i] * sexp( EFe_q[nion][i]/phycon.te));
		}

		/*RateCoefficient = pow(phycon.te, -1.5) * sum;*/
		RateCoefficient = sum / phycon.te32;
		strcpy(chDRDataSource[nAtomicNumberCScale][nAtomicNumberCScale-n_core_e_before_recomb] ,
				"Bad06D");

		return RateCoefficient;
	}

	/*Invalid entries returns '-2':more electrons than protons */
	else if( nAtomicNumberCScale < n_core_e_before_recomb )
	{
		RateCoefficient = -2;
	}
	/*Invalid entries returns '-2' if nAtomicNumberCScale and n_core_e_before_recomb are out of the range*/
	else if( nAtomicNumberCScale >= LIMELM )
	{
		RateCoefficient = -2;
	}
	/*undefined z and n returns '-1'*/
	else if( !lgDRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{
		RateCoefficient = -1;
	}
	else if( lgDRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{
		/* this branch, recombination coefficient has been defined */
		sum = 0;
		/* loop over all non-zero terms */
		for(i=0; i<nDRFitPar[nAtomicNumberCScale][n_core_e_before_recomb]; ++i  )
		{
			sum += (DRFitParPart1[nAtomicNumberCScale][n_core_e_before_recomb][i] *
				sexp( DRFitParPart2[nAtomicNumberCScale][n_core_e_before_recomb][i]/phycon.te));
		}

		strcpy(chDRDataSource[nAtomicNumberCScale][nAtomicNumberCScale-n_core_e_before_recomb] ,
				"BadWeb");

		/*RateCoefficient = pow(phycon.te, -1.5) * sum;*/
		RateCoefficient = sum / phycon.te32;
	}
	/*unknown invalid entries returns '-99'*/
	else
	{
		RateCoefficient = -99;
	}

	ASSERT( RateCoefficient < 1e-6 );

	return RateCoefficient;
}

/**Badnell_RR_rate_eval
\param z_val atomic number on C scale - He - 1
\param n_val number of core electrons before capture of free electron
*/
STATIC double Badnell_RR_rate_eval(
			/* atomic number on C scale - He - 1 */
			int nAtomicNumberCScale, 
			/* number of core electrons before capture of free electron */
			int n_core_e_before_recomb )
{
	double RateCoefficient;
	double B, D, F;

	DEBUG_ENTRY( "Badnell_RR_rate_eval()" );

	ASSERT( nAtomicNumberCScale>=0 && nAtomicNumberCScale<LIMELM );

	if( nAtomicNumberCScale==ipIRON && 
		n_core_e_before_recomb>=12 && n_core_e_before_recomb<=18 )
	{
		/* RR rate coefficients from Table 3 of
		*>>refer	Fe	RR	Badnell, N. 2006, ApJ, 651, L73 
		* Fe 8+ to Fe 12+, but also include Fe13+ and Fe 14+,
		* so these are 26-8=18 to 26-14=12 
		* increasing number of bound electrons, 0 is 14 elec, 1 is 15 elec 
		* Fe 3p^q, q=2-6
		* this is DR fit coefficients given in table 1 of Badnell 06 */
		static const double parFeq[7][6] ={
			{1.179e-9 , 0.7096, 4.508e2, 3.393e7, 0.0154, 3.977e6},
			{1.050e-9 , 0.6939, 4.568e2, 3.987e7, 0.0066, 5.451e5},
			{9.832e-10, 0.7146, 3.597e2, 3.808e7, 0.0045, 3.952e5},
			{8.303e-10, 0.7156, 3.531e2, 3.554e7, 0.0132, 2.951e5},
			{1.052e-9 , 0.7370, 1.639e2, 2.924e7, 0.0224, 4.291e5},
			{1.338e-9 , 0.7495, 7.242e1, 2.453e7, 0.0404, 4.199e5},
			{1.263e-9 , 0.7532, 5.209e1, 2.169e7, 0.0421, 2.917e5}
		};

		double temp;
		/* nion is for the above block of numbers */
		long int nion = n_core_e_before_recomb - 12;
		ASSERT( nion>=0 && nion <=6 );

		temp = -parFeq[nion][5]/phycon.te; /* temp = (-T2/T) */
		B = parFeq[nion][1] + parFeq[nion][4]*exp(temp);
		D = sqrt(phycon.te/parFeq[nion][2]); /* D = (T/T0)^1/2 */
		F = sqrt(phycon.te/parFeq[nion][3]); /* F = (T/T1)^1/2 */
		RateCoefficient = parFeq[nion][0]/(D*pow((1.+D),(1.-B))*pow((1.+F),(1.+B)));
		strcpy(chRRDataSource[nAtomicNumberCScale][nAtomicNumberCScale-n_core_e_before_recomb] ,"Bad06");

		return RateCoefficient;
	}

	/*Invalid entries returns '-2':if the z_values are smaller than equal to the n_values */
	else if( nAtomicNumberCScale < n_core_e_before_recomb )
	{
		RateCoefficient = -2;
	}
	/*Invalid entries returns '-2' if nAtomicNumberCScale and n_core_e_before_recomb are out of the range*/
	else if( nAtomicNumberCScale >= LIMELM )
	{
		RateCoefficient = -2;
	}
	/*undefined z and n returns '-1'*/
	else if( !lgRRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{
		RateCoefficient = -1;
	}
	/* coefficients:A=RRFitPar[0], B=RRFitPar[1], T0=RRFitPar[2], T1=RRFitPar[3], DRFitParPart1=RRFitPar[4], T2=RRFitPar[5] */
	else if( lgRRBadnellDefined[nAtomicNumberCScale][n_core_e_before_recomb] )
	{

		/* RateCoefficient=A*[(T/T0)^1/2*(1+(T/T0)^1/2)^1-B*(1+(T/T1)^1/2)^1+B]^-1 
		where B = B + DRFitParPart1*exp(-T2/T) */
		double temp;

		temp = -RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][5]/phycon.te; /* temp = (-T2/T) */
		B = RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][1] + 
			RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][4]*exp(temp);
		D = sqrt(phycon.te/RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][2]); /* D = (T/T0)^1/2 */
		F = sqrt(phycon.te/RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][3]); /* F = (T/T1)^1/2 */
		RateCoefficient = RRFitPar[nAtomicNumberCScale][n_core_e_before_recomb][0]/(D*pow((1.+D),(1.-B))*pow((1.+F),(1.+B)));
		strcpy(chRRDataSource[nAtomicNumberCScale][nAtomicNumberCScale-n_core_e_before_recomb] ,"Bad06");
	}

	/*unknown invalid entries returns '-99'*/
	else
		RateCoefficient = -99;

	return RateCoefficient;
}


/*Badnell_rec_init This code is written by Terry Yun, 2005 *
 * It reads rate coefficient fits into 3D arrays and output array.out for testing *
 * The testing can be commented out */
void Badnell_rec_init( void )
{
	double par_C[MAX_FIT_PAR_DR];
	double par_E[MAX_FIT_PAR_DR];
	int NuclearCharge=-1, NumberElectrons=-1;
	int count, number;
	double temp_par[MAX_FIT_PAR_RR];
	int M_state, W_state;

	const int NBLOCK = 2;
	int data_begin_line[NBLOCK];/*it tells you where the data set begins(begins with 'Z')*/
	int length_of_line;	/*this variable checks for a blank line*/
	FILE *ioDATA;
	const char* chFilename;
	int yr, mo, dy;
	size_t chs;

	const int BIGGEST_INDEX_TO_USE = 103;

	/* Declaration of data file name array - done by Kausalya */
	long TheirIndexToOurIndex[BIGGEST_INDEX_TO_USE];
	string chLine;
	double value;
	bool lgEOL;
	long int i1;
	long INDX=0,INDP=0,N=0,S=0,L=0,J=0,maxINDX=0,loopindex=0,max_N_of_data=-1;
	bool lgFlag = true;

	static int nCalled = 0;

	static const char* cdDATAFILE[] = 
	{
		/* the list of filenames for Badnell DR, one to two electron */
		"",
		"UTA/nrb00_h_he1ic12.dat",
		"UTA/nrb00_h_li2ic12.dat",
		"UTA/nrb00_h_be3ic12.dat",
		"UTA/nrb00_h_b4ic12.dat",
		"UTA/nrb00_h_c5ic12.dat",
		"UTA/nrb00_h_n6ic12.dat",
		"UTA/nrb00_h_o7ic12.dat",
		"UTA/nrb00_h_f8ic12.dat",
		"UTA/nrb00_h_ne9ic12.dat",
		"UTA/nrb00_h_na10ic12.dat",
		"UTA/nrb00_h_mg11ic12.dat",
		"UTA/nrb00_h_al12ic12.dat",
		"UTA/nrb00_h_si13ic12.dat",
		"UTA/nrb00_h_p14ic12.dat",
		"UTA/nrb00_h_s15ic12.dat",
		"UTA/nrb00_h_cl16ic12.dat",
		"UTA/nrb00_h_ar17ic12.dat",
		"UTA/nrb00_h_k18ic12.dat",
		"UTA/nrb00_h_ca19ic12.dat",
		"UTA/nrb00_h_sc20ic12.dat",
		"UTA/nrb00_h_ti21ic12.dat",
		"UTA/nrb00_h_v22ic12.dat",
		"UTA/nrb00_h_cr23ic12.dat",
		"UTA/nrb00_h_mn24ic12.dat",
		"UTA/nrb00_h_fe25ic12.dat",
		"UTA/nrb00_h_co26ic12.dat",
		"UTA/nrb00_h_ni27ic12.dat",
		"UTA/nrb00_h_cu28ic12.dat",
		"UTA/nrb00_h_zn29ic12.dat"
	};
	//End of modification

	DEBUG_ENTRY( "Badnell_rec_init()" );

	/* must only do this once */
	if( nCalled > 0 )
	{
		return;
	}
	++nCalled;

#	if defined(PRINT_DR) || defined(PRINT_RR)
	FILE *ofp = open_data( FILE_NAME_OUT, "w" );
#	endif

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					for( long k=0; k<NUM_DR_TEMPS; ++k )
						iso_sp[ipISO][nelem].fb[ipHi].DielecRecombVsTemp[k] = 0.;
				}
			}
		}
	}

	/* Modification  done by Kausalya 
	 * Start - Try to open all the 29 data files.*/
	for( long nelem=ipHELIUM; nelem<LIMELM; nelem++)
	{
		if( dense.lgElmtOn[nelem] )
		{
			ioDATA = open_data( cdDATAFILE[nelem], "r" );

			lgFlag = true;
			ASSERT(ioDATA);

			for( long i=0; i<BIGGEST_INDEX_TO_USE; i++ )
				TheirIndexToOurIndex[i] = -1;

			/* Reading lines */
			while(lgFlag)
			{
				if( read_whole_line(chLine,ioDATA) )
				{
					if( chLine.find("INDX  INDP ") != string::npos )
					{
						/* ignore next line of data */
						if( !read_whole_line( chLine, ioDATA ) )
						{
							fprintf( ioQQQ, " Badnell data file appears to be corrupted.\n");
							cdEXIT(EXIT_FAILURE);
						}

						/* This one should be real data */
						while( read_whole_line(chLine, ioDATA) )
						{
							if( chLine == "\n" )
							{
								lgFlag = false;
								break;
							}

							i1=3;
							INDX=(long)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
							if( INDX >= BIGGEST_INDEX_TO_USE )
							{
								INDX--;
								lgFlag = false;
								break;
							}

							ASSERT( INDX < BIGGEST_INDEX_TO_USE );									 

							INDP=(long)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
							ASSERT( INDP >= 1 );									 

							if(INDP==1)
							{
								size_t p;
								if( (p=chLine.find("1S1 ")) != string::npos )
								{
									i1 = p+4+1;
									N=(long)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
									ASSERT( N>=1 );
								}
								else
								{
									TotalInsanity();
								}

								if( (p=chLine.find("     (")) != string::npos )
								{
									i1 = p+6+1;
									S=(long)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
									/* S in file is 3 or 1, we need 1 or 0 */
									ASSERT( S==1 || S==3 );
								}
								else
								{
									TotalInsanity();
								}

								/* move i1 one further to get L */
								i1++;
								L=(long)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
								ASSERT( L >= 0 && L < N );

								/* move i1 two further to get J */
								i1 += 2;
								J=(long)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
								ASSERT( Triangle2(2*L, S-1, 2*J) );

								/* if line in data file is higher N than highest considered, stop reading.  */
								if( N <= iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max + iso_sp[ipHE_LIKE][nelem].nCollapsed_max )
								{
									TheirIndexToOurIndex[INDX] = iso_sp[ipHE_LIKE][nelem].QN2Index(N, L, S, 2*J+1);
									ASSERT( TheirIndexToOurIndex[INDX] >= 0 );
								}
								else
								{
									/* Current line is not being used, 
									 * decrement INDX so maxINDX is set correctly below. */
									INDX--;
									lgFlag = false;
									break;
								}

								max_N_of_data = MAX2( max_N_of_data, N );
							}    
							else                                                                                      
							{
								// Stop parsing the tuple since INDP!=1
								lgFlag = false;   
							}     
						}                              
					}                                                                           
				}
				else
				{    
					// End of file is reached.
					lgFlag = false;
				}   
			}

			maxINDX =INDX;
			ASSERT( maxINDX > 0 );
			ASSERT( maxINDX < BIGGEST_INDEX_TO_USE );
			/* reset INDX */
			INDX = 0;
			lgFlag = true;
			while(lgFlag)
			{
				if( read_whole_line(chLine,ioDATA) )
				{
					/* to access the first table whose columns are INDX ,INDP */
					if( chLine.find("INDX TE= ") != string::npos )
					{
						lgFlag = false;
						/* we found the beginning of the data array */
						/* ignore next line of data */
						if( !read_whole_line( chLine, ioDATA ) )
						{
							fprintf( ioQQQ, " Badnell data file appears to be corrupted.\n");
							cdEXIT(EXIT_FAILURE);
						}

						/* This one should be real data */
						while( read_whole_line(chLine, ioDATA) )
						{
							/* If we find this string, we have reached the end of the table. */
							if( chLine.find("PRTF") != string::npos || INDX >= maxINDX || INDX<0 )
								break;

							i1=3;
							INDX=(long)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
							if( INDX>maxINDX )
								break;
							
							freeBound *fb;
							
							if( TheirIndexToOurIndex[INDX] < iso_sp[ipHE_LIKE][nelem].numLevels_max && 
								TheirIndexToOurIndex[INDX] > 0 )
								fb = &iso_sp[ipHE_LIKE][nelem].fb[TheirIndexToOurIndex[INDX]];
							else
								continue;

							for(loopindex=0;loopindex<10;loopindex++)
							{
								value=(double)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
								fb->DielecRecombVsTemp[loopindex] += value;
							}

							/* data are broken into two lines, read second line here */
							if( !read_whole_line( chLine, ioDATA ) )
							{
								fprintf( ioQQQ, " Badnell data file appears to be corrupted.\n");
								cdEXIT(EXIT_FAILURE);
							}

							/* start of data for second line */
							i1 = 13;
							for(loopindex=10;loopindex<19;loopindex++)
							{
								value=(double)FFmtRead(chLine.c_str(),&i1,chLine.length(),&lgEOL);
								fb->DielecRecombVsTemp[loopindex] += value;
							}
						}
					}
				}
				else
					lgFlag = false;
			}  
			fclose(ioDATA);
			ASSERT( maxINDX > 0 );
			ASSERT( maxINDX < BIGGEST_INDEX_TO_USE );
			ASSERT( max_N_of_data > 0 );

			if( max_N_of_data < iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max + iso_sp[ipHE_LIKE][nelem].nCollapsed_max )
			{
				long indexOfMaxN;
				L = -1;
				S = -1; 

				/* This loop extrapolates nLS data to nLS states */
				for( long i=TheirIndexToOurIndex[maxINDX]+1; 
					i<iso_sp[ipHE_LIKE][nelem].numLevels_max-iso_sp[ipHE_LIKE][nelem].nCollapsed_max; i++ )
				{
					long ipISO = ipHE_LIKE;
					L = L_(i);
					S = S_(i);

					if( L > 4 )
						continue;

					indexOfMaxN = iso_sp[ipHE_LIKE][nelem].QN2Index(max_N_of_data, L, S);
					for(loopindex=0;loopindex<19;loopindex++)
					{
						iso_sp[ipHE_LIKE][nelem].fb[i].DielecRecombVsTemp[loopindex] = 
							iso_sp[ipHE_LIKE][nelem].fb[indexOfMaxN].DielecRecombVsTemp[loopindex] *
							pow3( (double)max_N_of_data/(double)iso_sp[ipHE_LIKE][nelem].st[i].n());
					}
				}

				/* Get the N of the highest resolved singlet P (in the model, not the data) */
				indexOfMaxN = 
					iso_sp[ipHE_LIKE][nelem].QN2Index(iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max, 1, 1);

				/* This loop extrapolates nLS data to collapsed n levels, just use highest singlet P data */
				for( long i=iso_sp[ipHE_LIKE][nelem].numLevels_max-iso_sp[ipHE_LIKE][nelem].nCollapsed_max;
					i<iso_sp[ipHE_LIKE][nelem].numLevels_max; i++ )
				{
					for(loopindex=0;loopindex<19;loopindex++)
					{
						iso_sp[ipHE_LIKE][nelem].fb[i].DielecRecombVsTemp[loopindex] = 
							iso_sp[ipHE_LIKE][nelem].fb[indexOfMaxN].DielecRecombVsTemp[loopindex] *
							pow3( (double)iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max/
							(double)iso_sp[ipHE_LIKE][nelem].st[i].n());
					}
				}
			}
		}
	}

	for( long i=0; i<NBLOCK; ++i )
	{
		/* set to really large negative number - crash if used before being redefined */
		data_begin_line[i] = INT_MIN;
	}

	chFilename = "badnell_dr.dat";
	ioDATA = open_data( chFilename, "r" );

	count = 0;
	number = 0;

	/*Find out the number line where the data starts 
	 * there are two main blocks of data and each starts with a Z in column 2 */
	while( read_whole_line(chLine, ioDATA) )
	{
		count++;

		if( chLine.find("Z ")  != string::npos )
		{
			/* number has to be 0 or 1, and indicates the first or second block of data
			 * count is the line number for the start of that block */
			data_begin_line[number] = count;
			ASSERT( number < NBLOCK );
			number++;
		}
	}

	/*set a flag for a undefined data*/
	if( lgMustAllocateRec )
	{
		nDRFitPar.reserve(LIMELM);
		DRFitParPart1.reserve(LIMELM);
		RRFitPar.reserve(LIMELM);
		for( long nelem=0; nelem < LIMELM; nelem++ )
		{
			nDRFitPar.reserve(nelem, nelem+1);
			DRFitParPart1.reserve(nelem, nelem+1);
			RRFitPar.reserve(nelem, nelem+1);
			for( long ion=0; ion < nelem+1; ++ion )
			{
				DRFitParPart1.reserve(nelem, ion, MAX_FIT_PAR_DR);
				RRFitPar.reserve(nelem, ion, MAX_FIT_PAR_RR);
			}
		}
		nDRFitPar.alloc();
		lgDRBadnellDefined.alloc(nDRFitPar.clone());
		lgDRBadnellDefinedPart2.alloc(nDRFitPar.clone());
		lgRRBadnellDefined.alloc(nDRFitPar.clone());
		lgDR_BadWeb_exist.alloc(nDRFitPar.clone());
		DRFitParPart1.alloc();
		DRFitParPart2.alloc(DRFitParPart1.clone());
		RRFitPar.alloc();
	
		lgMustAllocateRec = false;
	}

	lgDRBadnellDefined = false;
	lgDRBadnellDefinedPart2 = false;
	lgRRBadnellDefined = false;
	/*set fitting coefficients to zero initially*/
	DRFitParPart1 = 0.;
	DRFitParPart2 = 0.;
	RRFitPar = 0.;

	count = 0;
	/*Start from beginning to read in again*/  
	fseek(ioDATA, 0, SEEK_SET);
	/* read magic number for DR data */
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init could not read first line of badnell_dr.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	count++;

	/* look for ')' on the line, magic number comes after it */
	if( (chs = chLine.find(")")) == string::npos )
	{
		/* format is incorrect */
		fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init data file incorrect format.\n");
		cdEXIT(EXIT_FAILURE);
	}

	++chs;
	sscanf(chLine.c_str()+chs, "%4i%2i%2i",&yr, &mo, &dy);
	/* check magic number - the date on the line */
	int dr_yr = 2023, dr_mo = 5, dr_dy = 12;
	if((yr != dr_yr) || (mo != dr_mo) || (dy != dr_dy))
	{
		fprintf(ioQQQ,
			"DISASTER PROBLEM Badnell_rec_init The version of %s I found (%i %i %i) is not the current version (%i %i %i).\n",
			chFilename, yr, mo, dy, dr_yr, dr_mo, dr_dy);
		fprintf(ioQQQ," The first line of the file is the following\n %s\n", chLine.c_str() );
		cdEXIT(EXIT_FAILURE);
	}

	while( read_whole_line(chLine, ioDATA) )
	{
		count++;
		length_of_line = (int)chLine.length();

		/*reading in coefficients DRFitParPart1 */
		if( count > data_begin_line[0] && count < data_begin_line[1] && length_of_line >3 )
		{
			/*set array par_C to zero */
			for( long i=0; i<MAX_FIT_PAR_DR; i++ )
			{
				par_C[i] = 0;
			}
			sscanf(chLine.c_str(), "%i%i%i%i%lf%lf%lf%lf%lf%lf%lf%lf%lf",
				&NuclearCharge, &NumberElectrons, &M_state, &W_state, &par_C[0], &par_C[1], &par_C[2],
				&par_C[3], &par_C[4], &par_C[5], &par_C[6], &par_C[7], &par_C[8]);
			/* data files have atomic number on physics scale, convert to C scale
			 * for following code */
			long int NuclearChargeM1 = NuclearCharge-1;

			if(M_state == 1 && NuclearChargeM1 < LIMELM )
			{
				/*Set a flag to '1' when the indices are defined */
				ASSERT( NumberElectrons < LIMELM );
				ASSERT( NuclearChargeM1 < LIMELM );
				lgDRBadnellDefined[NuclearChargeM1][NumberElectrons] = true;

				/*counting the number of coefficients */
				nDRFitPar[NuclearChargeM1][NumberElectrons] = 9;
				for( long i=8; i>=0; i-- )
				{
					if( par_C[i] == 0 )
						--nDRFitPar[NuclearChargeM1][NumberElectrons];
					else
						break;
				}

				/*assign the values into array */
				for( long i=0; i<9; i++ )
					DRFitParPart1[NuclearChargeM1][NumberElectrons][i] = par_C[i];
			}
		}
	}

	/*starting again to read in E values */
	fseek(ioDATA, 0, SEEK_SET);
	count = 0; 
	while( read_whole_line(chLine, ioDATA) )
	{
		count++;
		length_of_line = (int)chLine.length();
		if( count > data_begin_line[1] && length_of_line > 3 )
		{

			/*set array par_E to zero*/
			for( long i=0; i<MAX_FIT_PAR_DR; i++ )
			{
				par_E[i] = 0;
			}
			sscanf(chLine.c_str(), "%i%i%i%i%lf%lf%lf%lf%lf%lf%lf%lf%lf",
				&NuclearCharge, &NumberElectrons, &M_state, &W_state, &par_E[0], &par_E[1], &par_E[2],
				&par_E[3], &par_E[4], &par_E[5], &par_E[6], &par_E[7], &par_E[8]);
			/* data file is on physics scale but we use C scale */
			long int NuclearChargeM1 = NuclearCharge-1;

			if(M_state == 1 && NuclearChargeM1<LIMELM)
			{
				ASSERT( NumberElectrons < LIMELM );
				ASSERT( NuclearChargeM1 < LIMELM );
				lgDRBadnellDefinedPart2[NuclearChargeM1][NumberElectrons] = true;

				/*counting the number of coefficients */
				nDRFitPar[NuclearChargeM1][NumberElectrons] = 9;
				for( long i=8; i>=0; i-- )
				{
					if( par_E[i] == 0 )
						--nDRFitPar[NuclearChargeM1][NumberElectrons];
					else
						break;
				}

				/*assign the values into array*/
				for( long i=0; i<nDRFitPar[NuclearChargeM1][NumberElectrons]; i++ )
					DRFitParPart2[NuclearChargeM1][NumberElectrons][i] = par_E[i];
			}
		}
	}

	fclose( ioDATA );

	/*output coefficients for defined values for testing */
#	ifdef PRINT_DR
	for( long nelem=0; nelem<LIMELM; nelem++ )
	{
		for( int ion=0; ion<nelem+1;++ion )
		{
			if( lgDRBadnellDefined[nelem][ion] )
			{
				fprintf(ofp, "%i %i %e %e %e %e %e %e %e %e %e\n",
					nelem, ion, DRFitParPart1[nelem][ion][0], 
					DRFitParPart1[nelem][ion][1], DRFitParPart1[nelem][ion][2], 
					DRFitParPart1[nelem][ion][3], DRFitParPart1[nelem][ion][4],
					DRFitParPart1[nelem][ion][5], DRFitParPart1[nelem][ion][6], 
					DRFitParPart1[nelem][ion][7], DRFitParPart1[nelem][ion][8]);
			}
		}
	}
	for( long nelem=0; nelem<LIMELM; nelem++ )
	{
		for( int ion=0; ion<nelem+1; ion++ )
		{
			if( lgDRBadnellDefinedPart2[nelem][ion] )
			{
				fprintf(ofp, "%i %i %e %e %e %e %e %e %e %e %e\n",
					nelem, ion, DRFitParPart2[nelem][ion][0], 
					DRFitParPart2[nelem][ion][1], DRFitParPart2[nelem][ion][2], 
					DRFitParPart2[nelem][ion][3], DRFitParPart2[nelem][ion][4],
					DRFitParPart2[nelem][ion][5], DRFitParPart2[nelem][ion][6], 
					DRFitParPart2[nelem][ion][7], DRFitParPart2[nelem][ion][8]);
			}
		}
	}
	fclose(ofp);
#	endif

	/*checking for the match of lgDRBadnellDefined and lgDRBadnellDefinedPart2 - 
	 * Both have to be defined*/
	bool lgDRBadnellBothDefined = true;
	for( int nelem=0; nelem<LIMELM; nelem++ )
	{
		for( int ion=0; ion<nelem+1; ion++ )
		{
			/* check that first and second half of DR fitting coefficients 
			 * are both defined */
			if( lgDRBadnellDefined[nelem][ion] != lgDRBadnellDefinedPart2[nelem][ion] )
			{
				fprintf( ioQQQ, "DR %i, RR %i: %c %c\n", nelem, ion, 
					 TorF(lgDRBadnellDefined[nelem][ion]), 
					 TorF(lgDRBadnellDefinedPart2[nelem][ion]));
				fprintf( ioQQQ, "PROBLEM ion_recomb_Badnell first and second half of Badnell DR not consistent.\n");
				lgDRBadnellBothDefined = false;
			}
		}
	}

	if( !lgDRBadnellBothDefined )
	{
		/* disaster - DR files are not consistent */
		fprintf(ioQQQ,
			"DISASTER PROBLEM The DR data files are corrupted - part 1 and 2 do not agree.\n");
		fprintf(ioQQQ," Start again with a fresh copy of the data directory\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* now do radiative recombination */
	chFilename = "badnell_rr.dat";
	ioDATA = open_data( chFilename, "r" );

	/* read magic number for RR data */
	{
		if( !read_whole_line( chLine, ioDATA ) )
		{
			fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init could not read first line of badnell_rr.dat.\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* this is just before date, which we use for magic number */
		if( (chs = chLine.find(")")) == string::npos )
		{
			/* format is incorrect */
			fprintf( ioQQQ, " DISASTER PROBLEM Badnell_rec_init data file incorrect format.\n");
			cdEXIT(EXIT_FAILURE);
		}
		++chs;
		sscanf(chLine.c_str()+chs, "%4i%2i%2i", &yr, &mo, &dy);
		int rr_yr = 2023, rr_mo = 5, rr_dy = 11;
		if((yr != rr_yr)||(mo != rr_mo)||(dy != rr_dy))
		{
			fprintf(ioQQQ,"DISASTER PROBLEM The version of %s I found (%i %i %i) is not the current version (%i %i %i).\n",
				chFilename, yr, mo, dy, rr_yr, rr_mo, rr_dy);
			fprintf(ioQQQ," The line was as follows:\n %s\n", chLine.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
	}

	while( read_whole_line(chLine, ioDATA) )
	{
		/*read in coefficients - first set array par to zero */
		for( long i=0; i<MAX_FIT_PAR_RR; i++ )
		{
			temp_par[i] = 0;
		}
		if(chLine[0] != '#')
		{
			sscanf(chLine.c_str(), "%i%i%i%i%lf%lf%lf%lf%lf%lf",
				&NuclearCharge, &NumberElectrons, &M_state, &W_state, &temp_par[0], &temp_par[1],
				&temp_par[2], &temp_par[3], &temp_par[4], &temp_par[5]);
			long NuclearChargeM1 = NuclearCharge-1;

			if(M_state == 1 && NuclearChargeM1<LIMELM)
			{
				ASSERT( NuclearChargeM1 < LIMELM );
				ASSERT( NumberElectrons <= LIMELM );
				/*Set a flag to '1' when the indices are defined */  
				lgRRBadnellDefined[NuclearChargeM1][NumberElectrons] = true;
				/*assign the values into array */
				for( long i=0; i<MAX_FIT_PAR_RR; i++ )
					RRFitPar[NuclearChargeM1][NumberElectrons][i] = temp_par[i];
			}
		}
	}

	/*output coefficients for defined values for testing */
#	ifdef PRINT_RR
	count = 0;
	for( long nelem=0; nelem<LIMELM; nelem++ )
	{
		for( long ion=0; ion<nelem+1; ion++ )
		{
			if( lgRRBadnellDefined[nelem][ion] )
			{
				fprintf(ofp, "%li %li %e %e %e %e %e %e\n",
					nelem, ion, RRFitPar[nelem][ion][0], 
					RRFitPar[nelem][ion][1], RRFitPar[nelem][ion][2], 
					RRFitPar[nelem][ion][3],
					RRFitPar[nelem][ion][4], RRFitPar[nelem][ion][5]);
				count++;
			}
		}
	}
	fprintf(ofp, "total lines are %i ", count);

	fclose(ofp);
#	endif

	fclose(ioDATA);

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			for( int nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
			{
				fprintf(ioQQQ,"\nDEBUG rr rec\t%i",nelem);
				for( int ion=0; ion<=nelem; ++ion )
				{
					fprintf(ioQQQ,"\t%.2e", Badnell_RR_rate_eval(nelem+1 , nelem-ion ) );
				}
				fprintf(ioQQQ,"\n");
				fprintf(ioQQQ,"DEBUG dr rec\t%i",nelem);
				for( int ion=0; ion<=nelem; ++ion )
				{
					fprintf(ioQQQ,"\t%.2e", Badnell_DR_rate_eval(nelem+1 , nelem-ion ) );
				}
				fprintf(ioQQQ,"\n");
			}
			cdEXIT(EXIT_FAILURE);
		}
	}

	// gaussian noise for dielectronic recombination coefficients guesses
	// set with SET DIELECTRONIC RECOMBINATION NOISE command
	if( ionbal.guess_noise !=0. )
	{
		for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
			/* log normal noise with dispersion entered on command line */
			/* NB the seed was set when the command was parsed */
			RecNoise[nelem] = exp10( ran.normal()*ionbal.guess_noise );
	}
	else
	{
		for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
				RecNoise[nelem] = 1.;
	}

	// initialize some products
	for( long nelem=0; nelem<LIMELM; ++nelem )
		DR_Badnell_rate_coef_mean_ion[nelem] = 0.;

	return;
}

/*ion_recom_calculate calculate radiative and dielectronic recombination rate coefficients */
void ion_recom_calculate( void )
{
	static double TeUsed = -1 , EdenUsed = -1.;

	DEBUG_ENTRY( "ion_recom_calculate()" );

	/* do not reevaluate if change in temperature is small */
	if( fp_equal(phycon.te,TeUsed) && fp_equal( dense.eden , EdenUsed ))
		return;

	TeUsed = phycon.te;
	EdenUsed = dense.eden;

	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{

		for( long ion=0; ion < nelem+1; ++ion )
		{
			long int n_bnd_elec_before_recom ,
				n_bnd_elec_after_recom;

			n_bnd_elec_before_recom = nelem-ion;
			n_bnd_elec_after_recom = nelem-ion+1;

			// will insure these are >=0 at end
			ionbal.DR_Badnell_rate_coef[nelem][ion] = -1.;
			ionbal.RR_rate_coef_used[nelem][ion] = 0.;
			strcpy( chDRDataSource[nelem][ion] , "none" );
			strcpy( chRRDataSource[nelem][ion] , "none" );

			/* Badnell dielectronic recombination rate coefficients */
			if( (ionbal.DR_Badnell_rate_coef[nelem][ion] = 
				Badnell_DR_rate_eval(
				/* atomic number on C scale */
				nelem, 
				/* number of core electrons before capture of free electron,
				 * for bare ion this is zero */
				n_bnd_elec_before_recom )) >= 0. )
			{
				lgDR_BadWeb_exist[nelem][ion] = true;
			}
			else
			{
				/* real rate does not exist, will use mean later */
				lgDR_BadWeb_exist[nelem][ion] = false;
			}

			/* save D. Verner's radiative recombination rate coefficient
			* needed for rec cooling, cm3 s-1 */
			ionbal.RR_Verner_rate_coef[nelem][ion] =
				t_ADfA::Inst().rad_rec(
				/* number number of physics scale */
				nelem+1 ,
				/* number of protons on physics scale */
				n_bnd_elec_after_recom ,
				phycon.te );

			/* Badnell radiative recombination rate coefficients */
			if( (ionbal.RR_Badnell_rate_coef[nelem][ion] = Badnell_RR_rate_eval(
				/* atomic number on C scale */
				nelem, 
				/* number of core electrons before capture of free electron */
				n_bnd_elec_before_recom )) >= 0. )
			{
				ionbal.RR_rate_coef_used[nelem][ion] = ionbal.RR_Badnell_rate_coef[nelem][ion];
			}
			else
			{
				strcpy( chRRDataSource[nelem][ion] , "Verner" );
				ionbal.RR_rate_coef_used[nelem][ion] = ionbal.RR_Verner_rate_coef[nelem][ion];
			}
		}
		// recombination to bare nuclei has no DR
		ionbal.DR_Badnell_rate_coef[nelem][nelem] = 0.;
		strcpy(chDRDataSource[nelem][nelem] , "NA");
	}

	/* this branch starts idiosyncratic single ions */
	static const double Fe_Gu_c[9][6] = {
		{ 2.50507e-11, 5.60226e-11, 1.85001e-10, 3.57495e-9, 1.66321e-7, 0. },/*fit params for Fe+6*/
		{ 9.19610e-11, 2.92460e-10, 1.02120e-9, 1.14852e-8, 3.25418e-7, 0. }, /* fitting params for Fe+7 */
		{ 9.02625e-11, 6.22962e-10, 5.77545e-9, 1.78847e-8, 3.40610e-7, 0. }, /* fitting params for Fe+8 */
		{ 9.04286e-12, 9.68148e-10, 4.83636e-9, 2.48159e-8, 3.96815e-7, 0. }, /* fitting params for Fe+9 */
		{ 6.77873e-10, 1.47252e-9, 5.31202e-9, 2.54793e-8, 3.47407e-7, 0. }, /* fitting params for Fe+10 */
		{ 1.29742e-9, 4.10172e-9, 1.23605e-8, 2.33615e-8, 2.97261e-7, 0. }, /* fitting params for Fe+11 */
		{ 8.78027e-10, 2.31680e-9, 3.49333e-9, 1.16927e-8, 8.18537e-8, 1.54740e-7 },/*fit params for Fe+12*/
		{ 2.23178e-10, 1.87313e-9, 2.86171e-9, 1.38575e-8, 1.17803e-7, 1.06251e-7 },/*fit params for Fe+13*/
		{ 2.17263e-10, 7.35929e-10, 2.81276e-9, 1.32411e-8, 1.15761e-7, 4.80389e-8 }/*fit params for Fe+14*/
	};

	static const double Fe_Gu_E[9][6] = {
		{ 8.30501e-2, 8.52897e-1, 3.40225e0, 2.23053e1, 6.80367e1, 0. }, /* fitting params for Fe+6 */
		{ 1.44392e-1, 9.23999e-1, 5.45498e0, 2.04301e1, 7.06112e1, 0. }, /* fitting params for Fe+7 */ 
		{ 5.79132e-2, 1.27852e0, 3.22439e0, 1.79602e1, 6.96277e1, 0. }, /* fitting params for Fe+8 */
		{ 1.02421e-1, 1.79393e0, 4.83226e0, 1.91117e1, 6.80858e1, 0. }, /* fitting params for Fe+9 */
		{ 1.24630e-1, 6.86045e-1, 3.09611e0, 1.44023e1, 6.42820e1, 0. }, /* fitting params for Fe+10 */
		{ 1.34459e-1, 6.63028e-1, 2.61753e0, 1.30392e1, 6.10222e1, 0. }, /* fitting params for Fe+11 */
		{ 7.79748e-2, 5.35522e-1, 1.88407e0, 8.38459e0, 3.38613e1, 7.89706e1 }, /*fitting params for Fe+12*/
		{ 8.83019e-2, 6.12756e-1, 2.36035e0, 9.61736e0, 3.64467e1, 8.72406e1 }, /*fitting params for Fe+13*/
		{ 1.51322e-1, 5.63155e-1, 2.57013e0, 9.08166e0, 3.69528e1, 1.08067e2 } /* fitting params for Fe+14*/
	};

	/* do a series of special cases for Fe DR  */
	double te_eV32 = powpq(phycon.te_eV,3,2);
		
	/* >>chng 06 jul 07 by Mitchell Martin, added DR rate coefficient 
	* calculations for Fe+6->Fe+5 through Fe+14->Fe+13
	* this is still for nelem = ipIRON from the previous calculation 
	* starts with Fe+6 -> Fe+5 and does the next ion with each iteration */
	for( long ion=0; ion<9; ion++ )
	{
		/* only do this rate if not already done by a previous approximation */
		if( ionbal.DR_Badnell_rate_coef[ipIRON][ion+5]<0. )
		{
			double fitSum = 0; /* resets the fitting parameter calculation */
			for( long i=0; i<6; i++ )
			{
				fitSum += Fe_Gu_c[ion][i] * sexp( Fe_Gu_E[ion][i]/phycon.te_eV );
			}
			strcpy(chDRDataSource[ipIRON][ion+5] , "GuPC");
			lgDR_BadWeb_exist[ipIRON][ion+5] = true;
			ionbal.DR_Badnell_rate_coef[ipIRON][ion+5] = fitSum / te_eV32;
		}
	}
	/* this is end of Fe DR rates */

	// use C08 mean for stability
	static const double BadnelDR_RateSave[LIMELM] =
	{
			3.78e-13, 1.70e-12, 8.14e-12, 1.60e-11, 2.38e-11,
			6.42e-11, 5.97e-11, 1.47e-10, 1.11e-10, 3.26e-10,
			1.88e-10, 2.06e-10, 4.14e-10, 3.97e-10, 2.07e-10,
			2.46e-10, 3.38e-10, 3.15e-10, 9.70e-11, 6.49e-11,
			6.93e-10, 3.70e-10, 3.29e-11, 4.96e-11, 5.03e-11,
			2.91e-12, 4.62e-14, 0.00e+00, 0.00e+00, 0.00e+00
	};
	for( long nelem=0; nelem < LIMELM; ++nelem )
	{
		DR_Badnell_rate_coef_mean_ion[nelem] =
			BadnelDR_RateSave[nelem] * RecNoise[nelem] *
			// default of unity, set with SET DIELECTRONIC RECOMBINATION KLUDGE SCALE command
			ionbal.DR_mean_scale[nelem];
	}

	// iron is special case with Arnaud & Raymond 1992
	// use mean which is low T dr and AR which is high temp
	for( long ion=0; ion < ipIRON+1; ++ion )
	{
		if( ionbal.DR_Badnell_rate_coef[ipIRON][ion] < 0. )
		{
			ionbal.DR_Badnell_rate_coef[ipIRON][ion] =
					DR_Badnell_rate_coef_mean_ion[ion]
					+ atmdat_dielrec_fe(ion+1, phycon.te );
			strcpy(chDRDataSource[ipIRON][ion] , "mean+");
		}
	}
	// this routine will return something for all ions - even if just a guess
	for( long nelem=0; nelem < LIMELM; ++nelem )
	{
		for( long ion=0; ion < nelem+1; ++ion )
		{
			if( ionbal.DR_Badnell_rate_coef[nelem][ion] < 0. )
			{
				strcpy(chDRDataSource[nelem][ion] , "mean");
				ionbal.DR_Badnell_rate_coef[nelem][ion] = DR_Badnell_rate_coef_mean_ion[ion];
			}
		}
	}

	/* collisional suppression of DR, disable with command
	 * set dielectronic recombination suppression off
	 */
	if( ionbal.lgDRsup )
	{
		for( long nelem=ipHELIUM; nelem < LIMELM; ++nelem )
		{
			for( long ion=0; ion < nelem; ++ion )
			{
				// ASSERT(DielSupprsFactor[ion]>=0 && DielSupprsFactor[ion]<=1. );
				// old very simple expression
				//ionbal.DR_Badnell_rate_coef[nelem][ion] *= DielSupprsFactor[ion];

				// DR collisional suppression based on Badnell rates
				ionbal.DR_Badnell_suppress_fact[nelem][ion] = CollisSuppres(
						/* This routine takes the following arguments:
						 * atomic_number = nuclear charge */
						nelem+1,
						/*ionic_charge = ionic charge*/
						ion+1,
						/*eden = electron density */
						dense.eden,
						/*T = temperature (K)*/
						phycon.te );
				ionbal.DR_Badnell_rate_coef[nelem][ion] *=
					ionbal.DR_Badnell_suppress_fact[nelem][ion]; 

				ASSERT(ionbal.DR_Badnell_rate_coef[nelem][ion] >= 0);
				ASSERT(ionbal.RR_rate_coef_used[nelem][ion] >= 0);
			}
		}
	}

	static bool lgRunOnce = true;
	/* this set true with PRINT RECOMBINATION and SAVE DATA SOURCES commands */
	if( lgRunOnce )
	{
		lgRunOnce = false;
		if( ionbal.lgRecom_Badnell_print || save.lgSDSOn )
		{
			FILE *ioREC = ioQQQ;
			if( save.lgSDSOn )
				ioREC = save.ipSDSFile;
			fprintf(ioREC, "\n\n##################################################\n");
			fprintf(ioREC,"radiative recombination data sources \n" );

			for( long loop=0;loop<30;loop+=10)
			{
				fprintf(ioREC,"\n         ");
				for(long  ion=loop; ion<loop+10; ++ion )
				{
					fprintf(ioREC,"&%7li",ion);
				}
				fprintf(ioREC,"\\\\\n" );
				for( long nelem=loop; nelem<LIMELM; ++nelem )
				{
					fprintf(ioREC,"%2li %5s ",nelem+1 , elementnames.chElementNameShort[nelem] );
					long limit = MIN2(nelem+1,loop+10);
					for( long ion=loop; ion<limit; ++ion )
					{
						fprintf(ioREC,"&%7s",chRRDataSource[nelem][ion] );
					}
					for( long ion=limit; ion<loop+10; ++ion )
					{
						fprintf(ioREC,"&%7s",chRRDataSource[nelem][ion] );
					}
					fprintf(ioREC,"\\\\\n" );
				}
			}
			fprintf(ioREC,"\nRadiative recombination data sources\n");
			fprintf(ioREC,"Bad06: Badnell, N., 2006, ApJ, 167, 334B\n");
			fprintf(ioREC,"Verner: Verner & Ferland, 1996, ApJS, 103, 467\n");

			fprintf(ioREC, "\n\n##################################################\n");
			fprintf(ioREC,"Dielectronic recombination data sources \n" );

			for( long loop=0;loop<30;loop+=10)
			{
				fprintf(ioREC,"\n         ");
				for(long  ion=loop; ion<loop+10; ++ion )
				{
					fprintf(ioREC,"&%7li",ion);
				}
				fprintf(ioREC,"\\\\\n" );
				for( long nelem=loop; nelem<LIMELM; ++nelem )
				{
					fprintf(ioREC,"%2li %5s ",
						nelem+1 , elementnames.chElementNameShort[nelem] );
					long limit = MIN2(nelem+1,loop+10);
					for( long ion=loop; ion<limit; ++ion )
					{
						fprintf(ioREC,"&%7s",chDRDataSource[nelem][ion] );
					}
					for( long ion=limit; ion<loop+10; ++ion )
					{
						fprintf(ioREC,"&%7s",chDRDataSource[nelem][ion] );
					}
					fprintf(ioREC,"\\\\\n" );
				}
			}
			fprintf(ioREC,"\nDielectronic data sources\nBadWeb: Badnell web site http://amdpp.phys.strath.ac.uk/tamoc/DR/\n");
			fprintf(ioREC,"Bad06D: Badnell, N., 2006, ApJ, 651, L73\n");
			fprintf(ioREC,"GuPC: Gu, M. private communication\n");
			fprintf(ioREC,"mean: DR recombination ion mean (see below)\n");
			fprintf(ioREC,"mean+: DR recombination ion mean (see below) + Arnaud & Raymond 1992, ApJ, 398, 394 \n");

			// option to print rates - PRINT RECOMBINATION
			if( ionbal.lgRecom_Badnell_print )
			{
				fprintf(ioQQQ,"\n\nBadnell recombination RR, then DR, T=%.3e\n", phycon.te );
				for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
				{
					fprintf(ioQQQ,"nelem=%li %s, RR then DR\n",
						nelem+1, elementnames.chElementNameShort[nelem] );
					for( long ion=0; ion<nelem+1; ++ion )
					{
						fprintf(ioQQQ," %.2e", ionbal.RR_rate_coef_used[nelem][ion] );
					}
					fprintf(ioQQQ,"\n" );
					for( long ion=0; ion<nelem+1; ++ion )
					{
						fprintf(ioQQQ," %.2e", ionbal.DR_Badnell_rate_coef[nelem][ion] );
					}
					fprintf(ioQQQ,"\n\n" );
				}
				/* now print mean recombination and standard deviation */
				fprintf(ioQQQ,"mean DR recombination ion mean \n" );
				for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
				{
					fprintf(ioQQQ," %2li %.2e \n",
						nelem + 1,
						DR_Badnell_rate_coef_mean_ion[nelem]  );
				}

				fprintf( ioQQQ, "\n\nCollisSuppres finds following dielectronic"
					" recom suppression factors, Te=%10.3e eden=%10.3e\n",
					phycon.te, dense.eden );
				fprintf( ioQQQ, "nelem ion    fac \n" );
				for( long nelem=ipHELIUM; nelem < LIMELM; ++nelem )
				{
					for( long ion=0; ion < nelem; ++ion )
					{
						fprintf( ioQQQ, "%3ld %4ld %10.3e\n", nelem+1, ion,
							ionbal.DR_Badnell_suppress_fact[nelem][ion] );

					}
					fprintf( ioQQQ, "\n");
				}
			}
		}
	}
	return;
}

