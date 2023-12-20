/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/************************************************************************************************/
/*H_photo_cs returns hydrogenic photoionization cross section in cm^2                           */
/*H_Einstein_A calculates Einstein A for any nlz in s^-1                                        */
/*hv calculates photon energy in ergs for n -> n' transitions for H and H-like ions             */
/************************************************************************************************/
#include "cddefines.h"
#include "hydro_bauman.h"
#include "thirdparty.h"
#include "dense.h"
#include "hydro_tbl.h"
#include "vectorize.h"
#include "iso.h"
#include "service.h"

struct mx
{
	double m;
	long int x;
	mx() : m(0.), x(0L) {}
};

// we keep m above bounded to |m| <= 10^maxe == maxm
static const int maxe = 100;
static const double maxm = 1.e100;

struct mxq
{
	struct mx mx;
	long int q;
	mxq() : q(0L) {}
};

inline double exp10i(long x)
{
	// distribution of arguments in test:
	// x == 0 : 706479395
	// x < 0  :   5020125
	// x > 0  :     34810

	if( LIKELY(x == 0) )
		return 1.;
	else if( x < 0 )
	{
		// do not support gradual underflow here...
		if( UNLIKELY(x < min_pow10) )
			return 0.;
		else
			return neg_pow10[-x];
	}
	else
	{
		if( UNLIKELY(x > max_pow10) )
			return pow2(BIGDOUBLE);
		else
			return pos_pow10[x];
	}
}

/************************************************************************************************/
/*    these routines were written by Robert Bauman                                              */
/*    The main reference for this section of code is                                            */
/*    M. Brocklehurst                                                                           */
/*    Mon. Note. R. astr. Soc. (1971) 153, 471-490                                              */
/*                                                                                              */
/*    The recombination coefficient is obtained from the                                        */
/*    photoionization cross-section (see Burgess 1965).                                         */
/*    We have:                                                                                  */
/*                                                                                              */
/*                  -                           -               l'=l+1                          */
/*                 |  2 pi^(1/2) alpha^4 a_o^2 c |  2 y^(1/2)     ---                           */
/* alpha(n,l,Z,Te)=|-----------------------------|  ---------  Z   >    I(n, l, l', t)          */
/*                 |            3                |     n^2        ---                           */
/*                  -                           -               l'=l-1                          */
/*                                                                                              */
/*      where                         OO                                                        */
/*                                    -                                                         */
/*                                   |                                                          */
/*      I(n, l, l', t) = max(l,l') y | (1 + n^2 K^2)^2 Theta(n,l; K, l') exp( -K^2 y ) d(K^2)   */
/*                                   |                                                          */
/*                                  -                                                           */
/*                                  0                                                           */
/*                                                                                              */
/*      Here K = k / Z                                                                          */
/*                                                                                              */
/*                                                                                              */
/*      and                                                                                     */
/*                                                                                              */
/*                                                                                              */
/*       y =   Z^2 Rhc/(k Te)= 15.778/t                                                         */
/*                                                                                              */
/*       where "t" is the scaled temperature, and "Te" is the electron Temperature              */
/*                                                                                              */
/*       t = Te/(10^4  Z^2)                                                                     */
/*           Te in kelvin                                                                       */
/*                                                                                              */
/*                                           |              |^2                                 */
/*         Theta(n,l; K, l') = (1 + n^2 K^2) | g(n,l; K,l') |                                   */
/*                                           |              |                                   */
/*                                                                                              */
/*                                                                                              */
/*                                                 ---- Not sure if this is K or k              */
/*                                  OO            /     but think it is K                       */
/*         where                    -            v                                              */
/*                         K^2     |                                                            */
/*         g(n,l; K,l') = -----    | R_nl(r) F_(K,l) r dr                                       */
/*                         n^2     |                                                            */
/*                                -                                                             */
/*                                0                                                             */
/*                                                                                              */
/*                                                                                              */
/*         -                           -                                                        */
/*        |  2 pi^(1/2) alpha^4 a_o^2 c |                                                       */
/*        |-----------------------------|                                                       */
/*        |            3                |                                                       */
/*         -                           -                                                        */
/*                                                                                              */
/*         = 2 * (3.141592654)^1/2 * (7.29735308e-3)^4                                          */
/*                      * (0.529177249e-10)^2 * (2.99792458e8) / 3                              */
/*                                                                                              */
/*         = 2.8129897e-21                                                                      */
/*         Mathematica gives 2.4764282710571237e-21                                             */
/*                                                                                              */
/*    The photoionization cross-section (see also Burgess 1965)                                 */
/*     is given by;                                                                             */
/*                   _               _        l'=l+1                                            */
/*                  |4 PI alpha a_o^2 |  n^2   ---     max(l,l')                                */
/*      a(Z';n,l,k)=|---------------- |  ---    >       --------- Theta(n,l; K, l')             */
/*                  |       3         |  Z^2   ---     (2l + 1)                                 */
/*                   _               _        l'=l-1                                            */
/*                                                                                              */
/*                                                                                              */
/*      where  Theta(n,l; K, l') is defined above                                               */
/************************************************************************************************/
/************************************************************************************************/
/*      For the transformation:                                                                 */
/*                              Z -> rZ = Z'                                                    */
/*                                                                                              */
/*                              k -> rk = k'                                                    */
/*      then                                                                                    */
/*                                                                                              */
/*                              K -> K = K'                                                     */
/*                                                                                              */
/*      and the cross-sections satisfy;                                                         */
/*                                                1                                             */
/*                               a(Z'; n,l,k') = --- a(Z; n,l,k)                                */
/*                                               r^2                                            */
/************************************************************************************************/

/** \verbatim
************************* for LOG version of the file  ***************************************
 In this version, quantities that would normal cause a 64-bit floating point processor        
 to underflowed or overflow on intermediate values (ones internal to the calculation)         
 are evaluated using logs. This allows us to use an upper principal quantum number `n'        
 greater than 50 which is where the other version begins to fail. The trade-off is,           
 of course, lower accuracy( or is it precision ) and perhaps speed.                           
      We use LOG_10 for convenience.                                                          
**********************************************************************************************
 The functions which are evaluated using logarithms are denoted with a trailing underscore.   
      example:   hri_() calculates the same thing as hri_log10()                              
      except it uses logs internally.                                                         
**********************************************************************************************
 these are the hydrogenic routines written by Robert Bauman                                   
      For references, see h_bauman.c                                                          
**********************************************************************************************
  IN THE FOLLOWING WE HAVE  n > n'                                                            
\endverbatim
\return returns hydrogenic photoionization cross section in cm^2
\param photon_energy incident photon energy
\param n principal quantum number, 1 for ground
\param l angular momentum, 0 for s
\param iz charge, 1 for H+, 2 for He++, etc
*/
double H_photo_cs_log10(
	double photon_energy, 
	long int n, 
	long int l, 
	long int iz
);

/**\verbatim
  Calculates the Einstein A's for hydrogen                           
   for the transition n,l --> n',l'                                   
  units of sec^(-1)                                                  

  In the following, we have n > n'
  \endverbatim                                
  \param n principal quantum number, 1 for ground, upper level
  \param l angular momentum, 0 for s
  \param np principal quantum number, 1 for ground, lower level
  \param lp angular momentum, 0 for s
  \param iz Nuclear charge, 1 for H+, 2 for He++, etc
  */
double H_Einstein_A_log10(
	long int n,
	long int l,
	long int np,
	long int lp,
	long int iz,
	double mass_nuc
);

/****************************************************************************************/
/* hv calculates photon energy in ergs for n -> n' transitions for H and H-like ions    */
/*  In the following, we have n > n'                                                    */
/****************************************************************************************/

inline double hv(
	/* returns energy in ergs */
	/* principal quantum number, 1 for ground, upper level     */
	long int n,
	/* principal quantum number, 1 for ground, lower level     */
	long int nprime,
	long int iz,
	double mass_nuc
);

/********************************************************************************/
/*  In the following, we have n > n'                                            */
/********************************************************************************/

STATIC double fsff_log10(
	/* principal quantum number, 1 for ground, upper level       */
	long int n,
	/* angular momentum, 0 for s                                 */
	long int l,
	/* principal quantum number, 1 for ground, lower level       */
	long int np
);

/********************************************************************************/
/*   F21_mx()                                                                   */
/*   Calculates the Hyper_Spherical_Function 2_F_1(a,b,c;y)                     */
/*   Here a,b, and c are (long int)                                             */
/*   y is of type (double)                                                      */
/*   Fa_plus holds 2_F_1(a+2,b,c;y) and is of type (mx&)                        */
/*   ab specifies whether recursion is over a or b, and is of type (char)       */
/********************************************************************************/

STATIC mx F21_mx(
	long int a,
	long int b,
	long int c,
	double y,
	mx& Fa_plus,
	char ab
);

STATIC mx F21i_log(
	long int a,
	long int b,
	long int c,
	double y,
	mx& Fap2
);

/**\verbatim
  This routine, hri_log10(), calculates the hydrogen radial integral,  
  for the transition n,l --> n',l'                                      
  It is, of course, dimensionless.                                      

  In the following, we have n > n'   \endverbatim
  \param n principal quantum number, 1 for ground, upper level
  \param l angular momentum, 0 for s
  \param np principal quantum number, 1 for ground, lower level
  \param lp angular momentum, 0 for s
  \param iz Nuclear charge, 1 for H+, 2 for He++, etc                             
  \param mass_nuc nuclear mass, in g
*/
inline double hri_log10(
	long int n,
	long int l,
	long int np,
	long int lp,
	long int iz,
	double mass_nuc
);

/******************************************************************************/
/*  In the following, we have n > n'                                          */
/******************************************************************************/
STATIC double hrii_log(
	/* principal quantum number, 1 for ground, upper level       */
	long int n,
	/* angular momentum, 0 for s                                 */
	long int l,
	/* principal quantum number, 1 for ground, lower level       */
	long int np,
	/* angular momentum, 0 for s                                 */
	long int lp
);

STATIC double  bh_log(
	double k,
	long int n,
	long int l,
	mxq *rcsvV_mxq
);

STATIC double bhintegrand_log(
	double k,
	long int n,
	long int l,
	long int lp,
	mxq *rcsvV_mxq
);

STATIC mx bhG_mx(
	double K,
	long int n,
	long int l,
	long int lp,
	mxq *rcsvV_mxq
);

STATIC mx bhGp_mx(
	long int q,
	double K,
	long int n,
	long int l,
	long int lp,
	mxq *rcsvV_mxq,
	const mx& GK_mx
);

STATIC mx bhGm_mx(
	long int q,
	double K,
	long int n,
	long int l,
	long int lp,
	mxq *rcsvV_mxq,
	const mx& GK_mx
);

STATIC double bhg_log(
	double K,
	long int n,
	long int l,
	long int lp,
	mxq *rcsvV_mxq
);

inline void normalize_mx( mx& target );
inline mx add_mx( const mx& a, const mx& b );
inline mx sub_mx( const mx& a, const mx& b );
inline mx mxify( double a );
inline double unmxify( const mx& a_mx );
inline mx mxify_log10( double log10_a );
inline mx mult_mx( const mx& a, const mx& b );
inline double prodxx_log10( long int lp, double Ksqrd );

/****************************************************************************************/
/*  64 pi^4 (e a_o)^2    64 pi^4 (a_o)^2    e^2     1                      1            */
/*  ----------------- = ----------------- -------- ----  = 7.5197711e-38 -----          */
/*     3 h c^3                3  c^2       hbar c  2 pi                   sec           */
/****************************************************************************************/

static const double CONST_ONE = 32.*pow3(PI)*pow2(BOHR_RADIUS_CM)*FINE_STRUCTURE/(3.*pow2(SPEEDLIGHT));

/************************************************************************/
/* (4.0/3.0) * PI * FINE_STRUCTURE_CONSTANT * BOHRRADIUS * BOHRRADIUS   */
/*                                                                      */
/*                                                                      */
/*      4 PI alpha a_o^2                                                */
/*      ----------------                                                */
/*             3                                                        */
/*                                                                      */
/*      where   alpha = Fine Structure Constant                         */
/*              a_o   = Bohr Radius                                     */
/*                                                                      */
/*      = 3.056708^-02 (au Length)^2                                    */
/*      = 8.56x10^-23 (meters)^2                                        */
/*      = 8.56x10^-19 (cm)^2                                            */
/*      = 8.56x10^+05 (barns)                                           */
/*      = 0.856 (MB or megabarns)                                       */
/*                                                                      */
/*                                                                      */
/*      1 barn = 10^-28 (meter)^2                                       */
/************************************************************************/

static const double PHYSICAL_CONSTANT_TWO = 4./3.*PI*FINE_STRUCTURE*pow2(BOHR_RADIUS_CM);

/************************Start of program***************************/

double H_photo_cs(
	/* incident photon energy relative to threshold       */
	double rel_photon_energy,
	/* principal quantum number, 1 for ground             */
	long int n,
	/* angular momentum, 0 for s                          */
	long int l,
	/* charge, 1 for H+, 2 for He++, etc                  */
	long int iz )
{
	DEBUG_ENTRY( "H_photo_cs()" );

	double result;
	if( n <= long(t_hydro_tbl::Inst().nmaxcs()) )
	{
		result = t_hydro_tbl::Inst().cs(rel_photon_energy, n, l, iz);
	}
	else
	{
		result = H_photo_cs_log10( rel_photon_energy , n , l , iz );
	}
	return result;
}

/*****************************************************************************/
/*H_photo_cs_log10 returns hydrogenic photoionization cross section in cm-2
 * this routine is called by H_photo_cs when n is large */
/*****************************************************************************/
double H_photo_cs_log10(
	/* photon energy relative to threshold energy         */
	double rel_photon_energy,
	/* principal quantum number, 1 for ground            */
	long int n,
	/* angular momentum, 0 for s                         */
	long int l,
	/* charge, 1 for H+, 2 for He++, etc                 */
	long int iz
)
{
	DEBUG_ENTRY( "H_photo_cs_log10()" );

	long int dim_rcsvV_mxq;

	double electron_energy;
	double t1;
	double result = 0.;
	double xn_sqrd = (double)(n*n);
	double z_sqrd = (double)(iz*iz);
	double Z = (double)iz;
	double K = 0.;   /* K = k / Z                            */
	double k = 0.;   /* k^2 = ejected-electron-energy (Ryd)  */

	/* expressions blow up at precisely threshold */
	if( rel_photon_energy < 1.+FLT_EPSILON )
	{
		/* below or very close to threshold, return zero */
		fprintf( ioQQQ,"PROBLEM IN HYDRO_BAUMAN: rel_photon_energy, n, l, iz: %e\t%li\t%li\t%li\n",
		         rel_photon_energy,
		         n,
		         l,
		         iz );
		cdEXIT(EXIT_FAILURE);
	}
	if( n < 1 || l >= n )
	{
		fprintf(ioQQQ," The quantum numbers are impossible.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* k^2 is the ejected photoelectron energy in ryd */
	/*electron_energy = SDIV( (photon_energy/ryd) - (z_sqrd/xn_sqrd) );*/
	electron_energy = (rel_photon_energy-1.) * (z_sqrd/xn_sqrd);

	k = sqrt( ( electron_energy ) );
	/* k^2 is the ejected photoelectron energy in ryd */
	/*k = sqrt( ( (photon_energy/ryd) - (z_sqrd/xn_sqrd) ) );*/

	K = (k/Z);

	dim_rcsvV_mxq = (((n * 2) - 1) + 1);

	/* create space */
	vector<mxq> rcsvV_mxq(dim_rcsvV_mxq);

	t1 = bh_log( K, n, l, get_ptr(rcsvV_mxq) );

	ASSERT( t1 > 0. );

	t1 = MAX2( t1, 1.0e-250 );

	result = PHYSICAL_CONSTANT_TWO * (xn_sqrd/z_sqrd) * t1;

	if( result <= 0. )
	{
		fprintf( ioQQQ, "PROBLEM: hydro_bauman...t1\t%e\n", t1 );
	}
	ASSERT( result > 0. );
	return result;
}

STATIC double bh_log(
	/* K = k / Z ::: k^2 = ejected-electron-energy (Ryd)    */
	double K,
	/* principal quantum number                             */
	long int n,
	/* angular  momentum quantum number                     */
	long int l,
	/* Temporary storage for intermediate                   */
	mxq *rcsvV_mxq
	/*   results of the recursive routine                   */
)
{
	DEBUG_ENTRY( "bh_log()" );

	long int lp = 0; /* l' */
	double sigma = 0.; /*      Sum in Brocklehurst eq. 3.13 */

	ASSERT( n > 0 );
	ASSERT( l >= 0 );
	ASSERT( n > l );

	if( l > 0 ) /* no lp=(l-1)  for  l=0 */
	{
		for( lp = l - 1; lp <= l + 1; lp = lp + 2 )
		{
			sigma += bhintegrand_log( K, n, l, lp, rcsvV_mxq );
		}
	}
	else
	{
		lp = l + 1;
		sigma = bhintegrand_log( K, n, l, lp, rcsvV_mxq );
	}
	ASSERT( sigma != 0. );
	return sigma;
}

/************************************************************************************************/
/*    The photoionization cross-section (see also Burgess 1965)                                 */
/*     is given by;                                                                             */
/*                   _               _        l'=l+1                                            */
/*                  |4 PI alpha a_o^2 |  n^2   ---     max(l,l')                                */
/*      a(Z';n,l,k)=|---------------- |  ---    >     ---------- Theta(n,l; K, l')              */
/*                  |       3         |  Z^2   ---     (2l + 1)                                 */
/*                   _               _        l'=l-1                                            */
/*                                                                                              */
/*                                                                                              */
/*      where  Theta(n,l; K, l') is defined                                                     */
/*                                                                                              */
/*                                           |              |^2                                 */
/*         Theta(n,l; K, l') = (1 + n^2 K^2) | g(n,l; K,l') |                                   */
/*                                           |              |                                   */
/*                                                                                              */
/*                                                                                              */
/*                                                 ---- Not sure if this is K or k              */
/*                                  OO            /     but think it is K                       */
/*         where                    -            v                                              */
/*                         K^2     |                                                            */
/*         g(n,l; K,l') = -----    | R_nl(r) F_(K,l) r dr                                       */
/*                         n^2     |                                                            */
/*                                -                                                             */
/*                                0                                                             */
/************************************************************************************************/

STATIC double bhintegrand_log(
	double K,     /* K = k / Z ::: k^2 = ejected-electron-energy (Ryd) */
	long int n,
	long int l,
	long int lp,
	/* Temporary storage for intermediate         */
	/* results of the recursive routine           */
	mxq *rcsvV_mxq
)
{
	DEBUG_ENTRY( "bhintegrand_log()" );

	double d2 = 0.;
	double d5 = 0.;
	double d7 = 0.;
	double Theta = 0.;
	double n2 = (double)(n*n);
	double Ksqrd = (K*K);
	double Two_L_Plus_One = (double)(2*l + 1);
	double lg = (double)max(l,lp);

	ASSERT( Ksqrd != 0. );
	ASSERT( K != 0. );
	ASSERT( lg != 0. );
	ASSERT( n2 != 0. );
	ASSERT( Two_L_Plus_One != 0. );

	ASSERT( n > 0);
	ASSERT( l >= 0);
	ASSERT( lp >= 0);

	/**********************************************/
	/*                                            */
	/*    l>                                      */
	/*  ------    Theta(nl,Kl')                   */
	/*   2l+2                                     */
	/*                                            */
	/*                                            */
	/* Theta(nl,Kl') =                            */
	/*   (1+n^2K^2) * | g(nl,Kl')|^2              */
	/*                                            */
	/**********************************************/

	d2 = ( 1. + n2 * (Ksqrd) );

	ASSERT( d2 != 0. );

	d5 = bhg_log( K, n, l, lp, rcsvV_mxq );

	d5 = MAX2( d5, 1.0E-150 );
	ASSERT( d5 != 0. );

	Theta = d2 * d5 * d5;
	ASSERT( Theta != 0. );

	d7 = (lg/Two_L_Plus_One) * Theta;

	ASSERT( d7 != 0. );
	return d7;
}

/****************************************************************************************/
/*  *** bhG ***                                                                         */
/*    Using various recursion relations                                                 */
/*    (for l'=l+1)                                                                      */
/*    equation: (3.23)                                                                  */
/*    G(n,l-2; K,l-1) = [ 4n^2-4l^2+l(2l-1)(1+(n K)^2) ] G(n,l-1; K,l)                  */
/*                      - 4n^2 (n^2-l^2)[1+(l+1)^2 K^2] G(n,l; K,l+1)                   */
/*                                                                                      */
/*    and (for l'=l-1)                                                                  */
/*    equation: (3.24)                                                                  */
/*    G(n,l-1; K,l-2) = [ 4n^2-4l^2 + l(2l-1)(1+(n K)^2) ] G(n,l; K,l-1)                */
/*                      - 4n^2 (n^2-(l+1)^2)[ 1+(lK)^2 ] G(n,l; K,l+1)                  */
/*                                                                                      */
/*    the starting point for the recursion relations are;                               */
/*    equation: (3.18)                                                                  */
/*                     | pi |(1/2)   8n                                                 */
/*     G(n,n-1; 0,n) = | -- |      ------- (4n)^n exp(-2n)                              */
/*                     | 2  |      (2n-1)!                                              */
/*                                                                                      */
/*     equation: (3.20)                                                                 */
/*                             exp(2n-2/K tan^(-1)(n K)                                 */
/*     G(n,n-1; K,n) =  -----------------------------------------  *  G(n,n-1; 0,n)     */
/*                      sqrt(1 - exp(-2 pi K)) * (1+(n K)^2)^(n+2)                      */
/*                                                                                      */
/*     equation: (3.20)                                                                 */
/*     G(n,n-2; K,n-1) = (2n-2)(1+(n K)^2) n G(n,n-1; K,n)                              */
/*                                                                                      */
/*     equation: (3.21)                                                                 */
/*                       (1+(n K)^2)                                                    */
/*     G(n,n-1; K,n-2) = ----------- G(n,n-1; K,n)                                      */
/*                           2n                                                         */
/*                                                                                      */
/*     equation: (3.22)                                                                 */
/*     G(n,n-2; K,n-3) = (2n-1) (4+(n-1)(1+n^2 K^2)) G(n,n-1; K,n-2)                    */
/****************************************************************************************/

/*************log version********************************/
STATIC mx bhG_mx(
	double K,
	long int n,
	long int l,
	long int lp,
	/* Temporary storage for intermediate        */
	/* results of the recursive routine          */
	mxq *rcsvV_mxq
)
{
	DEBUG_ENTRY( "bhG_mx()" );

	double log10_GK = 0.;
	double log10_G0 = 0.;

	double d1 = 0., d2 = 0., d3 = 0., d4 = 0., d5 = 0., d6 = 0.;
	double ld1 = 0., ld2 = 0., ld3 = 0., ld4 = 0., ld5 = 0., ld6 = 0.;

	double n1 = (double)n;
	double n2 = n1 * n1;
	double Ksqrd = K * K;

	mx GK_mx;

	/* l=l'-1 or l=l'+1 */
	ASSERT( (l == lp - 1) ||  (l == lp + 1) );
	ASSERT( K != 0. );
	ASSERT( n1 != 0. );
	ASSERT( n2 != 0. );
	ASSERT( Ksqrd != 0. );
	ASSERT( ((2*n) - 1) >= 0   );

	/******************************/
	/*                 n          */
	/*                ---         */
	/*    log( n! ) =  >  log(j)  */
	/*                ---         */
	/*                j=1         */
	/******************************/

	/*************************************************************/
	/*                     | pi |(1/2)   8n                      */
	/*     G(n,n-1; 0,n) = | -- |      ------- (4n)^n exp(-2n)   */
	/*                     | 2  |      (2n-1)!                   */
	/*************************************************************/

	/******************************/
	/*                            */
	/*                            */
	/*      log10( (2n-1)! )      */
	/*                            */
	/*                            */
	/******************************/

	ld1 = lfactorial( 2*n - 1 );
	ASSERT( ld1 >= 0. );

	/**********************************************/
	/*            (4n)^n                          */
	/**********************************************/
	/*    log10( 4n^n ) = n log10( 4n )           */
	/**********************************************/

	ld2 = n1 * log10( 4. * n1 );
	ASSERT( ld2 >= 0. );

	/**********************************************/
	/*            exp(-2n)                        */
	/**********************************************/
	/*    log10( exp( -2n ) ) = (-2n) * log10(e)  */
	/**********************************************/
	ld3 = (-(2. * n1)) * (LOG10_E);
	ASSERT( ld3 <= 0. );

	/******************************************************************************/
	/*    ********G0*******                                                       */
	/*                                                                            */
	/*            | pi |(1/2)   8n                                                */
	/*      G0 =  | -- |      ------- (4n)^n exp(-2n)                             */
	/*            | 2  |      (2n-1)!                                             */
	/******************************************************************************/

	log10_G0 = log10(SQRTPIBY2 * 8. * n1) + ( (ld2 + ld3) - ld1);

	/******************************************************************************/
	/*    *****GK*****                                                            */
	/*                                                                            */
	/*                            exp(2n- (2/K) tan^(-1)(n K)  )                  */
	/*      G(n,n-1; K,n) =  ----------------------------------------- *  G0      */
	/*                       sqrt(1 - exp(-2 pi/ K)) * (1+(n K))^(n+2)            */
	/******************************************************************************/

	ASSERT( K != 0. );

	/**********************************************/
	/*  sqrt(1 - exp(-2 pi/ K))                   */
	/**********************************************/
	/*    log10(sqrt(1 - exp(-2 pi/ K))) =        */
	/*            (1/2) log10(1 - exp(-2 pi/ K))  */
	/**********************************************/

	d1 = (1. - exp(-(2. *  PI )/ K ));
	ld4 = (1./2.) * log10( d1 );
	ASSERT( K != 0. );
	ASSERT( d1 != 0. );

	/**************************************/
	/*         (1+(n K)^2)^(n+2)          */
	/**************************************/
	/* log10( (1+(n K)^2)^(n+2) ) =       */
	/*    (n+2) log10( (1 + (n K)^2 ) )   */
	/**************************************/

	d2 = ( 1. + (n2 * Ksqrd));
	ld5 = (n1 + 2.) * log10( d2 );
	ASSERT( d2 != 0. );

	ASSERT( ld5 >= 0. );

	/**********************************************/
	/*   exp(2n- (2/K)*tan^(-1)(n K)  )           */
	/**********************************************/
	/* log10( exp(2n- (2/K) tan^(-1)(n K)  ) =    */
	/*     (2n- (2/K)*tan^(-1)(n K) ) * Log10(e)  */
	/**********************************************/

	/*    tan^(-1)(n K) )         */
	d3 = atan( n1 * K );
	ASSERT( d3 != 0. );

	/*    (2/K)*tan^(-1)(n K) )   */
	d4 = (2. / K) * d3;
	ASSERT( d4 != 0. );

	/*            2n              */
	d5 = (double) ( 2. * n1 );
	ASSERT( d5 != 0. );

	/*  (2n-2/K tan^(-1)(n K))    */
	d6 = d5 - d4;
	ASSERT( d6 != 0. );

	/* log10( exp(2n- (2/K) tan^(-1)(n K)  )      */
	ld6 = LOG10_E * d6;
	ASSERT( ld6 != 0. );

	/******************************************************************************/
	/*    *****GK*****                                                            */
	/*                                                                            */
	/*                            exp(2n- (2/K) tan^(-1)(n K)  )                  */
	/*      G(n,n-1; K,n) =  ----------------------------------------- *  G0      */
	/*                       sqrt(1 - exp(-2 pi/ K)) * (1+(n K))^(n+2)            */
	/******************************************************************************/

	log10_GK = (ld6 -(ld4 + ld5)) + log10_G0;
	ASSERT( log10_GK != 0. );

	GK_mx = mxify_log10( log10_GK );

	/*  GENERAL CASE: l = l'-1  */
	if( l == lp - 1 )
	{
		mx result_mx = bhGm_mx( l, K, n, l, lp, rcsvV_mxq , GK_mx );
		/* Here the m in bhGm() refers           */
		/* to the minus sign(-) in l=l'-1        */
		return result_mx;
	}
	/*  GENERAL CASE: l = l'+1  */
	else if( l == lp + 1 )
	{
		mx result_mx = bhGp_mx( l, K, n, l, lp, rcsvV_mxq , GK_mx );
		/* Here the p in bhGp() refers           */
		/* to the plus sign(+) in l=l'+1         */
		return result_mx;
	}
	else
	{
		printf( "BadMagic: l and l' do NOT satisfy dipole requirements.\n\n" );
		cdEXIT(EXIT_FAILURE);
	}
	printf( "This code should be inaccessible\n\n" );
	cdEXIT(EXIT_FAILURE);
}

/************************************************************************************************/
/*  ***   bhGp.c   ***                                                                          */
/*                                                                                              */
/*    Here we calculate G(n,l; K,l') with the  recursive formula                                */
/*    equation: (3.24)                                                                          */
/*                                                                                              */
/*    G(n,l-1; K,l-2) = [ 4n^2-4l^2 + l(2l+1)(1+(n K)^2) ] G(n,l; K,l-1)                        */
/*                                                                                              */
/*                      - 4n^2 (n^2-(l+1)^2)[ 1+(lK)^2 ] G(n,l+1; K,l)                          */
/*                                                                                              */
/*    Under the transformation l -> l + 1 this gives                                            */
/*                                                                                              */
/*    G(n,l+1-1; K,l+1-2) = [ 4n^2-4(l+1)^2 + (l+1)(2(l+1)+1)(1+(n K)^2) ] G(n,l+1; K,l+1-1)    */
/*                                                                                              */
/*                      - 4n^2 (n^2-((l+1)+1)^2)[ 1+((l+1)K)^2 ] G(n,l+1+1; K,l+1)              */
/*                                                                                              */
/*    or                                                                                        */
/*                                                                                              */
/*     G(n,l; K,l-1) = [ 4n^2-4(l+1)^2 + (l+1)(2l+3)(1+(n K)^2) ] G(n,l+1; K,l)                 */
/*                                                                                              */
/*                      - 4n^2 (n^2-(l+2)^2)[ 1+((l+1)K)^2 ] G(n,l+2; K,l+1)                    */
/*                                                                                              */
/*    from the reference                                                                        */
/*    M. Brocklehurst                                                                           */
/*    Mon. Note. R. astr. Soc. (1971) 153, 471-490                                              */
/*                                                                                              */
/*                                                                                              */
/*  * This is valid for the case l=l'+1  *                                                      */
/*  * CASE:     l = l'+1                 *                                                      */
/*  * Here the p in bhGp() refers        *                                                      */
/*  * to the Plus sign(+) in l=l'+1      *                                                      */
/************************************************************************************************/

/***********************log version*******************************/
STATIC mx bhGp_mx(
	long int q,
	double K,
	long int n,
	long int l,
	long int lp,
	/* Temporary storage for intermediate       */
	/* results of the recursive routine         */
	mxq *rcsvV_mxq,
	const mx& GK_mx
)
{
	DEBUG_ENTRY( "bhGp_mx()" );

	/* static long int rcsv_Level = 1;                                    */
	/*    printf( "bhGp(): recursion level:\t%li\n",rcsv_Level++ );       */

	ASSERT( l == lp + 1 );

	long int rindx = 2*q;

	if( rcsvV_mxq[rindx].q == 0 )
	{
		/*  SPECIAL CASE:   n = l+1 = l'+2 */
		if( q == n - 1 )
		{
			/******************************************************/
			/*                    (1+(n K)^2)                     */
			/*  G(n,n-1; K,n-2) = ----------- G(n,n-1; K,n)       */
			/*                        2n                          */
			/******************************************************/

			double Ksqrd = (K * K);
			double n2 = (double)(n*n);

			double dd1 = (double) (2 * n);
			double dd2 = ( 1. + ( n2 * Ksqrd));
			double dd3 = dd2/dd1;

			mx dd3_mx = mxify( dd3 );
			mx G1_mx = mult_mx( dd3_mx, GK_mx);

			normalize_mx( G1_mx );

			ASSERT( l == lp + 1 );
			ASSERT( Ksqrd != 0. );
			ASSERT( n2 != 0. );
			ASSERT( dd1 != 0. );
			ASSERT( dd2 != 0. );

			rcsvV_mxq[rindx].q = 1;
			rcsvV_mxq[rindx].mx = G1_mx;
			return G1_mx;
		}
		/*  SPECIAL CASE:   n = l+2 = l'+3 */
		else if( q == (n - 2) )
		{
			/****************************************************************/
			/*                                                              */
			/* G(n,n-2; K,n-3) = (2n-1) (4+(n-1)(1+(n K)^2)) G(n,n-1; K,n-2)*/
			/*                                                              */
			/****************************************************************/
			/*                    (1+(n K)^2)                               */
			/*  G(n,n-1; K,n-2) = ----------- G(n,n-1; K,n)                 */
			/*                        2n                                    */
			/****************************************************************/

			double Ksqrd = (K*K);
			double n2 = (double)(n*n);
			double dd1 = (double)(2 * n);
			double dd2 = ( 1. + ( n2 * Ksqrd) );
			double dd3 = (dd2/dd1);
			double dd4 = (double) ((2 * n) - 1);
			double dd5 = (double) (n - 1);
			double dd6 = (4. + (dd5 * dd2));
			double dd7 = dd4 * dd6;

			/****************************************************************/
			/*                                                              */
			/* G(n,n-2; K,n-3) = (2n-1) (4+(n-1)(1+(n K)^2)) G(n,n-1; K,n-2)*/
			/*                                                              */
			/****************************************************************/

			mx dd3_mx = mxify( dd3 );
			mx dd7_mx = mxify( dd7 );
			mx G1_mx = mult_mx( dd3_mx, GK_mx );
			mx G2_mx = mult_mx( dd7_mx, G1_mx );

			normalize_mx( G2_mx );

			ASSERT( l == lp + 1 );
			ASSERT( Ksqrd != 0. );
			ASSERT( n2 != 0. );
			ASSERT( dd1 != 0. );
			ASSERT( dd2 != 0. );
			ASSERT( dd3 != 0. );
			ASSERT( dd4 != 0. );
			ASSERT( dd5 != 0. );
			ASSERT( dd6 != 0. );
			ASSERT( dd7 != 0. );

			rcsvV_mxq[rindx].q = 1;
			rcsvV_mxq[rindx].mx = G2_mx;
			return G2_mx;
		}
		/* The GENERAL CASE n > l + 2*/
		else
		{
			/**************************************************************************************/
			/*  G(n,l; K,l-1) = [ 4n^2-4(l+1)^2 + (l+1)(2l+3)(1+(n K)^2) ] G(n,l+1; K,l)          */
			/*                                                                                    */
			/*                      - 4n^2 (n^2-(l+2)^2)[ 1+((l+1)K)^2 ] G(n,l+2; K,l+1)          */
			/*                                                                                    */
			/*            FROM   Eq. 3.24                                                         */
			/*                                                                                    */
			/*  G(n,l-1; K,l-2) = [ 4n^2-4l+^2 + l(2l+1)(1+(n K)^2) ] G(n,l; K,l-1)               */
			/*                                                                                    */
			/*                      - 4n^2 (n^2-(l+1)^2)[ 1+((lK)^2 ] G(n,l+1; K,l)               */
			/**************************************************************************************/

			long int lp1 = (q + 1);
			long int lp2 = (q + 2);

			double Ksqrd = (K * K);
			double n2 = (double)(n * n);
			double lp1s = (double)(lp1 * lp1);
			double lp2s = (double)(lp2 * lp2);

			double d1 = (4. * n2);
			double d2 = (4. * lp1s);
			double d3 = (double)((lp1)*((2 * q) + 3));
			double d4 = (1. + (n2 * Ksqrd));
			/* [ 4n^2 - 4(l+1)^2 + (l+1)(2l+3)(1+(n K)^2) ]       */
			double d5 = d1 - d2 + (d3 * d4);

			/* (n^2-(l+2)^2)      */
			double d6 = (n2 - lp2s);

			/* [ 1+((l+1)K)^2 ]   */
			double d7 = (1. + (lp1s * Ksqrd));

			/* { 4n^2 (n^2-(l+1)^2)[ 1+((l+1) K)^2 ] }    */
			double d8 = (d1 * d6 * d7);

			/**************************************************************************************/
			/*  G(n,l; K,l-1) = [ 4n^2 - 4(l+1)^2 + (l+1)(2l+3)(1+(n K)^2) ] G(n,l+1; K,l)        */
			/*                                                                                    */
			/*                      - 4n^2 (n^2-(l+2)^2)[ 1+((l+1)K)^2 ] G(n,l+2; K,l+1)          */
			/**************************************************************************************/

			mx d5_mx=mxify( d5 );
			mx d8_mx=mxify( d8 );

			mx t0_mx = bhGp_mx( (q+1), K, n, l, lp, rcsvV_mxq, GK_mx );
			mx t1_mx = bhGp_mx( (q+2), K, n, l, lp, rcsvV_mxq, GK_mx );

			mx d9_mx = mult_mx( d5_mx, t0_mx );
			mx d10_mx = mult_mx( d8_mx, t1_mx );

			mx result_mx = sub_mx( d9_mx, d10_mx );
			normalize_mx( result_mx );

			ASSERT( d1 != 0. );
			ASSERT( d2 != 0. );
			ASSERT( d3 != 0. );
			ASSERT( d4 != 0. );
			ASSERT( d5 != 0. );
			ASSERT( d6 != 0. );
			ASSERT( d7 != 0. );
			ASSERT( d8 != 0. );

			ASSERT( l == lp + 1 );
			ASSERT( Ksqrd != 0. );
			ASSERT( n2 != 0. );
			ASSERT( lp1s != 0. );
			ASSERT( lp2s != 0. );

			rcsvV_mxq[rindx].q = 1;
			rcsvV_mxq[rindx].mx = result_mx;
			return result_mx;
		}
	}
	else
	{
		ASSERT( rcsvV_mxq[rindx].q != 0 );
		rcsvV_mxq[rindx].q = 1;
		return rcsvV_mxq[rindx].mx;
	}
}

/************************************************************************************************/
/*  ***   bhGm.c  *** */
/*                                                                                              */
/*    Here we calculate G(n,l; K,l') with the  recursive formula                                */
/*    equation: (3.23)                                                                          */
/*                                                                                              */
/*    G(n,l-2; K,l-1) = [ 4n^2-4l^2 + l(2l-1)(1+(n K)^2) ] G(n,l-1; K,l)                        */
/*                                                                                              */
/*                      - 4n^2 (n^2-l^2)[ 1 + (l+1)^2 K^2 ] G(n,l; K,l+1)                       */
/*                                                                                              */
/*    Under the transformation l -> l + 2 this gives                                            */
/*                                                                                              */
/*    G(n,l+2-2; K,l+2-1) = [ 4n^2-4(l+2)^2 + (l+2)(2(l+2)-1)(1+(n K)^2) ] G(n,l+2-1; K,l+2)    */
/*                                                                                              */
/*                      - 4n^2 (n^2-(l+2)^2)[ 1 + (l+2+1)^2 K^2 ] G(n,l+2; K,l+2+1)             */
/*                                                                                              */
/*                                                                                              */
/*    or                                                                                        */
/*                                                                                              */
/*    G(n,l; K,l+1) = [ 4n^2-4(l+2)^2 + (l+2)(2l+3)(1+(n K)^2) ] G(n,l+1; K,l+2)                */
/*                                                                                              */
/*                      - 4n^2 (n^2-(l+2)^2)[ 1 + (l+3)^2 K^2 ] G(n,l+2; K,l+3)                 */
/*                                                                                              */
/*                                                                                              */
/*    from the reference                                                                        */
/*    M. Brocklehurst                                                                           */
/*    Mon. Note. R. astr. Soc. (1971) 153, 471-490                                              */
/*                                                                                              */
/*                                                                                              */
/*  * This is valid for the case l=l'-1 *                                                       */
/*  * CASE:     l = l'-1                *                                                       */
/*  * Here the p in bhGm() refers       *                                                       */
/*  * to the Minus sign(-) in l=l'-1    *                                                       */
/************************************************************************************************/

/************************log version***********************************/
STATIC mx bhGm_mx(
	long int q,
	double K,
	long int n,
	long int l,
	long int lp,
	mxq *rcsvV_mxq,
	const mx& GK_mx
)
{
	DEBUG_ENTRY( "bhGm_mx()" );

	/*static long int rcsv_Level = 1;                                     */
	/*printf( "bhGm(): recursion level:\t%li\n",rcsv_Level++ );           */

	ASSERT( l == lp - 1 );
	ASSERT( l < n );

	long int rindx = 2*q+1;

	if( rcsvV_mxq[rindx].q == 0 )
	{
		/*  CASE:     l = n - 1      */
		if( q == n - 1 )
		{
			mx result_mx = GK_mx;
			normalize_mx( result_mx );

			rcsvV_mxq[rindx].q = 1;
			rcsvV_mxq[rindx].mx = result_mx;

			ASSERT(l == lp - 1);
			return result_mx;
		}
		/*  CASE:     l = n - 2      */
		else if( q == n - 2 )
		{
			double Ksqrd = (K * K);
			double n1 = (double)n;
			double n2 = (double) (n*n);
			double dd1 = (double)((2 * n) - 1);
			double dd2 = (1. + (n2 * Ksqrd));
			/*(2n-1)(1+(n K)^2) n*/
			double dd3 = (dd1*dd2*n1); 

			/******************************************************/
			/*     G(n,n-2; K,n-1) =                              */
			/*            (2n-1)(1+(n K)^2) n G(n,n-1; K,n)       */
			/******************************************************/

			mx dd3_mx = mxify( dd3 );
			mx G2_mx = mult_mx( dd3_mx, GK_mx );

			normalize_mx( G2_mx );

			ASSERT( l == lp - 1);
			ASSERT( n1 != 0. );
			ASSERT( n2 != 0. );
			ASSERT( dd1 != 0. );
			ASSERT( dd2 != 0. );
			ASSERT( dd3 != 0. );
			ASSERT( Ksqrd != 0. );

			rcsvV_mxq[rindx].q = 1;
			rcsvV_mxq[rindx].mx = G2_mx;
			return G2_mx;
		}
		else
		{
			/******************************************************************************/
			/* G(n,l; K,l+1) = [ 4n^2-4(l+2)^2 + (l+2)(2l+3)(1+(n K)^2) ] G(n,l+1; K,l+2) */
			/*                                                                            */
			/*                    - 4n^2 (n^2-(l+2)^2)[ 1+((l+3)^2 K^2) ] G(n,l+2; K,l+3) */
			/*                                                                            */
			/*                                                                            */
			/*            FROM   Eq. 3.23                                                 */
			/*                                                                            */
			/* G(n,l-2; K,l-1) = [ 4n^2-4l^2 + (l+2)(2l-1)(1+(n K)^2) ] G(n,l-1; K,l)     */
			/*                                                                            */
			/*                   - 4n^2 (n^2-l^2)[ 1 + (l+1)^2 K^2 ] G(n,l; K,l+1)        */
			/******************************************************************************/

			long int lp2 = (q + 2);
			long int lp3 = (q + 3);

			double lp2s = (double)(lp2 * lp2);
			double lp3s = (double)(lp3 * lp3);
			double n2 = (double)(n*n);
			double Ksqrd = (K * K);

			/******************************************************/
			/*    [ 4n^2-4(l+2)^2 + (l+2)(2l+3)(1+(n K)^2) ]      */
			/******************************************************/

			double d1 = (4. * n2);
			double d2 = (4. * lp2s);
			double d3 = (double)(lp2)*((2*q)+3);
			double d4 = (1. + (n2 * Ksqrd));
			/* 4n^2-4(l+2)^2 + (l+2)(2l+3)(1+(n K)^2)           */
			double d5 = d1 - d2 + (d3 * d4);

			mx d5_mx=mxify(d5);

			/******************************************************/
			/*        4n^2 (n^2-(l+2)^2)[ 1+((l+3)^2 K^2) ]       */
			/******************************************************/

			double d6 = (n2 - lp2s);
			double d7 = (1. + (lp3s * Ksqrd));
			/* 4n^2 (n^2-(l+2)^2)[ 1 + (l+3)^2 K^2 ]            */
			double d8 = d1 * d6 * d7;

			mx d8_mx = mxify(d8);

			/******************************************************************************/
			/* G(n,l; K,l+1) = [ 4n^2-4(l+2)^2 + (l+2)(2l+3)(1+(n K)^2) ] G(n,l+1; K,l+2) */
			/*                                                                            */
			/*                   - 4n^2 (n^2-(l+2)^2)[ 1 + (l+3)^2 K^2 ] G(n,l+2; K,l+3)  */
			/******************************************************************************/

			mx d9_mx = bhGm_mx( (q+1), K, n, l, lp, rcsvV_mxq, GK_mx );
			mx d10_mx = bhGm_mx( (q+2), K, n, l, lp, rcsvV_mxq, GK_mx );
			mx d11_mx = mult_mx( d5_mx, d9_mx );
			mx d12_mx = mult_mx( d8_mx, d10_mx);
			mx result_mx = sub_mx( d11_mx , d12_mx );
			rcsvV_mxq[rindx].q = 1;
			rcsvV_mxq[rindx].mx = result_mx;

			ASSERT(l == lp - 1);
			ASSERT(n2 != 0. );
			ASSERT(lp2s != 0.);
			ASSERT( lp3s != 0.);
			ASSERT(Ksqrd != 0.);

			ASSERT(d1 != 0.);
			ASSERT(d2 != 0.);
			ASSERT(d3 != 0.);
			ASSERT(d4 != 0.);
			ASSERT(d5 != 0.);
			ASSERT(d6 != 0.);
			ASSERT(d7 != 0.);
			ASSERT(d8 != 0.);
			return result_mx;
		}
	}
	else
	{
		ASSERT(  rcsvV_mxq[rindx].q != 0 );
		return rcsvV_mxq[rindx].mx;
	}
}

/****************************************************************************************/
/*                                                                                      */
/*      bhg.c                                                                           */
/*                                                                                      */
/*      From reference;                                                                 */
/*      M. Brocklehurst                                                                 */
/*      Mon. Note. R. astr. Soc. (1971) 153, 471-490                                    */
/*                                                                                      */
/*                                                                                      */
/*      We wish to compute the following function,                                      */
/*                                                                                      */
/*    equation: (3.17)                                                                  */
/*                 -           s=l'              - (1/2)                                */
/*                |  (n+l)!   -----               |                                     */
/*    g(nl, Kl) = | --------   | |  (1 + s^2 K^2) |    *  (2n)^(l-n) G(n,l; K,l')       */
/*                | (n-l-1)!   | |                |                                     */
/*                 -           s=0               -                                      */
/*                                                                                      */
/*    Using various recursion relations (for l'=l+1)                                    */
/*                                                                                      */
/*    equation: (3.23)                                                                  */
/*    G(n,l-2; K,l-1) = [ 4n^2-4l^2+l(2l-1)(1+(n K)^2) ] G(n,l-1; K,l)                  */
/*                                                                                      */
/*                      - 4n^2 (n^2-l^2)[1+(l+1)^2 K^2] G(n,l; K,l+1)                   */
/*                                                                                      */
/*    and (for l'=l-1)                                                                  */
/*                                                                                      */
/*    equation: (3.24)                                                                  */
/*    G(n,l-1; K,l-2) = [ 4n^2-4l^2 + l(2l-1)(1+(n K)^2) ] G(n,l; K,l-1)                */
/*                                                                                      */
/*                      - 4n^2 (n^2-(l+1)^2)[ 1+(lK)^2 ] G(n,l; K,l+1)                  */
/*                                                                                      */
/*                                                                                      */
/*    the starting point for the recursion relations are:                               */
/*                                                                                      */
/*                                                                                      */
/*    equation (3.18):                                                                  */
/*                                                                                      */
/*                     | pi |(1/2)   8n                                                 */
/*     G(n,n-1; 0,n) = | -- |      ------- (4n)^2 exp(-2n)                              */
/*                     | 2  |      (2n-1)!                                              */
/*                                                                                      */
/*     equation (3.20):                                                                 */
/*                                                                                      */
/*                          exp(2n-2/K tan^(-1)(n K)                                    */
/*     G(n,n-1; K,n) =  ---------------------------------------                         */
/*                      sqrt(1 - exp(-2 pi/ K)) * (1+(n K)^(n+2)                        */
/*                                                                                      */
/*                                                                                      */
/*                                                                                      */
/*     equation (3.20):                                                                 */
/*     G(n,n-2; K,n-1) = (2n-2)(1+(n K)^2) n G(n,n-1; K,n)                              */
/*                                                                                      */
/*                                                                                      */
/*     equation (3.21):                                                                 */
/*                                                                                      */
/*                       (1+(n K)^2)                                                    */
/*     G(n,n-1; K,n-2) = ----------- G(n,n-1; K,n)                                      */
/*                           2n                                                         */
/****************************************************************************************/

/********************log version**************************/
STATIC double bhg_log(
	double K,
	long int n,
	long int l,
	long int lp,
	/* Temporary storage for intermediate      */
	/*   results of the recursive routine      */
	mxq *rcsvV_mxq
)
{
	/**************************************************************************************/
	/*    equation: (3.17)                                                                */
	/*                 -           s=l'              - (1/2)                              */
	/*                |  (n+l)!   -----               |                                   */
	/*    g(nl, Kl) = | --------   | |  (1 + s^2 K^2) |    *  (2n)^(l-n) G(n,l; K,l')     */
	/*                | (n-l-1)!   | |                |                                   */
	/*                 -           s=0               -                                    */
	/**************************************************************************************/

	DEBUG_ENTRY( "bhg_log()" );

	double d1 = (double)(2*n);
	double d2 = (double)(l-n);
	double Ksqrd = (K*K);

	/**************************************************************************************/
	/*                                                                                    */
	/*          |  (n+l)!  |                                                              */
	/*    log10 | -------- | = log10((n+1)!) - log10((n-l-1)!)                            */
	/*          | (n-l-1)! |                                                              */
	/*                                                                                    */
	/**************************************************************************************/

	double ld1 = lfactorial( n + l );
	double ld2 = lfactorial( n - l - 1 );

	/**********************************************************************/
	/*        |  s=l'                |     s=l'                           */
	/*        | -----                |     ---                            */
	/*  log10 |  | |  (1 + s^2 K^2)  | =    >    log10((1 + s^2 K^2))     */
	/*        |  | |                 |     ---                            */
	/*        |  s=0                 |     s=0                            */
	/**********************************************************************/

	double ld3 = prodxx_log10( lp, Ksqrd );

	/**********************************************/
	/*      -           s=l'              - (1/2) */
	/*     |  (n+l)!   -----               |      */
	/*     | --------   | |  (1 + s^2 K^2) |      */
	/*     | (n-l-1)!   | |                |      */
	/*      -           s=0               -       */
	/**********************************************/

	/***********************************************************************/
	/*                                                                     */
	/*            |  -           s=l'              - (1/2) |               */
	/*            | |  (n+l)!   -----               |      |               */
	/*       log10| | --------   | |  (1 + s^2 K^2) |      | ==            */
	/*            | | (n-l-1)!   | |                |      |               */
	/*            |  -           s=0               -       |               */
	/*                                                                     */
	/*              |                           |  s=l'               |  | */
	/*              |      |  (n+l)!  |         | -----               |  | */
	/*       (1/2)* |log10 | -------- | + log10 |  | |  (1 + s^2 K^2) |  | */
	/*              |      | (n-l-1)! |         |  | |                |  | */
	/*              |                           |  s=0                |  | */
	/*                                                                     */
	/***********************************************************************/

	double ld4 = (1./2.) * ( ld3 + ld1 - ld2 );

	/**********************************************/
	/*                   (2n)^(l-n)               */
	/**********************************************/
	/*    log10( 2n^(L-n) ) = (L-n) log10( 2n )   */
	/**********************************************/

	double ld5 = d2 * log10( d1 );

	/**************************************************************************************/
	/*    equation: (3.17)                                                                */
	/*                 -           s=l'              - (1/2)                              */
	/*                |  (n+l)!   -----               |                                   */
	/*    g(nl, Kl) = | --------   | |  (1 + s^2 K^2) |    *  (2n)^(l-n) * G(n,l; K,l')   */
	/*                | (n-l-1)!   | |                |                                   */
	/*                 -           s=0               -                                    */
	/**************************************************************************************/

	/****************************************************/
	/*                                                  */
	/*  -           s=l'              - (1/2)           */
	/* |  (n+l)!   -----               |                */
	/* | --------   | |  (1 + s^2 K^2) |  * (2n)^(L-n)  */
	/* | (n-l-1)!   | |                |                */
	/*  -           s=0               -                 */
	/****************************************************/

	double ld6 = (ld5+ld4);

	mx d6_mx = mxify_log10( ld6 );
	mx dd1_mx = bhG_mx( K, n, l, lp, rcsvV_mxq );
	mx dd2_mx = mult_mx( d6_mx, dd1_mx );
	normalize_mx( dd2_mx );
	double result = unmxify( dd2_mx );

	// This routine can return zero for very high angular momentum at high energies.
	// K = 1 corresponds to photon energy equal to ground state IP.
	// This was the old upper limit for excited states before changes to force 
	// consistency between recombination coefficients, RRC, and cooling coefficients.  
	// Those changes required higher energies in opacity stacks than the ground state IP.
	ASSERT( result > 0. || (result==0. && l > 50 && K > 1.) );

	ASSERT( Ksqrd != 0. );
	ASSERT( ld3 >= 0. );

	ASSERT( d1 > 0. );
	ASSERT( d2 < 0. );
	return result;
}

/********************************************************************************/
/*         m_e m_n       m_e                                                    */
/*    u = --------- = -----------                                               */
/*        m_e + m_n   1 + m_e/m_n                                               */
/*                                                                              */
/*        m_e                                                                   */
/*  Now  ----- = 0.000544617013                                                 */
/*        m_p                                                                   */
/*            u                                                                 */
/*  so that  --- =  0.999455679                                                 */
/*           m_e                                                                */
/*                                                                              */
/*  for the hydrogen atom                                                       */
/********************************************************************************/
STATIC inline double reduced_mass_rel( double mass_nuc )
{
	DEBUG_ENTRY( "reduced_mass_rel()" );

	ASSERT( mass_nuc > 0. );
	return 1./(1. + ELECTRON_MASS/mass_nuc); // u/m_e
}

/************************************************************************/
/*  Find the Einstein A's for hydrogen for a                            */
/*  transition n,l -->  n',l'                                           */
/*                                                                      */
/*  In the following, we will assume n > n'                             */
/************************************************************************/
/*******************************************************************************/
/*                                                                             */
/*   Einstein A() for the transition from the                                  */
/*   initial state n,l to the finial state n',l'                               */
/*   is given by oscillator f()                                                */
/*                                                                             */
/*                     hbar w    max(l,l')  |              | 2                 */
/*   f(n,l;n',l') = - -------- ------------ | R(n,l;n',l') |                   */
/*                     3 R_oo   ( 2l + 1 )  |              |                   */
/*                                                                             */
/*                                                                             */
/*                     E(n,l;n',l')     max(l,l')  |              | 2          */
/*   f(n,l;n',l') = -  ------------   ------------ | R(n,l;n',l') |            */
/*                      3 R_oo         ( 2l + 1 )  |              |            */
/*                                                                             */
/*                                                                             */
/*   See for example Gordan Drake's                                            */
/*     Atomic, Molecular, & Optical Physics Handbook pg.638                    */
/*                                                                             */
/*   Here R_oo is the infinite mass Rydberg length                             */
/*                                                                             */
/*                                                                             */
/*        h c                                                                  */
/*   R_oo --- = 13.605698 eV                                                   */
/*        {e}                                                                  */
/*                                                                             */
/*                                                                             */
/*   R_oo =  2.179874e-11 ergs                                                 */
/*                                                                             */
/*   w = omega                                                                 */
/*     = frequency of transition from n,l to n',l'                             */
/*                                                                             */
/*                                                                             */
/*                                                                             */
/*     here g_k are statistical weights obtained from                          */
/*      the appropriate angular momentum quantum numbers                       */
/*                                                                             */
/*                                                                             */
/*                                                          -            -  2  */
/*                   64 pi^4 (e a_o)^2    max(l,l')        |              |    */
/*   A(n,l;n',l') = -------------------  -----------  v^3  | R(n,l;n',l') |    */
/*                         3 h c^3         2*l + 1         |              |    */
/*                                                          -            -     */
/*                                                                             */
/*                                                                             */
/*  pi             3.141592654                                                 */
/*  plank_hbar     6.5821220         eV sec                                    */
/*  R_oo           2.179874e-11      ergs                                      */
/*  plank_h        6.6260755e-34     J sec                                     */
/*  e_charge       1.60217733e-19    C                                         */
/*  a_o            0.529177249e-10   m                                         */
/*  vel_light_c    299792458L        m sec^-1                                  */
/*                                                                             */
/*                                                                             */
/*                                                                             */
/*  64 pi^4 (e a_o)^2    64 pi^4 (a_o)^2    e^2     1                      1   */
/*  ----------------- = ----------------- -------- ----  = 7.5197711e-38 ----- */
/*     3 h c^3                3  c^2       hbar c  2 pi                   sec  */
/*                                                                             */
/*                                                                             */
/*            e^2               1                                              */
/*  using ---------- = alpha = ----                                            */
/*          hbar c             137                                             */
/*******************************************************************************/

realnum H_Einstein_A(/*  IN THE FOLLOWING WE HAVE  n > n'                        */
	/* principal quantum number, 1 for ground, upper level      */
	long int n,
	/* angular momentum, 0 for s                                */
	long int l,
	/* principal quantum number, 1 for ground, lower level      */
	long int np,
	/* angular momentum, 0 for s                                */
	long int lp,
	/* Nuclear charge, 1 for H+, 2 for He++, etc                */
	long int iz
)
{
	DEBUG_ENTRY( "H_Einstein_A()" );

	ASSERT( iz >= 1 && iz <= LIMELM );

	realnum result;
	if( size_t(n) > t_hydro_tbl::Inst().nmaxnl_u() || size_t(np) > t_hydro_tbl::Inst().nmaxnl_l() )
	{
		QNPair inPair(n, l, 2, -1, np, lp, 2, -1);

		auto p = iso_sp[ipH_LIKE][iz-1].CachedAs.find(inPair);
		if( p != iso_sp[ipH_LIKE][iz-1].CachedAs.end() )
			return p->second;

		double mass_nuc;
		if( iz == 1 )
		{
			mass_nuc = PROTON_MASS; 
		}
		else
		{
			// this is an approximation for the mass of the nucleus...
			mass_nuc = ATOMIC_MASS_UNIT * dense.AtomicWeight[iz-1];
		}
		double Aul = H_Einstein_A_log10( n,l,np,lp,iz, mass_nuc );
		iso_sp[ipH_LIKE][iz-1].CachedAs[inPair] = Aul;
		result = (realnum)Aul;
	}
	else
	{
		result = t_hydro_tbl::Inst().tp( np, lp, n, l, iz );
	}
	return result;
}

/**********************log version****************************/
double H_Einstein_A_log10(/* returns Einstein A in units of (sec)^-1                      */
	long int n,
	long int l,
	long int np,
	long int lp,
	long int iz,
	double mass_nuc
)
{
	DEBUG_ENTRY( "H_Einstein_A_log10()" );

	/* validate the incoming data                         */
	if( iz < 1 )
	{
		fprintf(ioQQQ," The charge is impossible.\n");
		cdEXIT(EXIT_FAILURE);
	}
	if( n < 1 || np < 1 || l < 0 || l >= n || lp < 0 || lp >= np || abs(l-lp) != 1 )
	{
		fprintf(ioQQQ," The quantum numbers are impossible.\n");
		cdEXIT(EXIT_FAILURE);
	}
	if( n <= np  )
	{
		fprintf(ioQQQ," The principal quantum numbers are such that n <= n'.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/*  hv calculates photon energy in ergs for n -> n' transitions for H and H-like ions */
	double d1 = hv( n, np , iz, mass_nuc );
	double d2 = d1 / HPLANCK; /* v = hv / h */
	double d3 = pow3(d2);
	double lg = (double)(l > lp ? l : lp);
	double Two_L_Plus_One = (double)(2*l + 1);
	double d6 = lg / Two_L_Plus_One;
	double d7 = hri_log10( n, l, np, lp , iz, mass_nuc );
	double d8 = d7 * d7;
	double result = CONST_ONE * d3 * d6 * d8;

	return result;
}

/********************************************************************************/
/* hv calculates photon energy for n -> n' transitions for H and H-like ions    */
/*              simplest case of no "l" or "m" dependence                       */
/*  epsilon_0 = 1 in vacu                                                       */
/*                                                                              */
/*                                                                              */
/*                         R_Z                                                  */
/*  Energy(n,Z) = -  Z^2 -------                                                */
/*                         n^2                                                  */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*  Friedrich -- Theoretical Atomic Physics, 2006,  pg. 87   eq. 2.14           */
/*                                                                              */
/*         u                                                                    */
/*  R_Z = --- R_inf  where                                                      */
/*        m_e                                                                   */
/*                                                                              */
/*          h c                                                                 */
/*  R_inf = --- =  2.179874e-11 ergs                                            */
/*           e                                                                  */
/*                                                                              */
/*  (Harmin Lecture Notes for course phy-651 Spring 1994)                       */
/*  where m_e (m_n) is the mass of and electron (nucleus)                       */
/*  and u is the reduced electron mass for the hydrogenic ion                   */
/*                                                                              */
/*                                                                              */
/*         m_e m_n       m_e                                                    */
/*    u = --------- = -----------                                               */
/*        m_e + m_n   1 + m_e/m_n                                               */
/*                                                                              */
/*        m_e                                                                   */
/*  Now  ----- = 0.000544617013                                                 */
/*        m_p                                                                   */
/*            u                                                                 */
/*  so that  --- =  0.999455679                                                 */
/*           m_e                                                                */
/*                                                                              */
/*  for the hydrogen atom                                                       */
/*                                                                              */
/*  the mass of the nucleus is in g                                             */
/*  returns energy of photon in ergs                                            */
/*                                                                              */
/*  hv (n,n',Z) is for transitions n -> n'                                      */
/*                                                                              */
/*  1 erg = 1e-07 J                                                             */
/********************************************************************************/

inline double hv( long int n, long int nprime, long int iz, double mass_nuc )
{
	DEBUG_ENTRY( "hv()" );

	ASSERT( nprime > 0 );
	ASSERT( n > nprime );
	ASSERT( iz > 0 );

	double n1 = (double)n;
	double n2 = n1*n1;
	double np1 = (double)nprime;
	double np2 = np1*np1;
	double rmr = reduced_mass_rel(mass_nuc);
	double izsqrd = (double)(iz*iz);

	double d1 = 1. / n2;
	double d2 = 1. / np2;
	double d3 = izsqrd * rmr * EN1RYD;
	double d4 = d2 - d1;
	double result = d3 * d4;

	return result;
}

/************************************************************************/
/*   hri_log10()                                                        */
/*   Calculate the hydrogen radial wavefunction integral                */
/*   for the dipole transition  l'=l-1  or  l'=l+1                      */
/*   for the higher energy state n,l  to the lower energy state n',l'   */
/*   no "m" dependence                                                  */
/************************************************************************/
/*      here we have a transition                                       */
/*      from the higher energy state n,l                                */
/*      to the lower energy state n',l'                                 */
/*      with a dipole selection rule on l and l'                        */
/************************************************************************/
/*                                                                      */
/*   hri_log10() test n,l,n',l'  for domain errors and                  */
/*              swaps n,l <--> n',l' for the case  l'=l+1               */
/*                                                                      */
/*   It then calls hrii_log()                                           */
/*                                                                      */
/*   Dec. 6, 1999                                                       */
/*   Robert Paul Bauman                                                 */
/************************************************************************/
/*   Account for the nucleus mass in the H-like atom.                   */
/*   Oct 8, 2016 - Peter van Hoof, Marios Chatzikos                     */
/************************************************************************/

inline double hri_log10( long int  n, long int l, long int np, long int lp , long int iz,
			 double mass_nuc )
{
	/**********************************************************************/
	/*    from higher energy -> lower energy                              */
	/*    Selection Rule for l and l'                                     */
	/*    dipole process only                                             */
	/**********************************************************************/

	DEBUG_ENTRY( "hri_log10()" );

	long int a;
	long int b;
	long int c;
	long int d;
	double ld1 = 0.;
	double Z = (double)iz;

	ASSERT( n > 0 );
	ASSERT( np > 0 );
	ASSERT( l >= 0 );
	ASSERT( lp >= 0 );
	ASSERT( n > l );
	ASSERT( np > lp );
	ASSERT( n > np );
	ASSERT( iz > 0 );
	ASSERT( lp == l + 1 || lp == l - 1 );
	ASSERT( mass_nuc > 0. );

	if( l == lp + 1 )
	{
		/*        Keep variable  the same                                 */
		a = n;
		b = l;
		c = np;
		d = lp;
	}
	else if( l == lp - 1 )
	{
		/* swap n,l with n',l'                                            */
		a = np;
		b = lp;
		c = n;
		d = l;
	}
	else
	{
		printf( "BadMagic: l and l' do NOT satisfy dipole requirements.\n\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/****************************************************************/
	/*    Take care of the Z- and reduced mass-dependence here.     */
	/*    Computed like this, hri_log10() is in atomic units (a0)   */
	/* The line strengths S of this transition is S = hri_log10()^2 */
	/****************************************************************/
	ld1 = hrii_log(a, b, c, d) / (Z * reduced_mass_rel(mass_nuc));

	return ld1;
}

STATIC double hrii_log( long int  n, long int  l, long int np, long int lp)
{
	/******************************************************************************/
	/*      this routine hrii_log() is internal to the parent routine hri_log10() */
	/*      this internal routine only considers the case l=l'+1                  */
	/*      the case l=l-1 is done in the parent routine hri_log10()              */
	/*      by the transformation n <--> n' and l <--> l'                         */
	/*      THUS WE TEST FOR                                                      */
	/*          l=l'-1                                                            */
	/******************************************************************************/
	/**************************************************************************************/
	/*   THIS HAS THE GENERAL FORM GIVEN BY (GORDON 1929):                                */
	/*                                                                                    */
	/*    R(n,l;n',l') = (-1)^(n'-l) [4(2l-1)!]^(-1) *                                    */
	/*                      [(n+l)! (n'+l-1)!/(n-l-1)! (n'-l)!]^(1/2) *                   */
	/*                      (4 n n')^(l+1) (n-n')^(n+n'-2l-2) (n+n')^(-n-n') *            */
	/*                      { F(-n+l+1,-n'+l;2l;-4nn'/[n-n']^2)                           */
	/*                    - (n-n')^2 (n+n')^(-2) F(-n+l-1,-n'+l;2l; -4nn'/[n-n']^2 ) }    */
	/**************************************************************************************/

	DEBUG_ENTRY( "hrii_log()" );

	double y = 0.;
	double log10_fsf = 0.;

	ASSERT( n != np && l == lp + 1 );

	if( lp == np-1 ) /* A Pair of Easy Special Cases          */
	{
		if( l == n-1 )
		{
			/**********************************************************************/
			/*   R(n,l;n',l') = R(n,n-l;n-1,n-2)                                  */
			/*                                                                    */
			/*                = [(2n-2)(2n-1)]^(1/2)   [4n(n-1)/(2n-1)^2]^n  *    */
			/*                            [(2n-1) - 1/(2n-1)]/4                   */
			/**********************************************************************/

			double d1 = (double)((2*n-2)*(2*n-1));
			double d2 = sqrt( d1 );
			double d3 = (double)(4*n*(n-1));
			double d4 = (double)(2*n-1);
			double d5 = d4*d4;
			double d7 = d3/d5;
			double d8 = powi(d7,n);
			double d9 = 1./d4;
			double d10 = d4 - d9;
			double d11 = 0.25*d10;
			double result = (d2 * d8 * d11);

			ASSERT( d1 >= 0. );
			ASSERT( d3 >= 0. );
			return result;
		}
		else
		{
			double result = 0.;
			double ld1 = 0., ld2 = 0., ld3 = 0., ld4 = 0., ld5 = 0., ld6 = 0., ld7 = 0.;

			/******************************************************************************/
			/*   R(n,l;n',l') = R(n,l;l,l-1)                                              */
			/*                                                                            */
			/*                = [(n-l) ... (n+l)/(2l-1)!]^(1/2) [4nl/(n-l)^2]^(l+1) *     */
			/*                            [(n-l)/(n+l)]^(n+l) {1-[(n-l)/(n+l)]^2}/4       */
			/******************************************************************************/
			/**************************************/
			/*    [(n-l) ... (n+l)]               */
			/**************************************/
			/*    log10[(n-l) ... (n+l)] =        */
			/*                                    */
			/*        n+l                         */
			/*        ---                         */
			/*         >  log10(j)                */
			/*        ---                         */
			/*      j=n-l                         */
			/**************************************/

			ld1 = 0.;
			for( long int i1 = (n-l); i1 <= (n+l); i1++ )  /* from n-l to n+l INCLUSIVE           */
			{
				double d1 = (double)(i1);
				ld1 += log10( d1 );
			}

			/**************************************/
			/*    (2l-1)!                         */
			/**************************************/
			/*    log10[ (2n-1)! ]                */
			/**************************************/

			ld2 = lfactorial( 2*l - 1 );

			ASSERT( ((2*l)+1) >= 0);

			/**********************************************/
			/* log10( [(n-l) ... (n+l)/(2l-1)!]^(1/2) ) = */
			/*    (1/2) log10[(n-l) ... (n+l)] -          */
			/*            (1/2) log10[ (2n-1)! ]          */
			/**********************************************/

			ld3 = 0.5 * (ld1 - ld2);

			/**********************************************/
			/*            [4nl/(n-l)^2]^(l+1)             */
			/**********************************************/
			/*  log10( [4nl/(n-l)^2]^(l+1) ) =            */
			/*       (l+1) * log10( [4nl/(n-l)^2] )       */
			/*                                            */
			/*    = (l+1)*[ log10(4nl) - 2 log10(n-l) ]   */
			/*                                            */
			/**********************************************/

			double d1 = (double)(l+1);
			double d2 = (double)(4*n*l);
			double d3 = (double)(n-l);
			double d4 = log10(d2);
			double d5 = log10(d3);

			ld4 = d1 * (d4 - 2.*d5);

			/**********************************************/
			/*             [(n-l)/(n+l)]^(n+l)            */
			/**********************************************/
			/*   log10( [ (n-l)/(n+l) ]^(n+l)  ) =        */
			/*                                            */
			/*    (n+l) * [ log10(n-l) - log10(n+l) ]     */
			/*                                            */
			/**********************************************/

			d1 = (double)(n-l);
			d2 = (double)(n+l);
			d3 = log10( d1 );
			d4 = log10( d2 );

			ld5 = d2 * (d3 - d4);

			/**********************************************/
			/*      {1-[(n-l)/(n+l)]^2}/4                 */
			/**********************************************/
			/*   log10[ {1-[(n-l)/(n+l)]^2}/4 ]           */
			/**********************************************/

			d1 = (double)(n-l);
			d2 = (double)(n+l);
			d3 = d1/d2;
			d4 = d3*d3;
			d5 = 1. - d4;
			double d6 = 0.25*d5;

			ld6 = log10(d6);

			/******************************************************************************/
			/*   R(n,l;n',l') = R(n,l;l,l-1)                                              */
			/*                                                                            */
			/*                = [(n-l) ... (n+l)/(2l-1)!]^(1/2) [4nl/(n-l)^2]^(l+1) *     */
			/*                            [(n-l)/(n+l)]^(n+l) {1-[(n-l)/(n+l)]^2}/4       */
			/******************************************************************************/

			ld7 = ld3 + ld4 + ld5 + ld6;

			result = exp10( ld7 );

			ASSERT( result > 0. );
			return result;
		}
	}
	else
	{
		double result = 0.;
		long int a = 0, b = 0, c = 0;
		double d1 = 0., d2 = 0., d3 = 0., d4 = 0., d5 = 0., d6 = 0., d7 = 0.;
		mx d00, d01, d02, d03;

		/* this decides whether to do recursion over a or b */
		char ab = ( n < np ) ? 'a' : 'b';

		/**************************************************************************************/
		/*   THIS HAS THE GENERAL FORM GIVEN BY (GORDON 1929):                                */
		/*                                                                                    */
		/*    R(n,l;n',l') = (-1)^(n'-l) [4(2l-1)!]^(-1) *                                    */
		/*                      [(n+l)! (n'+l-1)!/(n-l-1)! (n'-l)!]^(1/2) *                   */
		/*                      (4 n n')^(l+1) (n-n')^(n+n'-2l-2) (n+n')^(-n-n') *            */
		/*                      { F(-n+l+1,-n'+l;2l;-4nn'/[n-n']^2)                           */
		/*                    - (n-n')^2 (n+n')^(-2) F(-n+l-1,-n'+l;2l; -4nn'/[n-n']^2 ) }    */
		/**************************************************************************************/

		/****************************************************************************************************/
		/* Calculate the whole shootin match                                                                */
		/*                            -                - (1/2)                                              */
		/*             (-1)^(n'-l)   | (n+l)! (n'+l-1)! |                                                   */
		/*  fsff() =   ----------- * | ---------------- | * (4 n n')^(l+1) (n-n')^(n+n'-2l-2)(n+n')^(-n-n') */
		/*             [4 (2l-1)!]   | (n-l-1)! (n'-l)! |                                                   */
		/*                            -                -                                                    */
		/* This is used in the calculation of hydrogen radial wave function integral for dipole transitions */
		/****************************************************************************************************/

		log10_fsf = fsff_log10( n, l, np );

		/******************************************************************************************/
		/*  2_F_1( a, b; c; y )                                                                   */
		/*                                                                                        */
		/*        F21_mx(-n+l+1, -n'+l; 2l; -4nn'/[n-n']^2)                                       */
		/*                                                                                        */
		/*                                                                                        */
		/*      Use a -> a' + 1                                                                   */
		/*                                 _                 _                                    */
		/*              (a' + 1) (1 - x)  |                   |                                   */
		/*      F(a') = ----------------- | F(a'+1) - F(a'+2) | + (a' + 1 + bx -c) F(a'+1)        */
		/*                 (a' + 1 - c)   |                   |                                   */
		/*                                 -                 -                                    */
		/*                                                                                        */
		/*      For the first F() in the solution of the radial integral                          */
		/*                                                                                        */
		/*            a = ( -n + l + 1 )                                                          */
		/*                                                                                        */
		/*      a = -n + l + 1                                                                    */
		/*      max(a)  = max(-n)  + max(l)     + 1                                               */
		/*              = -n       + max(n-1)   + 1                                               */
		/*              = -n       + n-1        + 1                                               */
		/*              = 0                                                                       */
		/*                                                                                        */
		/*      similarly                                                                         */
		/*                                                                                        */
		/*      min(a) = min(-n)   + min(l)   + 1                                                 */
		/*             = min(-n)   + 0        + 1                                                 */
		/*             = -n        + 1                                                            */
		/*                                                                                        */
		/*      a -> a' + 1   implies                                                             */
		/*                                                                                        */
		/*      max(a') = -1                                                                      */
		/*      min(a') = -n                                                                      */
		/******************************************************************************************/

		/**********************************************************************************/
		/*  2_F_1( a, b; c; y )                                                           */
		/*                                                                                */
		/*        F21_mx(-n+l+1, -n'+l; 2l; -4nn'/[n-n']^2)                               */
		/*                                                                                */
		/*             -4 nn'                                                             */
		/*  where Y = -------- .                                                          */
		/*            (n-n')^2                                                            */
		/*                                                                                */
		/**********************************************************************************/

		d2 = (double)(n - np);
		d3 = d2 * d2;

		d4 = 1. / d3;
		d5 = (double)(n * np);
		d6 = d5 * 4.;

		d7 = -d6;
		y = d7 * d4;

		/**************************************************************************************************/
		/*                                   THE GENERAL CASE                                             */
		/*                   USE RECURSION RELATION FOR HYPERGEOMETRIC FUNCTIONS                          */
		/*     for F(a,b;c;x) we have from eq.4   D. Hoang-Bing Astron. Astrophys. 238: 449-451 (1990)    */
		/*                                                                                                */
		/*     (a-c) F(a-1) = a (1-x) [ F(a) - F(a+1) ] + (a + bx - c) F(a)                               */
		/*                                                                                                */
		/*                     a (1-x)                       (a + bx - c)                                 */
		/*           F(a-1) = --------- [ F(a) - F(a+1) ] + -------------- F(a)                           */
		/*                     (a - c)                          (a - c)                                   */
		/*                                                                                                */
		/*                                                                                                */
		/*       A similar recusion relation holds for b with a <--> b.                                   */
		/*                                                                                                */
		/*                                                                                                */
		/*   we have initial conditions                                                                   */
		/*                                                                                                */
		/*                                                                                                */
		/*    F(0) = 1              with a = -1                                                           */
		/*                                                                                                */
		/*                 b                                                                              */
		/*   F(-1) = 1 - (---) x    with a = -1                                                           */
		/*                 c                                                                              */
		/**************************************************************************************************/

		/* d00 = F(-n+l+1,-n'+l;2l;-4nn'/[n-n']^2) is calculated below */

		/**************************************************************/
		/*  For the second F() in the solution of the radial integral */
		/*                                                            */
		/*        a = ( -n + l - 1 )                                  */
		/*                                                            */
		/*  a = -n + l - 1                                            */
		/*  max(a) = max(-n) + max(l)   - 1                           */
		/*         = -n      + (n - 1)  - 1                           */
		/*         = -2                                               */
		/*                                                            */
		/*  similarly                                                 */
		/*                                                            */
		/*  min(a) = min(-n) + min(l)   - 1                           */
		/*         = (-n)    + 0        - 1                           */
		/*         = -n - 1                                           */
		/*                                                            */
		/*  a -> a' + 1   implies                                     */
		/*                                                            */
		/*  max(a') = -3                                              */
		/*  min(a') = -n - 2                                          */
		/**************************************************************/

		a = -n + l - 1;
		b = -np + l;
		c = 2 * l;

		d01 = F21_mx(a, b, c, y, d00, ab );

		/**************************************************************************************/
		/*   THIS HAS THE GENERAL FORM GIVEN BY (GORDON 1929):                                */
		/*                                                                                    */
		/*    R(n,l;n',l') = (-1)^(n'-l) [4(2l-1)!]^(-1) *                                    */
		/*                      [(n+l)! (n'+l-1)!/(n-l-1)! (n'-l)!]^(1/2) *                   */
		/*                      (4 n n')^(l+1) (n-n')^(n+n'-2l-2) (n+n')^(-n-n') *            */
		/*                      { F(-n+l+1,-n'+l;2l;-4nn'/[n-n']^2)                           */
		/*                    - (n-n')^2 (n+n')^(-2) F(-n+l-1,-n'+l;2l; -4nn'/[n-n']^2 ) }    */
		/*                                                                                    */
		/*                = fsf * ( F(a,b,c;y)  - d3 * F(a',b',c';y) )                        */
		/*                                                                                    */
		/*        where d3 = (n-n')^2 (n+n')^2                                                */
		/*                                                                                    */
		/**************************************************************************************/

		/**************************************************************/
		/*  Calculate                                                 */
		/*                                                            */
		/*  (n-n')^2                                                  */
		/*  --------                                                  */
		/*  (n+n')^2                                                  */
		/**************************************************************/

		d1 = (double)((n - np)*(n -np));
		d2 = (double)((n + np)*(n + np));
		d3 = d1 / d2;

		d02.x = d01.x;
		d02.m = d01.m * d3;

		while ( fabs(d02.m) > maxm )
		{
			d02.m /= maxm;
			d02.x += maxe;
		}

		d03.x = d00.x;
		d03.m = d00.m - d02.m * exp10i(d02.x - d00.x);

		result = exp10( log10_fsf + d03.x ) * d03.m;

		ASSERT( result != 0. );

		/* we don't care about the correct sign of result... */
		return fabs(result);
	}
}

STATIC double fsff_log10( long int n, long int l, long int np )
{
	DEBUG_ENTRY( "fsff_log10()" );

	double d0 = 0., d1 = 0.;
	double ld0 = 0., ld1 = 0., ld2 = 0., ld3 = 0., ld4 = 0.;
	double result = 0.;

	/******************************************************************************************************/
	/*    Calculate the log10 of the whole shootin match                                                  */
	/*                    -                - (1/2)                                                        */
	/*          1        | (n+l)! (n'+l-1)! |                                                             */
	/*     ----------- * | ---------------- |      * (4 n n')^(l+1) |(n-n')^(n+n'-2l-2)| (n+n')^(-n-n')   */
	/*     [4 (2l-1)!]   | (n-l-1)! (n'-l)! |                                                             */
	/*                    -                -                                                              */
	/******************************************************************************************************/

	/**********************/
	/* Calculate (2l-1)!  */
	/**********************/

	d0 = (double)(2*l - 1);

	ld0 = lfactorial( 2*l - 1 );
	ld1 = log10(4.);
	result = -(ld0 + ld1);

	/**********************************************************************/
	/* Calculate (4 n n')^(l+1)                                           */
	/* powi( m , n) calcs m^n                                             */
	/* returns long double with m,n ints                                  */
	/**********************************************************************/

	d0 = (double)(4 * n * np);
	d1 = (double)(l + 1);
	result += d1 * log10(d0);

	/**********************************************************************/
	/* Calculate |(n-n')^(n+n'-2l-2)|                                     */
	/*    NOTE: Here we are interested only                               */
	/*             magnitude of (n-n')^(n+n'-2l-2)                        */
	/**********************************************************************/

	d0 = (double)(n - np);
	d1 = (double)(n + np - 2*l - 2);
	result += d1 * log10(fabs(d0));

	/* Calculate (n+n')^(-n-n') */
	d0 = (double)(n + np);
	d1 = (double)(-n - np);
	result += d1 * log10(d0);

	/**********************************************************************/
	/*                  -                -  (1/2)                         */
	/*                 | (n+l)! (n'+l-1)! |                               */
	/*     Calculate   | ---------------- |                               */
	/*                 | (n-l-1)! (n'-l)! |                               */
	/*                  -                -                                */
	/**********************************************************************/

	ld0 = lfactorial( n + l );
	ld1 = lfactorial( np + l - 1 );
	ld2 = lfactorial( n - l - 1 );
	ld3 = lfactorial( np - l );
	ld4 = 0.5*((ld0+ld1)-(ld2+ld3));

	result += ld4;
	ASSERT( result != 0. );
	return result;
}

STATIC mx F21_mx( long int a , long int b, long int c, double y, mx& Fa_plus, char ab )
{
	DEBUG_ENTRY( "F21_mx()" );

	/**********************************************************************************************/
	/*  begin sanity check, check order, and that none negative                                   */
	/*                                  THE GENERAL CASE                                          */
	/*                  USE RECURSION RELATION FOR HYPERGEOMETRIC FUNCTIONS                       */
	/*    for F(a,b;c;x) we have from eq.4   D. Hoang-Bing Astron. Astrophys. 238: 449-451 (1990) */
	/*                                                                                            */
	/*    (a-c) F(a-1) = a (1-x) [ F(a) - F(a+1) ] + (a + bx - c) F(a)                            */
	/*                                                                                            */
	/*                    a (1-x)                        (a + bx - c)                             */
	/*           F(a-1) = --------- [ F(a) - F(a+1) ] + -------------- F(a)                       */
	/*                     (a-c)                            (a-c)                                 */
	/*                                                                                            */
	/*                                                                                            */
	/*      A similar recusion relation holds for b with a <--> b.                                */
	/*                                                                                            */
	/*                                                                                            */
	/*  we have initial conditions                                                                */
	/*                                                                                            */
	/*                                                                                            */
	/*   F(0) = 1              with a = -1                                                        */
	/*                                                                                            */
	/*                b                                                                           */
	/*  F(-1) = 1 - (---) x    with a = -1                                                        */
	/*                c                                                                           */
	/*                                                                                            */
	/*  For the first F() in the solution of the radial integral                                  */
	/*                                                                                            */
	/*        a = ( -n + l + 1 )                                                                  */
	/*                                                                                            */
	/*  a = -n + l + 1                                                                            */
	/*  max(a) = max(-n) + max(l)     + 1                                                         */
	/*         = max(-n) + max(n - 1) + 1                                                         */
	/*         = max(-n + n - 1)      + 1                                                         */
	/*         = max(-1)              + 1                                                         */
	/*         = 0                                                                                */
	/*                                                                                            */
	/*  similarly                                                                                 */
	/*                                                                                            */
	/*  min(a) = min(-n) + min(l)   + 1                                                           */
	/*         = min(-n) + 0        + 1                                                           */
	/*         = (-n)    + 0        + 1                                                           */
	/*         = -n + 1                                                                           */
	/*                                                                                            */
	/*  a -> a' + 1   implies                                                                     */
	/*                                                                                            */
	/*  max(a') = -1                                                                              */
	/*  min(a') = -n                                                                              */
	/*                                                                                            */
	/*  For the second F() in the solution of the radial integral                                 */
	/*                                                                                            */
	/*        a = ( -n + l - 1 )                                                                  */
	/*                                                                                            */
	/*  a = -n + l + 1                                                                            */
	/*  max(a) = max(-n) + max(l)   - 1                                                           */
	/*         = -n      + (n - 1)  - 1                                                           */
	/*         = -2                                                                               */
	/*                                                                                            */
	/*  similarly                                                                                 */
	/*                                                                                            */
	/*  min(a) = min(-n) + min(l)   - 1                                                           */
	/*         = (-n)    + 0        - 1                                                           */
	/*         = -n      - 1                                                                      */
	/*                                                                                            */
	/*  a -> a' + 1   implies                                                                     */
	/*                                                                                            */
	/*  max(a') = -3                                                                              */
	/*  min(a') = -n - 2                                                                          */
	/**********************************************************************************************/

	ASSERT( a <= -2 );
	ASSERT( b <= 0 );
	ASSERT( c > 0 );

	/**********************************************************************************************/
	/* For the 2F1 function the following holds: 2F1( a, b, c, y ) == 2F1( b, a, c, y ).          */
	/* So we can do recursion either over a or b. The parameter ab decides which is stable.       */
	/*                                                                                            */
	/* This routine returns 2F1( a, b, c, y ) as the function result and 2F1( a+2, b, c, y ) in   */
	/* Fa_plus. If we do the recurion over a, then Fa_plus can be obtained for free since it      */
	/* was already calculated in the recursion to get 2F1( a, b, c, y ).                          */
	/**********************************************************************************************/
	mx result_mx;
	if( ab == 'a' )
	{
		result_mx = F21i_log(a, b, c, y, Fa_plus);
	}
	else if( ab == 'b' )
	{
		mx dum;
		result_mx = F21i_log(b, a, c, y, dum);
		Fa_plus = F21i_log(b, a+2, c, y, dum);
	}
	else
	{
		TotalInsanity();
	}
	return result_mx;
}

STATIC mx F21i_log( long int a, long int b, long int c, double y, mx& Fap2 )
{
	DEBUG_ENTRY( "F21i_log()" );

	mx Fa, Fap1;
	Fa.m = 1.;
	if( a == 0 )
		return Fa;

	Fap1 = Fa;
	// |Fa.m| < 2*max(n,np)^3, so no normalization needed here for realistic n, n' 
	Fa.m = 1. - y*(double(b)/double(c));
	if( a == -1 )
		return Fa;

	for( long ap = -2; ap >= a; --ap )
	{
		Fap2 = Fap1;
		Fap1 = Fa;

		/******************************************************************************************/
		/*                             -             -                                            */
		/*             (a)(1 - y)     |               |       (a + bx - c)                        */
		/*   F(a-1) = --------------  | F(a) - F(a+1) | +   --------------- F(a)                  */
		/*               (a - c)      |               |         (a - c)                           */
		/*                             -             -                                            */
		/*                                                                                        */
		/*                                                                                        */
		/*   with     F(0) = 1                                                                    */
		/*                         b                                                              */
		/*   and     F(-1) = 1 - (---) y                                                          */
		/*                         c                                                              */
		/*                                                                                        */
		/*                                                                                        */
		/*                                                                                        */
		/*   Use a -> a' + 1                                                                      */
		/*                               _                 _                                      */
		/*           (a' + 1) (1  - x)  |                   |    (a' + 1 + bx - c)                */
		/*   F(a') = -----------------  | F(a'+1) - F(a'+2) | +  ----------------- F(a'+1)        */
		/*              (a' + 1 - c)    |                   |      (a' + 1 - c)                   */
		/*                               -                 -                                      */
		/*                                                                                        */
		/*   For the first F() in the solution of the radial integral                             */
		/*                                                                                        */
		/*         a = ( -n + l + 1 )                                                             */
		/*                                                                                        */
		/*   a = -n + l + 1                                                                       */
		/*   max(a) = max(-n) + max(l)     + 1                                                    */
		/*          = -n      + max(n-1)   + 1                                                    */
		/*          = -n      + n-1        + 1                                                    */
		/*          = 0                                                                           */
		/*                                                                                        */
		/*   similarly                                                                            */
		/*                                                                                        */
		/*   min(a) = min(-n)   + min(l)   + 1                                                    */
		/*          = min(-n)   + 0        + 1                                                    */
		/*          = -n + 1                                                                      */
		/*                                                                                        */
		/*   a -> a' + 1   implies                                                                */
		/*                                                                                        */
		/*   max(a') = -1                                                                         */
		/*   min(a') = -n                                                                         */
		/******************************************************************************************/

		double d0 = 1./double(ap + 1 - c);
		double d3 = double(ap+1)*(1.-y)*d0;
		double d6 = 1. + double(b)*y*d0;

		mx d10, d11;
		if( LIKELY(Fap1.m != 0.) )
		{
			d10.m = Fap1.m  - Fap2.m * exp10i(Fap2.x - Fap1.x);
			d10.x = Fap1.x;
		}
		else
		{
			d10.m = -Fap2.m;
			d10.x = Fap2.x;
		}

		d10.m *= d3;

		d11.x = Fap1.x;
		d11.m = d6 * Fap1.m;

		if( LIKELY(d11.m != 0.) )
		{
			Fa.m = d11.m + d10.m * exp10i(d10.x - d11.x);
			Fa.x = d11.x;
		}
		else
		{
			Fa = d10;
		}

		while( UNLIKELY(fabs(Fa.m) > maxm) )
		{
			Fa.m /= maxm;
			Fa.x += maxe;
		}
	}
	return Fa;
}

/********************************************************************************/

inline void normalize_mx( mx& target )
{
	while( fabs(target.m) > maxm )
	{
		target.m /= maxm;
		target.x += maxe;
	}
	while( fabs(target.m) < 1.0/maxm )
	{
		target.m *= maxm;
		target.x -= maxe;
	}
	return;
}

inline mx add_mx( const mx& a, const mx& b )
{
	mx result;

	if( a.m != 0. )
	{
		result.x = a.x;
		result.m = a.m + b.m * exp10i(b.x - a.x);
	}
	else
	{
		result = b;
	}
	normalize_mx( result );
	return result;
}

inline mx sub_mx( const mx& a, const mx& b )
{
	mx result;
	mx minusb = b;
	minusb.m = -minusb.m;

	result = add_mx( a, minusb );
	normalize_mx( result );

	return result;
}

inline mx mxify( double a )
{
	mx result_mx;

	result_mx.x = 0;
	result_mx.m = a;
	normalize_mx( result_mx );

	return result_mx;
}

inline double unmxify( const mx& a_mx )
{
	return a_mx.m * exp10i( a_mx.x );
}

inline mx mxify_log10( double log10_a )
{
	mx result_mx;

	while ( log10_a > double(maxe) )
	{
		log10_a -= double(maxe);
		result_mx.x += maxe;
	}

	while ( log10_a < double(-maxe) )
	{
		log10_a += double(maxe);
		result_mx.x -= maxe;
	}

	result_mx.m = exp10(log10_a);

	return result_mx;
}

inline mx mult_mx( const mx& a, const mx& b )
{
	mx result;

	result.m = (a.m * b.m);
	result.x = (a.x + b.x);
	normalize_mx( result );

	return result;
}

inline double prodxx_log10( long int lp, double Ksqrd )
{
	/**********************************************************************/
	/*        |  s=l'                |     s=l'                           */
	/*        | -----                |     ---                            */
	/*  log10 |  | |  (1 + s^2 K^2)  | =    >    log10((1 + s^2 K^2))     */
	/*        |  | |                 |     ---                            */
	/*        |  s=0                 |     s=0                            */
	/**********************************************************************/

	if( lp == 0 )
		return 0.;

	avx_ptr<double> arg(lp), res(lp);
	for( long int s = 1; s <= lp; s++ )
		arg[s-1] = 1. + double(s*s)*Ksqrd;
	vlog10( arg.ptr0(), res.ptr0(), 0, lp );
	return reduce_a( res.ptr0(), 0, lp );
}
