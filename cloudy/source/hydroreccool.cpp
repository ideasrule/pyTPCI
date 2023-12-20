/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HydroRecCool hydrogen recombination cooling, called by iso_cool */
#include "cddefines.h"
#include "hydrogenic.h"
#include "phycon.h"
#include "iso.h"
#include "freebound.h"

double HydroRecCool(
	/* n is the prin quantum number on the physical scale */
	long int n , 
	/* nelem is the charge on the C scale, 0 is hydrogen */
	long int nelem)
{
	long int nm1; /* save n - 1*/

	double fac, 
	  hclf_v, 
	  x;
	static const double  a[15]={-26.6446988,-26.9066506,-27.0619832,-27.1826903,
	  -27.2783527,-27.3595949,-27.569406,-27.611159,-27.654748,-27.70479,
	  -27.745398,-27.776126,-27.807261,-27.833093,-27.866134};
	static const double  b[15]={-0.40511045,-0.41644707,-0.45834359,-0.49137946,
	  -0.51931762,-0.54971231,-0.18555807,-0.29204736,-0.36741085,
	  -0.45843009,-0.49753713,-0.51418739,-0.52287028,-0.52445456,
	  -0.52292803};
	static const double  c[15]={11.29232731,11.71035693,12.89309608,13.85569087,
	  14.67354775,15.56090323,6.147461,9.0304953,11.094731,13.630431,
	  14.721959,15.172335,15.413946,15.458123,15.428761};
	static const double  d[15]={.067257375,.07638384,.089925637,.102252192,
	  .111016071,.119518918,0.0093832482,0.028119606,0.039357697,0.050378417,
	  0.051644049,0.051367182,0.04938724,0.050139066,0.043085968};
	static const double  e[15]={-1.99108378,-2.26898352,-2.65163846,-3.02333001,
	  -3.29462338,-3.56633674,-1.0019228,-1.5128672,-1.8327058,-2.1866371,
	  -2.2286257,-2.1932699,-2.1205891,-2.1317169,-1.9175186};
	static const double  f[15]={-0.0050802618,-0.005849291,-0.0074854405,-0.0085677543,
	  -0.0093067267,-0.0098455637,0.040903604,0.037491802,0.035618861,
	  0.034132954,0.032418252,0.02947883,0.027393564,0.027607009,0.02433868};
	static const double  g[15]={.166267838,.196780541,.242675042,.282237211,
	  .310204623,.335160025,-0.81087376,-0.73435108,-0.69164333,-0.64907209,
	  -0.61216299,-0.55239109,-0.51048669,-0.51963194,-0.4504203};
	static const double  h[15]={.00020528663,.00027588492,.00033980652,.00041445515,
	  .00046423276,.0005121808,-0.0011986559,-0.0011333973,-0.0010992935,
	  -0.0010878727,-0.0010412678,-0.00095539899,-0.00089141547,-0.00089294364,
	  -0.00079179756};
	static const double  i[15]={-0.0071357493,-0.0099630604,-0.01178647,-0.014696455,
	  -0.01670318,-0.01887373,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.};

	DEBUG_ENTRY( "HydroRecCool()" );

	/* confirm that n is 1 or greater. */
	ASSERT( n > 0 );

	/* this is log of (temperature divided by charge squared) since sqlogz is log10 Z^2 */
	x = phycon.telogn[0] - phycon.sqlogz[nelem];

	/* this will be called for very silly low temperatures, due to high charge.
	 * evaluate at lowest fitted temp, and then extrapolate using known form.
	 * this is ok since these very high charge species do not contribute much
	 * recombination cooling */

	if( n > 15 || x < 0.2 )
	{
		/*double oldfac;*/
		/* use scale factor from actual recombination rate for this level, this ion,
		 * fac accounts for decreasing efficiency high up
		 * fac is 0.38 at n=15, and 0.32 at 25 */
		/*oldfac = 0.219 + (2.58 - 2.586/POW2((double)n));
		oldfac = 0.35*0.69* pow( (15./(double)n ) , 0.67 );
		fprintf(ioQQQ,"%e %e %e\n", oldfac , fac , fac/oldfac );*/
		/* >>chng 00 dec 07 use LaMothe & Ferland rates */
		/* >>refer	H	recom cool	LaMothe, J., & Ferland, G.J., 2001, PASP, in press */
		fac = HCoolRatio( phycon.te*POW2((double)n) /POW2(nelem+1.) );
		hclf_v = iso_sp[ipH_LIKE][nelem].fb[n].RadRecomb[ipRecRad]*phycon.te*BOLTZMANN* fac;
		return( hclf_v );
	}

	/* bail if te too high (if this ever happens, change logic so that HCoolRatio is 
	 * used in this limit - the process must be small in this case and routine is
	 * well bounded at high-energy end)*/
	if( x > phycon.TEMP_LIMIT_HIGH_LOG )
	{
		fprintf( ioQQQ, " HydroRecCool called with invalid temperature=%e nelem=%li\n", 
		  phycon.te , nelem );
		cdEXIT(EXIT_FAILURE);
	}

	/* convert onto c array for n*/
	nm1 = n - 1;

	if( nelem == 0 )
	{
		/* simple case, most often called, for hydrogen itself */
		fac = (a[nm1] + 
		  c[nm1]*phycon.telogn[0] + 
		  e[nm1]*phycon.telogn[1] + 
		  g[nm1]*phycon.telogn[2] + 
		  i[nm1]*phycon.telogn[3])/
		  (1. + b[nm1]*phycon.telogn[0] + 
		  d[nm1]*phycon.telogn[1] + 
		  f[nm1]*phycon.telogn[2] + 
		  h[nm1]*phycon.telogn[3]);
	}
	else
	{
		/* hydrogenic ions, expand as powers in t-z2 */
		fac = (a[nm1] + 
		  c[nm1]*x + 
		  e[nm1]*POW2( x) + 
		  g[nm1]*POW3( x) + 
		  i[nm1]*powi( x,4) ) /
		  (1. + b[nm1]*x  + 
		  d[nm1]*POW2( x ) + 
		  f[nm1]*POW3( x ) + 
		  h[nm1]*powi( x ,4) );
	}

	hclf_v = exp10(fac)*POW3(nelem+1.);
	return( hclf_v );
}

/* this function returns the ratio of cooling to recombination as
 * derived in 
 * >>refer	H	rec cooling	LaMothe, J., & Ferland, G.J., 2001, PASP, 113, 165 */
double HCoolRatio( 
	/* the scaled temperature, Tn^2/Z^2 */
	double t )
{
	double gamma;

	DEBUG_ENTRY( "HCoolRatio()" );

	if( t< 1e0 )
	{
		gamma = 1.;
	}
	else if( t < 7.4e5 )
	{
		double y;
		double x1,x2,x3,x4;
		x1=t;
		x2=t*sqrt(t);
		x3=t*t;
		x4=t*t*log(t);
		y=1.000285197084355-7.569939287228937E-06*x1
		+2.791888685624040E-08*x2-1.289820289839189E-10*x3
		+7.829204293134294E-12*x4;
		gamma = y;
	}
	else if( t < 5e10 )
	{
		double y;
		double x1,x2,x3,x4,xl;
		xl = log(t);
		x1=t;
		x2=xl*xl;
		x3=1.0/sqrt(t);
		x4=xl/(t*t);
		y=0.2731170438382388+6.086879204730784E-14*x1
		-0.0003748988159766978*x2+270.2454763661910*x3
		-1982634355.349780*x4;
		gamma = y;
	}
	else if( t < 3e14 )
	{
		double y;
		double x1,x2;
		x1=sqrt(t);
		x2=log(t);
		y=-17.02819709397900+4.516090033327356E-05*x1
		+1.088324678258230*x2;
		gamma = 1/y;
	}
	else
	{
		/*gamma = 3.85e11 * pow(t , -1. );*/
		gamma = 1.289e11 * pow(t , -0.9705 );
	}

	gamma = MIN2( 1.0, gamma );
	gamma = MAX2( 0.0, gamma );

	return( gamma );
}
