/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*hypho - create hydrogenic photoionization cross sections */
#include "cddefines.h"
#include "hypho.h"

/* some numbers used below */
static const int NCM = 3000;
static const int NFREQ = NCM;
/*exp1 routine from ucl group to compute 1-exp(-x)
 * mod implicit none, and f77 control, gfj '96 jul 11 */

STATIC double exp1(double x)
{
	double dx, 
	  exp1_v;

	DEBUG_ENTRY( "exp1()" );

	dx = fabs(x);
	if( dx < 1.0e-9 )
	{
		exp1_v = -x;
	}
	else if( dx < 1.0e-5 )
	{
		exp1_v = ((-x*0.5) - 1.0)*x;
	}
	else if( dx < 1.0e-3 )
	{
		exp1_v = (((-x*0.1666666667) - 0.5)*x - 1.0)*x;
	}
	else
	{
		exp1_v = 1.0 - exp(x);
	}

	return exp1_v;
}

/*hypho create hydrogenic photoionization cross sections */
void hypho(
	/* Z-Nelec, the residual charge, 1 for hydrogen, 2 for helium */
	double zed, 
	/* principal quantum number */
  long int n, 
  /* lowest angular momentum */
  long int lmin, 
  /* highest angular momentum */
  long int lmax, 
  /* scaled lowest energy, in units of zed^2/n^2, at which the cs will be done */
  double en, 
  /* number of points to do */
  long int ncell, 
  /* array of energies (in units given above) */
  realnum anu[], 
  /* calculated photoionization cross section in cm^-2 */
  realnum bfnu[])
{
	long int i, 
	  ifp, 
	  il, 
	  ilmax, 
	  j, 
	  k, 
	  l, 
	  lg, 
	  ll, 
	  llk, 
	  lll, 
	  llm, 
	  lm, 
	  lmax1, 
	  m, 
	  mulp, 
	  mulr;
	double a, 
	  ai, 
	  al, 
	  alfac, 
	  con2, 
	  e, 
	  ee, 
	  fll, 
	  flm, 
	  fn, 
	  fth, 
	  gl, 
	  gn0, 
	  gn1e, 
	  gne, 
	  gnt, 
	  p, 
	  p1, 
	  p2, 
	  se, 
	  sl, 
	  sl4, 
	  sm, 
	  sm4, 
	  sn, 
	  sn4, 
	  sum, 
	  zed2;

	static double ab[NFREQ], 
	  alo[7000], 
	  fal[7000], 
	  freq[NFREQ], 
	  g2[2][NFREQ], 
	  g3[2][NFREQ];

	double anl, 
	  con3, 
	  fac, 
	  x;
	static int first = true;
	static const double zero = 0.;
	static const double one = 1.;
	static const double two = 2.;
	static const double four = 4.;
	static const double con1 = 0.8559;

	DEBUG_ENTRY( "hypho()" );

	/*generate hydrogenic photoionization cross sections */
	/*     *************************************************************
	 *          GENERATE HYDROGENIC PHOTOIONIZATION CROSS SECTIONS
	 *
	 *     Hypho calculates Hydrogenic bound-free (BF) xsectns for each n,l .
	 *     For Opacity work we calculate l-weighted average for each n .
	 *     BF xsectns for Hydorgenic ion are tabulated at E/z**2. E is
	 *     the energy in Ry corresponding to the frequencies on the Opacity
	 *     mesh. It can changed accoding to need.
	 *
	 * zed=Z-Nelc, 
	 * n=principal quantum number
	 * lmin=lowest angular momentum considered of level n
	 * lmx=highest angular momentum considered of level n
	 * en=zed*zed/(n*n) lowest energy at which the hydrogenic cross sections are
	 * to be calculated
	 * anu = photon energy
	 * bfnu=photoionization cross section in cm^-2
	 *
	 *     local variables */
	/* CON1=conversion factor from a0**2 to Mb, x-sectns in megabarns
	 *     CON1=8.5594e-19 *  1.e+18 */

	vector<long> mm(NFREQ);

	if( ncell > NCM )
	{
		fprintf( stderr, " increase ncm in hypho to at least%5ld\n", 
		  ncell );
		cdEXIT(EXIT_FAILURE);
	}

	if( first )
	{
		fal[0] = zero;
		for( i=1; i < 7000; i++ )
		{
			ai = (double)(i);
			alo[i-1] = log10(ai);
			fal[i] = alo[i-1] + fal[i-1];
		}
		first = false;
	}

	/* these may not be redefined, and so will crash */
	ll = -INT_MAX;
	lm = -INT_MAX;
	anl = -DBL_MAX;

	/*          CALCULATION OF HYDROGENIC PHOTOIONIZATION CROSS SECTION
	 *
	 *         INITIALIZATION FOR N
	 * * con2=(a0**2*1.e+18)*(n/zed)**2, zed=z-nelc
	 * * for hydrogen, the threshold, fth=1/n**2, and for hydrogenic atom,it is
	 * zed**2/n**2. */

	zed2 = zed*zed;
	fn = (double)(n);
	sn = fn*fn;
	sn4 = four*sn;
	con2 = con1*POW2(fn/ zed);
	fth = one/sn;

	gn0 = 2.3052328943 - 2.302585093*fal[n+n-1] - fn*0.61370563888 + 
	  alo[n-1]*(fn + one)*2.30258093;
	lmax1 = MIN2(lmax,n-1);
	ilmax = n - lmin;

	/* * calculate top-up st.wt for given n and lmin=lmax+1 (R-matrix). subtract
	 * the sum from 2n**2. */
	gl = 0.;
	if( lmin != 0 )
	{
		for( i=0; i < lmin; i++ )
		{
			lg = i;
			gl += 2.*lg + 1;
		}
	}

	gnt = 2.*(POW2( (double)n ) - gl);

	/* define photon energies epnu corresponding to hydrogenic atom with charge
	 * z and freq to hydrogen atom;  INITIALIZE G'S */
	for( i=0; i < ncell; i++ )
	{
		mm[i] = 0;
		bfnu[i] = (realnum)zero;
		freq[i] = anu[i]*zed2;
		for( j=0; j < 2; j++ )
		{
			g2[j][i] = zero;
			g3[j][i] = zero;
		}
	}

	/* gjf Aug 95 change
	 * original code returned total &(*(*$# if freq(1)<en */
	freq[0] = MAX2(freq[0],en);

	/* calculate cross section:  L LOOP */

	for( il=1; il <= ilmax; il++ )
	{
		l = n - il;
		m = 0;
		al = (double)(l);
		k = n - l - 2;
		con3 = con2/(two*al + one);

		/* FREQUENCY LOOP (FREQ UNITS ARE RYD*ZED**2) */
		for( ifp=0; ifp < ncell; ifp++ )
		{
			if( freq[ifp] < fth )
			{
				if( l <= lmax1 )
					anl = 0.;
			}
			else
			{
				/* >>chng 99 jun 24, move m from below to end, change m-1 to m */
				/*m += 1;*/
				se = freq[ifp] - fth;
				if( se < 0. )
				{
					fprintf( stderr, " %4ld%12.4e%12.4e\n", ifp, 
					  freq[ifp], fth );
				}

				/*            if ( se .lt. -1.e-30) se = 1.e-06 */
				e = sqrt(se);
				x = one + sn*se;
				if( k < 0 )
				{
					if( e >= 0.314 )
					{
						ee = 6.2831853/e;
						p = 0.5*log10(exp1(-ee));
					}
					else
					{
						p = zero;
					}

					if( e > 1.0e-6 )
					{
						a = two*(fn - atan(fn*e)/e);
					}
					else
					{
						a = zero;
					}

					ab[m] = (gn0 + a)/2.302585 - p - (fn + two)*
					  log10(x);
					gne = 0.1;
					gn1e = x*gne/(fn + fn);
					g3[1][m] = gne;
					g3[0][m] = gn1e;
					g2[1][m] = gne*fn*x*(fn + fn - one);
					g2[0][m] = gn1e*(fn + fn - one)*(four + 
					  (fn - one)*x);
				}

				double g11 = zero;
				double g12 = zero;
				double g22 = g2[1][m];
				double g32 = g3[1][m];
				double g21 = g2[0][m];
				double g31 = g3[0][m];
				double muls = mm[m];

				if( k < 0 )
				{
					/*                    l.eq.n-1
					 * */
					g11 = g31;
					if( l == 0 )
						g11 = zero;
					g12 = g32;
				}

				else if( k == 0 )
				{

					/*                    l.eq.n-2
					 * */
					g11 = g21;
					if( l == 0 )
						g11 = zero;
					g12 = g22;
				}

				else
				{

					/*                     l.lt.n-2
					 * */
					if( k <= 1 )
					{
						ll = n - 1;
						lm = n - 2;
					}
					sl = POW2( (double)ll );
					sl4 = four*sl;
					fll = (double)(ll);
					g12 = (sn4 - sl4 + (two*sl - fll)*x)*g22 - 
					  sn4*(sn - sl)*(one + POW2(fll + one)*se)*
					  g32;

					if( l != 0 )
					{
						sm = POW2( (double)lm );
						sm4 = four*sm;
						flm = (double)(lm);
						g11 = (sn4 - sm4 + (two*sm + flm)*x)*g21 - 
						  sn4*(sn - POW2(flm + one))*(one + 
						  sm*se)*g31;
						g31 = g21;
						g21 = g11;
					}

					g32 = g22;
					g22 = g12;

					if( (m+1) == ncell )
					{
						ll -= 1;
						if( l != 0 )
							lm -= 1;
					}

					if( g12 >= 1.e20 )
					{
						muls += 35;
						g22 *= 1.e-35;
						g32 *= 1.e-35;
						g12 *= 1.e-35;

						if( l != 0 )
						{
							g11 *= 1.e-35;
							g21 *= 1.e-35;
							g31 *= 1.e-35;
						}
					}
				}

				mm[m] = muls;
				g2[1][m] = g22;
				g3[1][m] = g32;
				g2[0][m] = g21;
				g3[0][m] = g31;

				alfac = fal[n+l] - fal[n-l-1] + two*(al - fn)*
				  alo[n*2-1];

				p1 = one;
				lll = l + 1;
				llm = l - 1;
				mulr = 0;

				if( llm >= 1 )
				{
					for( i=1; i <= llm; i++ )
					{
						ai = (double)(i);
						p1 *= one + ai*ai*se;
						if( p1 >= 1.e20 )
						{
							p1 *= 1.e-10;
							mulr += 10;
						}
					}
				}

				p2 = p1;
				llk = llm + 1;
				if( llk < 1 )
					llk = 1;

				for( i=llk; i <= lll; i++ )
				{
					ai = (double)(i);
					p2 *= one + ai*ai*se;
				}

				mulp = 0;
				while( g12 >= one )
				{
					mulp -= 10;
					g12 *= 1.e-10;
					if( l != 0 )
						g11 *= 1.e-10;
				}

				sum = alfac + (realnum)(mulr) + two*(ab[m] + (realnum)(muls-
				  mulp+1));

				fac = zero;

				/*     >>chng 96 apr 19, 35 had been 50 */
				if( fabs(sum) < 35. )
					fac = exp10(sum);
				if( l != 0 )
					g11 *= g11*p1;
				g12 *= g12*p2;

				/* anl=con3*(g11 *al + g12 *(al+one))
				 * *x*fac */

				if( l <= lmax1 )
					anl = fac*x*con3*(g11*al + g12*(al + 1.));
				anl *= 2.*(2.*al + 1.);

				bfnu[ifp] += (realnum)(anl*1e-18);
				/* >>chng 99 jun 24, move m inc to here */
				++m;
			}
		}
	}

	/*  bfmin = 1.e-03 * bfnu(1) / gnt */
	for( i=0; i < ncell; i++ )
	{
		bfnu[i] /= (realnum)gnt;
	}

	return;
}
