/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*he_1trans compute Aul for given line	*/
/*ritoa - converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
/*ForbiddenAuls calculates transition probabilities for forbidden transitions.	*/
/*scqdri - stands for Semi-Classical Quantum Defect Radial Integral	*/
/*Jint - used by scqdri	*/
/*AngerJ - used by scqdri */

#include "cddefines.h" 
#include "dense.h"
#include "trace.h"
#include "hydro_bauman.h"
#include "iso.h"
#include "helike.h"
#include "helike_einsta.h"
#include "freebound.h"
#include "lines_service.h"
#include "integrate.h"
#include "container_classes.h"
#include "parser.h"
#include "thirdparty.h"

/* the array of transitions probabilities read from data file.  */
static map<QNPair, array<double,LIMELM>> TransProbs;

STATIC double helike_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo );
STATIC double helike_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo, long sLo, long jLo );

/*ritoa converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
STATIC double ritoa( long li, long lf, long nelem, double k, double RI2 );

/* handles all forbidden transition probabilities. */
STATIC double ForbiddenAuls(long nelem, long nHi, long lHi, long sHi, long jHi,
			    long nLo, long lLo, long sLo, long jLo);

/* used as parameters in qg32 integration */
static double vJint , zJint;

/* the integrand used in the AngerJ function below. */
STATIC double Jint( double theta )
{
	/*	[ cos[vx - z sin[x]] ]  */
	double 
		d1 = vJint * theta,
		d2 = zJint * sin(theta),
		d3 = (d1 - d2),
		d4 = cos(d3);

	return( d4 );
}

/* AngerJ function. */
// See also http://dlmf.nist.gov/11.10 for series expansion and recurrences
STATIC double AngerJ_mac( double vv, double zz )
{
	// Maclaurin series for J_nu(z)
	// http://dlmf.nist.gov/11.10#iii
	// Series is absolutely convergent, but can take
	// numerous terms for large z
	// Problems when vv=2 -- should fall back to Bessel functions?
	// Worth trying the acceleration techinques of http://dlmf.nist.gov/3.9
	double g1=tgamma(1+0.5*vv), g2=tgamma(1-0.5*vv),
		g3=tgamma(1.5+0.5*vv), g4=tgamma(1.5-0.5*vv);
	double y1 = 0., y2=0., zp = 1.;
	for (long k=0; k<100; ++k)
	{
		double y10 = y1;
		y1 += zp/(g1*g2);
		if (y1 == y10)
			break;
		zp *= 0.5*zz;
		// Gamma(t+1) = t*Gamma(t)
		g1 *= k+0.5*vv+1;
		g2 *= k-0.5*vv+1;
		y2 += zp/(g3*g4);
		zp *= -0.5*zz;
		g3 *= k+0.5*vv+1.5;
		g4 *= k-0.5*vv+1.5;
	}
	double cosfac=cos(0.5*PI*vv), sinfac=sin(0.5*PI*vv);
	return y1*cosfac+y2*sinfac;
}

STATIC double bessjnu( double vv, double zz )
{
	// DLMF 10.17.3
	// Asymptotic series for Bessel function for non-integer index
	double rz=1./zz, akrz=1., bff=0., bgg=0., bfp=0.;
	for (long k=0; k<100; ++k)
	{
		if (k > 0 && fabs(akrz) > fabs(bfp))
			break;
		bfp = akrz;
		bff += akrz;
		akrz *= (4*vv*vv-(4*k+1)*(4*k+1))/(8*(2*k+1)*zz);
		bgg += akrz;
		akrz *= -(4*vv*vv-(4*k+3)*(4*k+3))/(8*(2*k+2)*zz);
	}
	double omega = zz-(2*vv+1)*PI/4;
	return sqrt(2/PI*rz)*(cos(omega)*bff-sin(omega)*bgg);
}
STATIC double AngerJ_asymp( double vv, double zz )
{
	double ff = 0., gg=0., fk=1., gk=1., rz=1./zz, rzsq = rz*rz, 
		fp=0.;
	for (long k=0; k<100; ++k)
	{
		if (k > 0)
		{
			fk *= (vv*vv-(2*k-1)*(2*k-1))*rzsq;
			gk *= (vv*vv-4*k*k)*rzsq;
			if (fk > fp)
				break;
		}
		fp = fk;
		ff += fk;
		gg += gk;
	}
	return bessjnu(vv,zz)+sin(PI*vv)/(PI*zz)*(ff-(vv*rz)*gg);
}
STATIC double AngerJ( double vv, double zz )
{
	DEBUG_ENTRY( "AngerJ()" );

	if( vv == int(vv) )
	{
		return bessel_jn(int(vv), zz);
	}
	else
	{
		/* Estimate number of peaks in integrand. */                            
		/* Integrate region four peaks at a time. */
		long divsor = max(nint(fabs(vv)/4.), 1);
		vJint = vv;
		zJint = zz;

		double y = 0.0;
		for( long rep = 0; rep < divsor; rep++ )
		{
			double x_low = PI * double(rep)/double(divsor);
			double x_up  = PI * double(rep+1)/double(divsor);
			y += qg32( x_low, x_up, Jint );
		}

		return y / PI;
	}
}

/******************************************************************************/
/******************************************************************************/
/*                                                                            */
/*    Semi-Classical Quantum Defect Radial Integral                           */      
/*                                                                            */
/*   See for example                                                          */
/*     Atomic, Molecular & Optical Physics Handbook                           */
/*     Gordon W. F. Drake; Editor                                             */
/*     AIP Press                                                              */
/*     Woddbury, New York.                                                    */
/*     1996                                                                   */
/*                                                                            */
/* NOTE:: we do not include the Bohr Radius a_o in the                        */
/*           definition of of  R(n,L;n'L') as per Drake.                      */
/*                                                                            */
/*                                                                            */
/*                   1  (n_c)^2 | {      D_l max(L,L') }                      */
/*    R(n,L;n'L') = --- ------- | { 1 -  ------------- } J( D_n-1; -x )  -    */
/*                   Z   2 D_n  | {          n_c       }                      */
/*                                                                            */
/*                                                                            */
/*                    {      D_L max(L,L') }                                  */
/*                -   { 1 +  ------------- } J( D_n+1; -x )                   */
/*                    {          n_c       }                                  */
/*                                                                            */
/*                                                                            */
/*                     2                    |                                 */
/*                  + ---  sin(Pi D_n)(1-e) |                                 */
/*                     Pi                   |                                 */
/*                                                                            */
/*  where                                                                     */
/*        n_c = (2n*'n*)/(n*'+n*)                                             */
/*                                                                            */
/*       Here is the quantity Drake gives...                                  */
/*            n_c = ( 2.0 * nstar * npstar ) / ( nstar + npstar );            */
/*                                                                            */
/*       while V.A. Davidkin uses                                             */
/*            n_c = sqrt(  nstar * npstar  );                                 */
/*                                                                            */
/*        D_n = n*' - n*                                                      */
/*                                                                            */      
/*        D_L = L*' - L*                                                      */
/*                                                                            */
/*        x = e D_n                                                           */
/*                                                                            */
/*        Lmx  = max(L',L)                                                    */
/*                                                                            */
/*        e = sqrt( 1 - {Lmx/n_c}^2 )                                         */
/*                                                                            */
/*                                                                            */
/*       Here n* = n - qd where qd is the quantum defect                      */
/*                                                                            */
/******************************************************************************/
/******************************************************************************/
STATIC double scqdri(/* upper and lower quantum numbers...n's are effective	*/
			  double nstar, long int l,
			  double npstar, long int lp,
              double iz
              )
{
	double n_c = ((2.0 * nstar * npstar ) / ( nstar + npstar ));

	double D_n = (nstar - npstar);
	double D_l = (double) ( l - lp );
	double lg  = (double) ( (lp > l) ? lp : l);

	double h = (lg/n_c);
	double g = h*h;
	double f = ( 1.0 - g );
	double e = (( f >= 0.0) ? sqrt( f ) : 0.0 );

	double x  = (e * D_n);
	double z  = (-1.0 * x);
	double v1 = (D_n + 1.0);
	double v2 = (D_n - 1.0); 

	double d1,d2,d7,d8,d9,d34,d56,d6_1; 

	DEBUG_ENTRY( "scqdri()" );

	if( iz == 0.0 )
		iz += 1.0;

	if( D_n == 0.0 )
	{
		return( -1.0 );
	}

	if( D_n < 0.0 )
	{
		return( -1.0 );
	}

	if( f < 0.0 )
	{
		/* This can happen for certain  quantum defects   */
		/* in the lower n=1:l=0 state. In which case you  */
		/* probably should be using some other alogrithm  */
		/* or theory to calculate the dipole moment.      */
		return( -1.0 );
	}

	d1 = ( 1.0 / iz );

	d2 = (n_c * n_c)/(2.0 * D_n);

	if (0)
	{
		// Confirm accuracy of result of AngerJ function

		// relative accuracy is ~1e-10 in the worst case, where this error
		// may well be in summing series representation
		double vv = v2, zz = z, y = AngerJ(vv,zz);
		double yt = AngerJ_mac(vv,zz);
		double vp = (zz > 0.) ? vv : -vv;
		double ya = fabs(zz) < fabs(vv) ? 0. : AngerJ_asymp(vp,fabs(zz));
		printf("ANGER %15.8g %15.8g: %15.8g %15.8g (%15.8g) asymp %15.8g\n",
				 zz,vv,yt,y,fabs(yt-y)/(1e-37+fabs(yt)),ya);
	}

	d34 = (1.0 - ((D_l * lg)/n_c)) * AngerJ( v1, z );

	d56 = (1.0 + ((D_l * lg)/n_c)) * AngerJ( v2, z );

	d6_1 = PI * D_n;

	d7 = (2./PI) * sin( d6_1 ) * (1.0 - e);

	d8 = d1 * d2 * ( (d34) - (d56) + d7 );

	d9 = d8 * d8;

	ASSERT( D_n  > 0.0 );
	ASSERT( l  >= 0  );
	ASSERT( lp >= 0 );
	ASSERT( (l == lp + 1) || ( l == lp - 1) );
	ASSERT( n_c != 0.0 );
	ASSERT( f >= 0.0 );
	ASSERT( d9  > 0.0 );

	return( d9 );
}

STATIC double ForbiddenAuls(long nelem, long nHi, long lHi, long sHi, long jHi,
			    long nLo, long lLo, long sLo, long jLo)
{
	double A;

	DEBUG_ENTRY( "ForbiddenAuls()" );

	QNPack inLo = QN2ind(nLo, lLo, sLo, 2*jLo+1);
	QNPack inHi = QN2ind(nHi, lHi, sHi, 2*jHi+1);

	QNPack inHe1s1S = QN2ind(1, 0, 1, 1);
	QNPack inHe2s3S = QN2ind(2, 0, 3, 3);
	QNPack inHe2s1S = QN2ind(2, 0, 1, 1);
	QNPack inHe2p3P0 = QN2ind(2, 1, 3, 1);
	QNPack inHe2p3P1 = QN2ind(2, 1, 3, 3);
	QNPack inHe2p3P2 = QN2ind(2, 1, 3, 5);
	QNPack inHe2p1P = QN2ind(2, 1, 1, 3);
	QNPack inHe3s3S = QN2ind(3, 0, 3, 3);
	QNPack inHe3s1S = QN2ind(3, 0, 1, 1);
	QNPack inHe3p3P = QN2ind(3, 1, 3, -1);
	QNPack inHe3p1P = QN2ind(3, 1, 1, 3);

	if( inLo == inHe1s1S && nHi == 2 )
	{
		if( inHi == inHe2s3S )
		{
			/* Parameters for 2^3S to ground transition.	*/
			if( nelem == ipHELIUM )
			{
				A = 1.272E-4;
				iso_put_error(ipHE_LIKE,nelem,inHe2s3S,inHe1s1S,IPRAD, 0.01f, 0.20f);
			}
			else
			{
				/* >>refer	Helike	As	Lin, C.D., Johnson, W.R., and Dalgarno, A. 1977,
				 * >>refercon	Phys. Rev. A 15, 1, 010015	*/
				// 2nu is treated elsewhere
				A = (3.9061E-7) * pow( (double)nelem+1., 10.419 );
				iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
			}
		}
		else if( inHi == inHe2s1S )
		{
			/* Parameters for 2^1S to ground transition.	*/
			if( nelem == ipHELIUM )
			{
				A = 1E-20;
				iso_put_error(ipHE_LIKE,nelem,inHe2s1S,inHe1s1S,IPRAD, 0.01f, 0.05f);
			}
			else
			{
				A = iso_ctrl.SmallA; // 2nu is treated elsewhere
				iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
			}
		}
		else if( inHi == inHe2p3P0 )
		{
			/* Parameters for 2^3P0 to ground transition.	*/
			if( nelem == ipHELIUM )
			{
				A = 1E-20;
				iso_put_error(ipHE_LIKE,nelem,inHe2p3P0,inHe1s1S,IPRAD, 0.00f, 0.00f);
			}
			else
			{
				A = iso_ctrl.SmallA;
				iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
			}
		}
		else if( inHi == inHe2p3P1 )
		{
			/* Parameters for 2^3P1 to ground transition.	*/
			if( nelem == ipHELIUM )
			{
				A = 177.58;
				iso_put_error(ipHE_LIKE,nelem,inHe2p3P1,inHe1s1S,IPRAD, 0.01f, 0.05f);
			}
			/* >>chng 06 jul 25, only use the fit to Johnson et al. values up to and
			 * including argon, where there calculation stops.  For higher Z use below */
			else if( nelem <= ipARGON )
			{
				A = ( 11.431 * pow((double)nelem, 9.091) );
				iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
			}
			else
			{
				/* a fit to the Lin et al. 1977. values, which go past zinc. */
				A = ( 383.42 * pow((double)nelem, 7.8901) );
				iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
			}
		}
		else if( inHi == inHe2p3P2 )
		{
			/* Parameters for 2^3P2 to ground transition.	*/
			if( nelem == ipHELIUM )
			{
				A = 0.327;
				iso_put_error(ipHE_LIKE,nelem,inHe2p3P2,inHe1s1S,IPRAD, 0.01f, 0.01f);
			}
			else
			{
				/* fit to Lin et al. 1977 values.  This fit is good
				 * to 7% for the range from carbon to iron. The Lin et al. values
				 * differ from the Hata and Grant (1984) values (only given from 
				 * calcium to copper) by less than 2%. */
				A = ( 0.11012 * pow((double)nelem, 7.6954) );
				iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
			} 
		}
		else
		{
			TotalInsanity();
		}
	}

	/* The next two cases are fits to probabilities given in 
	 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
	 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
	/* originally astro.ph. 0201454 */
	/* The involve Triplet P and Singlet S.  Rates for Triplet S to Singlet P 
	 * do not seem to be available.	*/

	/* Triplet P to Singlet S...Delta n not equal to zero!	*/
	else if( nelem > ipHELIUM && lHi == 1 && sHi == 3 && 
		lLo == 0 && sLo == 1 && nLo < nHi )
	{
		A = 8.0E-3 * exp(9.283/sqrt((double)nLo)) * pow((double)nelem,9.091) /
			pow((double)nHi,2.877);
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	/* Singlet S to Triplet P...Delta n not equal to zero!	*/
	else if( nelem > ipHELIUM && lHi == 0 && sHi == 1 && 
		lLo == 1 && sLo == 3 && nLo < nHi )
	{
		A = 2.416 * exp(-0.256*nLo) * pow((double)nelem,9.159) / pow((double)nHi,3.111);

		if( inLo == inHe2p3P0 || inLo == inHe2p3P2 )
		{
			/* This is divided by 3 instead of 9, because value calculated is specifically for 2^3P1.
			 * Here we assume statistical population of the other two.	*/
			A *= (2.*jLo+1.0)/3.0;
		}
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	else if( inLo == inHe2s3S && inHi == inHe3s3S )
	{
		double As_3TripS_to_2TripS[9] = {
			7.86E-9, 4.59E-6, 1.90E-4, 2.76E-3, 2.25E-2,
			1.27E-1, 5.56E-1, 2.01E0, 6.26E0 };

		/* These M1 transitions given by 
		 * >>refer He-like As Savukov, I.M., Labzowsky, and Johnson, W.R. 2005, PhysRevA, 72, 012504 */
		if( nelem <= ipNEON )
		{
			A = As_3TripS_to_2TripS[nelem-1];
			/* 20% error is based on difference between Savukov, Labzowsky, and Johnson (2005)
			 * and Lach and Pachucki (2001) for the helium transition. */
			iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.2f,0.2f);
		}
		else
		{
			/* This is an extrapolation to the data given above.  The fit reproduces
			 * the above values to 10% or better. */
			A= 7.22E-9*pow((double)nelem, 9.33);
			iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.3f,0.3f);
		}
	}

	else if( inLo == inHe2s3S && inHi == inHe2p1P )
	{
		/* This transition,1.549 , given by Lach and Pachucki, 2001 for the atom */
		if( nelem == ipHELIUM )
		{
			A = 1.549;
			iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
		}
		else
		{
			/* This is a fit to data given in
			 * >>refer	He-like	As	Savukov, I.M., Johnson, W.R., & Safronova, U.I. 
			 * >>refercon	astro-ph 0205163	*/
			A= 0.1834*pow((double)nelem, 6.5735);
			iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
		}
	}

	else if( nelem == ipHELIUM && inHi == inHe2p3P1 && inLo == inHe2p3P0 )
	{
		/* >>refer	He	As	Bulatov, N.N. 1976, Soviet Astronomy, 20, 315 */
		fixit("find a transition probability for this 2^3P0 - 2^3P1 transition.");
		/* This is the 29.6 GHz line that can be seen in radio galaxies. */
		/** \todo	2	find a transition probability for this 2^3P0 - 2^3P1 transition.
		 * It will require a bit of trickery to insert into the rate matrix, 
		 * because of the fact that the lower level has a higher index.  
		 * See discussion "Energy order within 2 3P" near the top of helike.c */
		A = 5.78E-12;
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);

	}

	else if( nelem == ipHELIUM && inHi == inHe2p3P2 && inLo == inHe2p3P1 )
	{
		fixit("find a transition probability for this 2^3P1 - 2^3P2 transition.");
		/* This is the 3 GHz line that can be seen in radio galaxies. */
		/** \todo	2	find a transition probability for this 2^3P1 - 2^3P2 transition.
		 * It will require a bit of trickery to insert into the rate matrix, 
		 * because of the fact that the lower level has a higher index.  
		 * See discussion "Energy order within 2 3P" near the top of helike.c */
		A = 3.61E-15;
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem == ipHELIUM && inHi == inHe3p3P && inLo == inHe1s1S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 44.326 * (1./3.) + 0.1146547 * (5./9.);
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem == ipHELIUM && inHi == inHe3p3P && inLo == inHe2s1S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 1.459495 * (1./3.) + 3.6558e-5 * (5./9.);
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem == ipHELIUM && inHi == inHe3p3P && inLo == inHe3s1S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 2.2297e-3 * (1./3.);
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem == ipHELIUM && inHi == inHe3p1P && inLo == inHe2s3S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 0.1815436;
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	else if( nelem == ipHELIUM && inHi == inHe3p1P && inLo == inHe3s3S )
	{
		// >>refer	He	As	Morton, D. \& Drake, G.W.F. 2011, PhysRevA, 83, 042503
		A = 0.14895852;
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	else if( ( inLo == inHe1s1S && lHi == 2 && sHi == 1 ) ||
		( inLo == inHe2s3S && lHi == 2 && sHi == 3 ) ||
		( inLo == inHe2s1S && lHi == 2 && sHi == 1 ) ||
		( nLo == 2 && lLo == 1 && sLo == 3 && nHi >= 3 && lHi == 1 && sHi == 3 ) ||
		( inLo == inHe2p1P && lHi == 1 && sHi == 1 ) )
	{
		// >>refer       Helike  QOS     Cann, N. M. \& Thakkar, A. J. 2002, J.Phys.B, 35, 421

		static const double f_params[5][4][3] = {
			{
		/* 1,0,3,2,1 */ {9.360591E-007, -3.1937E-006,  3.5186E-006},
		/* 1,0,4,2,1 */ {4.631435E-007, -1.4973E-006,  1.4848E-006},
		/* 1,0,5,2,1 */ {2.493912E-007, -7.8658E-007,  7.3994E-007},
		/* 1,0,6,2,1 */ {1.476742E-007, -4.5953E-007,  4.1932E-007}},
			{
		/* 2,0,3,2,3 */ {1.646733E-006, -2.0028E-006, -1.3552E-006},
		/* 2,0,4,2,3 */ {9.120593E-008,  3.1301E-007, -3.2190E-007},
		/* 2,0,5,2,3 */ {1.360965E-008,  1.1467E-007,  8.6977E-008},
		/* 2,0,6,2,3 */ {3.199421E-009,  4.5485E-008,  1.1016E-007}},
			{
		/* 2,0,3,2,1 */ {1.646733E-006, -2.9720E-006,  1.5367E-006},
		/* 2,0,4,2,1 */ {9.120593E-008, -3.9037E-008,  3.9156E-008},
		/* 2,0,5,2,1 */ {1.360965E-008,  1.4671E-008,  1.5626E-008},
		/* 2,0,6,2,1 */ {3.199421E-009,  8.9806E-009,  1.2436E-008}},
			{
		/* 2,1,3,1,3 */ {1.543812E-007, -2.8220E-007,  3.0261E-008},
		/* 2,1,4,1,3 */ {3.648237E-008, -6.6824E-008,  4.5879E-009},
		/* 2,1,5,1,3 */ {1.488556E-008, -2.7304E-008,  1.6628E-009},
		/* 2,1,6,1,3 */ {7.678610E-009, -1.4112E-008,  6.8453E-010}},
			{
		/* 2,1,3,1,1 */ {1.543812E-007, -2.8546E-007,  4.6605E-008},
		/* 2,1,4,1,1 */ {3.648237E-008, -6.8422E-008,  1.7038E-008},
		/* 2,1,5,1,1 */ {1.488556E-008, -2.8170E-008,  8.5012E-009},
		/* 2,1,6,1,1 */ {7.678610E-009, -1.4578E-008,  4.6686E-009}}
				};

		ASSERT( nLo <= 2 );
		long index_lo;
		if( inLo == inHe2p1P )
			index_lo = 4;
		else if( inLo == inHe2p3P0 || inLo == inHe2p3P1 || inLo == inHe2p3P2 )
			index_lo = 3;
		else if( inLo == inHe2s1S )
			index_lo = 2;
		else if( inLo == inHe2s3S )
			index_lo = 1;
		else if( inLo == inHe1s1S )
			index_lo = 0;
		else
			TotalInsanity();

		ASSERT( nHi >= 3 );
		long index_hi = MIN2( nHi, 6 ) - 3;
		double f_lu = POW2(nelem+1) * ( 
			f_params[index_lo][index_hi][0] + 
			f_params[index_lo][index_hi][1]/(nelem+1) + 
			f_params[index_lo][index_hi][2]/POW2(nelem+1) );

		// extrapolate for higher upper n
		if( nHi > 6 )
			f_lu *= pow3( 6. / nHi );

		double gLo = ( jLo < 0 ) ? double(sLo*(2*lLo+1)) : double(2*jLo+1);
		double gHi = ( jHi < 0 ) ? double(sHi*(2*lHi+1)) : double(2*jHi+1);
		// helike_energy() gives energy needed to ionize state
		double Enerwn = fabs(helike_energy(nelem, nLo, lLo, sLo, jLo) -
				     helike_energy(nelem, nHi, lHi, sHi, jHi));

		A = eina( gLo * f_lu, Enerwn, gHi );
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}
	
	else
	{
		/* Current transition is not supported.	*/
		A = iso_ctrl.SmallA;
		iso_put_error(ipHE_LIKE,nelem,inHi,inLo,IPRAD,0.1f,0.1f);
	}

	ASSERT( A > 0.);

	return A;   
}

/* Calculates Einstein A for a given transition.	*/
double he_1trans( 
		/* charge on c scale, Energy is wavenumbers, Einstein A	*/
		long nelem,
		/* quantum numbers of upper level:	*/
		double Eff_nupper, long lHi, long sHi, long jHi,
		/* and of lower level: */
		double Eff_nlower, long lLo, long sLo, long jLo)
{
	DEBUG_ENTRY( "he_1trans()" );

	ASSERT(nelem > ipHYDROGEN);

	/* Since 0.4 is bigger than any defect, adding that to the effective principle quantum number,
	 * and truncating to an integer will produce the principal quantum number.	*/
	long nHi = (int)(Eff_nupper + 0.4);
	long nLo = (int)(Eff_nlower + 0.4);

	/* Make sure this worked correctly.	*/
	ASSERT( fabs(Eff_nupper-(double)nHi) < 0.4 );
	ASSERT( fabs(Eff_nlower-(double)nLo) < 0.4 );

	/* check the validity of the J value */
	long gHi = -1, gLo = -1;
	if( jHi >= 0 )
	{
		gHi = 2*jHi + 1;
		ASSERT( Triangle2(sHi-1, 2*lHi, gHi-1) );
	}
	if( jLo >= 0 )
	{
		gLo = 2*jLo + 1;
		ASSERT( Triangle2(sLo-1, 2*lLo, gLo-1) );
	}

	QNPair inPair(nHi, lHi, sHi, gHi, nLo, lLo, sLo, gLo);

	auto p = iso_sp[ipHE_LIKE][nelem].CachedAs.find(inPair);
	if( p != iso_sp[ipHE_LIKE][nelem].CachedAs.end() )
		return p->second;

	/* First do allowed transitions	*/
	double Aul;
	if( (sHi == sLo) && (abs((int)(lHi - lLo)) == 1) )
	{
		Aul = -2.;

		/* For clarity, let's do this in two separate chunks...one for helium, one for everything else.	*/
		if( nelem == ipHELIUM )
		{
			/* Retrieve transition probabilities for Helium.	*/
			/* >>refer He	As	Drake, G.W.F., Atomic, Molecular, and Optical Physics Handbook */
			auto p = TransProbs.find(inPair);
			if( p != TransProbs.end() )
				Aul = p->second[nelem];

			if( Aul < 0. )
			{
				/* Here are the Lyman transitions.	*/
				if( nLo==1 )
				{
					ASSERT( (lHi == 1) && (sHi == 1) );

					/* these fits calculated from Drake A's (1996) */
					Aul = 1.59208e10 / pow3(Eff_nupper);
					ASSERT( Aul > 0.);
				}

				/* last resort for transitions involving significant defects, 
				 * except that highest lLo are excluded */
				else if( lHi>=2 && lLo>=2 && nHi>nLo )
				{
					Aul = H_Einstein_A(nHi, lHi, nLo, lLo, nelem);
					ASSERT( Aul >= 0. );
				}
				/* These fits come from extrapolations of Drake's oscillator strengths
				 * out to the series limit.  We also use this method to obtain threshold
				 * photoionization cross-sections for the lower level of each transition here. */
				/* See these two references :
				 * >>refer	He	As	Hummer, D. G. \& Storey, P. J. 1998, MNRAS 297, 1073 
				 * >>refer	Seaton's Theorem	Seaton, M. J. 1958, MNRAS 118, 504 */
				else if( nHi>10 && nLo<=5 && lHi<=2 && lLo<=2 )
				{
					int paramSet=0;
					double emisOscStr, x, a, b, c;
					static const double extrapol_Params[2][4][4][3] = {
						/* these are for singlets */
						{
							{ /* these are P to S */
								{ 0.8267396, 1.4837624, -0.4615955 },
								{ 1.2738405, 1.5841806, -0.3022984 },
								{ 1.6128996, 1.6842538, -0.2393057 },
								{ 1.8855491, 1.7709125, -0.2115213 }},
							{ /* these are S to P */
								{ -1.4293664, 2.3294080, -0.0890470 },
								{ -0.3608082, 2.3337636, -0.0712380 },
								{  0.3027974, 2.3326252, -0.0579008 },
								{  0.7841193, 2.3320138, -0.0497094 }},
							{ /* these are D to P */
								{ 1.1341403, 3.1702435, -0.2085843 },
								{ 1.7915926, 2.4942946, -0.2266493 },
								{ 2.1979400, 2.2785377, -0.1518743 },
								{ 2.5018229, 2.1925720, -0.1081966 }},
							{ /* these are P to D */
								{  0.0000000, 0.0000000,  0.0000000 },
								{ -2.6737396, 2.9379143, -0.3805367 },
								{ -1.4380124, 2.7756396, -0.2754625 },
								{ -0.6630196, 2.6887253, -0.2216493 }},
						},
						/* these are for triplets */
						{
							{ /* these are P to S */
								{ 0.3075287, 0.9087130, -1.0387207 },
								{ 0.6870690, 1.1485864, -0.6627317 },
								{ 0.9776064, 1.3382024, -0.5331906 },
								{ 1.2107725, 1.4943721, -0.4779232 }},
							{ /* these are S to P */
								{ -1.3659605, 2.3262253, -0.0306439 },
								{ -0.2899490, 2.3279391, -0.0298695 },
								{  0.3678878, 2.3266603, -0.0240021 },
								{  0.8427457, 2.3249540, -0.0194091 }},
							{ /* these are D to P */
								{ 1.3108281, 2.8446367, -0.1649923 },
								{ 1.8437692, 2.2399326, -0.2583398 },
								{ 2.1820792, 2.0693762, -0.1864091 },
								{ 2.4414052, 2.0168255, -0.1426083 }},
							{ /* these are P to D */
								{  0.0000000, 0.0000000,  0.0000000 },
								{ -1.9219877, 2.7689624, -0.2536072 },
								{ -0.7818065, 2.6595150, -0.1895313 },
								{ -0.0665624, 2.5955623, -0.1522616 }},
						}
					};

					if( lLo==0 )
					{
						paramSet = 0;
					}
					else if( lLo==1 && lHi==0 )
					{
						paramSet = 1;
					}
					else if( lLo==1 && lHi==2 )
					{
						paramSet = 2;
					}
					else if( lLo==2 )
					{
						paramSet = 3;
						ASSERT( lHi==1 );
					}

					ASSERT( (int)((sHi-1)/2) == 0 || (int)((sHi-1)/2) == 1 );
					a = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][0];
					b = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][1];
					c = extrapol_Params[(int)((sHi-1)/2)][paramSet][nLo-2][2];
					double LevelLoWN = helike_energy(nelem, nLo, lLo, sLo, jLo);
					double LevelHiWN = helike_energy(nelem, nHi, lHi, sHi, jHi);
					// helike_energy() gives energy needed to ionize state,
					double Enerwn = fabs(LevelLoWN - LevelHiWN);
					x = log( LevelLoWN/Enerwn );

					emisOscStr = exp(a+b*x+c*x*x)/pow3(Eff_nupper)*
						(2.*lLo+1)/(2.*lHi+1);

					Aul = TRANS_PROB_CONST*Enerwn*Enerwn*emisOscStr;

					if( nLo==2 && lLo==1 && sLo==3 )
						Aul *= double(gLo)/9.0;
					ASSERT( Aul > 0. );
				}
				else
				{
					// helike_energy() gives energy needed to ionize state,
					double Enerwn = fabs(helike_energy(nelem, nLo, lLo, sLo, jLo) -
							     helike_energy(nelem, nHi, lHi, sHi, jHi));
					/* Calculate the radial integral from the quantum defects.	*/
					double RI2 = scqdri(Eff_nupper,lHi,Eff_nlower,lLo,(double)(ipHELIUM));
					ASSERT( RI2 > 0. );
					/* Convert radial integral to Aul.	*/
					Aul = ritoa(lHi,lLo,ipHELIUM,Enerwn,RI2);
					/* radial integral routine does not recognize fine structure.
					 * Here we split 2^3P.	*/
					if( nLo==2 && lLo==1 && sLo==3 )
						Aul *= double(gLo)/9.0;
					ASSERT( Aul > 0. );
				}
			}
		}

		/* Heavier species	*/
		else
		{
			/* Retrieve transition probabilities for Helike ions.	*/
			/* >>refer He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
			 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J, originally astro.ph. 0201454 */
			auto p = TransProbs.find(inPair);
			if( p != TransProbs.end() )
				Aul = p->second[nelem];

			if( Aul < 0. )
			{
				/* Do same-n transitions. */
				if( nLo==nHi && lHi<=2 && lLo<=2 )
				{
					/* These are 2p3Pj to 2s3S fits to (low Z) Porquet & Dubau (2000) & 
					 * (high Z) NIST Atomic Spectra Database.	*/
					if( nLo==2 && lLo==0 && sLo==3 )
					{
						ASSERT( nHi==2 && lHi==1 && sHi==3 );

						if( jHi==0 )
							Aul = 3.31E7 + 1.13E6 * pow((double)nelem+1.,1.76);
						else if( jHi==1 )
							Aul = 2.73E7 + 1.31E6 * pow((double)nelem+1.,1.76);
						else if( jHi==2 )
							Aul = 3.68E7 + 1.04E7 * exp(((double)nelem+1.)/5.29);
						else
						{
							fprintf(ioQQQ," PROBLEM Impossible value for jHi in he_1trans.\n");
							TotalInsanity();
						}
					}

					/* These are 2p1P to 2s1S fits to data from TOPbase.	*/
					else if( ( nLo==2 && lLo==0 && sLo==1 ) && ( nHi==2 && lHi==1 && sHi==1 ) )
					{
						Aul = 5.53E6 * exp( 0.171*(nelem+1.) );
					}

					else 
					{
						/* This case should only be entered if n > 2.  Those cases were done above.	*/
						ASSERT( nLo > 2 );

						/* Triplet P to triplet S, delta n = 0	*/
						if( (lHi == 1) && (sHi == 3) && (lLo == 0) && (sLo == 3))
						{
							Aul = 0.4 * 3.85E8 * pow((double)nelem,1.6)/pow((double)nHi,5.328);
						}
						/* Singlet P to singlet D, delta n = 0	*/
						else if( (lHi == 1) && (sHi == 1) && (lLo == 2) && (sLo == 1))
						{
							Aul = 1.95E4 * pow((double)nelem,1.6) / pow((double)nHi, 4.269);
						}
						/* Singlet P to singlet S, delta n = 0	*/
						else if( (lHi == 1) && (sHi == 1) && (lLo == 0) )
						{
							Aul = 6.646E7 * powpq((double)nelem,3,2) / pow((double)nHi, 5.077);
						}
						else 
						{
							ASSERT( (lHi == 2) && (sHi == 3)  && (lLo == 1) );
							Aul = 3.9E6 * pow((double)nelem,1.6) / pow((double)nHi, 4.9);
							if( (lHi >2) || (lLo > 2) )
								Aul *= (lHi/2.);
							if(lLo > 2)
								Aul *= (1./9.);
						}
					}
					ASSERT( Aul > 0.);
				}

				/* assume transitions involving F and higher orbitals are hydrogenic.	*/
				else if( (nHi > nLo) && ((lHi > 2) || (lLo > 2)) )
				{
					Aul = H_Einstein_A(nHi ,lHi , nLo , lLo , nelem);
					ASSERT( Aul > 0.);
				}

				/* These transitions are of great importance, but the below radial integral 
				 * routine fails to achieve desirable accuracy, so these are fits as produced 
				 * from He A's for nupper through 9.  They are transitions to ground and 
				 * 2, 3, and 4 triplet S.	*/
				else if( ( nLo==1 ) || ( nLo==2 && lLo==0 ) ||
					( nLo==3 && lLo==0 && sLo==3 ) ||
					( nLo==4 && lLo==0 && sLo==3 ) )
				{
					/* Here are the Lyman transitions.	*/
					if( nLo==1 )
					{
						ASSERT( (lHi == 1) && (sHi == 1) );

						/* In theory, this Z dependence should be Z^4, but values from TOPbase 
						 * suggest 3.9 is a more accurate exponent.	Values from 
						 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
						 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
						/* originally astro.ph. 0201454  */
						Aul = 1.375E10 * pow((double)nelem, 3.9) / pow((double)nHi,3.1);
					}

					/* Here are the Balmer transitions.	*/
					else if( nLo==2 && lLo==0 && sLo==1 )
					{
						ASSERT( (lHi == 1) && (sHi == 1) );

						/* More fits, as above. */
						Aul = 5.0e8 * pow4((double)nelem) / pow((double)nHi, 2.889);
					}

					/* Here are transitions down to triplet S	*/
					else
					{
						ASSERT( lLo == 0 && sLo == 3 );
						ASSERT( lHi == 1 && sHi == 3 );

						/* These fits also as above. */
						if( nLo == 2 )
							Aul = 1.5 * 3.405E8 * pow4((double)nelem) / pow((double)nHi, 2.883);
						else if( nLo == 3 )
							Aul = 2.5 * 4.613E7 * pow4((double)nelem) / pow((double)nHi, 2.672);
						else 
							Aul = 3.0 * 1.436E7 * pow4((double)nelem) / pow((double)nHi, 2.617);
					}
					ASSERT( Aul > 0.);
				}

				/* Every other allowed transition is calculated as follows.	*/
				else
				{
					// helike_energy() gives energy needed to ionize state,
					double Enerwn = fabs(helike_energy(nelem, nLo, lLo, sLo, jLo) -
							     helike_energy(nelem, nHi, lHi, sHi, jHi));
					/* Calculate the radial integral from the quantum defects.	*/
					double RI2 = scqdri(Eff_nupper,lHi,Eff_nlower,lLo,(double)(nelem));
					/* Convert radial integral to Aul.	*/
					Aul = ritoa(lHi,lLo,nelem,Enerwn,RI2);
					/* radial integral routine does not recognize fine structure.
					 * Here we split 2^3P.	*/
					if( nLo==2 && lLo==1 && sLo==3 && Aul > iso_ctrl.SmallA )
						Aul *= double(gLo)/9.0;
				}
			}
		}
	}

	/* Now do forbidden transitions from 2-1 ... */
	/* and those given by  
	 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
	 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
	/* originally astro.ph. 0201454  
	 * for heavy elements. These are triplet P to singlet S, 
	 * going either up or down...Triplet S to Singlet P are not included, as they are far weaker.	*/
	else 
	{
		ASSERT( (sHi != sLo) || (abs((int)(lHi - lLo)) != 1) );
		Aul = ForbiddenAuls(nelem, nHi, lHi, sHi, jHi, nLo, lLo, sLo, jLo);
		ASSERT( Aul > 0. );
	}

	Aul = MAX2( Aul, iso_ctrl.SmallA );

	iso_sp[ipHE_LIKE][nelem].CachedAs[inPair] = Aul;
	return Aul;
}

STATIC double helike_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo )
{
	DEBUG_ENTRY( "helike_transprob_collapsed_to_collapsed()" );

	double Aul = 0.;
	ASSERT( nHi > nLo );

	for( long lLo=0; lLo < nLo; ++lLo )
	{
		double Aul_sing = helike_transprob_collapsed_to_resolved( nelem, nHi, nLo, lLo, 1, -1 );
		double Aul_trip = helike_transprob_collapsed_to_resolved( nelem, nHi, nLo, lLo, 3, -1 );
		Aul += Aul_sing + Aul_trip;
	}

	return Aul;
}

STATIC double helike_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo, long sLo, long jLo )
{
	DEBUG_ENTRY( "helike_transprob_collapsed_to_resolved()" );

	double n_eff_hi = nHi - helike_quantum_defect(nelem,nHi,-1,-1,-1);
	double n_eff_lo = nLo - helike_quantum_defect(nelem,nLo,lLo,sLo,jLo);

	// make sure 2^3Pj terms are fully specified
	if( nLo==2 && lLo==1 && sLo==3 )
		ASSERT( jLo>=0 && jLo<=2 );
	
	/* Lower level resolved, upper not. First calculate Aul
	 * from upper level with ang mom one higher.	*/
	double Aul = he_1trans( nelem, n_eff_hi, lLo+1, sLo, -1, n_eff_lo, lLo, sLo, jLo);

	Aul *= (2.*(lLo+1.)+1.) * sLo / (4.*(double)nHi*(double)nHi);

	if( lLo != 0 )
	{
		/* For all l>0, add in transitions from upper level with ang mom one lower.	*/
		double Aul1 = he_1trans( nelem, n_eff_hi, lLo-1, sLo, -1, n_eff_lo, lLo, sLo, jLo);

		/* now add in this part after multiplying by stat weight for lHi = lLo-1.	*/
		Aul += Aul1*(2.*(lLo-1.)+1.) * sLo / (4.*(double)nHi*(double)nHi);
	}

	ASSERT( Aul > 0.);

	return Aul;
}

/*ritoa converts the square of the radial integral for a transition 
 * (calculated by scqdri) to the transition probability, Aul.	*/
STATIC double ritoa(long li, long lf, long nelem, double k, double RI2)
{
	/*  Variables are as follows:               */
	/*  lg = larger of li and lf                */
	/*  fmean = mean oscillator strength        */
	/*    for a given level.                    */
	/*  mu = reduced mass of optical electron.  */
	/*  EinsteinA = Einstein emission coef.     */
	/*  w = angular frequency of transition.    */
	/*  RI2_cm = square of rad. int. in cm^2.   */
	long lg;
	double fmean,mu,EinsteinA,w,RI2_cm;

	DEBUG_ENTRY( "ritoa()" );

	mu = ELECTRON_MASS/(1+ELECTRON_MASS/(dense.AtomicWeight[nelem]*ATOMIC_MASS_UNIT));

	w = 2. * PI * k * SPEEDLIGHT;

	RI2_cm = RI2 * BOHR_RADIUS_CM * BOHR_RADIUS_CM;

	lg = max(lf, li);

	fmean = 2.0*mu*w*lg*RI2_cm/((3.0*H_BAR) * (2.0*li+1.0));

	EinsteinA = TRANS_PROB_CONST*k*k*fmean;

	/* ASSERT( EinsteinA > SMALLFLOAT ); */

	return EinsteinA;
}

realnum helike_transprob( long nelem, long ipHi, long ipLo )
{
	DEBUG_ENTRY( "helike_transprob()" );

	long nHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].n();
	long lHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].l();
	long sHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].S();
	long jHi = iso_sp[ipHE_LIKE][nelem].st[ipHi].j();
	long nLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].n();
	long lLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].l();
	long sLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].S();
	long jLo = iso_sp[ipHE_LIKE][nelem].st[ipLo].j();

	double n_eff_hi = nHi - helike_quantum_defect(nelem, nHi, lHi, sHi, jHi);
	double n_eff_lo = nLo - helike_quantum_defect(nelem, nLo, lLo, sLo, jLo);

	realnum Aul;
	if( nHi > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
	{
		if( nLo > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
		{
			/* Neither upper nor lower is resolved...Aul()'s are easy.	*/
			Aul = helike_transprob_collapsed_to_collapsed(nelem, nHi, nLo);
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.001f,0.001f);
		}
		else 
		{
			Aul = helike_transprob_collapsed_to_resolved(nelem, nHi, nLo, lLo, sLo, jLo);
			iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.01f,0.01f);
		}
	}
	else
	{
		/* Both levels are resolved...do the normal bit.	*/
		Aul = he_1trans(nelem, n_eff_hi, lHi, sHi, jHi, n_eff_lo, lLo, sLo, jLo);
		iso_put_error(ipHE_LIKE,nelem,ipHi,ipLo,IPRAD,0.01f,0.01f);
	}

	ASSERT( Aul > 0._r );
	return Aul;
}

/* This routine is pure infrastructure and bookkeeping, no physics. */
void HelikeTransProbSetup()
{
	DEBUG_ENTRY( "HelikeTransProbSetup()" );

	/********************************************************************/
	/*************** Read in data from he_transprob.dat	*****************/

	DataParser d("he_transprob.dat", ES_NONE);

	d.getline();
	d.checkMagic(TRANSPROBMAGIC, N_HE1_TRANS_PROB);

	for( long i=0; i < N_HE1_TRANS_PROB; i++ )
	{
		d.getline();
		long i1, i2;
		d.getToken(i1);
		d.getToken(i2);
		if( i1 < 0 || i2 <= i1 )
			d.errorAbort("invalid level indices");

		long nLo, lLo, sLo, jLo, nHi, lHi, sHi, jHi;
		d.getToken(nLo);
		d.getToken(lLo);
		d.getToken(sLo);
		d.getToken(jLo);
		d.getToken(nHi);
		d.getToken(lHi);
		d.getToken(sHi);
		d.getToken(jHi);
		QNPair inPair(nHi, lHi, sHi, 2*jHi+1, nLo, lLo, sLo, 2*jLo+1);

		double* p = TransProbs[inPair].data();
		d.getToken(p+1, LIMELM-1); // start reading data for helium
		d.checkEOL();
	}
	d.checkEOD();
}
