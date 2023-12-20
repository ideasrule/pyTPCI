/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HCSAR_interp interpolate on collision strengths */
/*C6cs123 line collision rates for lower levels of hydrogenic carbon, n=1,2,3 */
/*Ca20cs123 line collision rates for lower levels of hydrogenic calcium, n=1,2,3 */
/*Hydcs123 Hydrogenic de-excitation collision rates n=1,2,3 */
/*H1cs123 hydrogen collision data levels involving 1s,2s,2p,3. */
/*Ne10cs123 line collision rates for lower levels of hydrogenic neon, n=1,2,3 */
/*He2cs123 line collision strengths for lower levels of helium ion, n=1,2,3, by K Korista */
/*Fe26cs123 line collision rates for lower levels of hydrogenic iron, n=1,2,3 */
#include "cddefines.h"
#include "atmdat.h"
#include "atmdat_adfa.h"
#include "helike_cs.h"
#include "hydrogenic.h"
#include "hydro_vs_rates.h"
#include "iso.h"
#include "opacity.h"
#include "phycon.h"
#include "thirdparty.h"
#include "integrate.h"
#include "freebound.h"
#include "hydroeinsta.h"
#include "dense.h"

STATIC double Fe26cs123(QNPack inLo, QNPack inHi);
STATIC double He2cs123(QNPack inLo, QNPack inHi);
STATIC double Hydcs123(long nelem, long nHi, long lHi, long sHi, long nLo, long lLo, long sLo, char chType);
STATIC double C6cs123(QNPack inLo, QNPack inHi);
STATIC double Ca20cs123(QNPack inLo, QNPack inHi);
STATIC double Ne10cs123(QNPack inLo, QNPack inHi);

STATIC realnum HCSAR_interp(long nHi, long lHi, long sHi, long nLo, long lLo, long sLo);
STATIC realnum HlikeCSInterp( long nelem, long Collider, long nHi, long lHi, long sHi, long nLo, long lLo, long sLo );

/* 1 is 1s, 2 is 2s, 3 is 2p, and 4 is 3s, 5 is 3p, and 6 is 3d */
static QNPack in1 = QN2ind(1, 0, 2, 2);
static QNPack in2 = QN2ind(2, 0, 2, 2);
static QNPack in3 = QN2ind(2, 1, 2, -1);
static QNPack in4 = QN2ind(3, 0, 2, 2);
static QNPack in5 = QN2ind(3, 1, 2, -1);
static QNPack in6 = QN2ind(3, 2, 2, -1);

static const realnum HCSTE[NHCSTE] = {5802.f,11604.f,34812.f,58020.f,116040.f,174060.f,232080.f,290100.f};

/*HCSAR_interp interpolate on collision strengths */
STATIC realnum HCSAR_interp(long nHi, long lHi, long sHi, long nLo, long lLo, long sLo)
{
	DEBUG_ENTRY( "HCSAR_interp()" );

	if( nLo == 2 && nHi == 2 )
	{
		fprintf(ioQQQ,"HCSAR_interp was called for the 2s-2p transition, which it cannot do\n");
		cdEXIT(EXIT_FAILURE);
	}
	QNPair inPair(nHi, lHi, sHi, -1, nLo, lLo, sLo, -1);
	const realnum* csarr = t_ADfA::Inst().h_coll_str(nHi, lHi, nLo, lLo);
	realnum cs;
	/* this is to linearly extrapolate to cs(T=0)=0 for low energies*/
	if (phycon.te < HCSTE[0])
	{
		realnum lowint[2] = {0.f,HCSTE[0]};
		realnum lowarr[2] = {0.f,csarr[0]};
		cs = linint(lowint,lowarr,2,realnum(phycon.te));
	}
	else
		cs = linint(HCSTE, csarr, NHCSTE, realnum(phycon.te));

	ASSERT( cs >= 0._r );
	return cs;
}

/*Hydcs123 Hydrogenic de-excitation collision strengths between levels n=1,2,3,
 * for any charge.  routine only called by HydroCSInterp to fill in hydroline arrays
 * with collision strengths */
STATIC double Hydcs123(
	 /* charge, 0 for hydrogen, 1 for helium, etc */
	 long nelem,
	 /* quantum numbers for upper state */
	 long nHi, long lHi, long sHi,
	 /* quantum numbers for lower state */
	 long nLo, long lLo, long sLo,
	 /* = 'e' for electron collisions, ='p' for proton */
	 char chType)
{
	long int i,
	  j, 
	  k;
	double C, 
	  D, 
	  EE, 
	  expq ,
	  Hydcs123_v, 
	  logx,
	  Ratehigh, 
	  Ratelow, 
	  TeUse, 
	  gLo, 
	  gHi, 
	  q, 
	  rate, 
	  slope, 
	  temp, 
	  temphigh, 
	  templow, 
	  tev, 
	  x, 
	  QuanNLo, 
	  QuanNUp, 
	  Charge, 
	  ChargeSquared, 
	  zhigh, 
	  zlow;

	static const double ap[5] = {-2113.113,729.0084,1055.397,854.632,938.9912};
	static const double bp[5] = {-6783.515,-377.7190,724.1936,493.1107,735.7466}; 
	static const double cp[5] = {-3049.719,226.2320,637.8630,388.5465,554.6369};
	static const double dp[5] = {3514.5153,88.60169,-470.4055,-329.4914,-450.8459};
	static const double ep[5] = {0.005251557,0.009059154,0.008725781,0.009952418,0.01098687};
	static const double ae[5] = {-767.5859,-643.1189,-461.6836,-429.0543,-406.5285};
	static const double be[5] = {-1731.9178,-1442.548,-1055.364,-980.3079,-930.9266};
	static const double ce[5] = {-939.1834,-789.9569,-569.1451,-530.1974,-502.0939};
	static const double de[5] = {927.4773,773.2008,564.3272,524.2944,497.7763};
	static const double ee[5] = {-0.002528027,-0.003793665,-0.002122103,-0.002234207,-0.002317720};
	static const double A[2] = {4.4394,0.0};
	static const double B[2] = {0.8949,0.8879};
	static const double C0[2] = {-0.6012,-0.2474};
	static const double C1[2] = {-3.9710,-3.7562};
	static const double C2[2] = {-4.2176,2.0491};
	static const double D0[2] = {2.930,0.0539};
	static const double D1[2] = {1.7990,3.4009};
	static const double D2[2] = {4.9347,-1.7770};

	QNPack inLo = QN2ind(nLo, lLo, sLo, -1);
	QNPack inHi = QN2ind(nHi, lHi, sHi, -1);

	DEBUG_ENTRY( "Hydcs123()" );
	/* Hydrogenic de-excitation collision rates n=1,2,3 
	 * >>refer	h1	cs	Callaway, J. 1983, Phys Let A, 96, 83
	 * >>refer	h1	cs	Zygelman, B., & Dalgarno, A. 1987, Phys Rev A, 35, 4085 
	 * for 2p-2s only.
	 * The fit from Callaway is in nuclear charge for 1s - 2s,2p only.
	 * For transtions involving level 3, interpolation in Z involving
	 * the functions He2cs123,C6cs123,Ne10cs123,Ca20cs123, Fe26cs123.
	 *
	 * The fits from ZD are for 2p-2s for Z=2,6,12,16,18 only other charges are
	 * interpolated, both electron and proton rates are included,
	 * the variable chType is either 'e' or 'p'.
	 *
	 * the lower level runs from in1 to in3 (1s, 2s, 2p)
	 * the upper level runs from in2 to in6 (2s, 2p, 3s, 3p, 3d) */

	/* for Callaway fit: */
	/* for Zygelman and Dalgarno: */

	/* first entry is 2p, then 2s */

	/* fit in nuclear charge Z for 2p-2s collisions in hydrogenic species
	 * equation is a+bx+cx^2ln(x)+dexp(x)+eln(x)/x^2, where x=te/Z^2 in a.u.
	 * first are the proton rates: */
	/* these are electron rates: */

	/* following is for charged species */
	/* charge is on scale with He+=1, Li++=2, etc */
	ASSERT( nelem > ipHYDROGEN && nelem < LIMELM );

	ASSERT( nLo > 0 && nLo <= 2 );
	ASSERT( nHi > 1 && nHi <= 3 );

	/* set quantum numbers and stat. weights of the transitions: */
	if( inHi == in6 )
	{
		/* upper is n=3 then set level, stat. weight */
		QuanNUp = 3.;
		gHi = 10.;
		/* following will be set here even though it is not used in this case,
		 * to prevent good compilers from falsing on i not set,
		 * there is assert when used to make sure it is ok */
		i = -1;
	}
	else if( inHi == in5 )
	{
		/* upper is n=3 then set level, stat. weight */
		QuanNUp = 3.;
		gHi = 6.;
		/* following will be set here even though it is not used in this case,
		 * to prevent good compilers from falsing on i not set,
		 * there is assert when used to make sure it is ok */
		i = -1;
	}
	else if( inHi == in4 )
	{
		/* upper is n=3 then set level, stat. weight */
		QuanNUp = 3.;
		gHi = 2.;
		/* following will be set here even though it is not used in this case,
		 * to prevent good compilers from falsing on i not set,
		 * there is assert when used to make sure it is ok */
		i = -1;
	}
	else if( inHi == in3 )
	{
		/* upper is nl=2p then set level, stat. weight */
		QuanNUp = 2.;
		gHi = 6.;
		/* used to point within vectors defined above */
		i = 0;
	}
	else if( inHi == in2 )
	{
		/* upper is nl=2s then set level, stat. weight */
		QuanNUp = 2.;
		gHi = 2.;
		/* used to point within vectors defined above */
		i = 1;
	}
	else
	{
		fprintf( ioQQQ, " Insane levels in Hydcs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* which lower level? */
	if( inLo == in1 )
	{
		/* lower is n=1 then set level, stat. weight */
		QuanNLo = 1.;
		gLo = 2.;
	}
	else if( inLo == in2 )
	{
		/* lower is nl=2s then set level, stat. weight */
		QuanNLo = 2.;
		gLo = 2.;
	}
	else if( inLo == in3 )
	{
		QuanNLo = 2.;
		gLo = 6.;
	}
	else
	{
		fprintf( ioQQQ, " Insane levels in Hydcs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* following is the physical charge */
	Charge = (double)(nelem + 1);
	/* square of charge */
	ChargeSquared = Charge*Charge;

	if( inLo == in2 && inHi == in3 )
	{
		/*************************************** this is 2p-2s:
		 * series of if statements determines which entries Charge is between. */
		if( nelem < 1 )
		{
			/* this can't happen since routine returned above when ip=1 and 
			 * special atomic hydrogen routine called */
			fprintf( ioQQQ, " insane charge given to Hydcs123\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else if( nelem == 1 )
		{
			zlow = 2.;
			j = 1;
			zhigh = 2.;
			k = 1;
		}
		/* Li through C */
		else if( nelem <= 5 )
		{
			zlow = 2.;
			j = 1;
			zhigh = 6.;
			k = 2;
		}
		else if( nelem <= 11 )
		{
			zlow = 6.;
			j = 2;
			zhigh = 12.;
			k = 3;
		}
		else if( nelem <= 15 )
		{
			zlow = 12.;
			j = 3;
			zhigh = 16.;
			k = 4;
		}
		else if( nelem <= 17 )
		{
			zlow = 16.;
			j = 4;
			zhigh = 18.;
			k = 5;
		}
		/* following changed to else from else if, 
		 * to prevent false comment in good compilers */
		/*else if( nelem > 18 )*/
		else 
		{
			zlow = 18.;
			j = 5;
			zhigh = 18.;
			k = 5;
		}

		/* convert Te to a.u./Z^2
		 * determine rate at the low Charge */
		x = EVRYD/TE1RYD*phycon.te/(27.211396*pow2(zlow));
		TeUse = MIN2(x,0.80);
		x = MAX2(0.025,TeUse);
		logx = log(x);

		/* what type of collision are we dealing with? */
		if( chType == 'e' )
		{
			/* electron collisions */
			Ratelow = ae[j-1] + be[j-1]*x + ce[j-1]*pow2(x)*logx + de[j-1]*
			  exp(x) + ee[j-1]*logx/pow2(x);
		}
		else if( chType == 'p' )
		{
			Ratelow = ap[j-1] + bp[j-1]*x + cp[j-1]*pow2(x)*logx + dp[j-1]*
			  exp(x) + ep[j-1]*logx/pow2(x);
		}
		else 
		{
			/* this can't happen */
			fprintf( ioQQQ, " insane collision species given to Hydcs123\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( fp_equal( zlow, zhigh ) )
		{
			rate = Ratelow;
		}
		else
		{
			/* determine rate at the high Charge */
			x = EVRYD/TE1RYD*phycon.te/(27.211396*pow2(zhigh));
			TeUse = MIN2(x,0.80);
			x = MAX2(0.025,TeUse);
			logx = log(x);
			if( chType == 'e' )
			{
				Ratehigh = ae[k-1] + be[k-1]*x + ce[k-1]*pow2(x)*logx + 
					de[k-1]*exp(x) + ee[k-1]*logx/pow2(x);
			}
			else
			{
				Ratehigh = ap[k-1] + bp[k-1]*x + cp[k-1]*pow2(x)*logx + 
					dp[k-1]*exp(x) + ep[k-1]*logx/pow2(x);
			}
			/* linearly interpolate in charge */
			slope = (Ratehigh - Ratelow)/(zhigh - zlow);
			rate = slope*(Charge - zlow) + Ratelow;
		}
		rate = rate/ChargeSquared/Charge*1.0e-7;
		/* must convert to cs and need to know the valid temp range */
		templow = 0.025*27.211396*TE1RYD/EVRYD*ChargeSquared;
		temphigh = 0.80*27.211396*TE1RYD/EVRYD*ChargeSquared;
		TeUse = MIN2((double)phycon.te,temphigh);
		temp = MAX2(TeUse,templow);
		Hydcs123_v = rate*gHi*sqrt(temp)/COLL_CONST;

		if( chType == 'p' )
		{
			/* COLL_CONST is incorrect for protons, correct here */
			Hydcs123_v *= powpq( PROTON_MASS/ELECTRON_MASS, 3, 2 );
		}
	}
	else if( inHi == in4 || inHi == in5 || inHi == in6 )
	{
		/* n = 3
		 * for the rates involving n = 3, must do something different. */
		if( nelem < 1 )
		{
			fprintf( ioQQQ, " insane charge given to Hydcs123\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else if( nelem == 1 )
		{
			zlow = 2.;
			Ratelow = He2cs123(inLo,inHi);
			zhigh = 2.;
			Ratehigh = Ratelow;
		}
		else if( nelem > 1 && nelem <= 5 )
		{
			zlow = 2.;
			Ratelow = He2cs123(inLo,inHi);
			zhigh = 6.;
			Ratehigh = C6cs123(inLo,inHi);
		}
		else if( nelem > 5 && nelem <= 9 )
		{
			zlow = 6.;
			Ratelow = C6cs123(inLo,inHi);
			zhigh = 10.;
			Ratehigh = Ne10cs123(inLo,inHi);
		}
		else if( nelem > 9 && nelem <= 19 )
		{
			zlow = 10.;
			Ratelow = Ne10cs123(inLo,inHi);
			zhigh = 20.;
			Ratehigh = Ca20cs123(inLo,inHi);
		}
		else if( nelem > 19 && nelem <= 25 )
		{
			zlow = 20.;
			Ratelow = Ca20cs123(inLo,inHi);
			zhigh = 26.;
			Ratehigh = Fe26cs123(inLo,inHi);
		}
		/*>>>chng 98 dec 17, to else to stop comment from good compilers*/
		/*else if( nelem > 26 )*/
		else
		{
			Charge = 26.;
			zlow = 26.;
			Ratelow = Fe26cs123(inLo,inHi);
			zhigh = 26.;
			Ratehigh = Ratelow;
		}

		/* linearly interpolate */
		if( fp_equal( zlow, zhigh ) )
		{
			rate = Ratelow;
		}
		else
		{
			slope = (Ratehigh - Ratelow)/(zhigh - zlow);
			rate = slope*(Charge - zlow) + Ratelow;
		}

		/** NB - all of these actually calculate EIE collision strengths */
		Hydcs123_v = rate;

		/* the routines C6cs123, Ne10cs123, etc... do not resolve L for n>2 */
		/* dividing by N should roughly recover the original l-resolved data */
		Hydcs123_v /= 3.;
	}
	else
	{
		/* this branch 1-2s, 1-2p */
		if( nelem == 1 )
		{
			/* this brance for helium, then return */
			Hydcs123_v = He2cs123(inLo,inHi);
			return Hydcs123_v;
		}

		/* electron temperature in eV */
		tev = phycon.te / EVDEGK;
		/* energy in eV for hydrogenic species and these quantum numbers */
		EE = ChargeSquared*EVRYD*(1./QuanNLo/QuanNLo - 1./QuanNUp/QuanNUp);
		/* EE/kT for this transion */
		q = EE/tev;
		TeUse = MIN2(q,10.);
		/* q is now EE/kT but between 1 and 10 */
		q = MAX2(1.,TeUse);
		expq = exp(q);

		/* i must be 0 or 1 */
		ASSERT( i==0 || i==1 );
		C = C0[i] + C1[i]/Charge + C2[i]/ChargeSquared;
		D = D0[i] + D1[i]/Charge + D2[i]/ChargeSquared;

		/* following code changed so that e1 always returns e1,
		 * orifinal version only returned e1 for x < 1 */
		/* use disabled e1: */
		/*if( q < 1. )*/
		/*{*/
			/*rate = (B[i-1] + D*q)*exp(-q) + (A[i-1] + C*q - D*q*q)**/
			  /*e1(q);*/
		/*}*/
		/*else*/
		/*{*/
			/*rate = (B[i-1] + D*q) + (A[i-1] + C*q - D*q*q)*e1(q);*/
		/*}*/
		/*rate *= 8.010e-8/2./ChargeSquared/tev*sqrt(tev);*/
		/* convert to de-excitation */
		/*if( q < 1. )*/
		/*{*/
			/*rate = rate*exp(q)*gLo/gHi;*/
		/*}*/
		/*else*/
		/*{*/
			/*rate = rate*gLo/gHi;*/
		/*}*/

		/*>>>chng 98 dec 17, e1 always returns e1 */
		rate = (B[i] + D*q)/expq + (A[i] + C*q - D*q*q)*
			  e1(q);
		// Equation (2) of Callaway 1983 -- can simplify by cancelling tev factors
		// with R_0(T)
		//rate *= 8.010e-8/2./ChargeSquared/tev*sqrt(tev);
		rate *= 8.010e-8/(2. * ChargeSquared * sqrt(tev));
		/* convert to de-excitation */
		rate *= expq*gLo/gHi;

		/* convert to cs */
		Hydcs123_v = rate*gHi*phycon.sqrte/COLL_CONST;
	}
	return Hydcs123_v;
}

/*C6cs123 line collision rates for lower levels of hydrogenic carbon, n=1,2,3 */
STATIC double C6cs123(QNPack inLo, QNPack inHi)
{
	double C6cs123_v, 
	  TeUse, 
	  t, 
	  logx, 
	  x;
	static const double a[3] = {-92.23774,-1631.3878,-6326.4947};
	static const double b[3] = {-11.93818,-218.3341,-849.8927};
	static const double c[3] = {0.07762914,1.50127,5.847452};
	static const double d[3] = {78.401154,1404.8475,5457.9291};
	static const double e[3] = {332.9531,5887.4263,22815.211};

	DEBUG_ENTRY( "C6cs123()" );

	/* These are fits to Table 5 of
	 * >>refer	c6	cs	Aggarwal, K.M., & Kingston, A.E. 1991, J Phys B, 24, 4583
	 * C VI collision rates for 1s-3l, 2s-3l, and 2p-3l, 
	 * principal quantum numbers n and l.
	 *
	 * i is the lower level and runs from 1 to 3 (1s, 2s, 2p)
	 * j is the upper level and runs from 2 to 6 (2s, 2p, 3s, 3p, 3d)
	 * 1s-2s,2p is not done here.
	 * check temperature: fits only good between 3.8 < log Te < 6.2
	 */
	/* arrays for fits of 3 transitions see the code below for key: */

	TeUse = MAX2(phycon.te,6310.);
	t = MIN2(TeUse,1.6e6);
	x = log10(t);
	logx = log(x);

	if( inLo == in1 && inHi == in2 )
	{
		/* 1s - 2s (first entry) */
		fprintf( ioQQQ, " Carbon VI 2s-1s not done in C6cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && inHi == in3 )
	{
		/* 1s - 2p (second entry) */
		fprintf( ioQQQ, " Carbon VI 2p-1s not done in C6cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 1s - 3 (first entry) */
		C6cs123_v = a[0] + b[0]*x + c[0]*pow2(x)*sqrt(x) + d[0]*logx + 
		  e[0]*logx/pow2(x);
	}
	else if( inLo == in2 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 2s - 3 (second entry) */
		C6cs123_v = a[1] + b[1]*x + c[1]*pow2(x)*sqrt(x) + d[1]*logx + 
		  e[1]*logx/pow2(x);
	}
	else if( inLo == in3 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 2p - 3s (third entry) */
		C6cs123_v = a[2] + b[2]*x + c[2]*pow2(x)*sqrt(x) + d[2]*logx + 
		  e[2]*logx/pow2(x);
	}
	else
	{
		fprintf( ioQQQ, "  insane levels for C VI n=1,2,3 !!!\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return C6cs123_v;
}

/*Ca20cs123 line collision rates for lower levels of hydrogenic calcium, n=1,2,3 */
STATIC double Ca20cs123(QNPack inLo, QNPack inHi)
{
	double Ca20cs123_v, 
	  TeUse, 
	  t, 
	  logx, 
	  x;
	static const double a[3] = {-12.5007,-187.2303,-880.18896};
	static const double b[3] = {-1.438749,-22.17467,-103.1259};
	static const double c[3] = {0.008219688,0.1318711,0.6043752};
	static const double d[3] = {10.116516,153.2650,717.4036};
	static const double e[3] = {45.905343,685.7049,3227.2836};

	DEBUG_ENTRY( "Ca20cs123()" );

	/* 
	 * These are fits to Table 5 of
	 * >>refer	ca20	cs	Aggarwal, K.M., & Kingston, A.E. 1992, J Phys B, 25, 751
	 * Ca XX collision rates for 1s-3l, 2s-3l, and 2p-3l, 
	 * principal quantum numbers n and l.
	 *	
	 * i is the lower level and runs from 1 to 3 (1s, 2s, 2p)
	 * j is the upper level and runs from 2 to 6 (2s, 2p, 3s, 3p, 3d)
	 * 1s-2s,2p is not done here.
	 * check temperature: fits only good between 5.0 < log Te < 7.2
	 */

	/* arrays for fits of 3 transitions see the code below for key: */

	TeUse = MAX2(phycon.te,1.0e5);
	t = MIN2(TeUse,1.585e7);
	x = log10(t);
	logx = log(x);

	if( inLo == in1 && inHi == in2 )
	{
		/* 1s - 2s (first entry) */
		fprintf( ioQQQ, " Ca XX 2s-1s not done in Ca20cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && inHi == in3 )
	{
		/* 1s - 2p (second entry) */
		fprintf( ioQQQ, " Ca XX 2p-1s not done in Ca20cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && ( inHi == in4 || inHi == in5 || inHi == in6 ))
	{
		/* 1s - 3 (first entry) */
		Ca20cs123_v = a[0] + b[0]*x + c[0]*pow2(x)*sqrt(x) + d[0]*logx + 
		  e[0]*logx/pow2(x);
	}
	else if( inLo == in2 && ( inHi == in4 || inHi == in5 || inHi == in6 ))
	{
		/* 2s - 3 (second entry) */
		Ca20cs123_v = a[1] + b[1]*x + c[1]*pow2(x)*sqrt(x) + d[1]*logx + 
		  e[1]*logx/pow2(x);
	}
	else if( inLo == in3 && ( inHi == in4 || inHi == in5 || inHi == in6 ))
	{
		/* 2p - 3s (third entry) */
		Ca20cs123_v = a[2] + b[2]*x + c[2]*pow2(x)*sqrt(x) + d[2]*logx + 
		  e[2]*logx/pow2(x);
	}
	else
	{
		fprintf( ioQQQ, "  insane levels for Ca XX n=1,2,3 !!!\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return Ca20cs123_v;
}

/*Ne10cs123 line collision rates for lower levels of hydrogenic neon, n=1,2,3 */
STATIC double Ne10cs123(QNPack inLo, QNPack inHi)
{
	double Ne10cs123_v, 
	  TeUse, 
	  logx, 
	  t, 
	  x;
	static const double a[3] = {3.346644,151.2435,71.7095};
	static const double b[3] = {0.5176036,20.05133,13.1543};
	static const double c[3] = {-0.00408072,-0.1311591,-0.1099238};
	static const double d[3] = {-3.064742,-129.8303,-71.0617};
	static const double e[3] = {-11.87587,-541.8599,-241.2520};

	DEBUG_ENTRY( "Ne10cs123()" );

	/*  These are fits to Table 5 of
	 *  >>refer	ne10	cs	Aggarwal, K.M., & Kingston, A.E. 1991, PhyS, 44, 517
	 *  Ne X collision rates for 1s-3, 2s-3l, and 2p-3l, 
	 *  principal quantum numbers n and l.
	 *
	 *  i is the lower level and runs from 1 to 3 (1s, 2s, 2p)
	 *  j is the upper level and runs from 2 to 6 (2s, 2p, 3s, 3p, 3d)
	 *  1s-2s,2p is not done here.
	 *  check temperature: fits only good between 3.8 < log Te < 6.2
	 * */
	/* arrays for fits of 3 transitions see the code below for key: */

	TeUse = MAX2(phycon.te,6310.);
	t = MIN2(TeUse,1.6e6);
	x = log10(t);
	logx = log(x);

	if( inLo == in1 && inHi == in2 )
	{
		/* 1s - 2s (first entry) */
		fprintf( ioQQQ, " Neon X 2s-1s not done in Ne10cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && inHi == in3 )
	{
		/* 1s - 2p (second entry) */
		fprintf( ioQQQ, " Neon X 2p-1s not done in Ne10cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 1s - 3 (first entry) */
		Ne10cs123_v = a[0] + b[0]*x + c[0]*pow2(x)*sqrt(x) + d[0]*logx + 
		  e[0]*logx/pow2(x);
	}
	else if( inLo == in2 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 2s - 3 (second entry) */
		Ne10cs123_v = a[1] + b[1]*x + c[1]*pow2(x)*sqrt(x) + d[1]*logx + 
		  e[1]*logx/pow2(x);
	}
	else if( inLo == in3 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 2p - 3s (third entry) */
		Ne10cs123_v = a[2] + b[2]*x + c[2]*pow2(x)*sqrt(x) + d[2]*logx + 
		  e[2]*logx/pow2(x);
	}
	else
	{
		fprintf( ioQQQ, "  insane levels for Ne X n=1,2,3 !!!\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return Ne10cs123_v;
}

/*He2cs123 line collision strengths for lower levels of helium ion, n=1,2,3, by K Korista */
STATIC double He2cs123(QNPack inLo, QNPack inHi)
{
	double He2cs123_v, 
	  t;
	static const double a[11]={0.12176209,0.32916723,0.46546497,0.044501688,
	  0.040523277,0.5234889,1.4903214,1.4215094,1.0295881,4.769306,9.7226127};
	static const double b[11]={0.039936166,2.9711166e-05,-0.020835863,3.0508137e-04,
	  -2.004485e-15,4.41475e-06,1.0622666e-05,2.0538877e-06,0.80638448,2.0967075e-06,
	  7.6089851e-05};
	static const double c[11]={143284.77,0.73158545,-2.159172,0.43254802,2.1338557,
	  8.9899702e-06,-2.9001451e-12,1.762076e-05,52741.735,-2153.1219,-3.3996921e-11};

	DEBUG_ENTRY( "He2cs123()" );

	/* These are fits to Table 2
	 * >>refer	he2	cs	Aggarwal, K.M., Callaway, J., Kingston, A.E., Unnikrishnan, K.
	 * >>refercon	1992, ApJS, 80, 473
	 * He II collision rates for 1s-2s, 1s-2p, 1s-3s, 1s-3p, 1s-3d, 2s-3s, 2s-3p, 2s-3d,
	 * 2p-3s, 2p-3p, and 2p-3d. 
	 * principal quantum numbers n and l.
	 *
	 * i is the lower level and runs from 1 to 3 (1s, 2s, 2p)
	 * j is the upper level and runs from 2 to 6 (2s, 2p, 3s, 3p, 3d)
	 * check temperature: fits only good between 5,000K and 500,000K
	 * */
	/* array for fits of 11 transitions see the code below for key: */

	t = phycon.te;
	if( t < 5000. )
	{
		t = 5000.;
	}
	else if( t > 5.0e05 )
	{
		t = 5.0e05;
	}

	/**************fits begin here**************
	 * */
	if( inLo == in1 && inHi == in2 )
	{
		/* 1s - 2s (first entry) */
		He2cs123_v = a[0] + b[0]*exp(-t/c[0]);
	}
	else if( inLo == in1 && inHi == in3 )
	{
		/* 1s - 2p (second entry) */
		He2cs123_v = a[1] + b[1]*pow(t,c[1]);
	}
	else if( inLo == in1 && inHi == in4 )
	{
		/* 1s - 3s (third entry) */
		double logt = log(t);
		He2cs123_v = a[2] + b[2]*logt + c[2]/logt;
	}
	else if( inLo == in1 && inHi == in5 )
	{
		/* 1s - 3p (fourth entry) */
		He2cs123_v = a[3] + b[3]*pow(t,c[3]);
	}
	else if( inLo == in1 && inHi == in6 )
	{
		/* 1s - 3d (fifth entry) */
		He2cs123_v = a[4] + b[4]*pow(t,c[4]);
	}
	else if( inLo == in2 && inHi == in4 )
	{
		/* 2s - 3s (sixth entry) */
		He2cs123_v = (a[5] + c[5]*t)/(1 + b[5]*t);
	}
	else if( inLo == in2 && inHi == in5 )
	{
		/* 2s - 3p (seventh entry) */
		He2cs123_v = a[6] + b[6]*t + c[6]*t*t;
	}
	else if( inLo == in2 && inHi == in6 )
	{
		/* 2s - 3d (eighth entry) */
		He2cs123_v = (a[7] + c[7]*t)/(1 + b[7]*t);
	}
	else if( inLo == in3 && inHi == in4 )
	{
		/* 2p - 3s (ninth entry) */
		He2cs123_v = a[8] + b[8]*exp(-t/c[8]);
	}
	else if( inLo == in3 && inHi == in5 )
	{
		/* 2p - 3p (tenth entry) */
		He2cs123_v = a[9] + b[9]*t + c[9]/t;
	}
	else if( inLo == in3 && inHi == in6 )
	{
		/* 2p - 3d (eleventh entry) */
		He2cs123_v = a[10] + b[10]*t + c[10]*t*t;
	}
	else
	{
		/**************fits end here************** */
		fprintf( ioQQQ, "  insane levels for He II n=1,2,3 !!!\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return He2cs123_v;
}

/*Fe26cs123 line collision rates for lower levels of hydrogenic iron, n=1,2,3 */
STATIC double Fe26cs123(QNPack inLo, QNPack inHi)
{
	double Fe26cs123_v, 
	  TeUse, 
	  logx, 
	  t, 
	  x;
	static const double a[3] = {-4.238398,-238.2599,-1211.5237};
	static const double b[3] = {-0.4448177,-27.06869,-136.7659};
	static const double c[3] = {0.0022861,0.153273,0.7677703};
	static const double d[3] = {3.303775,191.7165,972.3731};
	static const double e[3] = {15.82689,878.1333,4468.696};

	DEBUG_ENTRY( "Fe26cs123()" );

	/*  These are fits to Table 5 of
	 *  >>refer	fe26	cs	Aggarwal, K.M., & Kingston, A.E. 1993, ApJS, 85, 187
	 *  Fe XXVI collision rates for 1s-3, 2s-3, and 2p-3, 
	 *  principal quantum numbers n and l.
	 *
	 *  i is the lower level and runs from 1 to 3 (1s, 2s, 2p)
	 *  j is the upper level and runs from 2 to 6 (2s, 2p, 3s, 3p, 3d)
	 *  1s-2s,2p is not done here.
	 *  check temperature: fits only good between 5.2 < log Te < 7.2
	 *  */
	/*  arrays for fits of 3 transitions see the code below for key: */

	TeUse = MAX2(phycon.te,1.585e5);
	t = MIN2(TeUse,1.585e7);
	x = log10(t);
	logx = log(x);

	if( inLo == in1 && inHi == in2 )
	{
		/* 1s - 2s (first entry) */
		fprintf( ioQQQ, " Fe XXVI 2s-1s not done in Fe26cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && inHi == in3 )
	{
		/* 1s - 2p (second entry) */
		fprintf( ioQQQ, " Fe XXVI 2p-1s not done in Fe26cs123\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( inLo == in1 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 1s - 3 (first entry) */
		Fe26cs123_v = a[0] + b[0]*x + c[0]*pow2(x)*sqrt(x) + d[0]*logx + 
		  e[0]*logx/pow2(x);
	}
	else if( inLo == in2 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 2s - 3 (second entry) */
		Fe26cs123_v = a[1] + b[1]*x + c[1]*pow2(x)*sqrt(x) + d[1]*logx + 
		  e[1]*logx/pow2(x);
	}
	else if( inLo == in3 && ( inHi == in4 || inHi == in5 || inHi == in6 ) )
	{
		/* 2p - 3s (third entry) */
		Fe26cs123_v = a[2] + b[2]*x + c[2]*pow2(x)*sqrt(x) + d[2]*logx + 
		  e[2]*logx/pow2(x);
	}
	else
	{
		fprintf( ioQQQ, "  insane levels for Ca XX n=1,2,3 !!!\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return Fe26cs123_v;
}


realnum HydroCSInterp(long nelem,
			long ipHi,
			long ipLo,
			long ipCollider )
{
	DEBUG_ENTRY( "HydroCSInterp()" );
	
	long nHi = iso_sp[ipH_LIKE][nelem].st[ipHi].n();
	long lHi = iso_sp[ipH_LIKE][nelem].st[ipHi].l();
	long sHi = iso_sp[ipH_LIKE][nelem].st[ipHi].S();
	long gHi = iso_sp[ipH_LIKE][nelem].st[ipHi].g();
	double IP_Ryd_Hi = iso_sp[ipH_LIKE][nelem].fb[ipHi].xIsoLevNIonRyd;
	long nLo = iso_sp[ipH_LIKE][nelem].st[ipLo].n();
	long lLo = iso_sp[ipH_LIKE][nelem].st[ipLo].l();
	long sLo = iso_sp[ipH_LIKE][nelem].st[ipLo].S();
	//	long gLo = iso_sp[ipH_LIKE][nelem].st[ipLo].g();
	double IP_Ryd_Lo = iso_sp[ipH_LIKE][nelem].fb[ipLo].xIsoLevNIonRyd;
	double Aul = iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().Aul();
	// collisions are from high to low level, then initial level lifetime is from higher level
	double tauLo = iso_sp[ipH_LIKE][nelem].st[ipHi].lifetime();
	double EnerWN = iso_sp[ipH_LIKE]->trans(ipHi,ipLo).EnergyWN();
	double EnerErg = iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).EnergyErg();
	const char *where="      ";

	realnum CStemp = GetHlikeCollisionStrength( nelem, ipCollider,
					nHi, lHi, sHi, gHi, IP_Ryd_Hi,
					nLo, lLo, sLo,/*gLo,*/ IP_Ryd_Lo,
					Aul, tauLo, EnerWN, EnerErg, &where );
	return CStemp;
}

realnum GetHlikeCollisionStrength( long nelem, long ipCollider,
					long nHi, long lHi, long sHi, long gHi, double IP_Ryd_Hi,
					long nLo, long lLo, long sLo, /*long gLo,*/ double IP_Ryd_Lo,
					double Aul, double tauLo, double EnerWN, double EnerErg, const char **where )
{
	DEBUG_ENTRY( "GetHlikeCollisionStrength()" );

	if( !iso_ctrl.lgColl_excite[ipH_LIKE] )
		return 0.;

	double CStemp;
	bool lgResolvedData = true;

	// Energy gap between levels in eV
	double deltaE_eV = EnerErg/EN1EV;

	*where = "      ";
	/* This set of collision strengths should only be used
	 * if the Storey and Hummer flag is set */
	if( opac.lgCaseB_HummerStorey )
	{
		if( nLo == nHi )
		{
			ASSERT( lLo >= 0);
			ASSERT( lHi >= 0);

			if( nHi <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max &&
				abs( lLo - lHi ) != 1 )
			{
				/* if delta L is not +/- 1, set collision strength to zero. */
				CStemp = 0.;
			}
			else if (ipCollider == ipELECTRON)
			{
				/* l-changing collisions are calculated for heavy slow projectile */
				CStemp = 0.;
			}
			else
			{
				//classical PS64 for Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165
				CStemp =  CS_l_mixing_PS64(
						nelem,
						ipH_LIKE,
						tauLo,
						nelem+1.-ipH_LIKE,
						nHi,
						lHi,
						gHi,
						lLo,
						deltaE_eV,
						ipCollider);
				if (0)
				{
				// New PS-M Refer F. Guzman MNRAS (2016)
					CStemp = CS_l_mixing_PS64_expI(
							nelem,
							ipH_LIKE,
							tauLo,
							nelem+1.-ipH_LIKE,
							nHi,
							lHi,
							gHi,
							lLo,
							deltaE_eV,
							ipCollider,
							iso_ctrl.lgCS_PSM20[ipH_LIKE]);
				}
			}
		}

		else 
		{

			if( nHi <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max &&
				abs( lLo - lHi ) != 1 )
			{
				/* if delta L is not +/- 1, set collision strength to zero. */
				CStemp = 0.;
			}
			else if( ipCollider == ipELECTRON )
			{
				CStemp = CS_ThermAve_PR78( ipH_LIKE, nelem, nHi, nLo,
								EnerErg / EN1RYD, phycon.te );

				/* The data retourned by the previous routine is not resolved in l */
				lgResolvedData = false;

			}
			else
				CStemp = 0.;
		}
	}
	else
	{
		/* HlikeCSInterp interpolates on a table to return R-matrix collision strengths
		 * for hydrogen only
		 * */
		if((CStemp = HlikeCSInterp( nelem, ipCollider, nHi, lHi, sHi, nLo, lLo, sLo )) >= 0.f )
		{
			 /* 14/02/2019 FG: The routine has been modified to provide collapsed to resolved and
			 * collapses to collapsed values
			 */
			*where = "table";
		}
		else if( nLo == nHi )
		{
			/*This is always resolved, hence the asserts */
			ASSERT( lLo >= 0);
			ASSERT( lHi >= 0);

			if ( ipCollider == ipELECTRON)
				/* l-changing collisions are calculated for heavy slow projectile */
				CStemp = 0.;
			else if( iso_ctrl.lgCS_Vrinceanu[ipH_LIKE] )
			{
				CStemp = CS_l_mixing_VF01(ipH_LIKE, nelem,
					nHi,
					lHi,
					lLo,
					sLo,
					gHi,
					tauLo,
					IP_Ryd_Hi,
					IP_Ryd_Lo,
					phycon.te, 
					ipCollider );
				*where ="VF01  ";
			}
			else if( iso_ctrl.lgCS_VOS12[ipH_LIKE] )
			{
				/* Use the method given in 
				 * >>refer He	CS	 Vrinceau etal ApJ 747, 56 (2012); semiclasical treatment: eq. (7)
				 * statistical weights included */
				CStemp = CS_l_mixing_VOS12(nHi,lHi,lLo,nelem,gHi,nelem+1.-ipH_LIKE,ipCollider,phycon.sqrte);
				*where = "VOS12 ";
			}
			else if( iso_ctrl.lgCS_VOS12QM[ipH_LIKE] )// && abs( lLo - lHi ) == 1)
			{
				/* Use the method given in
				 * >>refer He	CS	 Vrinceau etal ApJ 747, 56 (2012); quantal treatment: eq. (2)
				 * statistical weights included */
				CStemp = CS_l_mixing_VOS12QM(ipH_LIKE, nelem,
						nLo,
						lHi,
						lLo,
						sLo,
						gHi,
						tauLo,
						IP_Ryd_Hi,
						IP_Ryd_Lo,
						phycon.te,
						ipCollider );
				*where = "VOS12Q";
			}
			else if ( iso_ctrl.lgCS_PS64[ipH_LIKE] && abs( lLo - lHi ) == 1 )
			{
				if (iso_ctrl.lgCS_PSClassic[ipH_LIKE])
				{
				//classical PS64 for Pengelly, R.M., & Seaton, M.J., 1964, MNRAS, 127, 165
					CStemp = CS_l_mixing_PS64(
							nelem,
							ipH_LIKE,
							tauLo,
							nelem+1.-ipH_LIKE,
							nHi,
							lHi,
							gHi,
							lLo,
							deltaE_eV,
							ipCollider);
					*where = "PS64  ";
				}
				else
				{
				// PS-M Refer F. Guzman MNRAS (2016)
					CStemp = CS_l_mixing_PS64_expI(
							nelem,
							ipH_LIKE,
							tauLo,
							nelem+1.-ipH_LIKE,
							nHi,
							lHi,
							gHi,
							lLo,
							//sHi,
							deltaE_eV,
							ipCollider,
							iso_ctrl.lgCS_PSM20[ipH_LIKE]);
				if(iso_ctrl.lgCS_PSM20[ipH_LIKE])
					*where = "PSM   ";
				else
					*where = "PSM20 ";
				}
			}
			else
				CStemp = 0.;
		}
		else if ( nHi != nLo )
		{
			/* 1) If both Upper and Lower levels are resolved, then, if routine gives unresolved,
			 * the result should be divided by the degeneration of upper level
			 *
			 * 2) If both are collapsed, if a resolved routine is used, results should be
			 * added up.
			 *
			 * I3) If Lower level is resolved and upper is collapsed we need to
			 * we need to consider only dipole transitions to lower level
			 */

			/* F. Guzman 07/18/19 Heavy ions impact are 2 to three orders of magnitude lower than electron impacts at nebular typical
			 * temperatures. The heavy ion impact approach using Vriens & Smeets 1980 used here is valid for high energies
			 * and is based in a mass scalling of the collision strength (Burgess and Tully (2005).
			 * Cloudy implementation of Percival & Richards 1978 does not include this treatment. Furthermore the Coulomb
			 * repulsion between ions will make these cross sections smaller only to be relevant at even higher energies.  */
			/* collisions using Vriens & Smeets only for heavy atoms if not specified */
			/* V&S do not have Z dependence, so is only used for neutrals */
			if( nelem==ipHYDROGEN && (iso_ctrl.lgCS_Vriens[ipH_LIKE] || ipCollider != ipELECTRON) )
			{
				CStemp = CS_VS80( ipH_LIKE, nHi, IP_Ryd_Hi, nLo, IP_Ryd_Lo, Aul, nelem, ipCollider, phycon.te );
				*where = "Vriens";
				/* This routine gives NO l-resolved data */
				lgResolvedData = false;
			}
			/* only electron impact collisions */
			else if (ipCollider == ipELECTRON)
			{
				if(iso_ctrl.lgCS_Lebedev[ipH_LIKE])
				{
					/* Lebedev and Beigman (1998) Phys. Highly excited atoms and ions p. 222
					 * D. Docenko 2013 Baltic Astr. 22, 363 */
					CStemp = hydro_Lebedev_deexcit(ipH_LIKE, nelem, nHi, nLo, IP_Ryd_Lo);
					*where = "lebed";
					lgResolvedData = false;
				}
				else if ( iso_ctrl.lgCS_Fujim[ipH_LIKE] )
				{


					CStemp = hydro_Fujimoto_deexcit(ipH_LIKE, nHi, nLo, Aul, IP_Ryd_Hi, IP_Ryd_Lo);
					*where = "Fuji ";
					lgResolvedData = false;
				}
				else if( iso_ctrl.lgCS_vrgm[ipH_LIKE])
				{
					CStemp = hydro_vanRegemorter_deexcit( ipH_LIKE, nelem, nHi, lHi, lLo,
										EnerWN, deltaE_eV, Aul, ipCollider );
					*where = "vrgm ";
					lgResolvedData = false;
				}
				else
				{
					/* highly excited levels */
					/* Vriens&Smeets (1980) give no Z dependence on their cross sections,
					 * Percival and Richards (1978) have got a Z dependence so their rates are preferred */
					CStemp = CS_ThermAve_PR78( ipH_LIKE, nelem, nHi, nLo,
									EnerErg / EN1RYD, phycon.te );
					*where = "PR78  ";
					/* This routine gives NO l-resolved data */
					lgResolvedData = false;
				}

			}
			else
				CStemp = 0.;
		}
		else
			TotalInsanity();

	}

	/*
	 * Resolve collision strengths for n-changing collisions
	 * from collapsed to resolved levels
	 */
	if (!lgResolvedData && nLo <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max)
	{
		double l_resolve =
			CSresolver(ipH_LIKE, nHi, lHi, sHi, nLo, lLo, sLo,
					iso_sp[ipH_LIKE][nelem].n_HighestResolved_max);

		enum {DEBUG_LOC = false};
		if( DEBUG_LOC )
		{
			if( nelem == ipHYDROGEN &&
				(nHi == iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + 15) &&
				 nLo == iso_sp[ipH_LIKE][nelem].n_HighestResolved_max )
			{
				fprintf( ioQQQ, "nelem: %ld"
						"  (nHi, lHi, sHi): (%ld, %ld, %ld) ->"
						"  (nLo, lLo, sLo): (%ld, %ld, %ld)"
						"\t CS: %.4e"
						"\t l-reslv: %.4e ->"
						"\t final: %.4e\n",
						nelem,
						nHi, lHi, sHi,
						nLo, lLo, sLo,
						CStemp, l_resolve, CStemp * l_resolve );
			}
		}

		CStemp *= l_resolve;
	}

	return (realnum)CStemp;
}

STATIC realnum HlikeCSInterp( long nelem, long Collider, long nHi, long lHi, long sHi, long nLo, long lLo, long sLo )
{
	DEBUG_ENTRY( "HlikeCSInterp()" );

	long topLo =1;
	long topHi =1;
	long l = lHi;
	long lp = lLo;
	long s = sHi;
	long sp = sLo;
	
	realnum cs = 0.;


	if (nLo > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max)
	{
		topLo = nLo;
		sp = 2;
	}

	if (nHi > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max)
	{
		topHi = nHi;
		s = 2;
	}

	ASSERT( s==2 && sp==2 );


	/* The loops ensure the addition of coll. strengths when the routine is called for collapsed levels */
	for (long Lo=0; Lo< topLo; Lo ++)
	{
		realnum csHi=0.;
		if (nLo > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max)
			lp = Lo;
		for (long Hi=0; Hi < topHi; Hi++)
		{
			if (nHi > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max)
				l = Hi;

			realnum cs1 = -1.f;
	
			if( nelem==ipHYDROGEN && Collider==ipELECTRON && nHi <= 5 && ( nHi != nLo ) )
			{
				cs1 = HCSAR_interp(nHi,l,s,nLo,lp,sp);
			}
			else if( nelem>ipHYDROGEN && Collider==ipELECTRON && nHi <= 3 && nLo < 3 )
			{
				cs1 = Hydcs123(nelem,nHi,l,s,nLo,lp,sp,'e');
			}
			else if( nelem>ipHYDROGEN && Collider==ipPROTON && nHi==2 && l==1 && nLo==2 && lp==0 )
			{
				cs1 = Hydcs123(nelem,nHi,l,s,nLo,lp,sp,'p');
			}

			if (nHi > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max && cs1 >= 0.)
				csHi += cs1;
			else
			{
				csHi = cs1;
				break;
			}
		}

		if (nLo > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max && csHi >= 0.)
			cs += csHi;
		else
		{
			cs = csHi;
			break;
		}

	}

	return cs;
}
