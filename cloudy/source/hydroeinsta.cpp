/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HydroEinstA calculates Einstein A's from  osillator strengths*/
#include "cddefines.h"
#include "hydroeinsta.h"
#include "hydro_bauman.h"
#include "hydrooscilstr.h"
#include "iso.h"
#include "hydro_tbl.h"
#include "lines_service.h"

STATIC realnum hydro_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo );
STATIC realnum hydro_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo );

double HydroEinstA(long int n1, 
	  long int n2)
{
	long int lower, iupper;
	double EinstA_v, 
	  ryd, 
	  xl, 
	  xmicron, 
	  xu;

	DEBUG_ENTRY( "HydroEinstA()" );
	/* (lower,upper) of Johnson 1972.  */

	/* strictly n -> n' transition probabilities
	 * no attempt to distribute according to l,l' */

	/* sort out the order of upper and lower, so can be called either way */
	lower = MIN2( n1 , n2 );
	iupper = MAX2( n1, n2 );
	if( lower < 1 || lower == iupper )
	{
		fprintf(ioQQQ," HydroEinstA called with impossible ns, =%li %li\n", lower, iupper);
		cdEXIT(EXIT_FAILURE);
	}

	xl = (double)lower;
	xu = (double)iupper;
	ryd = 1./POW2(xl) - 1./POW2(xu);
	xmicron = 1.E4/(ryd*RYD_INF);
	EinstA_v = HydroOscilStr(xl,xu)*TRANS_PROB_CONST*1e8f/(POW2(xmicron))*xl*xl/xu/xu;
	return EinstA_v;
}

realnum hydro_transprob( long nelem, long ipHi, long ipLo )
{
	DEBUG_ENTRY( "hydro_transprob()" );

	realnum Aul;
	long ipISO = ipH_LIKE;
	realnum error = 0.001f;
	if( ipHi >= iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max )
	{
		if( ipLo >= iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max )
		{
			Aul = hydro_transprob_collapsed_to_collapsed( nelem, N_(ipHi), N_(ipLo) );
			ASSERT( Aul > 0.f );
		}
		else 
		{
			Aul = hydro_transprob_collapsed_to_resolved( nelem, N_(ipHi), N_(ipLo), L_(ipLo) );
			error = 0.01f;
		}
	}
	else
	{
		if( N_(ipHi) == N_(ipLo) )
		{	
			if( false && L_(ipHi) == L_(ipLo)+1 && size_t(N_(ipHi)) <= t_hydro_tbl::Inst().nmaxnn() )
				Aul = t_hydro_tbl::Inst().tp( N_(ipLo), L_(ipLo), L_(ipHi), nelem+1 );
			else
				Aul = SMALLFLOAT;
		}
		else if( ipLo == 0 && ipHi == 1 )
		{
			// M1 transition 2s_1/2 -> 1s_1/2
			// for low Z the 2E1 transition will dominate, but for Z > 40 the M1 transition is stronger
			// >> refer	H-like	As	Marrus, E. \& Mohr, P. J. Advances in Atomic and Molecular Physics, Vol. 14, Academic, New York, 1978, p. 181
			Aul = 2.46e-6*powi((double)(nelem+1.),10);
		}
		else if( abs( L_(ipLo) - L_(ipHi) ) == 1 )
		{
			Aul = H_Einstein_A( N_(ipHi), L_(ipHi), N_(ipLo), L_(ipLo), nelem+1 );
		}
		else
		{
			ASSERT( N_(ipHi) > N_(ipLo) );
			ASSERT( L_(ipHi) == L_(ipLo) || 
				abs(L_(ipHi)-L_(ipLo)) > 1 );
			Aul = SMALLFLOAT;
		}
	}
	iso_put_error(ipISO,nelem,ipHi,ipLo,IPRAD,error,error);

	return Aul;
}

STATIC realnum hydro_transprob_collapsed_to_collapsed( long nelem, long nHi, long nLo )
{
	DEBUG_ENTRY( "hydro_transprob_collapsed_to_collapsed()" );

	ASSERT( nHi > nLo );

	if( size_t(nHi) <= t_hydro_tbl::Inst().nmaxn() )
		return t_hydro_tbl::Inst().tp( nLo, nHi, nelem+1 );
	else
		return powi(double(nelem+1),4)*HydroEinstA( nLo, nHi );
}

STATIC realnum hydro_transprob_collapsed_to_resolved( long nelem, long nHi, long nLo, long lLo )
{
	DEBUG_ENTRY( "hydro_transprob_collapsed_to_resolved()" );
		
	/* Lower level resolved, upper not. First calculate Aul
	 * from upper level with ang mom one higher.	*/
	realnum Aul = H_Einstein_A( nHi, lLo+1, nLo, lLo, nelem+1 );

	Aul *= (2.*(lLo+1.)+1.) * 2. / (2.*(double)nHi*(double)nHi);

	if( lLo != 0 )
	{
		/* For all l>0, add in transitions from upper level with ang mom one lower.	*/
		double Aul1 = H_Einstein_A( nHi, lLo-1, nLo, lLo, nelem+1 );
	
		/* now add in this part after multiplying by stat weight for lHi = lLo-1.	*/
		Aul += Aul1*(2.*(lLo-1.)+1.) * 2. / (2.*(double)nHi*(double)nHi);
	}

	ASSERT( Aul > 0.);

	return Aul;
}
