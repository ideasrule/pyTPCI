/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HeLikeError fills uncertainty arrays */
#include "cddefines.h" 
#include "ran.h"
#include "iso.h"

/* TruncatedNormal takes as input a percent uncertainty less than 33.3%
 * (expressed as 0.333). The routine then assumes this input variable represents one
 * standard deviation about a mean of unity, and returns a random number within
 * that range.  A hard cutoff is imposed at three standard deviations, which 
 * eliminates roughly 0.3% of the normal distribution.  In other words, the routine
 * returns a number in a pseudo-normal distribution with standard deviation equal to
 * the input.  The random deviate will be between 1-3*stdev and 1+3*stdev. */
inline realnum TruncatedNormal( double PctUncertainty )
{
	ASSERT( PctUncertainty < 0.333 );

	double result;
	do
	{
		/* We want this "percent uncertainty" to represent one standard deviation */
		result = 1. + ran.normal()*PctUncertainty;
	}
	/* only allow values that are within 3 standard deviations */
	while( (result < 1.-3.*PctUncertainty) || (result > 1.+3.*PctUncertainty) );
	return (realnum)result;
}

/* This routine handles errors when that option is turned on (via the command
 * "atom he-like error generation" */
void iso_put_error(long int ipISO,
			  long int nelem,
			  long int ipHi,
			  long int ipLo,
			  long int whichData,
			  realnum errorOpt,
			  realnum errorPess)
{
	DEBUG_ENTRY( "iso_put_error()" );

	if( iso_ctrl.lgRandErrGen[ipISO] )
	{
		/* whichData is either IPRAD, IPCOLLIS, or IPENERGY */
		ASSERT( whichData <= 2 );
		ASSERT( ipISO < NISO );
		ASSERT( nelem < LIMELM );
		ASSERT( ipHi <= iso_sp[ipISO][nelem].numLevels_max );
		ASSERT( ipLo <= iso_sp[ipISO][nelem].numLevels_max );
		ASSERT( errorOpt >= 0. );
		ASSERT( errorPess >= 0. );
		
		if( !iso_ctrl.lgPessimisticErrors )
			iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[whichData] = errorOpt;
		else
			iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[whichData] = errorPess;
	}
	return;
}

void iso_put_error(long int ipISO,
		   long int nelem,
		   QNPack inHi,
		   QNPack inLo,
		   long int whichData,
		   realnum errorOpt,
		   realnum errorPess)
{
	DEBUG_ENTRY( "iso_put_error()" );

	if( iso_ctrl.lgRandErrGen[ipISO] )
	{
		long ipHi = iso_sp[ipISO][nelem].QN2Index(inHi);
		long ipLo = iso_sp[ipISO][nelem].QN2Index(inLo);

		if( ipHi >= 0 && ipLo >= 0 )
			iso_put_error(ipISO, nelem, ipHi, ipLo, whichData, errorOpt, errorPess);
	}
}

void iso_error_generation( long ipISO, long nelem )
{
	long ipHi, ipLo, typeOfRate;

	DEBUG_ENTRY( "iso_error_generation()" );

	for( ipHi=1; ipHi<= iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
	{
		/* >>chng 06 mar 15, the upper limit incorrectly went to numLevels_max */
		for( ipLo=0; ipLo< ipHi; ipLo++ )
		{
			for( typeOfRate=0; typeOfRate<=1; typeOfRate++ )
			{
				if( iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[typeOfRate] >= 0. )
				{
					iso_sp[ipISO][nelem].ex[ipHi][ipLo].ErrorFactor[typeOfRate] =
						TruncatedNormal(iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[typeOfRate]);
					ASSERT( iso_sp[ipISO][nelem].ex[ipHi][ipLo].ErrorFactor[typeOfRate] > 0. );
				}
				else
				{
					iso_sp[ipISO][nelem].ex[ipHi][ipLo].ErrorFactor[typeOfRate] = 1.0f;
				}
			}
		}
	}

	/* set flag saying that error generation has been done.  */
	iso_sp[ipISO][nelem].lgErrGenDone = true;
	return;
}
