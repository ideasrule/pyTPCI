/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ContNegative sanity check for negative continuum intensities */
#include "cddefines.h"
#include "continuum.h"
#include "rfield.h"
#include "dense.h"
#include "phycon.h"

void ContNegative(void)
{
	DEBUG_ENTRY( "ContNegative()" );

	/* look for negative continuum points */
	bool lgFNeg = false;
	for( long i=0; i < rfield.nflux; i++ )
	{
		if( rfield.flux[0][i] < 0. )
		{
			fprintf( ioQQQ, " FLUX negative, value, freq, pointer=%10.2e%10.2e%5ld %4.4s %4.4s\n", 
			  rfield.flux[0][i], rfield.anu(i), i, rfield.chLineLabel[i].c_str()
			  , rfield.chContLabel[i].c_str() );
			lgFNeg = true;
		}
		else if( rfield.otslin[i] < 0. )
		{
			fprintf( ioQQQ, " OTSLIN negative, value, freq, pointer=%10.2e%10.2e%5ld\n", 
			  rfield.otslin[i], rfield.anu(i), i );
			lgFNeg = true;
		}
		else if( rfield.otscon[i] < 0. )
		{
			fprintf( ioQQQ, " OTSCON negative, value, freq, pointer=%10.2e%10.2e%5ld\n", 
			  rfield.otscon[i], rfield.anu(i), i );
			lgFNeg = true;
		}
		else if( rfield.SummedCon[i] < 0. )
		{
			fprintf( ioQQQ, " OUTCON negative, value, freq, pointer=%10.2e%10.2e%5ld\n", 
			  rfield.ConInterOut[i], rfield.anu(i), i );
			lgFNeg = true;
		}
		else if( rfield.ConInterOut[i] < 0. )
		{
			fprintf( ioQQQ, " OUTCON negative, value, freq, pointer=%10.2e%10.2e%5ld\n", 
			  rfield.ConInterOut[i], rfield.anu(i), i );
			lgFNeg = true;
		}
		else if( rfield.outlin[0][i] < 0. )
		{
			fprintf( ioQQQ, " OUTLIN negative, value, freq, pointer=%10.2e%10.2e%5ld\n", 
			  rfield.outlin[0][i], rfield.anu(i), i );
			lgFNeg = true;
		}
	}

	if( !lgFNeg )
	{
		fprintf( ioQQQ, " No parts of the continuum were negative, the electron density was%10.2e te=%10.2e\n", 
		  dense.eden, phycon.te );
		fprintf( ioQQQ, " This is zone number%4ld\n", nzone );
	}
	return;
}
