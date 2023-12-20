/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "ran.h"
#include "random_crash.h"

// generate a random crash, used to test the infrastructure for handling such crashes
// this implements the CRASH GRID comand
// the fake return argument is there to prevent the statements below from being optimized away
double RandomCrash()
{
	DEBUG_ENTRY( "RandomCrash()" );

	vector<int> signals;
	signals.push_back( -5 );
	signals.push_back( -4 );
	signals.push_back( -3 );
	signals.push_back( -2 );
	signals.push_back( -1 );
#	ifdef SIGNAL_HANDLER
	signals.push_back( SIGFPE );
	signals.push_back( SIGSEGV );
#	endif
	// 1 in 4 models will crash
	if( (ran.u8()&0x3) == 0 )
	{
		unsigned int i;
		if( signals.size() > 1 )
			i = ran.u8()%signals.size();
		else
			i = 0;
		int s = signals[i];
		if( s == -5 )
			DOMAIN_ERROR( "failed due to CRASH GRID" );
		else if( s == -4 )
			cdEXIT(EXIT_FAILURE);
		else if( s == -3 )
			throw cloudy_abort( "failed due to CRASH GRID" );
		else if( s == -2 )
			throw bad_assert( __FILE__, __LINE__, "failed due to CRASH GRID" );
		else if( s == -1 )
			OUT_OF_RANGE( "failed due to CRASH GRID" );
		else if( s == SIGFPE )
		{
			double x = 2./ZeroNum;
			return x;
		}
		else if( s == SIGSEGV )
		{
			double* p = (double*)ZeroPtr;
			return *p;
		}
	}
	return 0.;
}
