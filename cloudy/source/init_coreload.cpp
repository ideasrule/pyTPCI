/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitCoreload one time initialization of core load, called from cdDrive, this sets
* minimum set of values needed for the code to start - called after
* input lines have been read in and checked for VARY or GRID - so
* known whether single or multiple sims will be run */
#include "cddefines.h"
#include "init.h"
#include "parse.h"
#include "h2.h"
#include "iso.h"

// //////////////////////////////////////////////////////////////////////////
//
//
// NB DO NOT ADD VARIABLES TO THIS FILE!  THE GOAL IS TO REMOVE THIS FILE
// initialization of variables done one time per coreload should be done in
// a constructor for the data
//
//
// //////////////////////////////////////////////////////////////////////////

/*InitCoreload one time initialization of core load, called from cdDrive, this sets
* minimum set of values needed for the code to start - called after
* input lines have been read in and checked for VARY or GRID - so
* known whether single or multiple sims will be run */
void InitCoreload( void )
{
	static int nCalled=0;

	DEBUG_ENTRY( "InitCoreload()" );

	/* return if already called */
	if( nCalled )
		return;

	++nCalled;

	iso_init();

	SaveFilesInit();

	diatoms_init();
}
