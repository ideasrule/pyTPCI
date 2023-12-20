/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdInit routine to initialize variables, called at start of calculation */
/*cdPrepareExit prepare termination of the code, but do not terminate yet */
/* unset EXTERN so that everything is defined here */
#include "cddefines.h"

/* used for saving map*/
FILE *ioMAP = NULL;

/* external ZeroNum used to div by zero 
 * ok here since never changed*/
const double ZeroNum = 0.;
/* external ZeroPtr used to force segfault */
const void* ZeroPtr = NULL;

/* this must go here since it defines NTA needed for other lines*/
#include "taulines.h"

/* following is true extern in taulines.h */
long nWindLine = NWINDDIM;

#include "monitor_results.h"
#include "called.h"
#include "cddrive.h"
/* this will be set true when cdInit is called.  The definition is in cdInit.
* Other routines will check that this is true when they are called, 
* to verify that cdInit was called first */
bool lgcdInitCalled=false;
#include "grid.h"
/* this is set true once space allocated, then never change
* number of levels again with hydrogenic command, 
* also to make sure allocation only happens one time  */
bool lgHydroAlloc = false;
/*  */
#include "input.h"
#include "parse.h"
/* */
#include "save.h"
/* set true when allocated, init to false */
bool lgRfieldAllocated=false;
bool lgOpacAllocated=false;
#include "init.h"
#include "ran.h"
#include "trace.h"


/* =================================================================== */
void cdInit(void)
{
	long i;
	double vtest;

	DEBUG_ENTRY( "cdInit()" );

	/* set flag saying that cdInit has been called */
	lgcdInitCalled = true;

	/* initialize the search path so that we can read data files */
	cpu.i().initPath();

	/*test if the following integer types have the correct width*/
	if( sizeof(int16) != 2 || sizeof(uint16) != 2 || sizeof(int32) != 4 || sizeof(uint32) != 4 )
		TotalInsanity();

	/*********************************************************
	 *  on a VAX compile with /G_FLOATING option on FORTRAN; *
	 *  following makes sure this happened.                  *
	 *********************************************************/
	vtest = 1e-35;
	vtest /= 1e35;
	if( vtest == 0. )
	{
		fprintf( ioQQQ, " Something is wrong with the double precision.  Use /g_floating on a VAX\n" );
	}

	// initialize random number generator -- also done in maincl.cpp, so this is normally ignored
	// but when calling Cloudy as a subroutine, this is the place that does the actual initialization
	ran.init();
	
	/* initialize some variables dealing with cloudy's interaction with machine environment */
	/* if TALK is true then do standard printout
	 * if false then never say anything */
	/* only the master rank produces output */
	called.lgTalk = cpu.i().lgMPI_talk();
	/* this flag is needed to turn print on to have effect */
	called.lgTalkIsOK = cpu.i().lgMPI_talk();
	/* means talk not forced off by call to cdTalk*/
	called.lgTalkForcedOff = false;

	optimize.lgNoVary = false;
	optimize.lgVaryOn = false;
	optimize.lgOptimr = false;
	grid.lgGrid = false;
	grid.nGridCommands = 0;

	for( i=0; i<NUM_OUTPUT_TYPES; i++ )
	{
		grid.lgOutputTypeOn[i] = false;
	}

	/* this is a global variable in monitor_results.h, and can be checked by
	 * other routines to see if asserts are ok - (most calculations will not use asserts,
	 * and this will be the only place values are set, although they will be checked in maincl) */
	lgMonitorsOK = true;
	lgBigBotch = false;
	lgPrtSciNot = false;

	/* nRead is the number of the command in the input stream - many optimize options
	 * point to it to refer to the original command.  it is incremented before
	 * it is used, so will become 0.  it is the array element within the stack
	 * of emission lines */
	input.nRead = -1;
	/* the call to clear() is needed when calling Cloudy as a subroutine in a loop
	 * it will wipe the input card deck from the previous iteration */
	input.clear();
	input.crd.reserve(500);

	/* start the timer to log execution time */
	cdSetExecTime();

	/* zero out lots of variables */
	zero();
	return;
}


/* =================================================================== */
/* cdPrepareExit prepare termination of the code, but do not terminate yet
 * this routine should only be called by exception handlers, never from the main code */
void cdPrepareExit(exit_type exit_status)
{
	enum {DEBUG_LOC=false};
	if( DEBUG_LOC )
		fprintf(ioQQQ," cdExit called\n");

	// turn trace output off now since we are about to close the output file
	trace.lgTrace = false;

	// make sure file descriptors are closed in case they were redirected
	cdInput( "", "" );
	cdOutput( "", "" );

	// make sure the error condition is logged in the SAVE GRID output
	// we do this here (and not SaveDo) to make sure that the output is complete
	if( grid.lgGrid && cpu.i().MPIMode() == MS_GRID && save.ipSaveGrid >= 0 )
		SaveGrid( grid.pnunit, exit_status );

	/* close any open units */
	CloseSaveFiles( true );
}

