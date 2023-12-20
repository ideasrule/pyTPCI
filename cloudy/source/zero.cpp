/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*zero zero out or initialize variables, called by cdInit, but also by optimize_func during optimization,
 * this is called before any commands are parsed, called one time per model, at very start */
#include "cddefines.h"
#include "init.h"
#include "iterations.h"
#include "ionbal.h"
#include "atmdat_adfa.h"
#include "yield.h"
#include "secondaries.h"
#include "called.h"
#include "mole.h"
#include "thermal.h"
#include "trace.h"

// //////////////////////////////////////////////////////////////////////////
//
//
// NB DO NOT ADD VARIABLES TO THIS FILE!  THE GOAL IS TO REMOVE THIS FILE
// initialization of variables should be done in one of the ini_*.cpp routines
//
//
// //////////////////////////////////////////////////////////////////////////

/* zero out or initialize variables, called by cdInit, but also by optimize_func 
 * during optimization, called before command parser, one time per model, 
 * in a grid one time per grid point (so called nGrid times), 
 * only one time in multi-iteration models */
void zero()
{
	DEBUG_ENTRY( "zero()" );

	/* this is used to signify the first call to this routine.  At that
	 * stage some memory has not been allocated so must not initialize,
	 * set false at very end of this routine */
	static bool lgFirstCall = true;
	/* >>chng 06 jun 27, only allocate on first call - memory leak */
	if (lgFirstCall)
	{
		iterations.alloc();
		secondaries.alloc();
		ionbal.alloc();
		mole.alloc();
		lgFirstCall = false;
	}

	/* this routine is called exactly one time at the start of
	 * the calculation of a single model.  When the code is used as a subroutine
	 * this routine is called one time for each model.  It is called before
	 * the boundary conditions are read in, and is never called again
	 * during that calculation of the one model.  
	 * All default variables that must be initialized before a calculation starts
	 * must appear in the routine.  In a grid they are reset for each model
	 */

	/* this option, use the new atmdat_rad_rec recombination rates */
	t_ADfA::Inst().set_version( PHFIT96 );

	/**********************************************************************
	 * these are options to print errors to special window,
	 * set with print errors command, 
	 * output will go to standard error
	 * defined in cdInit 
	 **********************************************************************/
	lgPrnErr = false;
	ioPrnErr = stderr;

	vector<module*>& mods = module_list::Inst().m_l;
	for (vector<module*>::iterator mi = mods.begin(); mi != mods.end(); ++mi)
	{
		if ( trace.lgTrace )
			fprintf(ioQQQ,"Zeroing %s\n",(*mi)->chName());
		(*mi)->zero();
	}

	/* now set physical conditions array 
	 * following will force updating all temperature - density variables */
	TempChange(1e4, true);

	fnzone = 0.;
	nzone = 0;
	/* save initial condition for talk in case PRINT LAST used */
	called.lgTalkSave = called.lgTalk;

	/* this tells the code to use standard Auger yields */
	t_yield::Inst().reset_yield();

	/* there was a call to TestCode */
	lgTestCodeCalled = false;
	/* test code enabled with set test command */
	lgTestCodeEnabled = false;
}
