/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "iterations.h"

t_iterations iterations;

void t_iterations::zero()
{
	for( long i=0; i < iter_alloc; i++ )
	{
		IterPrnt[i] = 10000;
	}
	itermx = 0;
	/* this implements set coverage command */
	lgConverge_set = false;
	iteration = 0;
	
	/* this is default number of zones
	 * >>chng 96 jun 5, from 400 to 500 for thickest corners4 grid */
	/* >>chng 04 jan 30, from 600 to 800, code uses finer zoning today */
	/* >>chng 04 dec 24, from 800 to 1400, so that HII region - molecular cloud
	 * sims do not need set nend - all sims in test suite will run ok without set nend */
	nEndDflt=1400;
	lgEndDflt = true;	

	for( long i=0; i < iter_alloc; i++ )
	{
		nend[i] = nEndDflt;
		/*>>chng 03 nov 13, from 1e30 to 1e31, because default inner radius raised to 1e30 */
		StopThickness[i] = 1e31;
		StopRadius[i] = -1.;
	}	
}
void t_iterations::alloc()
{
	/* this is number of iterations that have been allocated - we could 
	 * increase this if more iterations are needed */
	iter_alloc = 200;
	IterPrnt.resize((size_t)iter_alloc);
	nend.resize( (size_t)iter_alloc);
	StopThickness.resize((size_t)iter_alloc);
	StopRadius.resize((size_t)iter_alloc);
}
