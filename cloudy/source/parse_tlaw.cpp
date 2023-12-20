/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTLaw parse parameters on the tlaw command to set temperature as function of depth,
 * currently only does Bertoldi & Draine simple T law */

/* Despite outward appearances, this routine has no connection whatsoever with
   Bassetlaw (http://www.bassetlawmuseum.org.uk/) */

#include "cddefines.h"
#include "thermal.h"
#include "parser.h"

void ParseTLaw( Parser &p )
{
	DEBUG_ENTRY( "ParseTLaw()" );

	/* this says that some type of temperature law has been specified */
	thermal.lgTLaw = true;
	thermal.lgTemperatureConstant = true;
	thermal.lgTemperatureConstantCommandParsed = true;

	if( p.nMatch("DB96") )
	{
		/* this is to simulate the temperature law given by equation 41 in 
		 * >>refer	H2	temperature law	Draine, B.T., & Bertoldi, Frank, 1996, ApJ, 468, 269-289 */
		thermal.lgTeBD96 = true;

		/* this is the initial temperature for the BD96 temperature law */
		thermal.T0BD96 = 500.f;

		/* the coefficient on column density for temp dropoff */
		thermal.SigmaBD96 = 6e-22f;
	}
	else if( p.nMatch("SN99") )
	{
		/* this is to simulate the temperature law given by equation 16 in 
		 * >>refer	H2	temperature law	Sternberg, A., & Neufeld, D.A. 1999, ApJ, 516, 371-380 */
		thermal.lgTeSN99 = true;

		/* this is the inital temperature for the BD96 temperature law */
		thermal.T0SN99 = 500.f;
	} 
	else if( p.nMatch("TABL") )
	{
		/* when called, read in temperatures from input stream */
		thermal.lgTeTLaw = true;

		p.readLaw(thermal.tlaw);
	}
	else
	{
		fprintf(ioQQQ," There must be a keyword on this command."
			"  The ones I know about are BD96, SN99, and TABLe\n");
		cdEXIT(EXIT_FAILURE);
	}

#if 0
#include "dense.h"
#include "optimize.h"
#include "input.h"
	/* all remainder is currently dead code, a copy of DLAW command,
	 * which could be activated if needs arose */
	/* call fcn dense_fabden(RADIUS) which uses the ten parameters
	 * N.B.; existing version of dense_fabden must be deleted */
	if (0)
	{
		/* this is usual case, call dense_fabden to get density */
		for( long j=0; j < 10; j++ )
		{
			dense.DensityLaw[j] = p.FFmtRead();
		}
		
		/* set flag so we know which law to use later */
		strcpy( dense.chDenseLaw, "DLW1" );
		
		/* vary option */
		if( optimize.lgVarOn )
		{
			/* NB - there are 5 = LIMEXT numbers on this line - if LIMEXT ever changes,
			 * chnage this too */
			strcpy( optimize.chVarFmt[optimize.nparm], "DLAW %f %f %f %f %f " );
			
			/* index for where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			for( log j=0; j<LIMEXT; ++j )
			{
				optimize.vparm[j][optimize.nparm] = (realnum)dense.DensityLaw[j];
			}
			optimize.vincr[optimize.nparm] = 0.5;
			optimize.nvarxt[optimize.nparm] = LIMEXT;
			++optimize.nparm;
		}
	}
#	endif
	return;
}
