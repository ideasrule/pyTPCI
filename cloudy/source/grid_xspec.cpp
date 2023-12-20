/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*gridXspec handles all grid calculations, called by griddo */
/*GridRetrieveXSPECData - obtain the correct spectrum for each grid point */
#include "cddefines.h"
#include "cddrive.h"
#include "continuum.h"
#include "rfield.h"
#include "grid.h"
#include "ipoint.h"
#include "called.h"
#include "prt.h"

/*gridXspec handles all grid calculations, called by grid_do */
void gridXspec(realnum xc[], long int nInterpVars)
{
	long int i;

	DEBUG_ENTRY( "gridXspec()" );

	if( nInterpVars > LIMPAR )
	{
		fprintf( ioQQQ, "grid_do: too many parameters are varied, increase LIMPAR\n" );
		cdEXIT(EXIT_FAILURE);
	}

	optimize.nOptimiz = 0;
	grid.nintparm = nInterpVars;

	/* if this is changed there must be some change made to actually
	 * stuff the additional parameter information. */
	grid.naddparm = 0;

	ASSERT( grid.nintparm + grid.naddparm <= LIMPAR );

	grid.totNumModels = 1;
	/* >>chng 06 aug 21, allow the number of parameter values to be different for different parameters. */
	for( i=0; i<nInterpVars; i++ )
	{
		/* >>chng 06 sep 4, use grid variable instead of passing to routine. */
		grid.totNumModels *= grid.numParamValues[i];
	
		if( grid.paramValuesFromList[i].size() != 0 && grid.lgSaveXspec )
		{
			fprintf( ioQQQ, " PROBLEM The 'save XSPEC' and 'grid list' options do not work together.\n" );
			fprintf( ioQQQ, " If any 'save XSPEC' command is given, all grid commands must follow this syntax:\n" );
			fprintf( ioQQQ, " 	grid <p1> <p2> <p3> \n" );
			fprintf( ioQQQ, " where p1 and p2 are the limits and p3 is a regular increment.\n Sorry.\n" );
			cdEXIT( EXIT_FAILURE );
		}
	}

	// option to cycle through grid multiple times, default is 1
	grid.totNumModels *= grid.nCycle;
	ASSERT( grid.totNumModels > 1 );

	grid.paramNames.resize(nInterpVars+grid.naddparm);
	grid.paramMethods.resize(nInterpVars+grid.naddparm);
	grid.paramRange.alloc(nInterpVars+grid.naddparm, 6);
	grid.paramData.reserve(nInterpVars+grid.naddparm);
	for( i=0; i < nInterpVars+grid.naddparm; i++ )
		grid.paramData.reserve(i, grid.numParamValues[i]);
	grid.paramData.alloc();
	grid.interpParameters.alloc(grid.totNumModels, nInterpVars);

	for( i=0; i < nInterpVars+grid.naddparm; i++ )
	{
		ostringstream oss;
		oss << "PARAM" << i+1;
		grid.paramNames[i] = oss.str();
		/* Method is 0 for linear, 1 for logarithmic */
		grid.paramMethods[i] = 0;
		/* Initial */
		grid.paramRange[i][0] = xc[i]+grid.paramIncrements[i]*(grid.numParamValues[i]-1.f)/2.f;
		/* Delta */
		grid.paramRange[i][1] = grid.paramIncrements[i]/10.f;
		/* Minimum */
		grid.paramRange[i][2] = xc[i];
		/* Bottom */
		grid.paramRange[i][3] = xc[i]+grid.paramIncrements[i]/10.f;
		/* Top */
		grid.paramRange[i][4] = xc[i]+grid.paramIncrements[i]*(grid.numParamValues[i]-1.f)-grid.paramIncrements[i]/10.f;
		/* Maximum */
		grid.paramRange[i][5] = xc[i]+grid.paramIncrements[i]*(grid.numParamValues[i]-1.f);

		for( long j=0; j<grid.numParamValues[i]; j++ )
		{
			grid.paramData[i][j] = xc[i]+grid.paramIncrements[i]*j;
		}
	}

	for( i=0; i< grid.totNumModels; i++ )
	{
		long j;
		realnum variableVector[LIMPAR];

		for( j=0; j<nInterpVars; j++ )
		{
			int index;
			long volumeOtherDimensions = 1;

			/* by "volume", we simply mean the product of the parameter dimensions 
			 * AFTER the present index, i.e.:
			 * first "volume" is product of grid.numParamValues[1]*grid.numParamValues[2]*....grid.numParamValues[n]
			 * second "volume" is product of grid.numParamValues[2]*grid.numParamValues[3]*....grid.numParamValues[n]
			 * last "volume" is unity.  */
			for( long k=j+1; k<nInterpVars; k++ )
			{
				volumeOtherDimensions *= grid.numParamValues[k];
			}

			/* For each successive parameter, the "volume" is less than the previous one.
			 * So the left-hand side of this modulus operation increases for each parameter,
			 * which causes the index of each parameter to be incremented more often than the
			 * index of the previous parameter.  Thus, the last dimension is incremented
			 * every time, then the second to last dimension is incremented with each repeat
			 * of the last dimension.  This repeats until, finally, the first dimension is
			 * incremented.  For example, the indices of the parameter vectors for a 2x2x3
			 * box would be ordered as such:
			 * [0][0][0]
			 * [0][0][1]
			 * [0][0][2]
			 * [0][1][0]
			 * [0][1][1]
			 * [0][1][2]
			 * [1][0][0]
			 * [1][0][1]
			 * [1][0][2]
			 * [1][1][0]
			 * [1][1][1]
			 * [1][1][2]
			 */
			index = (int)( (i/volumeOtherDimensions)%(grid.numParamValues[j]) );

			/* this prevents parameter incrementation for debugging purposes.  */
			if( grid.lgStrictRepeat )
				variableVector[j] = xc[j];
			else if( grid.paramValuesFromList[j].size() != 0U )
				variableVector[j] = grid.paramValuesFromList[j][index];
			else
				variableVector[j] = xc[j] + grid.paramIncrements[j]*index;

			grid.interpParameters[i][j] = variableVector[j];

			if( grid.lgLinearIncrements[j] && !optimize.lgOptimizeAsLinear[j] )
				variableVector[j] = log10(variableVector[j]);
		}

		for( j=nInterpVars; j<LIMPAR; j++ )
		{
			variableVector[j] = xc[j];
		}

		if( i == grid.totNumModels - 1 )
		{
			fixit("is this needed ??");
			called.lgTalk = cpu.i().lgMPI_talk();
			called.lgTalkIsOK = cpu.i().lgMPI_talk();
			prt.lgFaintOn = true;
			grid.lgGridDone = true;
		}

		(void)optimize_func(variableVector,i);
	}
	return;
}

/*GridRetrieveXSPECData - obtain the correct spectrum for each grid point */
void GridRetrieveXSPECData(int option)
{
	DEBUG_ENTRY( "GridRetrieveXSPECData()" );

	if( !grid.lgGrid )
	{
		fprintf( ioQQQ," Cannot save xspec files unless doing a grid.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* allocate some arrays if first call and save continuum energies. */
	if( grid.Spectra.empty() )
	{
		grid.Spectra.reserve(NUM_OUTPUT_TYPES);
		for( long i1=0; i1 < NUM_OUTPUT_TYPES; i1++ )
		{
			if( grid.lgOutputTypeOn[i1] )
			{
				grid.Spectra.reserve(i1,grid.totNumModels);
				for( long i2=0; i2 < grid.totNumModels; i2++ )
				{
					grid.Spectra.reserve(i1,i2,rfield.nflux);
				}
			}
		}
		grid.Spectra.alloc();
	}

	ASSERT( optimize.nOptimiz >= 0 && optimize.nOptimiz < grid.totNumModels );
	ASSERT( grid.lgOutputTypeOn[option] );

	/* grab spectrum for xspec and store it in grid.Spectra */
	cdSPEC2( option, &grid.Spectra[option][optimize.nOptimiz][0]);
}
