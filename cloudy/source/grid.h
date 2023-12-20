/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GRID_H_
#define GRID_H_

#include "optimize.h"
#include "container_classes.h"

/**called by cdDrive */
void grid_do();

/**gridXspec
\param xc[]
\param nInterpVars
*/
void gridXspec(realnum *, long);

/**GridRetrieveXSPECData - obtain the correct spectrum for each grid point */
void GridRetrieveXSPECData(int option);

const int NUM_OUTPUT_TYPES = 11;

struct t_grid
{
	multi_arr<realnum,3> Spectra;
	vector<string> paramNames;
	vector<long> paramMethods;
	multi_arr<realnum,2> paramRange;
	multi_arr<realnum,2> paramData;
	multi_arr<realnum,2> interpParameters;

	realnum paramLimits[LIMPAR][2];
	realnum paramIncrements[LIMPAR];
	vector<realnum> paramValuesFromList[LIMPAR];
	bool lgLinearIncrements[LIMPAR];
	bool lgNegativeIncrements;
	bool lgSaveXspec;

	/** set true when grid command entered */
	bool lgGrid;
	bool lgGridDone;
	bool lgInsideGrid;
	bool lgStrictRepeat;
	/** option to run the grid in parallel mode */
	bool lgParallel;
	/** number of CPUs to be used in parallel grid
	    runs (ignored in MPI mode) */
	unsigned int useCPU;
	/** option to not gather the main output
	    for each grid point into a single file */
	bool lgKeepMainOutputSeparate;
	/** randomly crash some grid points, used for testing */
	bool lgCrash;
	/** has crashing this grid point already been attempted before? */
	bool lgCrashEval;

	/** number of grid commands entered */
	long int nGridCommands;

	long nintparm;
	long naddparm;
	long numParamValues[LIMPAR];
	long totNumModels;

	/** number of times the grid is executed, default is 1 */
	long nCycle;

	bool lgOutputTypeOn[NUM_OUTPUT_TYPES];

	FILE* pnunit;
	long seqNum;

	t_grid()
	{
		lgGridDone = false;
		lgInsideGrid = false;
		lgStrictRepeat = false;
		seqNum = 0;
	}
};
extern t_grid grid;

#endif /* GRID_H_ */
