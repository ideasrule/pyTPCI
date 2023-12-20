/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*vary_input sets input lines to feed into cloudy in optimization runs */
#include "cddefines.h"
#include "input.h"
#include "save.h"
#include "grid.h"

void vary_input(bool *lgLimOK,
		int grid_index)
{
	long int i, 
	  np;

	DEBUG_ENTRY( "vary_input()" );

	// this would indicate int overflow, but is mainly there to keep the compiler from
	// whining about an unused variable...
	if( grid_index < -1 )
		TotalInsanity();

	/* set up chCardSav(n) array like Gary's input file, using input
	 * variable parameters p(i), and format information held in
	 * the common block /parmv/. Results written to common /kardsv/.
	 */

	/* will be set false if limit to a variable exceeded
	 * this is returned to calling code as problem indication*/
	*lgLimOK = true;

	if( cpu.i().lgMaster() || !grid.lgGrid )
		fprintf( ioQQQ, " **************************************************\n" );

	/* echo the variable input lines for this run */
	for( i=0; i < optimize.nvary; i++ )
	{
		bool lgLimitHit = false;

		np = optimize.nvfpnt[i];

		// check if the keyword _LOG is present; the optimizer may not work
		// correctly if it is not optimizing logarithmic quantities.
		//
		// exceptions are the commands ILLUMINATE and RATIO since they vary
		// quantities of order unity anyway, and the commands DLAW and FUDGE
		// since they are entirely defined by the user.
		//
		// it is ok not to convert to upper case first since the command line
		// image is completely under our own control.
		if( !optimize.lgOptimizeAsLinear[i] )
		{
			if( !nMatch( " LOG", optimize.chVarFmt[i] ) )
			{
				fprintf( ioQQQ, " vary_input: internal error - keyword _LOG not found!\n" );
				TotalInsanity();
			}
		}

		/* write formatted to the character string chCardSav(np),
		 * using the format held in chVarFmt(np) */

		/* >>chng 05 aug 09, by RP, both were == change to > and < */
		if( grid.paramIncrements[i] >= 0. &&
		    ( optimize.vparm[0][i] < optimize.varang[i][0] ||
		      optimize.vparm[0][i] > optimize.varang[i][1] ) )
		{
			*lgLimOK = false;
			lgLimitHit = true;
		}
		if( grid.paramIncrements[i] < 0. &&
		    ( optimize.vparm[0][i] > optimize.varang[i][0] ||
		      optimize.vparm[0][i] < optimize.varang[i][1] ) )
		{
			*lgLimOK = false;
			lgLimitHit = true;
		}

		/* now generate the actual command with parameters */
		input.crd[np]->chCardSav = MakeInputLine(i);

		if( cpu.i().lgMaster() || !grid.lgGrid )
		{
			fprintf( ioQQQ, " %s\n", input.crd[np]->chCardSav.c_str() );
			if( lgLimitHit )
				fprintf( ioQQQ, " >>> Limit to variable exceeded.\n" );
		}
	}

	if( cpu.i().lgMaster() && grid.lgGrid )
	{
		// write the line images to an input script, one file for each grid point
		fstream io;
		string fnam = GridPointPrefix(grid_index) + save.chRedirectPrefix + ".in";
		open_data( io, fnam, mode_w );
		for( size_t i=0; i < input.crd.size(); ++i )
			if( input.crd[i]->InclLevel == 0 )
				io << input.crd[i]->chCardSav << endl;
	}

	return;
}
