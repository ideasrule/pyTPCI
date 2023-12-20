/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseDLaw parse parameters on the dlaw command */
#include "cddefines.h"
#include "dense.h"
#include "optimize.h"
#include "input.h"
#include "parser.h"

void ParseDLaw(Parser &p)
{
	long int j;

	DEBUG_ENTRY( "ParseDLaw()" );

	if( dense.gas_phase[ipHYDROGEN] > 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER More than one density command was entered.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* call fcn dense_fabden(RADIUS) which uses the ten parameters
	 * N.B.; existing version of dense_fabden must be deleted
	 * >>chng 96 nov 29, added table option */
	if( p.nMatch("TABL") )
	{
		/* when called, read in densities from input stream */
		strcpy( dense.chDenseLaw, "DLW2" );

		p.readLaw(dense.DLW);
	}
	else if( p.nMatch("WIND") )
	{
		strcpy( dense.chDenseLaw, "DLW3" );
		/* This sets up a steady-state "wind" profile parametrized as in Springmann (1994):
		 *
		 * v(r) = v_star + (v_inf - v_0) * sqrt( Beta1 x + (1-Beta1) x^Beta2 )
		 *
		 * A mass loss rate into 4pi sterradians Mdot then allows the density via continuity:
		 *
		 * n(r) = Mdot / ( 4Pi m_H * mu * r^2 * v(r) )
		 */

		/* The parameters must be specified in this order:
		 *
		 * Mdot, v_inf, Beta2, Beta1, v_0, v_star.
		 *
		 * Only the first three are required.  The final three may be omitted right to left and
		 * take default values Beta1 = v_0 = v_star = 0. */

		for( j=0; j < 6; j++ )
		{
			dense.DensityLaw[j] = p.FFmtRead();
			if( j <= 2 && p.lgEOL() )
				p.NoNumb("density law element");
		}
	}
	else
	{
		/* this is usual case, call dense_fabden to get density */
		for( j=0; j < 10; j++ )
		{
			dense.DensityLaw[j] = p.FFmtRead();
			if( j == 0 && p.lgEOL() )
				p.NoNumb("density law element");
		}

		/* set flag so we know which law to use later */
		strcpy( dense.chDenseLaw, "DLW1" );

		/* vary option */
		if( optimize.lgVarOn )
		{
			ostringstream oss;
			oss << "DLAW %f" << setprecision(7);
			for( j=1; j < 10; j++ )
				oss << " " << dense.DensityLaw[j];
			strcpy( optimize.chVarFmt[optimize.nparm], oss.str().c_str() );
			optimize.lgOptimizeAsLinear[optimize.nparm] = true;

			/* index for where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			optimize.vparm[0][optimize.nparm] = (realnum)dense.DensityLaw[0];
			optimize.vincr[optimize.nparm] = 0.5f;
			optimize.nvarxt[optimize.nparm] = 1;
			++optimize.nparm;
		}
	}

	/* set fake density to signal that density command was entered */
	/* real density will be set once all input commands have been read */
	/* this is necessary since density may depend on subsequent commands */
	dense.SetGasPhaseDensity( ipHYDROGEN, 1.f );

	return;
}
