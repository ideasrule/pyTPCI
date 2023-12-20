/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseConstant parse parameters from the 'constant ...' command */
#include "cddefines.h"
#include "optimize.h"
#include "thermal.h"
#include "dense.h"
#include "pressure.h"
#include "phycon.h"
#include "input.h"
#include "parser.h"
#include "physconst.h"

void ParseConstant(Parser &p )
{
	DEBUG_ENTRY( "ParseConstant()" );

	if( p.nMatch("GRAI") && p.nMatch("TEMP") )
	{
		/* constant grain temperature command */
		thermal.ConstGrainTemp = (realnum)p.FFmtRead();

		/* if linear option is not on the line, convert to exponent if <= 10 */
		if( !p.nMatch("LINE") )
		{
			if( thermal.ConstGrainTemp <= 10. )
				thermal.ConstGrainTemp = exp10(thermal.ConstGrainTemp);
		}

		if( p.lgEOL() )
		{
			p.NoNumb("grain temperature");
		}
	}

	else if( p.nMatch("TEMP") )
	{
		/* a constant temperature model */
		thermal.lgTemperatureConstant = true;
		thermal.lgTemperatureConstantCommandParsed = true;

		/* this is an option to specify the temperature in different units
		 * keV, eV now supported */
		realnum convert_to_Kelvin = 1;
		if( p.nMatch(" EV ") )
		{
			convert_to_Kelvin = (realnum)EVDEGK;
		}
		else if( p.nMatch(" KEV") )
		{
			convert_to_Kelvin = (realnum)(EVDEGK * 1000.);
		}

		/* this is the "force" temperature.  same var used in force temp
		 * command, but lgTSetOn is not set so then allowed to vary 
		 * so constant temperature requires both lgTSetOn true and ConstTemp > 0 */
		thermal.ConstTemp = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("temperature");

		/* if linear option is not on the line and T<=10, assume number is log */
		if( p.nMatch(" LOG") || (thermal.ConstTemp <= 10. && !p.nMatch("LINE")) )
		{
			if( thermal.ConstTemp > log10(BIGFLOAT) )
			{
				fprintf(ioQQQ," PROBLEM temperature entered as a log but is too large "
					"for this processor.  I am interpreting it as the linear temperature.\n");
			}
			else
				thermal.ConstTemp = exp10(thermal.ConstTemp);
		}
		/* do units conversion here */
		thermal.ConstTemp *= convert_to_Kelvin;

		/* check temperature bounds */
		if( thermal.ConstTemp < phycon.TEMP_LIMIT_LOW )
		{
			thermal.ConstTemp = (realnum)(1.0001*phycon.TEMP_LIMIT_LOW);
			fprintf( ioQQQ, " PROBLEM Te too low, reset to %g K.\n",
				 thermal.ConstTemp );
		}
		if( thermal.ConstTemp > phycon.TEMP_LIMIT_HIGH )
		{
			thermal.ConstTemp = (realnum)(0.9999*phycon.TEMP_LIMIT_HIGH);
			fprintf( ioQQQ, " PROBLEM Te too high, reset to %g K.\n",
				 thermal.ConstTemp );
		}

		/* set the real electron temperature to the forced value */
		TempChange( thermal.ConstTemp, false );

		/* vary option */
		if( optimize.lgVarOn )
		{
			/*  no luminosity options on vary */
			optimize.nvarxt[optimize.nparm] = 1;
			// the keyword LOG is not used above, but is checked elsewhere
			strcpy( optimize.chVarFmt[optimize.nparm], "CONSTANT TEMP %f LOG" );

			/*  pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;

			/*  log of temp will be pointer */
			optimize.vparm[0][optimize.nparm] = (realnum)log10(thermal.ConstTemp);
			optimize.vincr[optimize.nparm] = 0.1f;
			optimize.varang[optimize.nparm][0] = (realnum)log10(1.00001*phycon.TEMP_LIMIT_LOW);
			optimize.varang[optimize.nparm][1] = (realnum)log10(0.99999*phycon.TEMP_LIMIT_HIGH);
			++optimize.nparm;
		}
	}

	else if( p.nMatch("DENS") )
	{
		/* constant density */
		strcpy( dense.chDenseLaw, "CDEN" );
		/* turn off radiation pressure */
		pressure.lgPres_radiation_ON = false;
		pressure.lgPres_magnetic_ON = false;
		pressure.lgPres_ram_ON = false;
	}

	else if( p.nMatch("PRES") )
	{
		/* constant pressure  */
		strcpy( dense.chDenseLaw, "CPRE" );

		/* >>chng 06 jun 20, add reset option, to reset density to keep 
		 * initial pressure itself constant from iteration to iteration, 
		 * rather than initial density */
		if( p.nMatch("RESE") )
		{
			/* this says not to keep initial density constant, 
			 * reset it to keep pressure const */
			dense.lgDenseInitConstant = false;
		}
		else
		{
			/* this is default, says keep initial density constant, 
			 * so pressure from iter to iter not really const */
			dense.lgDenseInitConstant = true;
		}

		if( p.nMatch("TIME") )
		{
			// pressure varies as a function of time
			/* this says not to keep initial density constant,
			 * reset it to keep pressure const */
			dense.lgDenseInitConstant = false;
			dense.lgPressureVaryTime = true;

			//  required number is timescale for time variation
			dense.PressureVaryTimeTimescale = (realnum)p.FFmtRead();
			if( dense.PressureVaryTimeTimescale<SMALLFLOAT )
			{
				// need two numbers
				fprintf(ioQQQ," PROBLEM the constant pressure time command requires"
					" a positive timescale.\n");
		  		cdEXIT(EXIT_FAILURE);
			}
			//  optional number is index for time variation
			dense.PressureVaryTimeIndex = (realnum)p.FFmtRead();
			// pressure will be initial pressure - pressure.PresTotlInit - multiplied by
			// (time / time scale ) ^ index
			if( p.lgEOL() )
			{
				// need two numbers
				fprintf(ioQQQ," PROBLEM the constant pressure time command requires"
					" two numbers, the timescale for the variation and an index.\n");
		  		cdEXIT(EXIT_FAILURE);
			}
		}

		if( p.nMatch(" GAS") )
		{
			/*  constant gas pressure (no radiation)
			 *  turn off radiation pressure */
			pressure.lgPres_radiation_ON = false;

			/*  turn off incident continuum */
			pressure.lgContRadPresOn = false;

			/* turn off magnetic and ram pressure */
			pressure.lgPres_magnetic_ON = false;
			pressure.lgPres_ram_ON = false;

			/*  optional number is power law index */
			pressure.PresPowerlaw = (realnum)p.FFmtRead();
		}

		else
		{
			/*  constant total pressure, gas+rad+incident continuum
			 *  turn on radiation pressure */
			pressure.lgPres_radiation_ON = true;
			pressure.lgPres_magnetic_ON = true;
			pressure.lgPres_ram_ON = true;

			/*  option to turn off continuum pressure */
			if( p.nMatch("NO CO") )
			{
				pressure.lgContRadPresOn = false;
			}
			else
			{
				pressure.lgContRadPresOn = true;
			}

			/*  option to not abort when too much radiation pressure, no abort */
			if( p.nMatch("NO AB") )
			{
				pressure.lgRadPresAbortOK = false;
			}
			else
			{
				pressure.lgRadPresAbortOK = true;
			}
			/*  there is no optional power law option for constant total pressure */
			pressure.PresPowerlaw = 0.;

			/* option to set pressure */
			if( p.nMatch(" SET") )
			{
				/* number on line is log of nT - option to specify initial pressure */
				pressure.lgPressureInitialSpecified = true;
				/* this is log of nT product - if not present then set zero */
				pressure.PressureInitialSpecified = p.FFmtRead();
				if( p.lgEOL() )
					p.NoNumb("initial pressure" );
				else
					/* pressure in nkT units */
					pressure.PressureInitialSpecified = exp10(pressure.PressureInitialSpecified)*
						BOLTZMANN;
			}
			else
				pressure.lgPressureInitialSpecified = false;
		}

		/* vary option */
		if( optimize.lgVarOn && pressure.lgPressureInitialSpecified )
		{
			/*  no options on vary */
			optimize.nvarxt[optimize.nparm] = 1;
			// the keyword LOG is not used above, but is checked elsewhere
			strcpy( optimize.chVarFmt[optimize.nparm], "CONSTANT PRESSURE SET %f LOG" );
			if( !dense.lgDenseInitConstant )
				strcat( optimize.chVarFmt[optimize.nparm], " RESET" );
			if( !pressure.lgContRadPresOn )
				strcat( optimize.chVarFmt[optimize.nparm], " NO CONTINUUM" );
			if( !pressure.lgRadPresAbortOK )
				strcat( optimize.chVarFmt[optimize.nparm], " NO ABORT" );

			/*  pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;

			/*  log of temp will be pointer */
			optimize.vparm[0][optimize.nparm] = (realnum)log10(pressure.PressureInitialSpecified/BOLTZMANN);
			optimize.vincr[optimize.nparm] = 0.1f;
			++optimize.nparm;
		}
	}

	else
	{
		/* no keys were recognized */
		fprintf( ioQQQ, " The keyword should be TEMPerature, DENSity, GAS or PRESsure, sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
