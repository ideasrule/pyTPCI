/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseStop parse the stop command */
#include "cddefines.h"
#include "optimize.h"
#include "phycon.h"
#include "predcont.h"
#include "rfield.h"
#include "geometry.h"
#include "iterations.h"
#include "stopcalc.h"
#include "input.h"
#include "parser.h"
#include "flux.h"
#include "service.h"

void ParseStop(Parser &p)
{
	long int j;

	double effcol, 
	  tread;

	DEBUG_ENTRY( "ParseStop()" );

	/* time option, for stopping time dependent calculations, used to stop
	 * iterations rather than zones.  Only some stop commands have this option */
	bool lgStopZone = true;
	if( p.nMatch("TIME") )
		lgStopZone = false;

	if( p.nMatch("TEMP") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("temperature");
		}

		/* off option disables this stopping criterion */
		if( p.lgEOL() && p.nMatch(" OFF") )
		{
			/* this is special case for ending temperature - do not use -
			 * but will still stop if Te falls below TeLowest, the lowest
			 * possible temperature */
			if( lgStopZone )
				StopCalc.TempLoStopZone = -1.f;
			else
				StopCalc.TempLoStopIteration = -1.f;
		}
		else
		{
			/* lowest electron temperature allowed before stopping
			 * assumed to be the log of the temperature if <10
			 * optional keyword LINEAR forces linear */
			if( a <= 10. && !p.nMatch("LINE") )
			{
				tread = exp10(a);
			}
			else
			{
				tread = a;
			}

			/* tread is linear temperature*/
			if( tread < phycon.TEMP_LIMIT_LOW )
			{
				fprintf( ioQQQ, 
					" Temperatures below %.2e K not allowed. Reset to lowest value."
					"   I am doing this myself.\n" ,
					phycon.TEMP_LIMIT_LOW );
				/* set slightly off extreme limit for safety */
				tread = phycon.TEMP_LIMIT_LOW*1.01;
			}
			else if( tread > phycon.TEMP_LIMIT_HIGH )
			{
				fprintf( ioQQQ, 
					" Temperatures is above %.2e K not allowed. Reset to highest value."
					"   I am doing this myself.\n" ,
					phycon.TEMP_LIMIT_HIGH);
				/* set slightly off extreme limit for safety */
				tread = phycon.TEMP_LIMIT_HIGH*0.99;
			}

			if( p.nMatch("EXCE") )
			{
				/* option for this to be highest allowed temperature,
				 * stop temperate exceeds */
				if( lgStopZone )
					StopCalc.TempHiStopZone = (realnum)tread;
				else
					StopCalc.TempHiStopIteration = (realnum)tread;
			}
			else
			{
				/* this is ending temperature - we stop if kinetic temperature
				 * falls below this */
				if( lgStopZone )
					StopCalc.TempLoStopZone = (realnum)tread;
				else
					StopCalc.TempLoStopIteration = (realnum)tread;
			}
		}
	}

	/* stop at 21cm line center optical depth  */
	else if( p.nMatch("OPTI") && p.nMatch("21CM") )
	{
		/* default is for number to be log of optical depth */
		bool lgLOG = true;
		if( p.nMatch("LINE") )
		{
			/* force it to be linear not log */
			lgLOG = false;
		}
		j = (long int)p.FFmtRead();
		if( j!=21 )
		{
			fprintf( ioQQQ, " First number on STOP 21CM OPTICAL DEPTH command must be 21\n" );
			cdEXIT(EXIT_FAILURE);
		}
		/* now get the next number, which is the optical depth */
		double a = (long int)p.FFmtRead();

		/* tau must be a log, and second number is energy where tau specified */
		if( lgLOG )
		{
			StopCalc.tauend = (realnum)exp10(a);
		}
		else
		{
			StopCalc.tauend = (realnum)a;
		}
		/* this flag says that 21cm line optical depth is the stop quantity */
		StopCalc.lgStop21cm = true;
	}
	/* stop optical depth at some energy */
	else if( p.nMatch("OPTI") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("optical depth");
		}

		/* default is for number to be log of optical depth */
		bool lgLOG = true;
		if( p.nMatch("LINE") )
		{
			/* force it to be linear not log */
			lgLOG = false;
		}

		/* tau entered as a log unless liner specified */
		if( lgLOG )
		{
			if( a > 37. )
			{
				fprintf( ioQQQ, " The log optical depth of %.2e is too big, the largest is log t = 37\n",
						a);
				cdEXIT(EXIT_FAILURE);
			}
			StopCalc.tauend = (realnum)exp10(a);
		}
		else
		{
			StopCalc.tauend = (realnum)a;
		}

		/* energy where tau specified */
		StopCalc.taunu = (realnum)p.FFmtRead();

		if( p.lgEOL() )
		{
			if( p.nMatch("LYMA") )
			{
				/* largest Lyman limit optical depth */
				StopCalc.taunu = 1.;
			}
			else if( p.nMatch("BALM") )
			{
				/* stop at this Balmer continuum optical depth */
				StopCalc.taunu = 0.25;
			}
			else
			{
				fprintf( ioQQQ, " There must be a second number, the energy in Ryd.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		else
		{
			/* if second number is negative then log of energy in rydbergs */
			if( StopCalc.taunu < 0. )
			{
				StopCalc.taunu = exp10(StopCalc.taunu);
			}

			/* check that energy is within bounds of code */
			if( StopCalc.taunu < rfield.emm() || StopCalc.taunu > rfield.egamry() )
			{
				fprintf( ioQQQ, " The energy must be in the range %10.2e to %10.2e.  It was %10.2e. Sorry.\n", 
				  rfield.emm(), rfield.egamry(), StopCalc.taunu );
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* vary option */
		if( optimize.lgVarOn )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "STOP OPTICAL DEPTH = %f LOG AT %f RYD" );
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.tauend);
			optimize.vparm[1][optimize.nparm] = (realnum)StopCalc.taunu;
			optimize.vincr[optimize.nparm] = 0.5;
			optimize.nvarxt[optimize.nparm] = 2;
			++optimize.nparm;
		}
	}

	/* stop optical depth at extinction in V filter */
	else if( p.nMatch(" AV ") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("optical depth in V");
		}
		/* default is for number to be A_V, log if negative */
		if( a<=0. )
		{
			a = exp10(a);
		}
		/* A_V can be for either point or extended source, difference is (1-g) multiplied by scat opacity
		 * if keyword point occurs then for point source, otherwise for extended source */
		if( p.nMatch("EXTE" ) )
		{
			StopCalc.AV_extended = (realnum)a;
		}
		else
		{
			/* default is point, as measured in ism work */
			StopCalc.AV_point = (realnum)a;
		}
	}

	/* stop when a fraction of molecules frozen out on grain surfaces is reached */
	else if( p.nMatch("MOLE") && p.nMatch("DEPL") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("molecular depletion");
		}
		if( a <= 0. )
		{
			StopCalc.StopDepleteFrac = (realnum)exp10(a);
		}
		else
		{
			StopCalc.StopDepleteFrac = (realnum)a;
		}
	}

	/* stop when absolute value of flow velocity falls below this value */
	else if( p.nMatch("VELO") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("flow velocity");
		}
		/* entered in km/s but stored as cm/s */
		StopCalc.StopVelocity = (realnum)(a*1e5);
	}

	/* stop at a given computed mass */
	else if( p.nMatch("MASS") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("mass");
		}
		/* number of log of mass in gm if inner radius is specified,
		 * mass per unit area, gm cm-2 if not
		 * leave it as a log since dare not deal with linear mass */
		StopCalc.xMass = (realnum)a;
		/* NB 0 is sentinel for not set, if a is zero we must reset it */
		if( StopCalc.xMass == 0 )
			StopCalc.xMass = SMALLFLOAT;

		/* vary option */
		if( optimize.lgVarOn )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "STOP MASS = %f LOG" );
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			optimize.vparm[0][optimize.nparm] = (realnum)StopCalc.xMass;
			optimize.vincr[optimize.nparm] = 0.5;
			optimize.nvarxt[optimize.nparm] = 1;
			++optimize.nparm;
		}
	}

	/* stop thickness command, also stop depth, this must come after stop
	 * optical depth, since do not want to trigger on depth in optical depth */
	else if( p.nMatch("THIC") || p.nMatch("DEPT") || p.nMatch("RADI") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("distance");
		}
		const double convl = p.nMatch("PARS") ? log10( PARSEC ) : 0.;
		bool lgStopRadius = p.nMatch("RADI") ? true : false ;
		const char* what = lgStopRadius ? "radius" : "thickness";

		if( p.nMatch("LINE") )
		{
			if( a > 0. )
			{
				a = log10(a) + convl;
			}
			else
			{
				fprintf(ioQQQ,"The %s is negative and linear is set - this is impossible.\n", what);
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			a += convl;
		}
		if( a > 37. )
		{
			fprintf( ioQQQ, "DISASTER %s too large\n", what );
			cdEXIT(EXIT_FAILURE);
		}
		if( lgStopRadius )
			iterations.StopRadius[0] = exp10(a);
		else
			iterations.StopThickness[0] = exp10(a);

		/* can stop at different thickness on each iteration */
		for( j=1; j < iterations.iter_alloc; j++ )
		{
			a = p.FFmtRead();
			if( p.lgEOL() )
			{
				if( lgStopRadius )
					iterations.StopRadius[j] = iterations.StopRadius[j-1];
				else
					iterations.StopThickness[j] = iterations.StopThickness[j-1];
			}
			else
			{
				if( p.nMatch("LINE") )
				{
					if( a > 0. )
					{
						a = log10(a) + convl;
					}
					else
					{
						fprintf(ioQQQ,"The %s is negative and linear is set -"
							" this is impossible.\n", what);
						cdEXIT(EXIT_FAILURE);
					}
				}
				else
				{
					a += convl;
				}
				if( a > 37. )
				{
					fprintf( ioQQQ, "DISASTER %s too large\n", what );
					cdEXIT(EXIT_FAILURE);
				}
				if( lgStopRadius )
					iterations.StopRadius[j] = exp10(a);
				else
					iterations.StopThickness[j] = exp10(a);
			}
		}

		/* vary option */
		if( optimize.lgVarOn )
		{
			optimize.nvarxt[optimize.nparm] = 1;
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			if( lgStopRadius )
			{
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP RADIUS %f LOG" );
				optimize.vparm[0][optimize.nparm] = (realnum)log10(iterations.StopRadius[0]);
			}
			else
			{
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP THICKNESS %f LOG" );
				optimize.vparm[0][optimize.nparm] = (realnum)log10(iterations.StopThickness[0]);
			}
			optimize.vincr[optimize.nparm] = 0.5f;
			++optimize.nparm;
		}
	}

	/* stop at a particular zone, for each iteration */
	else if( p.nMatch("ZONE") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("zone number");
		}
		/* stop after computing this zone */
		/* >>chng 03 jun 06, do not let fall below 1, stop zone 0 has same effect
		 * as stop zone 1, bug caught by Joop Schaye */
		iterations.nend[0] = (long)MAX2(1.,a);
		geometry.lgZoneSet = true;

		/* this tells code that we intend to stop at this zone, so caution not generated*/
		iterations.lgEndDflt = false;

		long int nZoneMax = iterations.nend[0];
		for( j=1; j < iterations.iter_alloc; j++ )
		{
			iterations.nend[j] = (long)p.FFmtRead();
			/* if eol on this iteration, set to previous.  In most cases
			 * all will be equal to the first */
			if( p.lgEOL() )
			{
				iterations.nend[j] = iterations.nend[j-1];
			}
			else
			{
				/* do not let fall below 1, stop zone 0 has same effect
				 * as stop zone 1, bug caught by Joop Schaye */
				iterations.nend[j] = MAX2( 1 , iterations.nend[j] );
			}
			nZoneMax = max( nZoneMax , iterations.nend[j] );
		}

		if( nZoneMax>2000 )
			fprintf(ioQQQ,"CAUTION - it will take a lot of memory to save"
			" results for %li zones.  Is this many zones really necessary?\n",
			nZoneMax );
	}

	/* stop when a prescribed continuum flux is reached */
	else if( p.nMatch("CONT") && p.nMatch("FLUX") )
	{
		/* first read the continuum energy and add this point to PredCont */
		double energy = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("energy");
		const char* unit = p.StandardEnergyUnit();
		long ind = t_PredCont::Inst().add( energy, unit );
		Energy E( energy, unit );

		double flux = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("flux");
		if( flux <= 0. || p.nMatch( " LOG") )
			flux = exp10(flux);
		Flux F( E, flux, p.StandardFluxUnit() );

		StopCalc.ContIndex.push_back( ind );
		StopCalc.ContNFnu.push_back( F );
	}

	/* stop at this electron fraction, relative to hydrogen */
	else if( p.nMatch("EFRA") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("electron fraction");
		}
		if( a <= 0. )
		{
			StopCalc.StopElecFrac = (realnum)exp10(a);
		}
		else
		{
			StopCalc.StopElecFrac = (realnum)a;
		}
	}

	/* stop at a hydrogen molecular fraction, relative to total hydrogen,
	 * this is 2H_2 / H_total*/
	else if( p.nMatch("MFRA") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("hydrogen molecular fraction");
		}
		if( a <= 0. )
		{
			StopCalc.StopH2MoleFrac = (realnum)exp10(a);
		}
		else
		{
			StopCalc.StopH2MoleFrac = (realnum)a;
		}
	}

	/* stop at a ionized hydrogen fraction, relative to total hydrogen,
	 * this is H+ / H_total */
	else if( p.nMatch("PFRA") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("ionized hydrogen fraction");
		}
		if( a <= 0. )
		{
			StopCalc.StopHPlusFrac = (realnum)exp10(a);
		}
		else
		{
			StopCalc.StopHPlusFrac = (realnum)a;
		}
	}

	/* stop at a particular column density */
	else if( p.nMatch("COLU") )
	{
		string chLabel;
		// species label may contain numbers, so needs to be parsed first...
		bool lgFoundLabel = !p.GetQuote( chLabel );
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("column density");
		}

		/* check for linear option, if present take log since a being
		 * log column density is default */
		if( p.nMatch( "LINE" ) )
			a = log10(a);

		if( lgFoundLabel )
		{
			/* species was given in double quotes */
			StopCalc.col_species = (realnum)exp10(a);
			StopCalc.lgStopSpeciesColumn = true;
			StopCalc.chSpeciesColumn = chLabel;
			trimTrailingWhiteSpace( StopCalc.chSpeciesColumn );
		}
		/*  stop at an effective column density */
		else if( p.nMatch("EFFE") )
		{
			/* actually stop at certain optical depth at 1keV */
			effcol = exp10(a);
			StopCalc.tauend = (realnum)(effcol*2.14e-22);
			StopCalc.taunu = (realnum)(1000./EVRYD);
			/* vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP EFFECTIVE COLUMN DENSITY %f LOG" );
				/* pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/* log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(effcol);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}

		else if( p.nMatch("IONI") )
		{
			/*  this is ionized column */
			if( a > 37. )
			{
				fprintf( ioQQQ, " column too big\n" );
				cdEXIT(EXIT_FAILURE);
			}

			StopCalc.colpls = (realnum)exp10(a);

			/*  vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP IONIZED COLUMN DENSITY %f LOG" );
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/*  log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.colpls);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}

		/*  stop at a neutral column */
		else if( p.nMatch("NEUT") )
		{
			StopCalc.colnut = (realnum)exp10(a);

			/*  vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP NEUTRAL COLUMN DENSITY %f LOG");
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/*  log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.colnut);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}

		/* >>chng 03 apr 15, add this option 
		 *  stop at a molecular hydrogen column density, input parameter
		 * is log of the H2 column density */
		else if( p.nMatch(" H2 ") )
		{
			/* this command has a 2 in the H2 label - must not parse the two by
			 * accident.  Get the first number off the line image, and confirm that
			 * it is a 2 */
			j = (long int)a;
			if( j != 2 )
			{
				fprintf( ioQQQ, " Something is wrong with the order of the numbers on this line.\n" );
				fprintf( ioQQQ, " The first number I encounter should be the 2 in H2.\n Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			a = p.FFmtRead();
			StopCalc.col_h2 = (realnum)exp10(a);

			/*  vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP H2 COLUMN DENSITY %f LOG");
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/*  log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.col_h2);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}

		// atomic hydrogen
		else if( p.nMatch("ATOM") )
		{
			StopCalc.col_h2_nut = (realnum)exp10(a);
			/*  vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP ATOMIC COLUMN DENSITY %f LOG");
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/*  log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.col_h2_nut);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}

		else if( p.nMatch("H/TS") )
		{
			/* >> 05 jan 09, add stop integrated n(H0) / Tspin */
			StopCalc.col_H0_ov_Tspin = (realnum)exp10(a);
			/*  vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP H/TSPIN COLUMN DENSITY %f LOG");
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/*  log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.col_H0_ov_Tspin);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}

		else if( p.nMatch(" CO ") )
		{
			/* chng << 03 Oct. 27--Nick Abel, add this option */
			/* stop at a carbon monoxide column density */
			StopCalc.col_monoxco = (realnum)exp10(a);
			/*  vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP CO COLUMN DENSITY %f LOG");
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/*  log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.col_monoxco);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}

		/* fall through default is total hydrogen column density */
		else
		{
			/*  both HII and HI */
			if( a > 37. )
			{
				fprintf( ioQQQ, " column too big\n" );
				cdEXIT(EXIT_FAILURE);
			}

			StopCalc.HColStop = (realnum)exp10(a);

			/*  vary option */
			if( optimize.lgVarOn )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], "STOP COLUMN DENSITY %f LOG" );
				/*  pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;
				/*  log of temp will be pointer */
				optimize.vparm[0][optimize.nparm] = (realnum)log10(StopCalc.HColStop);
				optimize.vincr[optimize.nparm] = 0.5f;
				++optimize.nparm;
			}
		}
	}

	/* stop when electron density falls below this value, linear or log */
	else if( p.nMatch("EDEN") )
	{
		double a = p.FFmtRead();
		
		if( p.lgEOL() && !p.nMatch(" OFF") )
		{
			p.NoNumb("electron density");
		}
		/* stop if electron density falls below this value
		 * LINEAR option */
		if( p.nMatch("LINE") )
		{
			StopCalc.StopElecDensity = (realnum)a;
		}
		else
		{
			StopCalc.StopElecDensity = (realnum)exp10(a);
		}
	}

	/* stop at a particular line ratio - this must come last since many commands
	 * have linear option - don't want to trigger on that */
	else if( p.nMatch("LINE") )
	{
		string chLabel;
		/* first line wavelength, then intensity relative to Hbeta for stop
		 * if third number is entered, it is wl of line in denominator */

		/* get label for the line - must do this first so we clear the label string before
		 * trying to read the wavelength */
		if (p.GetQuote( chLabel ))
			p.StringError();

		/* copy first four char of label into caps, and null terminate*/
		strncpy( StopCalc.chStopLabel1[StopCalc.nstpl], chLabel.c_str() , NCHLAB-1 );
		StopCalc.chStopLabel1[StopCalc.nstpl][NCHLAB-1] = 0;
		trimTrailingWhiteSpace( StopCalc.chStopLabel1[StopCalc.nstpl] );

		// default intrinsic intensity, accept emergent
		StopCalc.nEmergent[StopCalc.nstpl] = 0;
		if( p.nMatch("EMER") )
			StopCalc.nEmergent[StopCalc.nstpl] = 1;

		/* get line wavelength */
		StopCalc.StopLineWl1[StopCalc.nstpl] = (realnum)p.getWave();

		/* get relative intensity */
		StopCalc.stpint[StopCalc.nstpl] = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There MUST be a relative intensity  entered "
				"for first line in STOP LINE command.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* check for second line - use Hbeta is not specified */

		/* get label for the line - must do this first so we clear the label string before
		 * trying to read the wavelength */
		
		j = p.GetQuote( chLabel );

		if( j != 0 )
		{
			/* no line label, normalization line will be H beta */
			strncpy( StopCalc.chStopLabel2[StopCalc.nstpl], "H  1" , NCHLAB-1 );
			StopCalc.chStopLabel2[StopCalc.nstpl][NCHLAB-1] = 0;
			trimTrailingWhiteSpace( StopCalc.chStopLabel2[StopCalc.nstpl] );
			StopCalc.StopLineWl2[StopCalc.nstpl] = Hbeta_WavLen;
		}
		else
		{
			/* copy first four char of label into caps, and null terminate*/
			strncpy( StopCalc.chStopLabel2[StopCalc.nstpl], chLabel.c_str() , NCHLAB-1 );
			StopCalc.chStopLabel2[StopCalc.nstpl][NCHLAB-1] = 0;
			trimTrailingWhiteSpace( StopCalc.chStopLabel2[StopCalc.nstpl] );

			/* wavelength of second line, may be absent and so zero -
			 * we will use Hbeta if not specified */
			StopCalc.StopLineWl2[StopCalc.nstpl] = (realnum)p.getWaveOpt();

		}
		/* increment number of stop lines commands entered */
		StopCalc.nstpl = MIN2(StopCalc.nstpl+1,MXSTPL-1);
	}

	/* oops! no keyword that we could find */
	else
	{
		fprintf( ioQQQ, " I did not recognize a keyword on this STOP line, line image follows;\n" );
		p.PrintLine(ioQQQ);
		fprintf( ioQQQ, "Sorry.\n");
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
