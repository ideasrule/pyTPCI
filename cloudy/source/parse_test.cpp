/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTest parse the test command */
#include "cddefines.h"
#include "iso.h"
#include "monitor_results.h"
#include "parser.h"
#include "lines.h"

inline void StuffCommand( const string& chCmd, Parser& p, bool lgPrint)
{
	if( lgPrint )
		fprintf(ioQQQ, "%s\n", chCmd.c_str());
	p.setline(chCmd);
	// is this needed?
	if( chCmd.substr(0,7) != "MONITOR" )
		p.set_point(4);
}

void ParseTest(Parser &p)
{
	DEBUG_ENTRY( "ParseTest()" );

	/* do a smoke test - 
	* generate test input stream for rapid testing of code */

	/* optional keyword PRINT will cause all the commands to be printed */
	bool lgPrintTest = p.nMatch("PRIN");

	bool lgH2 = p.nMatch(" H2 ");
	bool lgLARG = p.nMatch("LARG");
	bool lgMOLE = p.nMatch("MOLE");

	/* option to turn on the H2 molecule */
	if( lgH2 )
	{
		/* this will both enable the molecule and also set the threshold for actually
		* computing it to a very low value, so that it is actually used in ionized gas */
		StuffCommand( "DATABASE H2 LIMIT -20", p, lgPrintTest );
		ParseDatabaseH2(p);
	}

	/* use largest possible hydrogen atom */
	if( lgLARG )
	{
		StuffCommand( "DATABASE H-LIKE ELEMENT HYDROGEN LEVELS LARGER", p, lgPrintTest );
		ParseDatabaseISO(ipH_LIKE,p);
	}

	if( lgMOLE )
	{
#		if 0
		{
			/* hydrogen density */
			StuffCommand( "TRACE TEMPERATURE CONVERGENCE", p, lgPrintTest );
			ParseTrace(chStuff);
		}
#		endif

		/* hydrogen density */
		StuffCommand( "HDEN 5", p, lgPrintTest );
		ParseHDEN(p);

		/* make a constant temperature model */
		StuffCommand( "CONSTANT TEMPER 50K", p, lgPrintTest );
		ParseConstant(p);

		/* continuum to include full energy range */
		StuffCommand( "TABLE ISM", p, lgPrintTest );
		ParseTable(p);

		/* extinguish this continuum */
		StuffCommand( "EXTINGUISH 23", p, lgPrintTest );
		ParseExtinguish(p);
		
		/* stop in second zone, so we do use the zone increment logic */
		StuffCommand( "STOP ZONE 2", p, lgPrintTest );
		ParseStop(p);

		/* set thickness */
		StuffCommand( "SET DR 0", p, lgPrintTest );
		ParseSet(p);
		
		/* do Case B so lyman line pumping of H is not important */
		StuffCommand( "CASE B", p, lgPrintTest );
		ParseCaseB(p);

		/* set cosmic rays */
		//>>chng 12 apr 09, this uses old Williams CR backbround for continuity
		StuffCommand( "COSMIC RAY BACKGROUND linear 0.1266", p, lgPrintTest );
		ParseCosmicRays(p);

		/* add grains and ISM abundances */
		StuffCommand( "ABUNDANCES ISM", p, lgPrintTest );
		ParseAbundances(p);

		/* add grains and ism abundances */
		StuffCommand( "CONSTANT GRAIN TEMPERATURE 20K", p, lgPrintTest );
		ParseConstant(p);

		/* create series of monitor commands */
		StuffCommand( "MONITOR EDEN 0.5597", p, lgPrintTest );
		ParseMonitorResults(p);

		/* create series of monitor commands */
		if( lgH2 )
		{
			/*>>chng 13 apr 13, -1.528 to -1.450 we are not monitoring this in test suite */
			StuffCommand( "MONITOR MOLECULAR FRACTION H2 -1.450 error 0.1", p, lgPrintTest );
		}
		else
		{
			/*>>chng 13 apr 13, -1.528 to -1.542 we are not monitoring this in test suite */
			StuffCommand( "MONITOR MOLECULAR FRACTION H2 -1.542 error 0.1", p, lgPrintTest );
		}
		ParseMonitorResults(p);

		/* create series of monitor commands */
		/*>>chng 13 apr 13, 0.987 to 0.954 we are not monitoring this in test suite */
		/*>>chng 21 nov 04, 0.954 to 0.9946 we are not monitoring this in test suite */
		StuffCommand( "MONITOR COLUMN CO 0.9946 error 0.1", p, lgPrintTest );
		ParseMonitorResults(p);			

		/* create series of monitor commands */
		StuffCommand( "MONITOR TEMPERATURE \"H\" 50K", p, lgPrintTest );
		ParseMonitorResults(p);
	}
	else
	{
		/* hydrogen density */
		StuffCommand( "HDEN 4", p, lgPrintTest );
		ParseHDEN(p);

		/* make a constant temperature model */
		StuffCommand( "CONSTANT TEMPER 4", p, lgPrintTest );
		ParseConstant(p);

		/* continuum to include full energy range */
		StuffCommand( "TABLE AGN", p, lgPrintTest );
		ParseTable(p);

		/* set ionization parameter */
		StuffCommand( "IONIZATION PARAMETER -2", p, lgPrintTest );
		ParseIonParI(p);

		/* use old abundances */
		StuffCommand( "ABUNDANCES OLD SOLAR 84", p, lgPrintTest );
		ParseAbundances(p);

		/* >>chng 02 apr 19, add this */
		/* set total Lyman continuum depth - this is to prevent caution
		 * that Lyman continuum was thin but expected to be thick */
		StuffCommand( "STOP LYMAN OPTICAL -4", p, lgPrintTest );
		ParseStop(p);
		
		/* stop in second zone, so we do use the zone increment logic */
		StuffCommand( "STOP ZONE 2", p, lgPrintTest );
		ParseStop(p);

		/* set thickness */
		StuffCommand( "SET DR 0", p, lgPrintTest );
		ParseSet(p);

		/* create series of monitor commands */
		StuffCommand( "MONITOR HYDROGEN 1 IONIZATION -3.052", p, lgPrintTest );
		ParseMonitorResults(p);

		StuffCommand( "MONITOR HELIUM 2 IONIZATION -1.076", p, lgPrintTest );
		ParseMonitorResults(p);

		/*>>chng 13 feb 01, from -2.377 to -2.319. undo r6703 DR suppression */
		/*>>chng 13 apr 13, from -2.319 to -2.359, redo DR suppression */
		StuffCommand( "MONITOR CARBON 2 IONIZATION -2.359", p, lgPrintTest );
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.653 to -0.560. Badnell DR by default */
		/*>>chng 13 feb 01, from -0.596 to -0.565. undo r6703 DR suppression */
		/*>>chng 13 apr 13, from -0.565 to -0.586, redo DR suppression */
		StuffCommand( "MONITOR CARBON 3 IONIZATION -0.586", p, lgPrintTest );
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.348 to -0.373. Badnell DR by default */
		StuffCommand( "MONITOR CARBON 4 IONIZATION -0.361", p, lgPrintTest );
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.490 to -0.530. Badnell DR by default */
		StuffCommand( "MONITOR CARBON 5 IONIZATION -0.514", p, lgPrintTest );
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.800 to -0.861. Badnell DR by default */
		/*>>chng 11 jul 12, from -0.861 to -0.935, DR coll suppression */
		/*>>chng 13 feb 01, from -0.935 to -0.863. undo r6703 DR suppression */
		/*>>chng 13 apr 13, from -0.863 to -0.894, redo DR suppression */
		/*>>chng 15 may 14, from -0.894 to -0.864, CollisSuppres, final version from Dragan Nikolic */
		/*>>chng 15 oct 14, from -0.864 to -0.894, CollisSuppres, corrected version from Nikolic */
		/*>>chng 16 dec 20, from -0.894 to -0.865, add secondary autoionization to DR suppression */
		StuffCommand( "MONITOR OXYGEN 3 IONIZATION  -0.865", p, lgPrintTest );
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.180 to -0.157. Badnell DR by default */
		StuffCommand( "MONITOR OXYGEN 4 IONIZATION -0.148", p, lgPrintTest );
		ParseMonitorResults(p);

		/*>>chng 06 dec 01, from -0.770 to -0.808. Badnell DR by default */
		/*>>chng 13 feb 01, from -0.790 to -0.807. undo r6703 DR suppression */
		/*>>chng 13 mar 13, from -0.807 to -0.798. undo DR suppression */
		StuffCommand( "MONITOR OXYGEN 5 IONIZATION -0.798", p, lgPrintTest );
		ParseMonitorResults(p);

		/* >>chng 02 apr 19, from 0.7258 to 0.946, due to adding Lyman cont depth */
		/* >>chng 07 oct 22, from 0.946  to 1.108, resolve l-levels of h-like sequence */
		ostringstream wavlen;
		wavlen << setprecision( LineSave.sig_figs ) << Hbeta_WavLen;	 

		StuffCommand( "MONITOR LINE \"CA B\"" + wavlen.str() + " 1.108", p, lgPrintTest );
		ParseMonitorResults(p);

		/* >>chng 02 apr 19, from 2.4603 to 3.25, due to adding Lyman cont depth 
		 * >>chng 06 nov 17, asserted value from 3.25 to 3.11, drift over last year */
		/* >>chng 06 dec 01, from 3.11 to 2.72. Badnell DR by default */
		/* >>chng 07 oct 22, from 2.72 to 3.18, resolve l-levels of h-like sequence */
		/* >>chng 11 jul 12, from 3.18 to 2.695, DR coll suppression */
		/* >>chng 13 feb 01, from 2.695 to 3.189. undo r6703 DR suppression */
		/* >>chng 13 apr 13, from 3.189 to 2.961, redo DR suppression */
		/* >>chng 15 may 14, from 2.961 to 3.275, CollisSuppres, final version from Dragan Nikolic */
		/* >>chng 15 oct 14, from 3.275 to 3.063, CollisSuppres, corrected version from Nikolic */
		/* >>chng 16 dec 20, from 3.063 to 3.207, add secondary autoionization to DR suppression */
		StuffCommand( "MONITOR LINE \"O  3\" 5006.84 3.207", p, lgPrintTest );
		ParseMonitorResults(p);

		StuffCommand( "MONITOR HTOT -15.019", p, lgPrintTest );
		ParseMonitorResults(p);
	}
	
	return;
}
