/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCosmicRays parse the cosmic rays command */
#include "cddefines.h"
#include "optimize.h"
#include "hextra.h"
#include "ionbal.h"
#include "input.h"
#include "parser.h"

/*ParseCosmicRays parse the cosmic rays command */
void ParseCosmicRays( Parser &p )
{
	int npar = 0;
	realnum a;
	double var;
	string ExtraPar;

	DEBUG_ENTRY( "ParseCosmicRays()" );

	/* cosmic ray density, log of rate relative to background, log of H0 rate in neutral gas,
	 * or density of rel. electrons,
	 * quantity is log unless keyword linear appears */
	/* if no number is present FFmtRead returns zero */
	a = (realnum)p.FFmtRead();
	if( p.lgEOL() )
		a = 0.;

	/* if keyword LINEAR not present, then log, and make linear */
	if( !p.nMatch("LINE") )
		a = exp10(a);
	/* a is now linear scale factor, or linear density, with default of 1 if no number  */

	/* default is cosmic ray ionization rate relative to galactic background, but can
	 * also give density, which was the only option originally */
	if( p.nMatch("DENS") )
	{
		if( p.lgEOL() )
		{
			p.NoNumb("cosmic ray density");
		}
		hextra.cryden = a;

		/*  optional power law density  */
		hextra.crpowr = (realnum)p.FFmtRead();

		/*  option to specify a temp for non-rel electrons - but only when a density */
		hextra.crtemp = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			/* relativistic limit (Balbus and McKee) */
			hextra.crtemp = 2.6e9;
		}
		else
		{
			var = exp10((double)hextra.crtemp);
			hextra.crtemp = (realnum)MIN2(var,2.6e9);
		}
		npar = 3;
		ExtraPar = "DENSITY";
	}
	else if( p.nMatch( "RATE"  ) )
	{
		/* this sets rate - use stored density and rate for background to set
		 * new density since code works with density */
		ASSERT( a > 0. );
		hextra.cryden = hextra.background_density * a / hextra.background_rate;
		hextra.crtemp = 2.6e9f;
		npar = 1;
		ExtraPar = "RATE";
	}
	else if( p.nMatch( "BACKGROU"  ) )
	{
		/* >>chng 06 may 28, require explicit BACKGROUnd to hit background for safety */
		/* cr relative to galactic background BACK - no check on string since default */
		/* >>chng 04 mar 10, background is now 
		 * >>refer	cr	ion	Williams, J.P., Bergin, E.A., Caseli, P., Myers, P.C., & Plume, R. 1998, ApJ, 503, 689 */
		/* galactic background cosmic ray density to produce
		 * secondary ionization rate quoted by Tielens and Hollenbach */
		/* hextra.cryden = 2e-9f;*/
		/* >>chng 99 jun 24, slight change to value
		 * quoted by 
		 * >>refer	cosmic ray	ionization rate	McKee, C.M., 1999, astro-ph 9901370
		 * this will produce a total
		 * secondary ionization rate of 2.5e-17 s^-1, as tested in 
		 * tsuite secondary.in.  If each ionization produces 2.4 eV of heat,
		 * the background heating rate should be 9.6e-29 * n*/
		/* >>chng 00 nov 28, changed density to 4.9e-9 to reproduce TH85a
		 * when photoionization is turned off. 
		 >>refer	cosmic ray	ionization rate	Tielens, A.G.G.M., & Hollenbach, D., 1998, ApJ, 291, 722
		 */
		/* hextra.cryden = 7.07e-9f;*/
		/* this value reproduces the TH cr ionization rate when the factor
		 * of 0.46 is included.  This will directly go onto the h ionization rate
		 * without the factor of 0.46 there.  this is necessary for the more
		 * general case where cr ionization is actually self-consistently determined
		 * from rate hot electrons injected into the plasma */
		/*hextra.cryden = 2.25e-9f;*/
		ASSERT( a > 0. );
		hextra.cryden = hextra.background_density * a;
		hextra.crtemp = 2.6e9f;
		npar = 1;
		ExtraPar = "BACKGROUND";
	}
	else if( p.nMatch( "EQUI"  ) )
	{
		/* equipartition cosmic rays, set from B */
		hextra.lg_CR_B_equipartition = true;
		/* this has to be positive for cr's to be on 
		 * it will be reevaluated when B is known */
		hextra.cryden = SMALLFLOAT;
		hextra.crtemp = 2.6e9f;
	}

	else
	{
		/* no keyword found */
		fprintf( ioQQQ, " There must be a keyword on this COSMIC RAY command.\n" );
		fprintf( ioQQQ, " The keywords are DENSITY, RATE, BACKGROUND, and EQUIPARTITION.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* this is current cosmic ray density divided by background - used in
	 * a few chemical reactions */
	hextra.cryden_ov_background = hextra.cryden / hextra.background_density;
	/* >>chng 05 jan 05, 
	 * set the cr ionization rate to very rough value, before we have enough
	 * information to evaluate it - may be needed in initial guess of H and He ionization*/
	ionbal.CosRayIonRate = hextra.cryden_ov_background * 2.5e-17;

	/* vary option */
	if( optimize.lgVarOn && ExtraPar.length() > 0 )
	{
		/* will be one parameter */
		optimize.nvarxt[optimize.nparm] = npar;
		sprintf( optimize.chVarFmt[optimize.nparm], "COSMic rays %s= %%f LOG", ExtraPar.c_str() );
		/* log of cosmic rays rates relative to background */
		optimize.vparm[0][optimize.nparm] = (realnum)log10(a);
		if( npar == 3 )
		{
			strcat( optimize.chVarFmt[optimize.nparm], " %f %f" );
			optimize.vparm[1][optimize.nparm] = hextra.crpowr;
			optimize.vparm[2][optimize.nparm] = realnum(log10(hextra.crtemp));
		}
		/* array index for where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* the increment in the first steps away from the original value */
		optimize.vincr[optimize.nparm] = 0.2f;
		++optimize.nparm;
	}

	return;
}
