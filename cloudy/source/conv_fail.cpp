/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvFail handle convergence failure */
#include "cddefines.h"
#include "prt.h"
#include "phycon.h"
#include "hextra.h"
#include "pressure.h"
#include "dense.h"
#include "thermal.h"
#include "called.h"
#include "hcmap.h"
#include "conv.h"

/*ConvFail handle convergence failure - aborts if too many failures occur */
void ConvFail(
	/* chMode is one of "pres", "chem", "eden", "ioni", "pops", "grai", "temp" */
	const char chMode[], /* chMode[5] */
	/* chDetail - string giving details about the convergence failure */
	const char chDetail[] )
{
	double relerror;

	DEBUG_ENTRY( "ConvFail()" );

	/* pressure failure */
	if( strcmp( chMode , "pres" )==0 )
	{
		/* record number of pressure failures */
		++conv.nPreFail;
		if( called.lgTalk )
		{
			fprintf( ioQQQ, 
				" PROBLEM  ConvFail %li, pressure not converged; itr %li, zone %.2f Te:%.3e Hden:%.4e curr Pres:%.4e Error:%.4e%% Pra/gas:%.3e\n", 
			  conv.nPreFail,
			  iteration,
			  fnzone, 
			  phycon.te, 
			  dense.gas_phase[ipHYDROGEN], 
			  pressure.PresTotlCurr, 
			  pressure.PresTotlError*100.,
			  pressure.pbeta);

			/* this identifies new dynamics that failed near the sonic point */
			if( fabs(pressure.PresGasCurr - pressure.PresRamCurr)/pressure.PresGasCurr < 0.1 &&
				 strcmp(dense.chDenseLaw,"DYNA") == 0 )
			{
				fprintf( ioQQQ, 
					"\n PROBLEM continued, pressure not converged; we are stuck at the sonic point.\n\n");
				pressure.lgSonicPoint = true;
			}
		}
	}

	/* electron density failure */
	else if( strcmp( chMode, "eden" ) == 0 )
	{
		/* record number of electron density failures */
		++conv.nNeFail;

		if( called.lgTalk )
		{
			fprintf( ioQQQ, 
				" PROBLEM  ConvFail %li, eden not converged itr %li zone %li fnzone %.2f correct=%.3e "
				"assumed=%.3e.", 
			  conv.nNeFail,
			  iteration ,
			  nzone ,
			  fnzone,
			  dense.EdenTrue, 
			  dense.eden
			  );

			/* some extra information that may be printed */
			/* heating cooling failure */
			if( !conv.lgConvTemp )
			{
				fprintf( ioQQQ, "  Temperature failure also." );
			}

			/* heating cooling failure */
			if( !conv.lgConvIoniz() )
			{
				fprintf( ioQQQ, "  Ionization failure also." );
			}
		}
		fprintf( ioQQQ, " \n");
	}

	else if( strcmp( chMode, "ioni" ) == 0 )
	{
		/* ionization failure */
		++conv.nIonFail;
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " PROBLEM  ConvFail %li, %s ionization not converged"
				" iteration %li zone %li fnzone %.2f reason %s BadConvIoniz0:%g [1]=%g\n", 
			  conv.nIonFail, 
			  chDetail ,
			  iteration ,
			  nzone,
			  fnzone ,
			  conv.chConvIoniz(),
			  conv.convIonizOldVal(),
			  conv.convIonizNewVal());
		}
	}

	else if( strcmp( chMode, "pops" ) == 0 )
	{
		/* populations failure */
		++conv.nPopFail;
		conv.lgConvPops = false;
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " PROBLEM  ConvFail %li, %s population not converged"
				" iteration %li zone %li fnzone %.2f %s %g %g\n", 
			  conv.nPopFail, 
			  chDetail ,
			  iteration,
			  nzone , 
			  fnzone  ,
			  conv.chConvIoniz(),
			  conv.convIonizOldVal(),
			  conv.convIonizNewVal());
		}
	}

	else if( strcmp( chMode, "grai" ) == 0 )
	{
		/* ionization failure */
		++conv.nGrainFail;
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " PROBLEM  ConvFail %ld, a grain failure occurred"
				" iteration %li zone %li fnzone  %.2f %s %g %g\n", 
			  conv.nGrainFail, 
			  iteration , 
			  nzone ,
			  fnzone ,
			  conv.chConvIoniz(),
			  conv.convIonizOldVal(),
			  conv.convIonizNewVal());
		}
	}

	/* rest of routine is temperature failure */
	else if( strcmp( chMode, "temp" ) == 0 )
	{
		ASSERT( fabs((thermal.htot - thermal.ctot)/thermal.htot ) > conv.HeatCoolRelErrorAllowed );
		++conv.nTeFail;
		if( called.lgTalk )
		{
			fprintf( ioQQQ, 
				" PROBLEM  ConvFail %ld, Temp not converged itr %li zone %li fnzone %.2f Te=%.4e"
				" Htot=%.3e Ctot=%.3e rel err=%.3e rel tol:%.3e\n", 
			  conv.nTeFail, 
			  iteration ,
			  nzone ,
			  fnzone, 
			  phycon.te, 
			  thermal.htot, 
			  thermal.ctot, 
			  (thermal.htot - thermal.ctot)/ thermal.htot,
			  conv.HeatCoolRelErrorAllowed );

			/* not really a temperature failure, but something else */
			if( !conv.lgConvIoniz() )
			{
				fprintf( ioQQQ, " Solution not converged due to %10.10s\n", 
							conv.chConvIoniz() );
			}
		}
	}
	else
	{
		fprintf( ioQQQ, " ConvFail called with insane mode %s detail %s\n", 
		  chMode , 
		  chDetail );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* increment total number of failures */
	++conv.nTotalFailures;

	/* now see how many total failures we have, and if it is time to abort */
	/* remember which zone this is */
	conv.ifailz[MIN2(conv.nTotalFailures,10)-1] = nzone;

	/* remember the relative error
	 * convert to single precision for following max, abs (vax failed here) */
	relerror = fabs(safe_div(thermal.htot-thermal.ctot, thermal.htot,0.0));

	conv.failmx = MAX2(conv.failmx,(realnum)MIN2(double(BIGFLOAT),relerror));

	/* this branch is non-abort exit - we have not exceeded the limit to the number of failures */
	if( conv.nTotalFailures < conv.LimFail )
	{ 
		return;
	}

	fprintf( ioQQQ, " Stop due to excessive convergence failures - there have been %ld so far. \n", 
		conv.LimFail );
	fprintf( ioQQQ, " This limit can be reset with the FAILURES command.\n" );

	/* check whether went into cold neutral gas without cosmic rays */
	if( phycon.te < 1e3 && dense.eden/dense.gas_phase[ipHYDROGEN] < 0.1 &&
		(hextra.cryden == 0.) )
	{
		fprintf( ioQQQ,"\n This problem may be solved by adding cosmic rays.\n");
		fprintf( ioQQQ,"\n The gas was cold and neutral.\n");
		fprintf( ioQQQ,"\n The chemistry is not designed to work without a source of ionization.\n");
		fprintf( ioQQQ, " >>> Add galactic background cosmic rays with the COSMIC RAYS BACKBOUND command and try again.\n\n" );
	}

	/* if due to pressure failures then recommend looking at pressure map */
	if( conv.nPreFail==conv.nTotalFailures )
	{
		fprintf( ioQQQ, " These were all pressure failures - we may be near an unstable point in the cooling curve. \n");
		fprintf( ioQQQ, " The PUNCH PRESSURE HISTORY command will show the n-T-P curve, and may be interesting.\n\n");
	}

	/* punt */
	if( conv.lgMap )
	{
		/* only do map if requested */
		/* adjust range of punting map */
		hcmap.RangeMap[0] = (realnum)(phycon.te/100.);
		hcmap.RangeMap[1] = (realnum)MIN2(phycon.te*100.,9e9);
		/* need to make printout out now, before disturbing solution with map */
		PrtZone();
		hcmap.lgMapBeingDone = true;
		map_do(ioQQQ,"punt");
	}

	if( called.lgTalk )
	{
		fprintf( ioQQQ, " ConvFail aborts since nTotalFailures=%ld is >= LimFail=%ld\n", 
		  conv.nTotalFailures, 
		  conv.LimFail );
		fprintf( ioQQQ, " This limit can be reset with the FAILURES command.\n");
		fflush( ioQQQ );
	}
	throw cloudy_abort("too many convergence failures");
}
