/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvTempEdenIoniz determine  temperature, called by ConPresTempEdenIoniz,
 * calls ConvEdenIoniz to get electron density and ionization */
/*lgConvTemp returns true if heating-cooling is converged */
/*CoolHeatError evaluate ionization, and difference in heating and cooling, for temperature temp */
/*DumpCoolStack helper routine to dump major coolants */
/*DumpHeatStack helper routine to dump major heating agents */
#include "cddefines.h"
#include "hmi.h"
#include "thermal.h"
#include "colden.h"
#include "pressure.h"
#include "dense.h"
#include "trace.h"
#include "phycon.h"
#include "conv.h"
#include "physconst.h"
#include "iter_track.h"
#include "radius.h"

/*lgConvTemp returns true if heating-cooling is converged */
STATIC bool lgConvTemp(const iter_track& TeTrack);
/*CoolHeatError evaluate ionization, and difference in heating and cooling, for temperature temp */
STATIC double CoolHeatError( double temp );

// debugging routines to print main sources of cooling and heating
STATIC void DumpCoolStack(double thres);
STATIC void DumpHeatStack(double thres);

/*ConvTempEdenIoniz determine  temperature, called by ConPresTempEdenIoniz,
 * calls ConvEdenIoniz to get electron density and ionization */
void ConvTempEdenIoniz()
{
	DEBUG_ENTRY( "ConvTempEdenIoniz()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "\n  ConvTempEdenIoniz called\n" );
	}
	if( trace.nTrConvg >= 2 )
	{
		fprintf( ioQQQ, "  ConvTempEdenIoniz called, entering temp loop using solver %s.\n",
			 conv.chSolverTemp );
	}

	// deal with special temperature laws first
	if( thermal.lgTemperatureConstant || thermal.lgTLaw )
	{
		if( thermal.lgTLaw )
		{
			double TeNew = phycon.te;
			if( thermal.lgTeBD96 )
			{
				/* Bertoldi & Drain 96 temp law specified by TLAW BD96 command */
				TeNew = thermal.T0BD96 / (1. + thermal.SigmaBD96 * colden.colden[ipCOL_HTOT]);
			}
			else if( thermal.lgTeSN99 )
			{
				/* Sternberg & Neufeld 99 temp law specified by TLAW SN99 command,
				 * this is equation 16 of 
				 * >>refer	H2	temp	Sternberg, A., & Neufeld, D. A. 1999, ApJ, 516, 371-380 */
				TeNew = thermal.T0SN99 / 
					(1. + 9.*POW4(2.*hmi.H2_total/dense.gas_phase[ipHYDROGEN]) );
			}
			else if( thermal.lgTeTLaw )
			{
				/* Tabulate temperature tlaw */
				TeNew = thermal.tlaw.tabval(radius.Radius, radius.depth);
			}
			else
				TotalInsanity();

			TempChange( TeNew, false );
		}

		ConvEdenIoniz();
		PresTotCurrent();
		
		// convergence is automatic...
		conv.lgConvTemp = true;
		
		if( trace.lgTrace || trace.nTrConvg >= 2 )
		{
			fprintf( ioQQQ, "  ConvTempEdenIoniz: Te %e C %.4e H %.4e\n",
						phycon.te, thermal.ctot, thermal.htot );
			fprintf( ioQQQ, "  ConvTempEdenIoniz returns ok.\n" );
		}
		return;
	}

	// this branch uses the van Wijngaarden-Dekker-Brent method
	if( strcmp( conv.chSolverTemp , "vWDB" ) == 0 )
	{
		conv.lgConvTemp = false;

		// here starts the standard solver for variable temperature
		iter_track TeTrack;
		double t1=0, error1=0, t2, error2;

		t2 = phycon.te;
		error2 = CoolHeatError( t2 );

		for( int n=0; n < 5; ++n )
		{
			const int DEF_ITER = 10;
			const double DEF_FACTOR = 0.2;
			double step, factor = DEF_FACTOR;

			TeTrack.clear();

			step = min( abs(safe_div( error2, conv.dCmHdT, 0. )), factor*t2 );

			// set up an initial guess for the bracket
			// t2 was already initialized outside the main loop, or is copied from the
			// previous iteration. don't record this evaluation, it may be poorly converged
			for( int i=0; i < 100; ++i )
			{
				t1 = t2;
				error1 = error2;

				double maxstep = factor*t1;
				// limited testing on the auto test suite shows that sqrt(2)
				// is close to the optimal value
				step = SQRT2*step;
				if( step == 0.0 || step > maxstep )
					step = maxstep;
				t2 = max( t1 + sign( step, -error1 ), phycon.TEMP_LIMIT_LOW );
				error2 = CoolHeatError( t2 );
				TeTrack.add( t2, error2 );

				// if n > 0, this indicates a previous failure to solve Te
				// this could be due to hysteresis (e.g. O-H charge transfer)
				// so ignore the first n steps, even if they seem to indicate
				// that a bracket is found, to allow the code some time to settle
				if( i >= n && error1*error2 <= 0. )
					break;

				// test for i >= n here to give the code a chance to declare
				// "bracket found" before aborting...
				if( i >= n && fp_equal( t2, phycon.TEMP_LIMIT_LOW ) )
				{
					/* temp is too low */
					fprintf(ioQQQ," PROBLEM DISASTER - the kinetic temperature appears to be below the lower limit of the code,"
							  " %.3eK.  It does not bracket thermal balance.\n",
							  phycon.TEMP_LIMIT_LOW );
					fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n");
					throw cloudy_abort("the kinetic temperature is below the lower limit");
				}
			}

			if( trace.nTrConvg >= 2 && error1*error2 > 0. )
			{
				fprintf( ioQQQ, "  ConvTempEdenIoniz: bracket1 fails t1: %e %e t2: %e %e\n",
					 t1, error1, t2, error2 );
				TeTrack.print_history();
			}

			// keeping the history up until now has a bad effect on convergence
			// so we wipe the slate clean....
			TeTrack.clear();

			// the bracket should have been found, now set up the Brent solver
			if( TeTrack.init_bracket( t1, error1, t2, error2 ) == 0 )
			{
				// The convergence criterion is based on the relative accuracy of Cool-Heat,
				// combined with a relative accuracy on the temperature itself. We need to
				// keep iterating until both accuracies are reached. Here we set tolerance on
				// Te to 2 ulp. If bracket gets narrower than 3 ulp we declare a convergence
				// failure to avoid changes getting lost in machine precision.
				TeTrack.set_tol(2.*DBL_EPSILON*t2);

				if( error1 != 0.0 || error2 != 0.0 )
					t2 = (t1*error2-t2*error1)/(error2-error1);
				else
					t2 = 0.5*(t1+t2);

				for( int i = 0; i < (1<<(n/2))*DEF_ITER; i++ )
				{
					// check for convergence, as well as a pathologically narrow bracket
					if( lgConvTemp(TeTrack) || TeTrack.bracket_width() < 3.*DBL_EPSILON*t2 )
						break;

					error2 = CoolHeatError( t2 );
					TeTrack.add( t2, error2 );
					t2 = TeTrack.next_val(factor);
				}

				if( conv.lgConvTemp )
					break;

				if( trace.nTrConvg >= 2 && !conv.lgConvTemp )
				{
					fprintf( ioQQQ, "  ConvTempEdenIoniz: brent fails\n" );
					TeTrack.print_history();
				}
			}
		}

		// only declare solution unstable if it is at least at the 2-sigma confidence level
		thermal.lgUnstable = ( conv.dCmHdT + 2.*conv.sigma_dCmHdT < 0. );

		if( trace.lgTrace || trace.nTrConvg >= 2 )
		{
			fprintf( ioQQQ, "  ConvTempEdenIoniz: Te %e C %.4e H %.4e (C-H)/H %.2f%%"
				 " d(C-H)/dT %.2e +/- %.2e\n",
				 phycon.te, thermal.ctot, thermal.htot,
				 (thermal.ctot/thermal.htot-1.)*100.,
				 conv.dCmHdT, conv.sigma_dCmHdT );
			fprintf( ioQQQ, "  ConvTempEdenIoniz returns converged=%c\n", TorF(conv.lgConvTemp) );
		}
	}
	else
	{
		fprintf( ioQQQ, "ConvTempEdenIoniz finds insane solver %s\n", conv.chSolverTemp );
		ShowMe();
	}
}


/* returns true if heating-cooling is converged */
STATIC bool lgConvTemp(const iter_track& TeTrack)
{
	DEBUG_ENTRY( "lgConvTemp()" );

	// The explicit test for H-C == 0. is needed since the requirement on the temperature
	// bracket width may not be simultaneously satisfied. Since we have found the zero point
	// exactly, we don't care about that... The temperature is as accurate as it is ever going
	// to be. Not doing the explicit test for H-C == 0. is a bug. If the temparture bracket is
	// too wide when we hit H-C == 0., this algorithm would never converge since vWDB requires
	// that the endpoints of the bracket have non-zero function values, so they cannot get
	// updated and the bracket width never gets smaller. So the requirement on the temperature
	// bracket width would never be satisfied...
	if( thermal.htot - thermal.ctot == 0.
	    || ( abs(thermal.htot - thermal.ctot)/thermal.htot <= conv.HeatCoolRelErrorAllowed &&
		 TeTrack.bracket_width()/phycon.te <= conv.HeatCoolRelErrorAllowed/3. )
	    || thermal.lgTemperatureConstant )
	{
		/* announce that temp is converged if relative heating - cooling mismatch
		 * is less than the relative heating cooling error allowed and the width of
		 * the temperature bracket is sufficiently small (this assures that the
		 * temperature is also well determined if H-C is a shallow function of T).
		 * If this is a constant temperature model, force convergence */
		conv.lgConvTemp = true;
		// remember numerical derivative to estimate initial stepsize on next call
		conv.dCmHdT = TeTrack.deriv(conv.sigma_dCmHdT);
	}
	else
	{
		/* big mismatch, this has not converged */
		conv.lgConvTemp = false;
	}

	if( trace.nTrConvg >= 2 )
		fprintf( ioQQQ, "  lgConvTemp: C-H rel err %.4e Te rel err %.4e converged=%c\n",
			 abs(thermal.htot - thermal.ctot)/thermal.htot,
			 TeTrack.bracket_width()/phycon.te,
			 TorF(conv.lgConvTemp) );

	return conv.lgConvTemp;
}

/*CoolHeatError evaluate ionization, and difference in heating and cooling, for temperature temp */
STATIC double CoolHeatError( double temp )
{
	DEBUG_ENTRY( "CoolHeatError()" );

	static ConvergenceCounter cctr=conv.register_("TEMP_CHANGES");
	++cctr;
	TempChange( temp, false );

	/* converge the ionization and electron density; 
	 * this calls ionize until lgIonDone is true */
	/* NB should NOT set insanity - but rather return error condition */
	ConvEdenIoniz();

	/* >>chng 01 mar 16, evaluate pressure here since changing and other values needed */
	/* reevaluate pressure */
	/* this sets values of pressure.PresTotlCurr */
	PresTotCurrent();

	/* keep track of temperature solver in this zone
	 * conv.hist_temp_nzone is reset in ConvInitSolution */
	if( nzone != conv.hist_temp_nzone )
	{
		/* first time in this zone - reset history */
		conv.hist_temp_nzone = nzone;
		conv.hist_temp_temp.clear();
		conv.hist_temp_heat.clear();
		conv.hist_temp_cool.clear();
	}

	conv.hist_temp_temp.push_back( phycon.te );
	conv.hist_temp_heat.push_back( thermal.htot );
	conv.hist_temp_cool.push_back( thermal.ctot );

	// dump major contributors to heating and cooling - for debugging purposes
	if( false )
	{
		DumpCoolStack( conv.HeatCoolRelErrorAllowed/5.*thermal.ctot );
		DumpHeatStack( conv.HeatCoolRelErrorAllowed/5.*thermal.htot );
	}

	if( trace.nTrConvg >= 2 )
	{
#if 0
		species *spCO = findspecies("CO")->local()->dbase;
		
		fprintf( ioQQQ, "DEBUGGG CO %ld %ld\n",nzone, 
					spCO ? spCO->numLevels_local : -1);
#endif

		fprintf( ioQQQ, "  CoolHeatError: Te: %.4e C: %.4e H: %.4e (C-H)/H: %.4e\n",
			 temp, thermal.ctot, thermal.htot, thermal.ctot/thermal.htot-1. );
	}

	double error = thermal.ctot - thermal.htot;

	// this can get set if temperature drops below floor temperature -> fake convergence
	if( thermal.lgTemperatureConstant )
		error = 0.;

	return error;
}

STATIC void DumpCoolStack(double thres)
{
	multimap<double,string> output;
	char line[200];

	for( int i=0; i < thermal.ncltot; ++i )
	{
		double fraction;
		if( abs(thermal.heatnt[i]) > thres )
		{
			fraction = thermal.heatnt[i]/thermal.ctot;
			sprintf( line, "heat %s %e: %e %e\n",
				 thermal.chClntLab[i], thermal.collam[i], thermal.heatnt[i], fraction );
			output.insert( pair<const double,string>( fraction, string(line) ) );
		}
		if( abs(thermal.cooling[i]) > thres )
		{
			fraction = thermal.cooling[i]/thermal.ctot;
			sprintf( line, "cool %s %e: %e %e\n",
				 thermal.chClntLab[i], thermal.collam[i], thermal.cooling[i], fraction );
			output.insert( pair<const double,string>( fraction, string(line) ) );
		}
	}

	dprintf( ioQQQ, " >>>>>>> STARTING COOLING DUMP <<<<<<\n" );
	dprintf( ioQQQ, "total cooling %e\n", thermal.ctot );
	// this will produce sorted output in reverse order (largest contributor first)
	for( multimap<double,string>::reverse_iterator i=output.rbegin(); i != output.rend(); ++i )
		dprintf( ioQQQ, "%s", i->second.c_str() );
	dprintf( ioQQQ, " >>>>>>> FINISHED COOLING DUMP <<<<<<\n" );
}

STATIC void DumpHeatStack(double thres)
{
	multimap<double,string> output;
	char line[200];

	for( int nelem=0; nelem < LIMELM; ++nelem )
	{
		for( int i=0; i < LIMELM; ++i )
		{
			double fraction = thermal.heating(nelem,i)/thermal.htot;
			if( abs(thermal.heating(nelem,i)) > thres )
			{
				sprintf( line, "heating(%i,%i): %e %e\n",
					 nelem, i, thermal.heating(nelem,i), fraction );
				output.insert( pair<const double,string>( fraction, string(line) ) );
			}
		}
	}

	dprintf( ioQQQ, " >>>>>>> STARTING HEATING DUMP <<<<<<\n" );
	dprintf( ioQQQ, "total heating %e\n", thermal.htot );
	// this will produce sorted output in reverse order (largest contributor first)
	for( multimap<double,string>::reverse_iterator i=output.rbegin(); i != output.rend(); ++i )
		dprintf( ioQQQ, "%s", i->second.c_str() );
	dprintf( ioQQQ, " >>>>>>> FINISHED HEATING DUMP <<<<<<\n" );
}
