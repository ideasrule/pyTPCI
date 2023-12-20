/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvEdenIoniz called by ConvTempIonz, calls ConvIoniz solving for eden */
/*lgConvEden returns true if electron density is converged */
/*EdenError evaluate ConvIoniz() until ionization has converged and return error on eden */
#include "cddefines.h"
#include "dense.h"
#include "trace.h"
#include "conv.h"
#include "cooling.h"
#include "thermal.h"
#include "iter_track.h"
#include "rt.h"

/*lgConvEden returns true if electron density is converged */
STATIC bool lgConvEden();
/*EdenError evaluate ConvIoniz() until ionization has converged and return error on eden */
STATIC double EdenError(double eden);

/*ConvEdenIoniz called by ConvTempEdenIoniz, calls ConvIoniz solving for eden */
void ConvEdenIoniz()
{
	DEBUG_ENTRY( "ConvEdenIoniz()" );

	/* this routine is called by ConvTempEdenIoniz, it calls ConvIoniz
	 * and changes the electron density until it converges */

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "\n" );
		fprintf( ioQQQ, "   ConvEdenIoniz entered\n" );
	}
	if( trace.nTrConvg>=3 )
	{
		fprintf( ioQQQ, 
			"   ConvEdenIoniz called, entering eden loop using solver %s.\n",
			conv.chSolverEden);
	}

	/* save entry value of eden */
	double EdenEntry = dense.eden;

	// this branch uses the van Wijngaarden-Dekker-Brent method
	if( strcmp( conv.chSolverEden , "vWDB" )== 0 )
	{
		conv.lgConvEden = false;

		iter_track NeTrack;
		double n1, error1, n2, error2;
		// this is the maximum relative step in eden
		double factor = 0.02;
		bool lgHysteresis = false;

		for( int n=0; n < 3; ++n )
		{
			const int DEF_ITER = 10;
			// if hysteresis is detected, we should lower the maximum
			// step to avoid upsetting the lower solvers by large steps
			if( lgHysteresis )
				factor /= 5.;

			NeTrack.clear();

			// when dense.EdenTrue becomes negative, error1 > n1 (since n1 > 0.)
			// a straight copy EdenTrue -> eden would then imply a step > 100% down
			// all of the code below will cap that to n1*(1.-factor), and eden stays > 0.
			// this is also asserted in EdenError, the ONLY place where dense.eden is set

			n1 = dense.eden;
			error1 = EdenError( n1 );
			NeTrack.add( n1, error1 );

			if( conv.lgSearch && dense.EdenTrue > SMALLFLOAT )
				n2 = sqrt(dense.eden*dense.EdenTrue);
			else if( abs(safe_div( error1, n1 )) < factor )
				n2 = dense.EdenTrue;
			else
				n2 = ( error1 > 0. ) ? n1*(1.-factor) : n1*(1.+factor);

			// n1 == n2 will occur if SET EDEN command was given
			if( !fp_equal( n1, n2 ) )
				error2 = EdenError( n2 );
			else
				error2 = error1;
			NeTrack.add( n2, error2 );

			int j = 0;

			// now hunt until we have bracketed the solution
			while( error1*error2 > 0. && j++ < DEF_ITER )
			{
				n1 = n2;
				error1 = error2;
				double deriv = NeTrack.deriv(5);
				// this can occur in fully stripped conditions where
				// dense.EdenTrue is essentially independent of dense.eden
				if( deriv == 0. )
					deriv = 1.;
				// the factor 1.2 creates 20% safety margin
				double step = safe_div( -1.2*error1, deriv, 0. );
				step = sign( min( abs(step), factor*n1 ), step );
				n2 = n1 + step;
				error2 = EdenError( n2 );
				NeTrack.add( n2, error2 );
			}

			if( error1*error2 > 0. && trace.nTrConvg >= 3 )
			{
				fprintf( ioQQQ, "   ConvEdenIoniz: bracket failure 1  n1: %e %e n2: %e %e\n",
					 n1, error1, n2, error2 );
				NeTrack.print_history();
			}

			// using the derivative failed, so simply start hunting up or downwards
			// we may need to take a big step, so max_iter should be big
			while( error1*error2 > 0. && j++ < 20*DEF_ITER )
			{
				n1 = n2;
				error1 = error2;
				n2 = ( error1 > 0. ) ? n1*(1.-factor) : n1*(1.+factor);
				error2 = EdenError( n2 );
				NeTrack.add( n2, error2 );
			}

			if( error1*error2 > 0. && trace.nTrConvg >= 3 )
			{
				fprintf( ioQQQ, "   ConvEdenIoniz: bracket failure 2  n1: %e %e n2: %e %e\n",
					 n1, error1, n2, error2 );
				NeTrack.print_history();
			}

			NeTrack.clear();

			// the bracket should have been found, now set up the Brent solver
			if( NeTrack.init_bracket( n1, error1, n2, error2 ) == 0 )
			{
				int nBound = 0;

				// set tolerance to 2 ulp; if bracket gets narrower than 3 ulp we declare
				// a convergence failure to avoid changes getting lost in machine precision
				NeTrack.set_tol(2.*DBL_EPSILON*n2);

				double NeNew = 0.5*(n1+n2);
				for( int i = 0; i < (1<<(n/2))*DEF_ITER; i++ )
				{
					// check for convergence, as well as a pathologically narrow bracket
					if( lgConvEden() || NeTrack.bracket_width() < 3.*DBL_EPSILON*n2 )
						break;

					NeTrack.add( NeNew, EdenError( NeNew ) );
					NeNew = NeTrack.next_val(factor);

					// this guards against hysteresis. the symptom of hysteresis is
					// that EdenTrue ends up being outside the bracket consistently.
					// if this happens several times in a row, break out of this loop
					int nVal = NeTrack.in_bounds(dense.EdenTrue);
					if( nVal == 0 )
						nBound = 0;
					else
						nBound += nVal;
					if( abs(nBound) >= 3 )
					{
						lgHysteresis = true;
						if( trace.nTrConvg >= 3 )
							fprintf( ioQQQ, "   ConvEdenIoniz: hysteresis detected\n" );
						break;
					}
				}
			}

			if( conv.lgConvEden )
				break;

			if( trace.nTrConvg >= 3 )
			{
				fprintf( ioQQQ, "   ConvEdenIoniz: brent fails\n" );
				NeTrack.print_history();
			}
		}

	}
	else if( strcmp( conv.chSolverEden , "SECA" )== 0 )
	{
		conv.lgConvEden = false;

		
		double n1=0., error1=0., n2=0., error2=0.;
		const int MAX_ITER = 20;
		for( int n=0; n < MAX_ITER; ++n )
		{
			if ( n == 0 )
			{
				n2 = dense.eden;
				if ( dense.EdenTrue > 0. )
					n2 = sqrt( n2*dense.EdenTrue );
			}
			else if ( n == 1 || (n1-n2)*(error1-error2) <= 0.0 )
			{
				n1 = n2;
				error1 = error2;
				if ( dense.EdenTrue > 0. )
					n2 = sqrt(dense.EdenTrue*dense.eden);
				else
					n2 = dense.eden/2.;
			}
			else
			{
				double nt = (n1*error2 - n2*error1)/(error2-error1);
				if (fabs(nt-n2) > 2*fabs(n2-n1))
					nt = 3*n2 - 2*n1;
				// fprintf(ioQQQ,"3 %lg %lg %lg\n",n1,n2,nt);
				if (nt < 0.1*n2)
					nt = 0.1*n2;
				n1 = n2;
				error1 = error2;
				n2 = nt;
			}				
			error2 = EdenError( n2 );
			
			if (0)
				fprintf(ioQQQ,"LONG Nzone %ld Loop %d density %15.8g true %15.8g error %15.8g\n",
						  nzone, n, n2, dense.EdenTrue, error2);

			if( lgConvEden() )
				break;

		}
	}
	else
	{
		fprintf( ioQQQ, "ConvEdenIoniz finds insane solver %s\n", conv.chSolverEden );
		ShowMe();
	}

	if( trace.lgTrace || trace.nTrConvg >= 3 )
	{
		fprintf( ioQQQ, "   ConvEdenIoniz: entry eden %.4e -> %.4e rel chng %.2f%% accuracy %.2f%%\n",
					EdenEntry, dense.eden, (safe_div(dense.eden,EdenEntry,1.)-1.)*100.,
					(safe_div(dense.eden,dense.EdenTrue,1.)-1.)*100. );
		fprintf( ioQQQ, "   ConvEdenIoniz returns converged=%c reason %s\n",
					TorF(conv.lgConvEden), conv.chConvIoniz() );
	}

	if (!lgConvBaseHeatTest)
	{
		//HeatZero();

		/* get total cooling, thermal.ctot = does not occur since passes as pointer.  This can add heat.
		 * it calls coolSum at end to sum up the total cooling */
		CoolEvaluate( &thermal.ctot );
		
		HeatSum();
	}
}

/* returns true if electron density is converged */
STATIC bool lgConvEden()
{
	conv.lgConvEden = 
		( abs(dense.eden-dense.EdenTrue) < abs(dense.eden)*conv.EdenErrorAllowed )
		|| ( dense.eden > 0 && 
			  abs(dense.eden-dense.EdenTrue) < 1e-15*dense.xNucleiTotal );
	if( !conv.lgConvEden )
	{
		conv.setConvIonizFail( "Ne big chg" , dense.EdenTrue, dense.eden);
	}
	return conv.lgConvEden;
}

/* evaluate ConvIoniz() until ionization has converged and return error on eden */
STATIC double EdenError(double eden)
{
	// this is the only place where the new electron density is set
	static ConvergenceCounter cctr=conv.register_("EDEN_CHANGES");
	++cctr;
	EdenChange( eden );

	//RT_OTS();
	double SumOTS;
	RT_OTS_Update(&SumOTS);
	RT_line_all_escape( NULL );

	for( int i=0; i < 5; ++i )
	{
		ConvIoniz();
		if( conv.lgConvIoniz() )
			break;
	}

	double error = dense.eden - dense.EdenTrue;

	if( trace.nTrConvg >= 3 )
		fprintf( ioQQQ, "   EdenError: eden %.4e EdenTrue %.4e rel. err. %.4e\n",
			 dense.eden, dense.EdenTrue, safe_div(dense.eden,dense.EdenTrue,1.)-1. );

	return error;
}

