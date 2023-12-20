/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "iter_track.h"

namespace {
	// test a function that would have gone into a limit cycle
	TEST(IterTrackBasicFloat)
	{
		iter_track_basic<sys_float> track;
		sys_float x = 1.f;
		sys_float xnew = 2.f;
		while( abs(x-xnew) > 2.f*FLT_EPSILON*abs(x) )
		{
			x = xnew;
			// derivative at root is -1
			xnew = track.next_val( x, 2.f/x );
		}
		CHECK( fp_equal( xnew, (sys_float)sqrt(2.) ) );
		// now clear the tracker and try again
		track.clear();
		x = 1.f;
		xnew = 2.f;
		while( abs(x-xnew) > 2.f*FLT_EPSILON*abs(x) )
		{
			x = xnew;
			// derivative at root is -1
			xnew = track.next_val( x, 3.f/x );
		}
		CHECK( fp_equal( xnew, (sys_float)sqrt(3.) ) );
	}

	// now test the same for doubles
	TEST(IterTrackBasicDouble)
	{
		iter_track_basic<double> track;
		double x = 1.;
		double xnew = 2.;
		while( abs(x-xnew) > 2.*DBL_EPSILON*abs(x) )
		{
			x = xnew;
			// derivative at root is -1
			xnew = track.next_val( x, 2./x );
		}
		CHECK( fp_equal( xnew, sqrt(2.) ) );
	}

	// test a function that would have diverged
	TEST(IterTrackBasicUnstable)
	{
		iter_track_basic<double> track;
		double x = 1.;
		double xnew = 2.;
		while( abs(x-xnew) > 2.*DBL_EPSILON*abs(x) )
		{
			x = xnew;
			// derivative at root is -5
			xnew = track.next_val( x, 1./x - 2.*x );
		}
		// possible solutions are +/-sqrt(1/3), either one may be found
		CHECK( fp_equal( abs(xnew), 1./sqrt(3.) ) );
	}

	// test a function that would have converged anyway
	TEST(IterTrackBasicStableNeg)
	{
		iter_track_basic<double> track;
		double x = 1.;
		double xnew = 2.;
		// this would have diverged without iter_tracking
		while( abs(x-xnew) > 2.*DBL_EPSILON*abs(x) )
		{
			x = xnew;
			// derivative at root is -1/3
			xnew = track.next_val( x, 1./x + x/3. );
		}
		CHECK( fp_equal( xnew, sqrt(1.5) ) );
	}

	TEST(IterTrackBasicStablePos)
	{
		iter_track_basic<double> track;
		double x = 1.;
		double xnew = 2.;
		// this would have diverged without iter_tracking
		while( abs(x-xnew) > 2.*DBL_EPSILON*abs(x) )
		{
			x = xnew;
			// derivative at root is 1/3
			xnew = track.next_val( x, 1./x + 2.*x/3. );
		}
		CHECK( fp_equal( xnew, sqrt(3.) ) );
	}

	double testfun(double x)
	{
		return sin(x)-0.5;
	}

	TEST(IterTrack)
	{
		double x1, fx1, x2, fx2, x3, fx3;
		iter_track track;
		x1 = 0.;
		fx1 = testfun(x1);
		x3 = 1.5;
		fx3 = testfun(x3);
		double tol = 1.e-12;
		track.set_tol(tol);
		CHECK_EQUAL( 0, track.init_bracket(x1,fx1,x3,fx3) );
		CHECK( fp_equal( track.bracket_width(), abs(x3-x1) ) );
		CHECK( !track.lgConverged() );
		x2 = 0.5*(x1+x3);
		fx2 = testfun(x2);
		track.add( x2, fx2 );
		double xnew = track.next_val(0.01);
		CHECK( fp_equal_tol( abs(xnew/x2-1.), 0.01, 1.e-12 ) );
		const int navg = 5;
		vector<double> xvals( navg ); // keep track of the last navg x-values
		for( int i=0; i < 100 && !track.lgConverged(); ++i )
		{
			x2 = track.next_val();
			fx2 = testfun(x2);
			track.add( x2, fx2 );
			// use xvals as circular buffer
			xvals[i%navg] = x2;
		}
		CHECK( track.lgConverged() );
		double exact_root = asin(0.5);
		CHECK( fp_equal_tol( track.root(), exact_root, tol ) );
		double sigma;
		double val = track.deriv( navg, sigma );
		double delta_lo = *min_element( xvals.begin(), xvals.end() ) - exact_root;
		double delta_hi = *max_element( xvals.begin(), xvals.end() ) - exact_root;
		CHECK( delta_lo < 0. );
		CHECK( delta_hi > 0. );
		// the exact derivative at the root is sqrt(3)/2 = 0.8660254...
		// the exact 2nd derivative at the root is -1/2
		double err_lo = -0.5*delta_lo;
		double err_hi = -0.5*delta_hi;
		CHECK( fp_bound( sqrt(3.)/2.+err_hi, val, sqrt(3.)/2.+err_lo ) );
		// the tangent at the exact root is given by asin(0.5) + sqrt(3)/2*(x-x0)
		// if we subtract that from the Taylor expansion of testfun we get:
		// residual = -1/4*(x-x0)^2 + O((x-x0)^3), hence sigma should be less
		// than the maximum absolute value of -1/4*(x-x0)^2 (the actual fit
		// should run slightly closer to the maximum deviant value than the tangent).
		CHECK( sigma < max( pow2(err_lo), pow2(err_hi) ) );
		// ask for more points than are available to see if that is handled correctly
		val = track.deriv( 200 );
		double val2 = track.deriv();
		CHECK( fp_equal( val, val2 ) );
		val = track.deriv( 200, sigma );
		double sigma2;
		val = track.deriv( sigma2 );
		CHECK( fp_equal( sigma, sigma2 ) );

		// now do the same thing for the zero_fit() methods...
		val = track.zero_fit( navg, sigma );
		// the exact root is asin(0.5) = 0.52359877...
		CHECK( fp_equal_tol( exact_root, val, 2.*sigma ) );
		val = track.zero_fit( 200 );
		val2 = track.zero_fit();
		CHECK( fp_equal( val, val2 ) );
		val = track.zero_fit( 200, sigma );
		val = track.zero_fit( sigma2 );
		CHECK( fp_equal( sigma, sigma2 ) );
	}

	// this is the short version of the above...
	TEST(AmsterdamMethod)
	{
		double x1, fx1, x2, fx2;
		x1 = 0.;
		fx1 = testfun(x1);
		x2 = 1.5;
		fx2 = testfun(x2);
		double tol = 1.e-12;
		int err = -1;
		double x = Amsterdam_Method( testfun, x1, fx1, x2, fx2, tol, 1000, err );
		CHECK_EQUAL( 0, err );
		CHECK( fp_equal_tol( x, asin(0.5), tol ) );
	}

	// now test some unstable functions
	double testfun2(double x)
	{
		return exp(x)-3.;
	}

	// the derivative at the root is 3
	TEST(AmsterdamMethod2)
	{
		double x1, fx1, x2, fx2;
		x1 = 0.;
		fx1 = testfun2(x1);
		x2 = 3.;
		fx2 = testfun2(x2);
		double tol = 1.e-12;
		int err = -1;
		double x = Amsterdam_Method( testfun2, x1, fx1, x2, fx2, tol, 1000, err );
		CHECK_EQUAL( 0, err );
		CHECK( fp_equal_tol( x, log(3.), tol ) );
	}

	double testfun3(double x)
	{
		return 1./x - 2.*x;
	}
	
	// the derivative at the root is -4
	TEST(AmsterdamMethod3)
	{
		double x1, fx1, x2, fx2;
		x1 = 0.1;
		fx1 = testfun3(x1);
		x2 = 3.;
		fx2 = testfun3(x2);
		double tol = 1.e-12;
		int err = -1;
		double x = Amsterdam_Method( testfun3, x1, fx1, x2, fx2, tol, 1000, err );
		CHECK_EQUAL( 0, err );
		CHECK( fp_equal_tol( x, sqrt(0.5), tol ) );
	}
}
