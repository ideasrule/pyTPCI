/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "iter_track.h"
#include "thirdparty.h"

////////////////////////////////////////////////////////////////////////////////
// The class iter_track and routine Amsterdam_Method were derived from this   //
// original source: http://www.mymathlib.com/roots/amsterdam.html             //
//                                                                            //
// The original code was heavily modified by: Peter van Hoof, ROB             //
//                                                                            //
// The Amsterdam method is more commonly known as the:                        //
//                  van Wijngaarden-Dekker-Brent method                       //
////////////////////////////////////////////////////////////////////////////////


realnum secant_track::next_val(realnum current, realnum next_est)
{
	realnum diff = next_est-current;
	realnum proposal = next_est;
	if (m_lastnext != 0.0 && diff != 0.0)
	{
		realnum lastdiff = m_lastnext-m_lastcurr;
		realnum numer = m_lastcurr*next_est-m_lastnext*current;
		realnum denom = diff-lastdiff;
		if (fabs(numer) >= 2.0*fabs(current*denom))
		{
			m_confidence *= 0.5;
		}
		else
		{
			proposal = numer/denom;
			//fprintf(ioQQQ,"!!!! %g %g %g %g\n",proposal,current,next_est,diff);
			//fflush(ioQQQ);
			if ((proposal-current)*(proposal-next_est) > 0.2*diff*diff)
			{
				proposal = next_est; 
				m_confidence *= 0.5;
			}
			else
			{
				m_confidence *= 2.0;
			}
		}
		m_confidence = min(1.0,max(0.25,m_confidence));

		proposal = m_confidence*proposal+(1.0-m_confidence)*current;
	}
	m_lastnext = next_est;
	m_lastcurr = current;
	return proposal;
}


double iter_track::next_val()
{
	// If the function at the left endpoint is positive, and the function //
	// at the right endpoint is negative.  Iterate reducing the length    //
	// of the interval by either bisection or quadratic inverse           //
	// interpolation.  Note that the function at the left endpoint        //
	// remains nonnegative and the function at the right endpoint remains //
	// nonpositive.                                                       //

	if( p_y(p_a) > 0.0 ) 
	{
		// Check that we are converging or that we have converged near //
		// the left endpoint of the inverval.  If it appears that the  //
		// interval is not decreasing fast enough, use bisection.      //
		if( (p_x(p_b)-p_x(p_a)) < p_tol )
		{
			if( p_y(p_b) > 0 )
				p_a = p_b;
			else
				p_set_root(p_x(p_b));
			return p_midpoint();
		}

		// Check that we are converging or that we have converged near //
		// the right endpoint of the inverval.  If it appears that the //
		// interval is not decreasing fast enough, use bisection.      //
		if( (p_x(p_c)-p_x(p_b)) < p_tol )
		{
			if( p_y(p_b) < 0 )
				p_c = p_b;
			else
				p_set_root(p_x(p_b));
			return p_midpoint();
		}

		// If quadratic inverse interpolation is feasible, try it. //

		if(  ( p_y(p_a) > p_y(p_b) ) && ( p_y(p_b) > p_y(p_c) ) )
		{
			double delta = p_denominator(p_y(p_a),p_y(p_b),p_y(p_c));
			if( delta != 0.0 )
			{
				double dab = p_x(p_a)-p_x(p_b);
				double dcb = p_x(p_c)-p_x(p_b);
				delta = safe_div( p_numerator(dab,dcb,p_y(p_a),p_y(p_b),p_y(p_c)), delta );

				// Will the new estimate of the root be within the   //
				// interval?  If yes, use it and update interval.    //
				// If no, use the bisection method.                  //

				if( delta > dab && delta < dcb )
				{
					if( p_y(p_b) > 0.0 )
						p_a = p_b;
					else if( p_y(p_b) < 0.0 )
						p_c = p_b;
					else
						p_set_root(p_x(p_b));
					return p_x(p_b) + delta;
				}
			}   
		}

		// If not, use the bisection method. //

		if( p_y(p_b) > 0.0 )
			p_a = p_b;
		else
			p_c = p_b;
		return p_midpoint();
	}
	else 
	{
		// If the function at the left endpoint is negative, and the function //
		// at the right endpoint is positive.  Iterate reducing the length    //
		// of the interval by either bisection or quadratic inverse           //
		// interpolation.  Note that the function at the left endpoint        //
		// remains nonpositive and the function at the right endpoint remains //
		// nonnegative.                                                       //

		if( (p_x(p_b)-p_x(p_a)) < p_tol )
		{ 
			if( p_y(p_b) < 0 )
				p_a = p_b;
			else
				p_set_root(p_x(p_b));
			return p_midpoint();
		}

		if( (p_x(p_c)-p_x(p_b)) < p_tol )
		{
			if( p_y(p_b) > 0 )
				p_c = p_b;
			else
				p_set_root(p_x(p_b));
			return p_midpoint();
		}

		if(  ( p_y(p_a) < p_y(p_b) ) && ( p_y(p_b) < p_y(p_c) ) )
		{
			double delta = p_denominator(p_y(p_a),p_y(p_b),p_y(p_c));
			if( delta != 0.0 )
			{
				double dab = p_x(p_a)-p_x(p_b);
				double dcb = p_x(p_c)-p_x(p_b);
				delta = safe_div( p_numerator(dab,dcb,p_y(p_a),p_y(p_b),p_y(p_c)), delta );
				if( delta > dab && delta < dcb )
				{
					if( p_y(p_b) < 0.0 )
						p_a = p_b;
					else if( p_y(p_b) > 0.0 )
						p_c = p_b;
					else
						p_set_root(p_x(p_b));
					return p_x(p_b) + delta;
				}
			}
		}

		if( p_y(p_b) < 0.0 )
			p_a = p_b;
		else
			p_c = p_b;
		return p_midpoint();
	}
}

double iter_track::deriv(int n, double& sigma) const
{
	n = min( n, p_history.size() );
	ASSERT( n >= 2 );
	double *x = new double[n];
	double *y = new double[n];
	for( int i=0; i < n; ++i )
	{
		x[i] = p_x(p_history.size() - n + i);
		y[i] = p_y(p_history.size() - n + i);
	}
	double a,b,siga,sigb;
	linfit( n, x, y, a, siga, b, sigb );
	delete[] y;
	delete[] x;
	sigma = sigb;
	return b;
}

double iter_track::zero_fit(int n, double& sigma) const
{
	n = min( n, p_history.size() );
	ASSERT( n >= 2 );
	double *x = new double[n];
	double *y = new double[n];
	for( int i=0; i < n; ++i )
	{
		x[i] = p_y(p_history.size() - n + i);
		y[i] = p_x(p_history.size() - n + i);
	}
	double a,b,siga,sigb;
	linfit( n, x, y, a, siga, b, sigb );
	delete[] y;
	delete[] x;
	sigma = siga;
	return a;
}

////////////////////////////////////////////////////////////////////////////////
//  double Amsterdam_Method( double (*f)(double), double a, double fa,        //
//                           double c, double fc, double tolerance,           //
//                           int max_iterations, int *err)                    //
//                                                                            //
//  Description:                                                              //
//     Estimate the root (zero) of f(x) using the Amsterdam method where      //
//     'a' and 'c' are initial estimates which bracket the root i.e. either   //
//     f(a) > 0 and f(c) < 0 or f(a) < 0 and f(c) > 0.  The iteration         //
//     terminates when the zero is constrained to be within an interval of    //
//     length < 'tolerance', in which case the value returned is the best     //
//     estimate that interval.                                                //
//                                                                            //
//     The Amsterdam method is an extension of Mueller's successive bisection //
//     and inverse quadratic interpolation.  Later extended by Van Wijnaarden,//
//     Dekker and still later by Brent.  Initially, the method uses the two   //
//     bracketing endpoints and the midpoint to estimate an inverse quadratic //
//     to interpolate for the root.  The interval is successively reduced by  //
//     keeping endpoints which bracket the root and an interior point used    //
//     to estimate a quadratic.  If the interior point becomes too close to   //
//     an endpoint and the function has the same sign at both points, the     //
//     interior point is chosen by the bisection method.  If inverse          //
//     quadratic interpolation is not feasible, the new interior point is     //
//     chosen by bisection, and if inverse quadratic interpolation results    //
//     in a point exterior to the bracketing interval, the new interior point //
//     is again chosen by bisection.                                          //
//                                                                            //
//  Arguments:                                                                //
//     double *f         Pointer to function of a single variable of type     //
//                       double.                                              //
//     double  a         Initial estimate.                                    //
//     double  fa        function value at a                                  //
//     double  c         Initial estimate.                                    //
//     double  fc        function value at c                                  //
//     double  tolerance Desired accuracy of the zero.                        //
//     int     max_iterations The maximum allowable number of iterations.     //
//     int     *err      0 if successful, -1 if not, i.e. if f(a)*f(c) > 0,   //
//                       -2 if the number of iterations > max_iterations.     //
//                                                                            //
//  Return Values:                                                            //
//     A zero contained within the interval (a,c).                            //
//                                                                            //
//  Example:                                                                  //
//     {                                                                      //
//        double f(double), zero, a, fa, c, fc, tol = 1.e-6;                  //
//        int err, max_iter = 20;                                             //
//                                                                            //
//        (determine lower bound, a, and upper bound, c, of a zero)           //
//        fa = f(a);  fc = f(c);                                              //
//        zero = Amsterdam_Method( f, a, fa, c, fc, tol, max_iter, &err);     //
//        ...                                                                 //
//     }                                                                      //
//     double f(double x) { define f }                                        //
////////////////////////////////////////////////////////////////////////////////

double Amsterdam_Method( double (*f)(double), double a, double fa, double c, double fc,
			 double tol, int max_iter, int& err )
{
	iter_track track;

	double result;
	set_NaN( result );

	track.set_tol(tol);

	// If the initial estimates do not bracket a root, set the err flag.  //
   	if( ( err = track.init_bracket( a, fa, c, fc ) ) < 0 )
		return result;

	double b = 0.5*(a + c);
	for( int i = 0; i < max_iter && !track.lgConverged(); i++ )
	{
		track.add( b, (*f)(b) );
		b = track.next_val();
	}

	if( track.lgConverged() )
	{
		err = 0;
		result = track.root();
	}
	else
	{
		err = -2;
	}
	return result;
}
