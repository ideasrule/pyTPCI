/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "monointerp.h"

/** Monotonic piecewise cubic interpolation, see
		>>refer Fritsch & Butland, 1984, SIAM J Sci Stat Comput, 5, 300 */

// Internal utility functions
namespace {
	// Hermite polynomial basis functions
	inline double h00( double t )
	{
		return (1.0+2.0*t)*(1.0-t)*(1.0-t);
	}
	inline double h10( double t )
	{
		return t*(1.0-t)*(1.0-t);
	}
	// Bisection search
	inline long bisect ( const std::vector<double> &x, double xval )
	{
		long n = x.size();
		long ilo = 0, ihi = n-1;
		if (xval <= x[0])
			return -1;
		if (xval >= x[n-1])
			return n-1;
		while (ihi-ilo!=1)
		{
			long imid = (ilo+ihi)/2;
			if (x[imid] > xval)
				ihi = imid;
			else
				ilo = imid;
		}
		return ilo;
	}
}

// Constructor for interpolation function object
Monointerp::Monointerp ( const double x[], const double y[], long n ) 
	: m_x(x,x+n), m_y(y,y+n), m_g(n)
{
	ASSERT(m_x.size() == m_y.size() && m_x.size() == m_g.size());
	std::vector<double> d(n-1),h(n-1);			
	for (long k=0;k<n-1;++k) 
	{
		h[k] = (x[k+1]-x[k]);
		d[k] = (y[k+1]-y[k])/h[k];
	}
	m_g[0] = d[0];
	for (long k=1;k<n-1;++k) 
	{
		m_g[k] = d[k]*d[k-1];
		if (m_g[k] > 0.0)
		{
			double a = (h[k-1]+2.0*h[k])/(3.0*(h[k-1]+h[k]));
			m_g[k] /= (a*d[k]+(1-a)*d[k-1]);
		}
		else
		{
			m_g[k] = 0.0;
		}
	}
	m_g[n-1] = d[n-2];
}

Monointerp::~Monointerp( void )
{
}

// Evaluate interpolant
double Monointerp::operator() ( double xval ) const
{
	double yval;
	
	if( xval <= m_x[0] )
	{
		yval = m_y[0];
	}
	else if( xval >= m_x[m_x.size()-1] )
	{
		yval = m_y[m_x.size()-1];
	}
	else
	{
		long k = bisect( m_x, xval );
		double h = m_x[k+1]-m_x[k], t = (xval-m_x[k])/h;
		yval = m_y[k]*h00(t) + h*m_g[k]*h10(t) 
			+ m_y[k+1]*h00(1.0-t) - h*m_g[k+1]*h10(1.0-t);
	}
	return yval;
}
