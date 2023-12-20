/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include "vectorize.h"

// Define integration methods
typedef enum { Gaussian32, VecGaussian32, Legendre } methods;

static const int numPoints = 32;

ALIGNED(CD_ALIGN) static const double qg32_w[numPoints] = {
	.35093050047350483e-2, .81371973654528350e-2, .12696032654631030e-1, .17136931456510717e-1,
	.21417949011113340e-1, .25499029631188088e-1, .29342046739267774e-1, .32911111388180923e-1,
	.36172897054424253e-1, .39096947893535153e-1, .41655962113473378e-1, .43826046502201906e-1,
	.45586939347881942e-1, .46922199540402283e-1, .47819360039637430e-1, .48270044257363900e-1,
	.48270044257363900e-1, .47819360039637430e-1, .46922199540402283e-1, .45586939347881942e-1,
	.43826046502201906e-1, .41655962113473378e-1, .39096947893535153e-1, .36172897054424253e-1,
	.32911111388180923e-1, .29342046739267774e-1, .25499029631188088e-1, .21417949011113340e-1,
	.17136931456510717e-1, .12696032654631030e-1, .81371973654528350e-2, .35093050047350483e-2
};

ALIGNED(CD_ALIGN) static const double qg32_a[numPoints] = {
	-.498631930924740780,  -.49280575577263417,  -.4823811277937532200, -.46745303796886984000,
	-.448160577883026060,  -.42468380686628499,  -.3972418979839712000, -.36609105937014484000,
	-.331522133465107600,  -.29385787862038116,  -.2534499544661147000, -.21067563806531767000,
	-.165934301141063820,  -.11964368112606854,  -.7223598079139825e-1, -.24153832843869158e-1,
	 .24153832843869158e-1, .7223598079139825e-1, .11964368112606854,    .165934301141063820,
	 .21067563806531767000, .2534499544661147000, .29385787862038116,    .331522133465107600,
	 .36609105937014484000, .3972418979839712000, .42468380686628499,    .448160577883026060,
	 .46745303796886984000, .4823811277937532200, .49280575577263417,    .498631930924740780
};

// define an integrator class.  Currently hard-wired to 32-point Gaussian
template<typename Integrand, methods Method>
class Integrator
{
	const double *p_w, *p_a;
public:
	Integrand func;
	explicit Integrator( const Integrand& fun ) : p_w(qg32_w), p_a(qg32_a), func(fun) {}
	double sum(double min, double max) const
	{
		if( Method == Gaussian32 )
		{
			double a = 0.5*(max+min);
			double b = max-min;
			double total = 0.;
			for( long i=0; i < numPoints; i++ )
				total += b * p_w[i] * func(a + b*p_a[i]);
			return total;
		}
		else
			TotalInsanity();
	}
};

// define a vectorized integrator class.  Currently hard-wired to 32-point Gaussian
template<typename Integrand, methods Method>
class VecIntegrator
{
	const double *p_w, *p_a;
public:
	Integrand func;
	explicit VecIntegrator( const Integrand& fun ) : p_w(qg32_w), p_a(qg32_a), func(fun) {}
	double sum(double min, double max) const
	{
		if( Method == VecGaussian32 )
		{
			double a = 0.5*(max+min);
			double b = max-min;
			ALIGNED(CD_ALIGN) double x[numPoints], y[numPoints];
			for( long i=0; i < numPoints; i++ )
				x[i] = a + b*p_a[i];
			func(x, y, numPoints);
			return b * reduce_ab(p_w, y, 0, numPoints);
		}
		else
			TotalInsanity();
	}
};

 /**			 
 32 point gaussian quadrature integration
 \param xl lower limit to integration
 \param xu - upper limit to integration
 \param (*fct) - pointer to routine to be integrated, arg is x val<BR>
 */ 
double qg32( double, double, double(*)(double) );
/* declar of optimize_func, the last arg, changed from double(*)() to above,
 * seemed to fix flags that were raised */

// Open and closed integration routines, developed on qromb etc. in
// Press et al.
namespace integrate
{
	template <class T>
	class Trapezium
	{
		T m_f;
		double m_a, m_b, m_sum;
		long m_n;
	public:
		static const int NREF=2;
		Trapezium(const T& f, double a, double b) : m_f(f), m_a(a), m_b(b), m_sum(0.0), m_n(0) {}
		void step()
		{
			if (m_n == 0)
			{
				m_sum = 0.5*(m_b-m_a)*(m_f(m_a)+m_f(m_b));
				m_n = 1;
			}
			else
			{
				double rn = 1.0/m_n;
				double si = 0.0;
				for (int i=0; i<m_n; ++i)
				{
					double x = (m_a*(m_n-i-0.5)+m_b*(i+0.5))*rn;
					si += m_f(x);
				}
				m_sum = (m_sum+si*(m_b-m_a)*rn)/NREF;
				m_n *= NREF;
			}
		}
		double sum() const
		{
			return m_sum;
		}
		long evals() const
		{
		return m_n;
		}
	};

	template <class T>
	class Midpoint
	{
		T m_f;
		double m_a, m_b, m_sum;
		long m_n;
	public:
		static const int NREF=3;
		Midpoint(const T& f, double a, double b) : m_f(f), m_a(a), m_b(b), m_sum(0.0), m_n(0) {}
		void step()
		{
			if (m_n == 0)
			{
				m_sum = (m_b-m_a)*m_f(0.5*(m_a+m_b));
				m_n = 1;
			}
			else
			{
				double rn = 1.0/m_n;
				const double sixth = 1./6.;
				double si = 0.0;
				for (int i=0; i<m_n; ++i)
				{
					double x1 = (m_a*(m_n-i-sixth)+m_b*(i+sixth))*rn;
					si += m_f(x1);
					double x2 = (m_a*(m_n-i-1+sixth)+m_b*(i+1-sixth))*rn;
					si += m_f(x2);
				}
				m_sum = (m_sum+si*(m_b-m_a)*rn)/NREF;
				m_n *= NREF;
			}
		}
		double sum() const
		{
			return m_sum;
		}
		long evals() const
		{
			return m_n;
		}
	};
	
	template <class T>
	class Romberg
	{
		T m_f;
		double m_sum, m_dy;
	public:
		explicit Romberg(const T& f) : m_f(f), m_sum(0.0), m_dy(0.0) {}
		void update(double eps);
		double sum() const
		{
			return m_sum;
		}
		double error() const
		{
			return m_dy;
		}
		long evals() const
		{
			return m_f.evals();
		}
	};
	
	template <class T>
	inline void Romberg<T>::update(double eps)
	{
		const int itmax=40/m_f.NREF, npt=5;
		double d1[npt], d0[npt-1] = {};
		double y = 0.0;
		const double w1 = m_f.NREF*m_f.NREF, w2 = 1.0/w1;
		for (int i=0; i<npt; ++i)
		{
			d1[i] = 0.;
		}
		for (int i=0; i<itmax; ++i)
		{
			m_f.step();
			y = m_f.sum();
			
			int l = (i<npt-1) ? i : npt-1;
			for (int m=0; m<l; ++m)
			{
				d0[m] = d1[m];
			}
			d1[0] = y;
			double fr = 1.0;
			for (int m=0; m<l; ++m)
			{
				d1[m+1] = (d1[m]-d0[m]*fr)/(w1-fr);
				y += d1[m+1];
				fr *= w2;
			}
			
			m_dy = fabs(m_sum-y);
			m_sum = y;
			//printf("Level %d value %g->%g error %g samples %ld\n",l,d1[0],m_sum,m_dy,m_f.evals());
			if ( i > 2 && m_dy <= eps*fabs(y) )
			{
				return;
			}
		}
		return;
	}

	template <class T>
	class Simple
	{
		T m_f;
		double m_sum, m_dy;
	public:
		explicit Simple(const T& f) : m_f(f), m_sum(0.0), m_dy(0.0) {}
		void update(double eps);
		double sum() const
		{
			return m_sum;
		}
		double error() const
		{
			return m_dy;
		}
		long evals() const
		{
			return m_f.evals();
		}
	};
	
	template <class T>
	inline void Simple<T>::update(double eps)
	{
		const int itmax=40/m_f.NREF;
		double coarse=0.,fine=0.;
		for (int i=0; i<itmax; ++i)
		{
			coarse = fine;
			m_f.step();
			fine = m_f.sum();
						
			m_dy = fabs(coarse-fine);
			m_sum = fine;
			// printf("Level %d value %g error %g samples %ld\n",i,m_sum,m_dy,m_f.evals());
			if ( i > 2 && m_dy <= eps*fabs(fine) )
			{
				return;
			}
		}
		return;
	}
}
#endif /* INTEGRATE_H_ */
