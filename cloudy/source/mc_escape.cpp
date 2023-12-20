/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "mc_escape.h"

#include "integrate.h"
#include "rt_escprob.h"
#include "thirdparty.h"
#include "ran.h"

namespace
{
	static void mc_table(bool lgHeader, double tau, double a, double pesc1, 
								double nbar, double lbar, double taul, double pdest)
	{
		if (lgHeader)
		{
			//       1234567890---1234567890---1234567890---1234567890--
			fprintf(ioQQQ,
					  "#   tau1/2        pesc1         <N>          <L>"
					  "        <taul>"
					  "        pdest       esca0k2      esctot\n");
		}
		else
		{
			fprintf(ioQQQ, 
					  "%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
					  tau, pesc1, nbar, lbar, taul, pdest,esca0k2(tau/SQRTPI),
					  esc_CRDwing_1side(tau/SQRTPI,a));
		}
		fflush(ioQQQ);
	}
	
	static double phigauss(double nu)
	{
		return exp(-nu*nu)/SQRTPI;
	}
	
	// Implement equation (157) of Avrett & Loeser 1966
	// 1966SAOSR.201.....A for comparison with escape probability with
	// damping.  Uses transformation x -> 1/y-1 to map domain of
	// integral to (0,1]
	class k2DampArg
	{
		realnum damp, tau;
	public:
		k2DampArg(realnum damp, realnum tau) : damp(damp), tau(tau) {}
		
		realnum operator()(realnum y) const
			{
				if (y <= 0.0)
					return 0.0;
				realnum x = 1./y-1.;
				realnum phi;
				VoigtU(damp,&x,&phi,1);
				if (phi <= 0.)
				{
					return 0.;
				}
				else
				{
					return phi*expn(2,phi*tau)/POW2(y);
				}
			}
	};

	class overlap
	{
		realnum damp, beta;
	public:
		overlap(realnum damp, realnum beta) : damp(damp), beta(beta) {}
		
		realnum operator()(realnum x) const
			{
				realnum phi;
				if (0)
				{
					phi = (fabs(x)<0.5);
				}
				else
				{
					VoigtU(damp,&x,&phi,1);
				}
				return beta*phi/(beta+phi);
			}
	};
}				  
	

// Simple complete redistribution Monte Carlo
// tau_in is half-depth in mean optical depth scale
void mc_escape(double tau_in, double a, double beta)
{
	// Parameters
	// beta = k_c/k_L -- Hummer & Kunasz 1980

	fprintf(ioQQQ, "# a = %.4e, beta = %.4e\n", a, beta);
	if (0)
	{
		// Confirm distributions
		const long NSAMP = 1e7;
		double totnu1 = 0.0, totnu2 = 0.0, totphi = 0.0;
		for (long i=0; i<NSAMP; ++i)
		{
			double nu = ran.normal()/sqrt(2.);
			totnu1 += nu*nu;
			double nu2 = (10.*(i+0.5))/NSAMP;
			double phi = phigauss(nu2);
			totphi += phi;
			totnu2 += nu2*nu2*phi;
		}
		// All these numbers should be 0.5 -- 1 and 3 are mean square of
		// distributions, 2 is integral over half of profile
		fprintf(ioQQQ,"Check %g %g %g\n",totnu1/NSAMP,totphi/(0.1*NSAMP),
				  totnu2/totphi);
	}

	realnum width = 10.;
	long npt = 1000;
	vector<realnum> v(npt), y(npt), cv(npt+1);
	if (a != 0.0)
	{
		// Construct cumulative pdf for generation of random deviates

		// Other distributions of ordinates will be better...
		for (long i = 0; i<npt; ++i)
		{
			v[i] = width*(i+0.5)/npt;
		}
		VoigtU(a,&v[0],&y[0],npt);
		cv[0] = 0.;
		for (long i = 0; i<npt; ++i)
		{
			cv[i+1] = cv[i] + (width/npt)*y[i];
		}
		for (long i = 0; i<npt; ++i)
		{
			cv[i+1] /= cv[npt];
		}
	}	

	const long NPART = 10000;
	mc_table(true,0.,a,0.,0.,0.,0.,0.);
	for (double tau0=0.01; tau0<tau_in; tau0 *= 1.1)
	{
		double nfirst = 0.0, ntot = 0.0, ltot=0.0, taulast=0.0,
			ndest = 0.0;
		for (long i=0; i<NPART; ++i)
		{
			double tau = tau0;
			double nscat = 0;
			for (;;)
			{
				// Save last scattering point
				double tauprev = tau;
				// CRD frequency, angle
				double nu,phix;
				if (a == 0.0)
				{
					nu = ran.normal()/sqrt(2.);
					phix = phigauss(nu); // Line profile function
				}
				else
				{
					realnum offx = 2.f*ran.rnm()-1.f;
					int sign = 1;
					if (offx < 0.0)
					{
						sign = -1;
						offx = -offx;
					}
					long ii = hunt_bisect(&cv[0],npt,offx)+1;
					ASSERT (ii > 0 && ii < npt);
					realnum frac = (offx-cv[ii])/SDIV(cv[ii+1]-cv[ii]);
					nu = sign*((1.0-frac)*v[ii]+frac*v[ii+1]);
					phix = (1.0-frac)*y[ii]+frac*y[ii+1];
				}
				double costheta = 2.0*ran.dbl()-1.0;

				// Random tau step at nu
				double dtaunu = -log(ran.dbl());
				// Convert tau step from nu to mean optical depth scale
				// by dividing by relative asorption at nu, including 
				// continuous absorption
				double dlength = dtaunu / SDIV(beta+phix); 
				// Convert to mean optical depth step perpendicular to slab
				// cf HK equation 2.5, dtau ~ mu/(beta+phix)
				double dtauperp = dlength * costheta;
				// Calculate length of step within slab
				if (tau-dtauperp < 0.0)
					ltot += tau/costheta;
				else if (tau-dtauperp > 2.0*tau0)
					ltot += (2.0*tau0-tau)/costheta;
				else
					ltot += dlength;
				// Calculate new position, reflecting across midplane
				tau -= dtauperp;
				if (tau > tau0)
					tau = 2.0*tau0-tau;
				if (tau < 0.0)
				{
					// Has escaped
					taulast += tauprev;
					if (nscat == 0)
						++nfirst;
					break;
				}
				if (beta > 0.0 && ran.dbl() < beta/(beta+phix) )
				{
					// beta = k_c/k_L
					// Destruction fraction is k_c/(k_c+k_L*phi)
					++ndest;
					break;
				}
				// Has not escaped
				++nscat;
			} 
			ntot += nscat;
		}
		// Print half-depth of slab, escape probability on midplane
		// emission, mean number of scatterings and mean path length
		// within slab
		mc_table(false,tau0,a,nfirst/NPART,ntot/NPART,ltot/NPART,
					taulast/(NPART-ndest),ndest/NPART);
	}

	fprintf(ioQQQ,"\n\n#Avrett tau a esca0k2 esctot \n");

	// try to compare the formulae from esc_CRDwing_1side with 
	// Avrett & Loeser exact expression
	double taustep = 2., taumin=1e-8,taumax=1e14;
	for (double tau = taumin; tau < taumax; tau *= taustep)
	{
		double esccom_v = esca0k2(tau);
		double esctot = esc_CRDwing_1side(tau,a);
		k2DampArg k2fun(a,SQRTPI*tau);
		double rmax1 = 1.0,rmax;
		for(int idiv=0;idiv<200;++idiv)
		{
			rmax = rmax1;
			rmax1 *= 0.75;
			if (k2fun(rmax1) != 0.0)
			{
				break;
			}
		}
		Integrator<k2DampArg,Gaussian32> IntDamp( k2fun );
		double intgral = IntDamp.sum(0.0,rmax);
		class integrate::Romberg<integrate::Midpoint<k2DampArg> >
			k2r(integrate::Midpoint<k2DampArg>(k2fun,0.0,rmax));
		k2r.update(1e-10);
		fprintf(ioQQQ,"%15.8g %15.8g %15.8g %15.8g %15.8g %15.8g\n",
				  tau,a,esccom_v,esctot,2.0*intgral,2.0*k2r.sum());
	}

	fprintf(ioQQQ,"\n\n#betaFbeta at a \n");
	double betastep = 2., betamin=1e-6,betamax=1e3;
	for (double beta = betamin; beta < betamax; beta *= betastep)
	{
		overlap ov(a,beta);
		class integrate::Romberg<integrate::Midpoint<overlap> >
			ovr(integrate::Midpoint<overlap>(ov,-100.0,100.0));
		ovr.update(1e-10);
		fprintf(ioQQQ,"%15.8g %15.8g %15.8g %15.8g\n",
				  beta,a,beta/(1.0+beta),ovr.sum());
	}
}
