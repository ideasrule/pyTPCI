/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "ipoint.h"
#include "rfield.h"
#include "transition.h"
#include "two_photon.h"

void TwoPhotonSetup( vector<two_photon> &tnu_vec, const long &ipHi, const long &ipLo, const double &Aul, const TransitionProxy &tr, const long ipISO, const long nelem )
{
	DEBUG_ENTRY( "TwoPhotonSetup()" );

	tnu_vec.resize( tnu_vec.size() + 1 );
	two_photon &tnu = tnu_vec.back();

	tnu.ipHi = ipHi;
	tnu.ipLo = ipLo;
	tnu.AulTotal = Aul;
	tnu.Pop = &(*tr.Hi()).Pop();
	tnu.E2nu = tr.EnergyRyd();
	tnu.induc_dn_max = 0.;
	
	/* pointer to the Lya energy */
	tnu.ipTwoPhoE = ipoint(tnu.E2nu);
	while( rfield.anu(tnu.ipTwoPhoE) > tnu.E2nu )
	{
		--tnu.ipTwoPhoE;
	}
	tnu.ipSym2nu.resize( tnu.ipTwoPhoE );
	tnu.As2nu.resize( tnu.ipTwoPhoE );
	tnu.local_emis.resize( tnu.ipTwoPhoE );

	/* >>chng 02 aug 14, change upper limit to full Lya energy */
	for( long i=0; i < tnu.ipTwoPhoE; i++ )
	{
		/* energy is symmetric energy, the other side of half E2nu */
		double energy = tnu.E2nu - rfield.anu(i);
		/* this is needed since mirror image of cell next to two-nu energy
		* may be slightly negative */
		energy = MAX2( energy, rfield.anumax(0) );
		/* find index for this symmetric energy */
		tnu.ipSym2nu[i] = ipoint(energy);
		while( rfield.anu(tnu.ipSym2nu[i]) > tnu.E2nu ||
			tnu.ipSym2nu[i] >= tnu.ipTwoPhoE)
		{
			--tnu.ipSym2nu[i];
		}
		ASSERT( tnu.ipSym2nu[i] >= 0 );
	}

	double SumShapeFunction = 0., Renorm= 0.;

	/* ipTwoPhoE is the cell holding the 2nu energy itself, and we do not want
	 * to include that in the following sum */
	ASSERT( rfield.anu(tnu.ipTwoPhoE-1)<=tnu.E2nu );
	for( long i=0; i < tnu.ipTwoPhoE; i++ )
	{
		double ShapeFunction;

		ShapeFunction = atmdat_2phot_shapefunction( rfield.anu(i)/tnu.E2nu, ipISO, nelem )*rfield.widflx(i)/tnu.E2nu;
		SumShapeFunction += ShapeFunction;

		/* >>refer	HI	2nu	Spitzer, L., & Greenstein, J., 1951, ApJ, 114, 407 */
		/* As2nu will add up to the A, so its units are s-1	*/ 
		tnu.As2nu[i] = (realnum)( tnu.AulTotal * ShapeFunction );
	}

	/* The spline function in twophoton.c causes a bit of an error in the integral of the
	 * shape function.  So we renormalize the integral to 1.	*/
	Renorm = 1./SumShapeFunction;

	for( long i=0; i < tnu.ipTwoPhoE; i++ )
	{
		tnu.As2nu[i] *= (realnum)Renorm;
	}

	/* The result should be VERY close to 1.	*/
	ASSERT( fabs( SumShapeFunction*Renorm - 1. ) < 0.00001 );

	return;
}

void CalcTwoPhotonRates( two_photon& tnu, bool lgDoInduced )
{
	DEBUG_ENTRY( "CalcTwoPhotonRates()" );
			
	/* this could fail when pops very low and Pop2Ion is zero */
	ASSERT( tnu.ipTwoPhoE>0 && tnu.E2nu>0. );

	double sum = 0.;
	tnu.induc_up = 0.;
	tnu.induc_dn = 0.;
	/* two photon emission, ipTwoPhoE is 
	* continuum array index for Lya energy */
	ASSERT( rfield.anu(tnu.ipTwoPhoE-1) < 1.01*tnu.E2nu || rfield.anu(tnu.ipTwoPhoE-2)<tnu.E2nu );
	for( long nu=0; nu < tnu.ipTwoPhoE; nu++ )
	{
		// We do not assert rfield.anu(nu)<=tnu.E2nu because the maximum 
		// index could be set to point to the next highest bin.

		// As2nu[nu] is transition probability A per bin
		// So sum is the total transition probability
		sum += tnu.As2nu[nu];
					
		// only include this if induced processes turned on,
		// otherwise inconsistent with rate solver treatment.
		if( lgDoInduced )
		{
			double rate_up = tnu.As2nu[nu] *
				rfield.SummedOcc[nu] * rfield.SummedOcc[tnu.ipSym2nu[nu]-1];
			tnu.induc_up += rate_up;
			tnu.induc_dn += rate_up + tnu.As2nu[nu] *
				(rfield.SummedOcc[nu] + rfield.SummedOcc[tnu.ipSym2nu[nu]-1]);
		}
	}

	/* a sanity check on the code, see Spitzer and Greenstein (1951), eqn 4.	*/
	/* >>refer	HI	2nu	Spitzer, L., & Greenstein, J., 1951, ApJ, 114, 407 */
	ASSERT( fabs( 1.f - (realnum)sum/tnu.AulTotal ) < 0.01f );

	return;
}

void CalcTwoPhotonEmission( two_photon& tnu, bool lgDoInduced )
{
	DEBUG_ENTRY( "CalcTwoPhotonEmission()" );
			
	/* this could fail when pops very low and Pop2Ion is zero */
	ASSERT( tnu.ipTwoPhoE>0 );

	/* two photon emission, ipTwoPhoE is 
	 * continuum array index for Lya energy */
	for( long nu=0; nu < tnu.ipTwoPhoE; nu++ )
	{
		// Pop has dimension cm^-3.  The factor of 2 is for two photons per 
		// transition. Thus two_photon_emiss has dimension photons cm-3 s-1.
		tnu.local_emis[nu] = 2.f * (realnum)(*tnu.Pop) * tnu.As2nu[nu];
	}

	// only include this if induced processes turned on,
	// otherwise inconsistent with rate solver treatment.
	if( lgDoInduced )
	{
		for( long nu=0; nu < tnu.ipTwoPhoE; nu++ )
		{
			// this is the total rate (in this energy bin)
			tnu.local_emis[nu] *= (1.f + rfield.SummedOcc[nu]) *
				(1.f+rfield.SummedOcc[tnu.ipSym2nu[nu]-1]);
		}
	}

	return;
}

/* option to print hydrogen and helium two-photon emission coefficients?	*/
void PrtTwoPhotonEmissCoef( const two_photon& tnu, const double& densityProduct )
{
	DEBUG_ENTRY( "PrtTwoPhotonEmissCoef()" );
				
	fprintf( ioQQQ, "\ny\tGammaNot(2q)\n");

	for( long yTimes20=1; yTimes20<=10; yTimes20++ )
	{
		double y = yTimes20/20.;

		fprintf( ioQQQ, "%.3e\t", (double)y );	

		long i = ipoint(y*tnu.E2nu);
		fprintf( ioQQQ, "%.3e\n", 
			8./3.*HPLANCK*(*tnu.Pop)/densityProduct*y*tnu.As2nu[i]*tnu.E2nu/rfield.widflx(i) );
	}

	return;
}

