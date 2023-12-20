/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef TWO_PHOTON_H_
#define TWO_PHOTON_H_

class TransitionProxy;

class two_photon
{
public:
	two_photon()
	{
		ipHi = -1;
		ipLo = -1;
		Pop = NULL;
		E2nu = 0.;
		AulTotal = 0.f;
		ipTwoPhoE = -1;
		ipSym2nu.clear();
		As2nu.clear();
		local_emis.clear();
		induc_up = 0.;
		induc_dn = 0.;
		induc_dn_max = 0.;	
	}
	
	void Reset()
	{
		induc_up = 0.;
		induc_dn = 0.;
		induc_dn_max = 0.;	
	}

	long ipHi, ipLo;
	double* Pop;
	double E2nu;
	realnum AulTotal;

	// pointer to the energy representing the two-photon gap,
	long ipTwoPhoE;

	// series of symmetric indices for two photon 
	vector<long> ipSym2nu;
	// two photon transition probabilities per energy bin
	vector<realnum> As2nu;
	// local emission per ( photons cm-3 s-1 bin-1 )
	vector<realnum> local_emis;

	// the induced upward two-photon rate 
	double induc_up;
	// the induced downward two-photon rate 
	double induc_dn;
	// the largest induced downward two photon rate 
	double induc_dn_max;
};

void atmdat_2phot_setSplineCoefs();

 /**
  atmdat_2phot_shapefunction two photon emission function for all atomic and ionic species 
  \param  EbyE2nu 
  \param  ipISO
  \param  nelem 
 */ 
double atmdat_2phot_shapefunction( double EbyE2nu, long ipISO, long nelem );

void CalcTwoPhotonRates( two_photon& tnu, bool lgDoInduced );
void CalcTwoPhotonEmission( two_photon& tnu, bool lgDoInduced );

void PrtTwoPhotonEmissCoef( const two_photon& tnu, const double& densityProduct );

// note the default values for the last two parameters -- the code uses HI 2s-1s shapefunctions by default
void TwoPhotonSetup( vector<two_photon> &tnu_vec, const long &ipHi, const long &ipLo, const double &Aul, const TransitionProxy &tr, const long ipISO, const long nelem );

#endif /* TWO_PHOTON_H_ */
