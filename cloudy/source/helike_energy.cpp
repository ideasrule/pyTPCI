/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "helike.h"
#include "iso.h"
#include "physconst.h"

/*
  Energy order within 2 3P

  The order of the levels within the 2 3P level of atomic helium is different
  from the order in other ions in the He-like iso sequence.  For simplicity the
  code assumes that the levels are always ordered in the sequence 2 3P0, 2 3P1, 2 3P2.
  So the J levels within 2 3P are not in increasing energy order for helium itself. 
*/

/* helike_quantum_defect - calculate quantum defect for a given level and nuclear charge. */
double helike_quantum_defect( long nelem, long n, long lqn, long S, long j )
{
	DEBUG_ENTRY( "helike_quantum_defect()" );

	/* The quantum defect, and parameters a,b, and c	*/
	double qd,a,b,c;	

	/* These are values of quantum defects of Helium levels at n=10. 
	 * First dimension is spin, second is angular momentum.  
	 * The defects are assumed to be constant for all n>10 and
	 * equal to these values. */
	static const double HeDefectAsymptotes[2][10] = {
		{1.40005E-01,-1.20673E-02,2.08056E-03,4.21484E-04,1.14868E-04,
			4.08648E-05,1.73548E-05,8.33891E-06,4.39680E-06,2.42075E-06},
		{2.97063E-01,6.81567E-02,2.82381E-03,4.27703E-04,1.17319E-04,
			4.25254E-05,1.85549E-05,9.24641E-06,5.30882E-06,3.02877E-06}
	};

	/* Parameters for fits to quantum defects for	*/
	/* P triplet and S orbitals. The dimensions are	*/
	/*	first: l									*/
	/*	second: n									*/
	/* 	third: parameters a,b,and c.				*/
	static const double param[3][4][3]=  
	{ 
		{{0.6451941,0.3119437,-1.2722842},		/* ^3S	*/
			{0.7664874,0.3455675,-1.3976462}, 
			{0.8247101,0.3603131,-1.4520500}, 
			{0.8878402,0.3714450,-1.4995732}}, 

		{{1.4203514,0.5311096,-2.6728087},		/* ^1S	*/
			{1.5733513,0.5997339,-2.9253834}, 
			{1.4531025,0.5924751,-2.8662756}, 
			{1.6038999,0.6342552,-3.0298071}}, 

		{{-2.2323488,0.0890840,-0.5166053},		/* ^3P	*/
			{-2.0463691,0.1222081,-0.6672983}, 
			{-1.9904104,0.1328918,-0.7150879}, 
			{-1.9500974,0.1452111,-0.7649031}} 
	};   

	/* Because they cannot be fit to a funtion of the same form as the other orbitals,
	 * the P singlets are fit to a different function, with these parameters and dimensions */
	/*	first: n									*/
	/*	second: parameters a and b.					*/
	static const double P1[4][2]=
	{
		{-56.65245,-3.661923},
		{-52.03411,-4.941075},
		{-50.43744,-5.525750},
		{-49.45137,-5.908615}
	};

	long int s;

	if( S==1 )
		s = 0;
	else if( S==3 )
		s = 1;
	else if( S < 0 )
	{
		ASSERT( n > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max );
		s = S;
	}
	else
		TotalInsanity();

	ASSERT(n >= 1L);
	ASSERT(lqn >= 0 || n > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max);
	ASSERT(n > lqn);
	/* Only Helium and up, and only those turned on.	*/
	ASSERT((nelem >= ipHELIUM) && (nelem < LIMELM));       

	if( n > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
	{
		/* collapsed levels are assumed to have zero quantum defect. */
		qd = 0.;
	}
	else if( nelem == ipHELIUM )
	{
		if( n<=10 && n<=iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
		{
			double Econt = iso_sp[ipHE_LIKE][nelem].energy_ioniz(n, lqn, S, 2*j+1);
			ASSERT( Econt >= 0. );
			qd = n - sqrt(HE_RYD_FACTOR*RYD_INF/Econt);
		} 
		else if( lqn<=9 )
		{
			ASSERT( s >= 0 && lqn >= 0 );
			/* defects are set equal to the values at n=10. */
			qd = HeDefectAsymptotes[s][lqn];
		}
		else if( s == 0 )
		{
			/* a simple fit for singlet high-l defects. */
			qd = 0.0497*pow((double)lqn, -4.4303);
		}
		else
		{
			/* a simple fit for triplet high-l defects. */
			qd = 0.0656*pow((double)lqn, -4.5606);
		}
	}
	else if( n==1 )
	{
		/* Quantum defects for ground state are found from the rydberg
		 * equation, and the ionization potential of the ion. */
		double EionRyd = iso_sp[ipHE_LIKE][nelem].IonPot/RYD_INF;
		qd = 1.0 - nelem * sqrt(1./EionRyd);
	}
	else
	{
		/* For levels with n > 5, the quantum defect	*/
		/* is approximately the same as if n equaled 5.	*/
		if( n > 5L )
		{
			n = 5L;
		}
		/* For P singlets	*/
		if( lqn==1L && s==0L )
		{
			qd = 1./(P1[n-2][0] + P1[n-2][1] * (nelem+1) * log((double)nelem+1.) );
		}
		/* Defects for orbitals with l>2 are approximately equal to zero.	*/
		else if( lqn < 2L )
		{
			a = param[2*(lqn+1)-s-1][n-2][0];  
			b = param[2*(lqn+1)-s-1][n-2][1];  
			c = param[2*(lqn+1)-s-1][n-2][2];  
			qd = exp((a+c*(nelem+1))/(1.0+b*(nelem+1)));  
		}
		/* This fit is a simplification of table 11.9 from 
		 * >>refer Helike	defects	Drake, G.W.F., editor.  Atomic, Molecular & Optical Physics Handbook.
		 * >>refercon	Chapter 11, "High Precision Calculations for Helium", G.W.F. Drake. 
		 * >>refercon	AIP Press: Woodbury, New York, 1996
		 * This will give quasi-real energies for all transitions, allowing a reasonable
		 * determination of which decays are zeroed due to being below the plasma frequency.
		 * The 1/nelem dependence is arbitray.  	*/
		else
		{
			ASSERT( lqn >= 2L );
			qd = ( ( 0.0612/(double)nelem ) / pow((double)lqn, 4.44) );
		}
	}
	return qd;
}

/* helike_energy calculates energy of a given level in cm^-1. */
double helike_energy(long nelem, long n, long l, long s, long j)
{
	DEBUG_ENTRY( "helike_energy()" );

	double Ef;
	if( n > iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max )
	{
		Ef = HE_RYD_FACTOR * RYD_INF * POW2((double)nelem/(double)n);  
	}
	else
	{
		Ef = iso_sp[ipHE_LIKE][nelem].energy_ioniz(n, l, s, 2*j+1);

		enum{ DEBUG_LOC = false };
		if( DEBUG_LOC )
		{
			fprintf( ioQQQ, "debug E_Helium n l s QN2ind %.2e %ld %ld %ld %ld\n",
					Ef, n, l, s, QN2ind(n, l, s, 2*j+1) );
			fprintf( ioQQQ, "debug E_Helium level_ion Ionpot %.2e %.2e \n",
					iso_sp[ipHE_LIKE][nelem].energy_ioniz(n, l, s, 2*j+1),
					iso_sp[ipHE_LIKE][nelem].IonPot);
		}

		if( Ef < 0. )
		{
			double Eff_n = n - helike_quantum_defect( nelem, n, l, s, j );
			/* quantum defect can only be negative for singlet P */
			ASSERT( ( l==1 && s==1 ) || ( n - Eff_n >= 0. ) );

			/* energies (in wavenumbers) that correspond to quantum defect */
			Ef = HE_RYD_FACTOR * RYD_INF * POW2((double)nelem/Eff_n);
		}
	} 

	ASSERT(Ef > 0.);

	return Ef;
}
