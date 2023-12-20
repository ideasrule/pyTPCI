/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "h2_priv.h"
#include "hmi.h"
#include "rfield.h"
#include "ipoint.h"
#include "mole.h"
	
STATIC void GetDissociationRateCoeff( diss_tran& tran );

void diatomics::Read_Mol_Diss_cross_sections()
{
	DEBUG_ENTRY( "diatomics::Read_Mol_Diss_cross_sections()" );

	/* these can be generalized later for other molecules */
	const int ipNUM_FILES = 1;
	string chFileNames[ipNUM_FILES] = { "cont_diss.dat" };

	// these cross sections are not meaningful without big H2 enabled. 
	ASSERT( lgEnabled );

	// for the summed-over-final-state rates, allocate space for all states in model. 
	// Will be zero-filled and overwritten with whatever data is available.
	Cont_Dissoc_Rate.reserve( n_elec_states );
	for( int iElecLo=0; iElecLo<n_elec_states; ++iElecLo )
	{
		Cont_Dissoc_Rate.reserve( iElecLo, nVib_hi[iElecLo]+1 );
		for( int iVibLo=0; iVibLo<=nVib_hi[iElecLo]; ++iVibLo )
		{
			Cont_Dissoc_Rate.reserve( iElecLo, iVibLo, nRot_hi[iElecLo][iVibLo]+1 );
		}
	}
	Cont_Dissoc_Rate.alloc();

	/* now open the data file */
	string chPath = "h2" + cpu.i().chDirSeparator() + chFileNames[0];
	FILE* ioDATA = open_data( chPath, "r" );

	Diss_Trans.clear();

	/* track number of datasets in file */
	string chLine;
	while( read_whole_line( chLine, ioDATA ) )
	{
		static bool skipData = false;
		long ini=0, inf, iv, ij;

		/* ignore all lines that begin with a '#', unless they are telling us the next quantum numbers.  If so,
		 * then read the quantum numbers and start a new dataset */
		if( chLine.length() >= 6 && chLine[0] == '#' && chLine[1] == '!' &&
		    chLine[4] == 'n'  && chLine[5] =='e' )
		{
			sscanf(chLine.c_str(), "#!  nei=%li, nef=%li, vi= %li, ji= %li", &ini, &inf, &iv, &ij);
			// skip data that is for levels beyond our model
			if( ini > n_elec_states )
				skipData = true;
			else if( iv > nVib_hi[ini] )
				skipData = true;
			else if( ij > nRot_hi[ini][iv] )
				skipData = true;
			else
			{
				skipData = false;
				diss_level Leveli = {ini, iv, ij};
				diss_level Levelf = {inf, -1, -1};
				diss_tran thisTran( Leveli, Levelf );
				Diss_Trans.push_back( thisTran );
			}
		}
		
		if( skipData )
			continue;
	
		/* if line does not begin with a '#', then they are numbers in the datset.  Read in energy in 
		 * wavenumbers and cross section in Angstroms ^2 */
		if( chLine[0] != '#' )
		{
			double energy, xsection;
			const double AngSquared = 1e-16;
			sscanf(chLine.c_str(), "%lf,%lf", &energy, &xsection);

			/* convert energy in wavenumbers to Ryd */
			Diss_Trans.back().energies.push_back( energy*WAVNRYD );

			/* convert cross section from Angstrom^2 to cm^2 by multiplying by 1e-16 */
			Diss_Trans.back().xsections.push_back( xsection*AngSquared );
		}
	}

	fclose(ioDATA);
	return;
}

double diatomics::MolDissocOpacity( const diss_tran& tran, const double& Mol_Ene )
{
	DEBUG_ENTRY( "diatomics::MolDissocOpacity()" );

	double cross_section = MolDissocCrossSection( tran, Mol_Ene ) * 
		states[ ipEnergySort[ tran.initial.n ][ tran.initial.v ][ tran.initial.j ] ].Pop();
	return cross_section;
}

double MolDissocCrossSection( const diss_tran& tran,  const double& Mol_Ene )
{
	DEBUG_ENTRY( "MolDissocCrossSection()" );

	double photodiss_cs = 0.;

	/* if energy is less than threshold, then set cross section is zero */
	if (Mol_Ene < tran.energies[0])
		photodiss_cs = 0.;
	/* if energy is greater than the highest energy defined in Philip's data, set 
	 * cross section to zero for now.  This may need to be changed later */
	else if(Mol_Ene > tran.energies.back())
		photodiss_cs = tran.xsections.back()/sqrt(powi(Mol_Ene/tran.energies.back(),7));
	/* If energy is in between high and low energy limits, do linear interpolation
	 * on Philip's data to compute cross section */
	else
	{
		ASSERT( Mol_Ene > tran.energies[0] && Mol_Ene < tran.energies.back() );
		photodiss_cs = linint(&tran.energies[0],
			&tran.xsections[0],
			tran.xsections.size(),
			Mol_Ene);
	}

	return photodiss_cs;
}

void diatomics::Mol_Photo_Diss_Rates( void )
{
	DEBUG_ENTRY( "diatomics::Mol_Photo_Diss_Rates()" );

	/* Start Calculation of H2 continuum dissociation rate here */

	// these cross sections are not meaningful without big H2 enabled. 
	ASSERT( lgEnabled && mole_global.lgStancil );

	/* Zero out dissociation rates */
	Cont_Dissoc_Rate.zero();
	Cont_Dissoc_Rate_H2g = 0.;
	Cont_Dissoc_Rate_H2s = 0.;
	
	for_each( Diss_Trans.begin(), Diss_Trans.end(), GetDissociationRateCoeff );
	for( vector< diss_tran >::const_iterator dt = Diss_Trans.begin(); dt != Diss_Trans.end(); ++dt )
	{
		double rate = (*this).GetDissociationRate( *dt );
		// this is summed over final electron state, tran.rate_coeff is not (because it's used to conserve energy)	
		Cont_Dissoc_Rate[dt->initial.n][dt->initial.v][dt->initial.j] += dt->rate_coeff;
		if( states[ ipEnergySort[dt->initial.n][dt->initial.v][dt->initial.j] ].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
		{
			/* Add up total rate due to levels which are above arbitrarily defined 
			 * H2* threshold.  Call this dissociation rate of H2* */
			Cont_Dissoc_Rate_H2s += rate;
		}
		else
		{
			/* Add up total rate due to levels which are below arbitrarily defined 
			 * H2* threshold.  Call this dissociation rate of H2g (H2 ground) */
			Cont_Dissoc_Rate_H2g += rate;
		}
	}

	/* This has current units of cm-3 s-1.  We want just the total rate coefficient for
	 * insertion into h_step.  Therefore, divide rate by density of H2* and H2g to get rate 
	 * coefficient for reaction */
	Cont_Dissoc_Rate_H2g = Cont_Dissoc_Rate_H2g / SDIV(H2_den_g);
	Cont_Dissoc_Rate_H2s = Cont_Dissoc_Rate_H2s / SDIV(H2_den_s);

	return;
}

STATIC void GetDissociationRateCoeff( diss_tran& tran )
{
	DEBUG_ENTRY( "GetDissociationRateCoeff()" );
	
	long index_min = ipoint( tran.energies[0] );
	long index_max = ipoint( tran.energies.back() ); 
	index_max = MIN2( index_max, rfield.nflux-1 );

	tran.rate_coeff = 0.;
	for(long i = index_min; i <= index_max; ++i)
	{
		/* Find cross section at given n, v, j, and photon energy */
		double photodiss_cs = MolDissocCrossSection( tran, rfield.anu(i) );
				
		/* Compute integral of cross section and photon energy -- rate coefficient */
		tran.rate_coeff += photodiss_cs*( rfield.flux[0][i] + rfield.ConInterOut[i]+
			rfield.outlin[0][i]+ rfield.outlin_noplot[i]  );
	}
}

double diatomics::GetDissociationRate( const diss_tran& tran )
{
	DEBUG_ENTRY( "diatomics::GetDissociationRate()" );

	long iElecLo = tran.initial.n;
	long iVibLo = tran.initial.v;
	long iRotLo = tran.initial.j;
					
	/* Compute rate, product of rate coefficient times density */
	return states[ ipEnergySort[iElecLo][iVibLo][iRotLo] ].Pop() * tran.rate_coeff;
}

double diatomics::Cont_Diss_Heat_Rate( void )
{
	DEBUG_ENTRY( "diatomics::Cont_Diss_Heat_Rate()" );

	/* 29 Jul 10 - NPA - this tells code to not compute detailed H2 dissociation heating rate
	 * if Stancil rate is not used.  However, we also do not want to compute this rate if Stancil 
	 * rate is on, but small H2 molecule is used.  Insert logic hear to not use rate if 
	 * small H2 atom is on */

	if( !mole_global.lgStancil || !lgEnabled )
		return 0.;
	
	Mol_Photo_Diss_Rates();

	double Cont_Dissoc_Heat_Rate = 0.0;
	for( vector< diss_tran >::const_iterator dt = Diss_Trans.begin(); dt != Diss_Trans.end(); ++dt )
		Cont_Dissoc_Heat_Rate += (*this).GetHeatRate( *dt );

	return Cont_Dissoc_Heat_Rate;
}

double diatomics::GetHeatRate( const diss_tran& tran )
{
	DEBUG_ENTRY( "diatomics::GetHeatRate()" );
					
	long index_min = ipoint( tran.energies[0] );
	long index_max = ipoint( tran.energies.back() ); 
	index_max = MIN2( index_max, rfield.nflux-1 );
 
	long iElecLo = tran.initial.n;
	long iVibLo = tran.initial.v;
	long iRotLo = tran.initial.j;
	double rate = 0.;

	for( long i = index_min; i<= index_max; ++i )
	{
		/* Photon energy minus threshold (Cloudy convention) */
		double energy = MAX2( 0., rfield.anu(i) - tran.energies[0] );

		/* The density of the current H2(v,J) level */
		double density = states[ ipEnergySort[iElecLo][iVibLo][iRotLo] ].Pop();

		double photodiss_cs = MolDissocCrossSection( tran, rfield.anu(i) );
		double Rate_Coeff = photodiss_cs*( rfield.flux[0][i] + rfield.ConInterOut[i]+
			rfield.outlin[0][i]+ rfield.outlin_noplot[i]  );

		/* The heating rate.  This is the product of the H2(v,J) density, the dissociation rate,
		 * and the energy of the photon (minus threshold).  Units are erg cm-3 s-1.  
		 * Sum over all possible H2 levels and dissociation energies. */
		rate += EN1RYD * energy * Rate_Coeff * density;
	}

	return rate;
}

