/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*MeanInc increment mean ionization fractions and temperatures over computed structure */
/*MeanZero zero mean of ionization fractions / temperatures arrays */
/*MeanIon derive mean ionization fractions or mean temperatures for some element */
#include "cddefines.h"
#include "physconst.h"
#include "radius.h"
#include "dense.h"
#include "hyperfine.h"
#include "magnetic.h"
#include "hmi.h"
#include "phycon.h"
#include "mean.h"
#include "mole.h"

t_mean mean;

t_mean::t_mean()
{
	DEBUG_ENTRY( "t_mean()" );

	// dim 3 for average over radius, area, and volume
	mean.xIonMean.reserve(3);
	for( int j=0; j < 3; ++j )
	{
		// compute average for every element
		mean.xIonMean.reserve(j,LIMELM);
		for( int nelem=0; nelem < LIMELM; ++nelem )
		{
			int limit = max(3,nelem+2);
			// compute average for every ionization stage (incl. H2 for hydrogen)
			mean.xIonMean.reserve(j,nelem,limit);
			for( int ion=0; ion < limit; ++ion )
				// dim 2 for storing Sum(quant*norm) and Sum(norm)
				mean.xIonMean.reserve(j,nelem,ion,2);
		}
	}
	mean.xIonMean.alloc();
	mean.xIonEdenMean.alloc( mean.xIonMean.clone() );
	mean.TempIonMean.alloc( mean.xIonMean.clone() );
	mean.TempIonEdenMean.alloc( mean.xIonMean.clone() );
	mean.TempB_HarMean.alloc(3,2);
	mean.TempHarMean.alloc(3,2);
	mean.TempH_21cmSpinMean.alloc(3,2);
	mean.TempMean.alloc(3,2);
	mean.TempEdenMean.alloc(3,2);
}

void t_mean::setup_molecules()
{
	DEBUG_ENTRY( "setup_molecules()" );

	vector<string> allMolecules;
	getMolecules( allMolecules );

	for( vector<string>::iterator it = allMolecules.begin();
		it != allMolecules.end(); it++ )
	{
		mean.molecules[ *it ].alloc(3,2);
	}

}

/*MeanZero zero mean of ionization fractions / temperatures arrays */
void t_mean::zero()
{
	DEBUG_ENTRY( "t_mean::zero()" );

	/* MeanZero is called at start of calculation by zero, and by
	 * startenditer to initialize the variables */

	mean.xIonMean.zero();
	mean.xIonEdenMean.zero();
	mean.TempIonMean.zero();
	mean.TempIonEdenMean.zero();
	mean.TempB_HarMean.zero();
	mean.TempHarMean.zero();
	mean.TempH_21cmSpinMean.zero();
	mean.TempMean.zero();
	mean.TempEdenMean.zero();

	for( auto it = mean.molecules.begin(); it != mean.molecules.end(); it++ )
	{
		it->second.zero();
	}

	return;
}

/*MeanInc increment mean ionization fractions and temperatures over computed structure */
void t_mean::MeanInc()
{
	/* this routine is called by radius_increment during the calculation to update the 
	 * sums that will become the rad, area, and vol weighted mean ionizations */

	DEBUG_ENTRY( "t_mean::MeanInc()" );

	/* Jacobian used in the integrals below */
	double Jac[3] = { radius.drad_x_fillfac, radius.darea_x_fillfac, radius.dVeffVol };

	for( int d=0; d < 3; ++d )
	{
		double dc;
		for( int nelem=0; nelem < LIMELM; nelem++ )
		{
			int limit = max(3,nelem+2);
			/* this normalizes xIonMean and xIonEdenMean,
			 * use gas_phase which includes molecular parts */
			double norm = dense.gas_phase[nelem]*Jac[d];
			for( int ion=0; ion < limit; ion++ )
			{
				if( nelem == ipHYDROGEN && ion == 2 )
					/* hydrogen is special case since must include H2,
					 * and must mult by 2 to conserve total H - will need to div
					 * by two if column density ever used */
					dc = 2.*hmi.H2_total*Jac[d];
				else
					dc = dense.xIonDense[nelem][ion]*Jac[d];
				mean.xIonMean[d][nelem][ion][0] += dc;
				mean.xIonMean[d][nelem][ion][1] += norm;
				mean.TempIonMean[d][nelem][ion][0] += dc*phycon.te;
				mean.TempIonMean[d][nelem][ion][1] += dc;

				dc *= dense.eden;
				mean.xIonEdenMean[d][nelem][ion][0] += dc;
				mean.xIonEdenMean[d][nelem][ion][1] += norm*dense.eden;
				mean.TempIonEdenMean[d][nelem][ion][0] += dc*phycon.te;
				mean.TempIonEdenMean[d][nelem][ion][1] += dc;
			}
		}

		/* =============================================================
		 *
		 * these are some special quantities for the mean 
		 */

		/* used to get magnetic field weighted wrt harmonic mean spin temperature, 
		 * * H0 - as measured by 21cm temperature */
		dc = ( hyperfine.Tspin21cm > SMALLFLOAT ) ? dense.xIonDense[ipHYDROGEN][0]*Jac[d]/phycon.te : 0.;
		/* mean magnetic field weighted wrt 21 cm opacity, n(H0)/T */
		mean.TempB_HarMean[d][0] += sqrt(fabs(magnetic.pressure)*PI8) * dc;
		/* this assumes that field is tangled */
		mean.TempB_HarMean[d][1] += dc;

		/* used to get harmonic mean temperature with respect to H,
		 * for comparison with 21cm temperature */
		dc = dense.xIonDense[ipHYDROGEN][0]*Jac[d];
		/* harmonic mean gas kinetic temp, for comparison with 21 cm obs */
		/*HEADS UP - next are the inverse of equation 3 of
		 * >>refer	Tspin	21 cm	Abel, N.P., Brogan, C.L., Ferland, G.J., O'Dell, C.R., 
		 * >>refercon	Shaw, G. & Troland, T.H. 2004, ApJ, 609, 247 */
		mean.TempHarMean[d][0] += dc;
		mean.TempHarMean[d][1] += dc/phycon.te;

		/* harmonic mean of computed 21 cm spin temperature - this is what 21 cm actually measures */
		mean.TempH_21cmSpinMean[d][0] += dc;
		mean.TempH_21cmSpinMean[d][1] += dc / SDIV( hyperfine.Tspin21cm );

		dc = Jac[d];
		mean.TempMean[d][0] += dc*phycon.te;
		mean.TempMean[d][1] += dc;

		dc = Jac[d]*dense.eden;
		mean.TempEdenMean[d][0] += dc*phycon.te;
		mean.TempEdenMean[d][1] += dc;

		for( auto it = mean.molecules.begin(); it != mean.molecules.end(); it++ )
		{
			dc = species_gasphase_density( it->first ) * Jac[d];
			it->second[d][0] += dc*phycon.te;
			it->second[d][1] += dc;
		}
	}
	return;
}

/*MeanIon derive mean ionization fractions or mean temperatures for some element */
void t_mean::MeanIon(
	/* either 'i' or 't' for ionization or temperature */
	char chType,
	/* atomic number on C scale */
	long int nelem, 
	/* type of average: 0=radius, 1=area, 2=volume */
	long int dim,
	/* this will say how many ion stages in arlog have non-zero values */
	long int *n, 
	/* array of values, log both cases, 
	 * hydrogen is special case since [2] will be H2  */
	realnum arlog[],
	/* true, include electron density, false do not*/
	bool lgDensity ) const
{
	DEBUG_ENTRY( "t_mean::MeanIon()" );

	/* limit is on C scale, such that ion < limit */
	int limit = max( 3, nelem+2 );

	/* fills in array arlog with log of ionization fractions
	 * n is set to number of non-zero abundances
	 * n set to 0 if element turned off
	 *
	 * first call MeanZero to zero out arrays
	 * next call MeanInc in zone calc to enter ionziation fractions or temperature
	 * finally this routine computes actual means
	 * */
	if( !dense.lgElmtOn[nelem] )
	{
		/* no ionization for this element */
		for( int ion=0; ion < limit; ion++ )
			arlog[ion] = -30.f;
		*n = 0;
		return;
	}

	/* set high ion stages, with zero abundances, to -30 */
	*n = limit;
	while( *n > 0 && mean.xIonMean[0][nelem][*n-1][0] <= 0. )
	{
		arlog[*n-1] = -30.f;
		--*n;
	}

	double meanv, normv;
	for( int ion=0; ion < *n; ion++ )
	{
		/* mean ionization */
		if( chType == 'i' )
		{
			if( lgDensity )
			{
				meanv = mean.xIonEdenMean[dim][nelem][ion][0];
				normv = mean.xIonEdenMean[dim][nelem][ion][1];
			}
			else
			{
				meanv = mean.xIonMean[dim][nelem][ion][0];
				normv = mean.xIonMean[dim][nelem][ion][1];
			}
			arlog[ion] = ( meanv > 0. ) ? (realnum)log10(max(1e-30,meanv/normv)) : -30.f;
		}
		/* mean temperature */
		else if( chType == 't' )
		{
			if( lgDensity )
			{
				meanv = mean.TempIonEdenMean[dim][nelem][ion][0];
				normv = mean.TempIonEdenMean[dim][nelem][ion][1];
			}
			else
			{
				meanv = mean.TempIonMean[dim][nelem][ion][0];
				normv = mean.TempIonMean[dim][nelem][ion][1];
			}
			arlog[ion] = ( normv > SMALLFLOAT ) ? (realnum)log10(max(1e-30,meanv/normv)) : -30.f;
		}
		else
		{
			fprintf(ioQQQ," MeanIon called with insane job: %c \n",chType);
			TotalInsanity();
		}
	}
	return;
}

bool t_mean::MeanMoleculeTemp( const string &chSpecies, long dim, double &TeMean )
{
	DEBUG_ENTRY( "MeanMoleculeTemp()" );

	auto mole = mean.molecules.find( chSpecies );

	if( mole == mean.molecules.end() )
	{
		TeMean = 0.;
		return 1;
	}

	double meanv = mole->second[dim][0];
	double normv = mole->second[dim][1];
	TeMean = ( normv > SMALLFLOAT ) ? realnum(max(1e-30,meanv/normv)) : 1.e-30f;
	return 0;
}
