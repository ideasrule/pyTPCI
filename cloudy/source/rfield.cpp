/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "rfield.h"

#include "opacity.h"
#include "continuum.h"

t_rfield rfield;

void t_rfield::zero()
{
	DEBUG_ENTRY( "t_rfield::zero()" );
	lgHabing = false;

	/* flag to turn off Lya ots */
	lgLyaOTS = true;
	/* HeII rec and Lya ots */
	lgHeIIOTS = true;
	lgKillOTSLine = false;
	lgKillOutLine = false;
	lgKillOutCont = false;

	/* DiffPumpOn is unity unless process disabled by setting to 1
	 * with no diffuse line pumping command */
	DiffPumpOn = 1.;

	/* >>chng 03 nov 28, add option to not do line transfer */
	lgDoLineTrans = true;

	/* flag saying whether to constantly reevaluated opacities -
	 * set false with no opacity reevaluate command */
	lgOpacityReevaluate = true;

	/* flag saying whether to constantly reevaluated ionization -
	 * set false with no ionization reevaluate command */
	lgIonizReevaluate = true;
	/* this element is default for choosing line width */
	fine_opac_nelem = ipIRON;
	/* there will be this many resolution elements in one FWHM for this element,
	 * at the lowest temperature to be considered */
	fine_opac_nresolv = 1;
	/* continuum scale factor for case of time varying continuum */
	time_continuum_scale = 1.;
	/* will fine optical depths be punched? */
	lgSaveOpacityFine = false;

	/* first is set true if one of the incident continua needs to have
	 * H-ionizing radiation blocked.  Second is set true is it is blocked
	 * with extinguish command - want both true if first is true */
	lgMustBlockHIon = false;
	lgBlockHIon = false;
	strncpy( rfield.chCumuType, "MASS", sizeof(rfield.chCumuType));
}

void Spectrum::resize( long int nbins )
{
	ASSERT( nbins > 0 );

	nflux = nbins;
	flux[0].resize( nflux );
	flux[1].resize( nflux );

	vzero( flux[1] );
}

const realnum *t_rfield::getCoarseTransCoef()
{
	// average opacity transmission coefficient fine to coarse
	if( opac.lgScatON && trans_coef_total_stale)
	{
		/* sum over coarse continuum */
		for( long i=0; i < nflux-1; i++ )
		{
			// find transmission coefficient if lower and upper bounds 
			// of coarse continuum is within boundaries of fine continuum 
			// unity is default
			if( ipnt_coarse_2_fine[i] && ipnt_coarse_2_fine[i+1] )
			{
				// first branch is normal case, where fine continuum is finer than
				// coarse continuum.  But, when end temp is very high, fine continuum is
				// very coarse, so may be just one cell, and following will not pass
				if( ipnt_coarse_2_fine[i+1]>ipnt_coarse_2_fine[i] )
				{
					trans_coef_total[i] = 0.;
					for( long j=ipnt_coarse_2_fine[i]; j<ipnt_coarse_2_fine[i+1]; ++j )
						trans_coef_total[i] += sexp(fine_opt_depth[j]);
					trans_coef_total[i] /= (ipnt_coarse_2_fine[i+1]-ipnt_coarse_2_fine[i]);
				}
				else
				{
					// in case where fine is coarser than coarse, 
					// just use first cell
					trans_coef_total[i] = sexp(fine_opt_depth[ipnt_coarse_2_fine[i]]);
				}
			}
		}
		trans_coef_total_stale = false;
	}
	return trans_coef_total.data();
}

void t_rfield::setTrimming()
{
	DEBUG_ENTRY( "t_rfield::setTrimming()" );

	double threshold = 0.;
	if( rfield.FluxFaint > 0. )
	{
		double peak = 0.;
		for( long i=0; i < rfield.nflux; ++ i )
			peak = max( peak, rfield.SummedCon[i]*rfield.anu(i)/rfield.widflx(i) );
		threshold = rfield.FluxFaint * peak;
	}
		
	// search for highest energy cell with positive flux, use 'i < nPositive' in loops
	nPositive = nflux;
	while( nPositive > 0 && SummedCon[nPositive-1] <= threshold )
		--nPositive;
}

/* report the total flux not including attenuated isotropic continua, if requested */
double flux_correct_isotropic( const bool lgSaveIsotr, const int nEmType, const int iflux )
{
	if( iflux < 0 || iflux >= rfield.nflux )
		return  0.;

	double this_flux = rfield.flux.get_flux( nEmType, iflux );
	if( nEmType == 0 &&
	    ( ! continuum.lgPrtIsotropicCont || ! lgSaveIsotr ) )
	{
		this_flux = (double)( rfield.flux_beam_const[ iflux ] +
					rfield.flux_beam_time[ iflux ] );
	}

	return this_flux;
}

double flux_correct_isotropic( const int nEmType, const int iflux )
{
	return flux_correct_isotropic( true, nEmType, iflux );
}
