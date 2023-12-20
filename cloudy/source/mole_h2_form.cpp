/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*mole_H2_form find state specific rates grains and H- form H2 */
#include "cddefines.h" 
#include "grainvar.h" 
#include "phycon.h" 
#include "hmi.h" 
#include "dense.h" 
#include "h2.h" 
#include "deuterium.h"
#include "mole.h"

/*mole_H2_form find state specific rates grains and H- form H2 */
void diatomics::mole_H2_form( void )
{
	DEBUG_ENTRY( "mole_H2_form()" );

	/* rate of entry into X from H- and formation on grain surfaces 
	 * will one of several distribution functions derived elsewhere
	 * first zero out formation rates and rates others collide into particular level */
	for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
	{
		for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
		{
			/* this will be the rate formation (s-1) of H2 due to
			 * both formation on grain surfaces and the H minus route,
			 * also H2+ + H => H2 + H+ into one vJ level */
			/* units cm-3 s-1 */
			H2_X_formation[iVib][iRot] = 0.;
			H2_X_Hmin_back[iVib][iRot] = 0.;
		}
	}

	/* loop over all grain types, finding total formation rate into each ro-vib level,
	 * also keeps trace of total formation into H2 ground and star, as defined by Tielens & Hollenbach,
	 * these are used in the H molecular network */
	hmi.H2_forms_grains = 0.;
	hmi.H2star_forms_grains = 0.;

	double sum_check = 0.;

	/* >>chng 02 oct 08, resolved grain types */
	for( size_t nd=0; nd < gv.bin.size(); ++nd )
	{
		int ipH2 = (int)gv.which_H2distr[gv.bin[nd].matType];
		for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
		{
			for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
			{
				/* >>chng 02 nov 14, changed indexing into H2_X_grain_formation_distribution and gv.bin, PvH */
				double one = 
					/* H2_X_grain_formation_distribution is normalized to a summed total of unity */
					H2_X_grain_formation_distribution[ipH2][iVib][iRot] * 
					/* units of following are s-1 */
					gv.bin[nd].rate_h2_form_grains_used;

				sum_check += one;

				/* final units are s-1 */
				/* units cm-3 s-1 */
				/* >>chng 04 may 05, added atomic hydrogen density, units cm-3 s-1 */
				H2_X_formation[iVib][iRot] += (realnum)one*dense.xIonDense[ipHYDROGEN][0];

				/* save rates for formation into "H2" and "H2*" in the chemical
				 * network - it resolves the H2 into two species, as in 
				 * Hollenbach / Tielens work - these rates will be used in the
				 * chemistry solver to get H2 and H2* densities */
				if( states[ ipEnergySort[0][iVib][iRot] ].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
				{
					hmi.H2star_forms_grains += one;
				}
				else
				{   
					/*  unit s-1, h2 means h2g*/ 
					hmi.H2_forms_grains += one;
				}
			}
		}
	}

	ASSERT( fp_equal_tol( sum_check, gv.rate_h2_form_grains_used_total, 1e-5*sum_check + DBL_MIN ) );

	/* convert to dimensionless factors that add to unity */
	/* >>chng 02 oct 17, use proper distribution function */
	hmi.H2star_forms_hminus = 0.;
	hmi.H2_forms_hminus = 0.;
	double frac_lo , frac_hi; 
	long ipT;

	/* which temperature point to use? */
	if( phycon.alogte<=1. )
	{
		ipT = 0;
		frac_lo = 1.;
		frac_hi = 0.;
	}
	else if( phycon.alogte>=4. )
	{
		ipT = nTE_HMINUS-2;
		frac_lo = 0.;
		frac_hi = 1.;
	}
	else
	{
		/* find the temp */
		for( ipT=0; ipT<nTE_HMINUS-1; ++ipT )
		{
			if( H2_logte_hminus[ipT+1]>phycon.alogte )
				break;
		}
		frac_hi = (phycon.alogte-H2_logte_hminus[ipT])/(H2_logte_hminus[ipT+1]-H2_logte_hminus[ipT]);
		frac_lo = 1.-frac_hi;
	}

	/* formation of H2 in excited states from H- H minus */
	/* >>chng 02 oct 17, temp dependent fit to rate, updated reference,
	 * about 40% larger than before */
	/* rate in has units cm-3 s-1 */
	double rate = (mole.findrk("H,H-=>H2,e-")+mole.findrk("H,H-=>H2*,e-")) * dense.xIonDense[ipHYDROGEN][0] * findspecieslocal("H-")->den;
	/* backward rate, s-1 */
	double rateback = (mole.findrk("H2,e-=>H,H-")+mole.findrk("H2*,e-=>H,H-"))*dense.eden;
	double oldrate = 0.;

	/* we now know how to interpolate, now fill in H- formation sites */
	for( long iVib=0; iVib<=nVib_hi[0]; ++iVib )
	{
		for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
		{
			/* the temperature-interpolated distribution function, adds up to unity, 
			 * dimensionless */
			double rate_interp =
				frac_lo*H2_X_hminus_formation_distribution[ipT][iVib][iRot] +
				frac_hi*H2_X_hminus_formation_distribution[ipT+1][iVib][iRot];

			/* above rate was set, had dimensions cm-3 s-1 
			 * rate is product of parent densities and total formation rate */
			double one = rate * rate_interp;

			/* save this rate [s-1] for back reaction in levels solver for v,J */
			H2_X_Hmin_back[iVib][iRot] = (realnum)(rateback * rate_interp);

			/* units cm-3 s-1 */
			H2_X_formation[iVib][iRot] += (realnum)one;

			oldrate += rate_interp;

			/* save rates to pass back into molecule network */
			if( states[ ipEnergySort[0][iVib][iRot] ].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
			{
				hmi.H2star_forms_hminus += one;
			}
			else
			{	
				/*  unit cm-3s-1, h2 means h2g*/
				hmi.H2_forms_hminus += one;
			}
		}
	}
	/* confirm that shape function is normalized correctly */
	ASSERT( fabs(1.-oldrate)<1e-4 );

	/* >>chng 03 feb 10, add this population process */
	/* H2+ + H => H2 + H+,
	 * >>refer	H2	population	Krstic, P.S., preprint 
	 * all goes into v=4 but no J information, assume into J = 0 */
	/* >>chng 04 may 05, add density at end */
	rate = mole.findrk("H,H2+=>H+,H2*") * dense.xIonDense[ipHYDROGEN][0] * findspecieslocal("H2+")->den;
	/* units cm-3 s-1 */
	H2_X_formation[4][0] += (realnum)rate;

	fixit("this code is still far too H2-centric. Kick the can a bit");
	ASSERT( this==&h2 || this==&hd );	
	if( this != &h2 )
	{
		for( long iVib=0; iVib<=nVib_hi[0]; ++iVib )
		{
			for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
			{
				// rescale everything in this routine by D/H 
				// THIS MUST NOT BE KEPT FOR OTHER DIATOMS!!!!
				H2_X_formation[iVib][iRot] *= deut.gas_phase/dense.gas_phase[ipHYDROGEN];
			}
		}
	}

	return;
}
