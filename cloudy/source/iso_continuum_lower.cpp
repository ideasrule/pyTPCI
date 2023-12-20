/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_continuum_lower - limit max prin. quan. no. due to continuum lowering processes	*/
#include "cddefines.h"
#include "phycon.h"
#include "dense.h"
#include "conv.h"
#include "iso.h"
#include "trace.h"

/* iso_continuum_lower - implement continuum lowering as per section 3.1 of
 * >>refer	ionization	 Bautista, M.A., & Kallman, T.R., 2000, ApJ, 544, 581
 * http://adsabs.harvard.edu/abs/2000ApJ...544..581B
 */
void iso_continuum_lower( long ipISO, long nelem )
{
	long nc;
	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/* size of rate matrices will be defined according to the n calculated here	*/

	ASSERT( dense.xNucleiTotal < MAX_DENSITY );
	ASSERT( nelem < LIMELM );
	/* this may change at a future date. */
	ASSERT( ipISO <= 1 );

	long eff_charge = nelem + 1 - ipISO;

	/* Particle packing - the density here should be density of all nuclei in the plasma */
	/* This one is just nuclear charge, which is independent of iso, and always nelem+1. */
	double np = sqrt( 1.8887E8 * (nelem+1.) / cbrt((double)dense.xNucleiTotal) );
	np = floor(np);
	ASSERT( np > 0. );

	/* Debye shielding - the density here is electron density	*/
	/* This one depends on effective charge. */
	double nd = 2.6E7 * eff_charge * eff_charge * powpq( phycon.te/dense.eden, 1, 4 );
	nd = floor(nd);
	ASSERT( nd > 0. );

	/* Stark broadening - this should be the density of singly charged ions, 
	 * both positive and negative.  The sum of protons, electrons, and HeII should be
	 * good enough.	*/
	/* This one depends on effective charge. */
	double ns = 3171. * pow( (double)eff_charge, 0.8 ) * pow( dense.eden + (double)dense.xIonDense[ipHYDROGEN][1]
		+ (double)dense.xIonDense[ipHELIUM][1], -0.1333);
	ns = floor(ns);
	ASSERT( ns > 0. );

	{
		double nc_doub = MIN3(np, nd, ns);
		if( nc_doub > INT16_MAX )
			nc = INT16_MAX;
		else
			nc = (long)nc_doub;
	}

	/* Don't allow continuum to be lowered below n=3. */
	nc = MAX2( nc, 3 );

	if( nc <= sp->n_HighestResolved_max)
	{
		sp->lgLevelsLowered = true;
		sp->lgLevelsEverLowered = true;
		sp->lgMustReeval = true;
		sp->n_HighestResolved_local = nc;
		sp->nCollapsed_local = 0;
		sp->numLevels_local = iso_get_total_num_levels( ipISO, nc, 0 );
	}
	/* Here is the case where the critical n lies among the one or more collapsed levels */
	/* we just get rid of any that are too high. */
	else if( nc <= sp->n_HighestResolved_max + sp->nCollapsed_max )
	{
		sp->lgLevelsLowered = true;
		sp->lgLevelsEverLowered = true;
		sp->lgMustReeval = true;
		sp->n_HighestResolved_local = sp->n_HighestResolved_max;
		sp->nCollapsed_local = nc - sp->n_HighestResolved_local;
		sp->numLevels_local = 
			iso_get_total_num_levels( ipISO, sp->n_HighestResolved_max, sp->nCollapsed_local );
	}
	/* This is usually where control will flow, because in most conditions the continuum will not be lowered.
	* Nothing changes in this case. */
	else
	{
		sp->numLevels_local = sp->numLevels_max;
		sp->nCollapsed_local = sp->nCollapsed_max;
		sp->n_HighestResolved_local = sp->n_HighestResolved_max;
		
		/* if levels were lowered on last pass but are not now, must reeval */
		if( sp->lgLevelsLowered )
		{
			sp->lgMustReeval = true;
		}
		else
		{
			sp->lgMustReeval = false;
		}

		sp->lgLevelsLowered = false;
	}
	
	if( !conv.nTotalIoniz )
		sp->lgMustReeval = true;

	/* None of these can be greater than that which was originally allocated. */
	ASSERT( sp->numLevels_local <= sp->numLevels_max );
	ASSERT( sp->nCollapsed_local <= sp->nCollapsed_max );
	ASSERT( sp->n_HighestResolved_local <= sp->n_HighestResolved_max );

	/* Lyman lines can not be greater than original allocation or critical pqn. */
	iso_ctrl.nLyman[ipISO] = MIN2( nc, iso_ctrl.nLyman_max[ipISO]);

	// zero out cooling and heating terms involving unused levels
	for( long ipHi=sp->numLevels_local; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( long ipLo=0; ipLo < ipHi; ++ipLo )
			CollisionZero( sp->trans(ipHi,ipLo).Coll() );
	}

	if( trace.lgTrace && (trace.lgHBug||trace.lgHeBug)  )
	{
		fprintf( ioQQQ,"     iso_continuum_lower: ipISO %li nelem %li nc %li (np:%.0f,nd:%.0f,ns:%.0f) numLevels %li nCollapsed %li n_HighestResolved %li \n",
			ipISO, 
			nelem,
			nc,
			np, nd, ns, 
			sp->numLevels_local,
			sp->nCollapsed_local,
			sp->n_HighestResolved_local
			);
	}

	return;
}
