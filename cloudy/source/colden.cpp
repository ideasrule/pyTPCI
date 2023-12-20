/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "colden.h"
t_colden colden;

void t_colden::zero()
{
	DEBUG_ENTRY( "t_colden::zero()" );
		/* zero out some column densities */
	for( long i=0; i < NCOLD; i++ )
	{
		colden[i] = 0.;
	}
	coldenH2_ov_vel = 0.;

	/* F=0 and F=1 column densities of H0*/
	H0_21cm_upper =0;
	H0_21cm_lower =0;

	/* variables to do with Jeans mass and radius */
	TotMassColl = 0.;
	tmas = 0.;
	wmas = 0.;
	rjnmin = FLT_MAX;
	ajmmin = FLT_MAX;

}
