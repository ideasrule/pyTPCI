/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "radius.h"
t_radius radius;

void t_radius::zero()
{
	DEBUG_ENTRY( "t_radius::zero()" );

	rinner = 0.;
	distance = 0.;
	Radius = 0.;
	Radius_mid_zone = 0.;
	depth = DEPTH_OFFSET;
	depth_mid_zone = DEPTH_OFFSET/2.;
	depth_x_fillfac = 0.;
	lgRadiusKnown = false;
	drad = 0.;
	drad_mid_zone = 0.;
	r1r0sq = 1.;
	PI4_Radius_sq = 0.;
	PI4_rinner_sq = 0.;
	/* this is changed with the roberto command, to go from out to in */
	dRadSign = 1.;

	/* RDFALT is default starting radius (cm) */
	/* >>chng 03 nov 12, from 25 to 30 for Lya clouds */
	rdfalt = 1e30;

	/* set default cylinder thickness */
	CylindHigh = 1e35;
	lgCylnOn = false;

	drad_x_fillfac = 1.;
	darea_x_fillfac = 1.;
	dVeffVol = 1.;
	dVeffAper = 1.;
	drNext = 1.;
	dRNeff = 1.;
	lgdR2Small = false;

	sdrmin = SMALLFLOAT;
	lgSdrminRel = false;
	sdrmax = 1e30;
	lgSdrmaxRel = false;
	lgSMinON = false;
	lgDrMnOn = true;
	lgFixed = false;
	sdrmin_rel_depth = 1e-5;

	lgDrMinUsed = false;
}
