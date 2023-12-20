/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "pressure.h"
#include "rfield.h"
#include "transition.h"
#include "rt_escprob.h"

t_pressure pressure;

void t_pressure::zero()
{
	DEBUG_ENTRY( "t_pressure::zero()" );
	/* pressure related variables */

	RhoGravity_dark = 0.;
	RhoGravity_self = 0.;
	RhoGravity_external = 0.;
	RhoGravity = 0.;
	IntegRhoGravity = 0.;
	gravity_symmetry = -1;
	self_mass_factor = 1.;

	PresRamCurr = 0.;
	pres_radiation_lines_curr = 0.;
	lgPradCap = false;
	lgPradDen = false;
	lgLineRadPresOn = true;
	/* normally abort when radiation pressure exceeds gas pressure in const pres mod,
	 * this is option to keep going, set with NO ABORT on constant pressure command */
	lgRadPresAbortOK = true;
	/* Ditto for whether to stop at sonic point, this gets set to false
	 * for some of the dynamics pressure modes (strongd, shock, antishock)*/
	lgSonicPointAbortOK = true;
	/* this flag will say we hit the sonic point */
	lgSonicPoint = false;
	/* True when no physical solution for desired pressure in strong D fronts */
	lgStrongDLimbo = false; 

	RadBetaMax = 0.;
	ipPradMax_nzone = -1;
	PresMax = 0.;

	/* initial and current pressures */
	PresTotlInit = 0.;
	PresTotlCurr = 0.;

}

double PressureRadiationLine( const TransitionProxy &t, realnum DopplerWidth )
{
	DEBUG_ENTRY( "PressureRadiationLine()" );

	/* return zero if below plasma frequency */	
	if( t.EnergyErg() / EN1RYD <= rfield.plsfrq )
		return 0.;

	/* RT_LineWidth gets line width in terms of RT effects */
	double width = RT_LineWidth(t, DopplerWidth);

	double PopOpc = t.Emis().PopOpc()/(*t.Lo()).g();
	/* return zero line radiation PressureReturned if line mases or has
	 * zero opacity */
	/* \todo 2 1e-22 is arbtrary but roughly 1/kpc.  Replace with a cloud width if available? */
	if( t.Emis().VoigtLineCen() * PopOpc*t.Emis().opacity()/ DopplerWidth <= 1.e-22 || width<=0. )
		return 0.;

	double PressureReturned = PI8 * HPLANCK / 3. * POW4(t.EnergyWN()) *
		((*t.Hi()).Pop()/(*t.Hi()).g())/PopOpc * width;

	/* this prevents line radiation PressureReturned from being very large when line 
	 * is not optically thick but total opacity at that energy is large 
	 * due to overlapping transitions */
	long int ipLineCenter = t.Emis().ipFine() + rfield.ipFineConVelShift;
	if( ipLineCenter > 0 && ipLineCenter < rfield.nfine && rfield.lgOpacityFine &&
		rfield.fine_opac_zone[ipLineCenter] > SMALLFLOAT )
	{
		double FractionThisLine = t.Emis().VoigtLineCen() * t.Emis().PopOpc() * t.Emis().opacity() /
			DopplerWidth / rfield.fine_opac_zone[ipLineCenter];
		if( FractionThisLine<1e-5 )
			FractionThisLine = 0.;
		/* fine opacities are only reevaluated one time per zone due
		 * to the expense - PopOpc is for the current solution - but the two
		 * may be out of step by a few percent, due to the variation in
		 * abundance from zone to zone.  This prevents the change
		 * in solution from increasing the radiation pressure.
		 * This correction is mainly an order of magnitude scaler to prevent
		 * optically thin lines from appearing to be optically thick due to
		 * overlapping lines */
		FractionThisLine = MIN2(1., FractionThisLine);
		ASSERT( FractionThisLine >= 0. && FractionThisLine <= 1.0 );
		PressureReturned *= FractionThisLine;
	}

	return PressureReturned;
}
