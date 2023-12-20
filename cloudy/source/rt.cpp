/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "rt.h"
t_rt rt;

void t_rt::zero()
{
	DEBUG_ENTRY( "t_rt::zero()" );
	dTauMase = 0.;
	lgMaserCapHit = false;
	lgMaserSetDR = false;

	DoubleTau = 1.;
	lgFstOn = true;
	lgElecScatEscape = true;
	lgElecScatInten = true;
	lgAbsorLineEscape = true;
}
