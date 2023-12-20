/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "oxy.h"
t_oxy oxy;

void t_oxy::zero()
{
	DEBUG_ENTRY( "t_oxy::zero()" );
	poiii2 = 0.;
	poiii3 = 0.;
	poiexc = 0.;

	d5007r = 0.;
	d6300 = 0.;
}
