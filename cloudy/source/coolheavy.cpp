/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "coolheavy.h"
t_CoolHeavy CoolHeavy;

void t_CoolHeavy::zero()
{
	DEBUG_ENTRY( "t_CoolHeavy::zero()" );
	/* free free heating, cooling, net */
	lgFreeOn = true;
	brems_cool_h = 0.;
	colmet = 0.;
	brems_cool_net = 0.;
	
}
