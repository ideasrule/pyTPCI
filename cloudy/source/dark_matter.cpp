/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "dark_matter.h"
t_dark_matter dark;

void t_dark_matter::zero()
{
	DEBUG_ENTRY( "t_dark_matter::zero()" );
	lgNFW_Set = false;
	r_200 = 0.;
	r_s = 0.;
}
