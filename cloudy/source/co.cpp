/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "co.h"
t_co co;

void t_co::zero()
{
	codfrc = 0.;
	codtot = 0.;
	CODissHeat = 0.;	
}
