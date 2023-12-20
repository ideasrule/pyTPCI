/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "ca.h"
t_ca ca;

void t_ca::zero()
{
	Ca2RmLya = 0.;
	popca2ex = 0.;
	Ca3d = 0.;
	Ca4p = 0.;
	dstCala = 0.;
}
