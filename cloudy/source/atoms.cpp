/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "atoms.h"
t_atoms atoms;

void t_atoms::zero()
{
	DEBUG_ENTRY( "t_atoms::zero()" );
	popMg2 = 0.;
	rateMg2 = 0.;
}
