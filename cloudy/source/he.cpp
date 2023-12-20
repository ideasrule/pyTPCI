/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "he.h"
t_he he;

void t_he::zero()
{
	DEBUG_ENTRY( "t_he::zero()" );
	/* zero fractions of He0 destruction due to 23S */
	this->nzone = 0;
	frac_he0dest_23S = 0.;
	frac_he0dest_23S_photo = 0.;
	
}
