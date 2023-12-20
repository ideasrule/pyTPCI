/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_zero zero out heating save arrays */
#include "cddefines.h"
#include "thermal.h"
#include "ionbal.h"

void ion_zero(long int nelem)
{
	long int i;

	DEBUG_ENTRY( "ion_zero()" );

	/* heating array, but only for direct photoionization */
	for( i=0; i < nelem; i++ )
	{
		thermal.setHeating(nelem,i,0.);
	}
	return;
}
