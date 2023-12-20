/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "abund.h"
t_abund abund;

void t_abund::zero()
{
	DEBUG_ENTRY( "t_abund::zero()" );
	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		/* depletion scale factors */
		DepletionScaleFactor[nelem] = 1.;
	}

	lgDepln = false;
	ScaleMetals = 1.;
}

double t_abund::SumDepletedAtoms()
{
	DEBUG_ENTRY( "t_abund::SumDepletedAtoms()" );

	double sum=0.;
	for( int nelem=0; nelem<LIMELM; ++nelem )
	{
		sum += (1.-DepletionScaleFactor[nelem]) * ReferenceAbun[nelem];
	}
	return sum;
}
