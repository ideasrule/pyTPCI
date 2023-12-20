/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atmdat_outer_shell determine outer shell, and statistical weights of that and higher ion, for any ion
 * written by Dima Verner */
#include "cddefines.h"
#include "atmdat.h"

void atmdat_outer_shell(long int iz, /* atomic number from 1 to 30 */
  long int in, /* number of electrons from 1 to iz */
  long int *imax, /* number of the outer shell */
  long int *ig0, /* statistical weight of (iz,in) ground state */
  long int *ig1) /* statistical weight of (iz,in-1) ground state */
{
	long int kg;

	static long iss[30]={1,1,2,2,3,3,3,3,3,3,4,4,5,5,5,5,5,5,6,6,6,
	  6,6,6,6,6,6,6,7,7};

	static long igl[30]={2,1,2,1,2,1,4,5,4,1,2,1,2,1,4,5,4,1,4,5,4,
	  1,6,9,10,9,6,1,2,1};

	static long iga[12]={2,1,4,5,4,7,6,9,10,9,2,1};

	DEBUG_ENTRY( "atmdat_outer_shell()" );
	/*determine outer shell for some species */
	/******************************************************************************
	 *** Input parameters:  iz - atomic number from 1 to 30 (integer) 
	 ***          in - number of electrons from 1 to iz (integer)
	 *** Output parameters: imax - number of the outer shell
	 ***          ig0  - statistical weight of (iz,in) ground state
	 ***          ig1  - statistical weight of (iz,in-1) ground state
	 ****************************************************************************** */

	if( iz < 1 || iz > 30 )
	{
		fprintf( ioQQQ, " ***ERROR: wrong atomic number\n" );
		return;
	}

	if( in < 1 || in > iz )
	{
		fprintf( ioQQQ, " ***ERROR: wrong number of electrons\n" );
		return;
	}

	/*** Number of the outer shell and statistical weight */
	*imax = iss[in-1];
	*ig0 = igl[in-1];

	/* in is 1 or greater - as verified above */
	if( in == 1 )
	{
		*ig1 = 1;
	}

	else if( in > 1 )
	{
		*ig1 = igl[in-2];
	}

	else
	{
		/* this is total insanity, cannot happen*/
		fprintf( ioQQQ, " ***ERROR: in insaniy in atmdat_outer_shell\n" );
		return;
	}

	if( in > 18 && iz == in )
	{
		*imax = 7;
		kg = iz - 18;
		*ig0 = iga[kg-1];
		if( iz == 20 )
			*ig1 = 2;
		if( iz == 21 )
			*ig1 = 3;
		if( iz == 22 )
			*ig1 = 4;
		if( iz == 25 )
			*ig1 = 7;
		if( iz == 26 )
			*ig1 = 10;
		if( iz == 30 )
			*ig1 = 2;
	}

	if( in > 18 && (iz - in) == 1 )
	{
		if( iz == 20 )
		{
			*imax = 7;
			*ig0 = 2;
		}

		if( iz == 21 )
		{
			*imax = 7;
			*ig0 = 3;
		}

		if( iz == 22 )
		{
			*imax = 7;
			*ig0 = 4;
		}

		if( iz == 25 )
		{
			*imax = 7;
			*ig0 = 7;
		}

		if( iz == 26 )
		{
			*imax = 7;
			*ig0 = 10;
		}

		if( iz == 30 )
		{
			*imax = 7;
			*ig0 = 2;
		}

	}

	return;
}
