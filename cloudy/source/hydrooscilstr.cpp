/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HydroOscilStr computes hydrogenic oscillator strengths, used in the function hdexct. */
#include "cddefines.h"
#include "hydrooscilstr.h"
#include "physconst.h"

double HydroOscilStr(double xLower, 
  double Upper)
{
	double fosc_v, 
	  gaunt, 
	  gnt0, 
	  gnt1, 
	  gnt2, 
	  x;

	DEBUG_ENTRY( "HydroOscilStr()" );

	/*	>>refer	H1	As	Johnson L.C., 1972 ApJ 174 227*/
	/* check order, and that none negative */
	ASSERT( xLower < Upper );
	ASSERT( xLower*Upper >0 );

	x = 1.0 - POW2(xLower/Upper);
	if( xLower >= 3 )
	{
		gnt0 = 0.9935 + 0.2328/xLower - 0.1296/xLower/xLower;
		gnt1 = -(0.6282 - 0.5598/xLower + 0.5299/xLower/xLower)/xLower;
		gnt2 = (0.3887 - 1.181/xLower + 1.470/xLower/xLower)/xLower/
		  xLower;
	}
	else if( xLower == 2 )
	{
		gnt0 = 1.0785;
		gnt1 = -.2319;
		gnt2 = 0.02947;
	}
	else
	{
		gnt0 = 1.1330;
		gnt1 = -.4059;
		gnt2 = 0.07014;
	}
	gaunt = gnt0 + gnt1/x + gnt2/x/x;
	fosc_v = 32./3./PI/sqrt(3.)*xLower/POW3(Upper)*gaunt/x/x/x;
	return( fosc_v );
}
