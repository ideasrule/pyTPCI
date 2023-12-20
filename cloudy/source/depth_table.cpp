/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"

#include "depth_table.h"
#include "thirdparty.h"

double DepthTable::tabval(double r0, 
  double depth) const
{
	double tabden_v, 
	  x;

	DEBUG_ENTRY( "DepthTable::tabval()" );
	/*interpolate on table of points for density with dlaw table command, by K Volk
	 *each line is log radius and H density per cc. */

	if( nzone == 0 ){ /* radius values are not set correctly yet */
	  tabden_v = val[0];
	}
	else
        {
	    /*begin sanity check */
	    if( r0 <= 0. || depth <= 0. )
	      {
		fprintf( ioQQQ, " dense_tabden called with insane depth, radius, =%10.2e%10.2e\n", 
			 depth, r0 );
	      }
	    /*end sanity check */

	    if (lgLinear) {
	      if (lgDepth) x=depth;
	      else x=r0;
	    }
	    else {
	    /* interpolate on radius or depth? */
	      if( lgDepth )
		{
		  /* depth key appeared = we want depth */
		  x = log10(depth);
		}
	      else
		{
		  /* use radius */
		  x = log10(r0);
		}
	    }
	    
	    /* set to impossible value, will crash if not reset */
	    tabden_v = -DBL_MAX;
	    
	    if( x < dist[0] || x >= dist[nvals-1] )
	      {
		fprintf( ioQQQ, " requested radius outside range of dense_tabden\n" );
		if (lgLinear)
		  {
		  fprintf( ioQQQ, " radius was%10.2e min, max=%10.2e%10.2e\n", 
			   x, dist[0], dist[nvals-1] );
		  }
		else
		  {
		    fprintf( ioQQQ, " log10(radius/depth) was %10.2f min, max=%10.2f %10.2f\n", 
			     x, dist[0], dist[nvals-1] ); 
		  }
		cdEXIT(EXIT_FAILURE);
	      }
	else
	  {
	    tabden_v = linint(get_ptr(dist),get_ptr(val),nvals,x);
	  }
	}
	/* got it, now return value, not log of density */
	if (!lgLinear)
	  tabden_v = exp10(tabden_v);
	return( tabden_v );
}
