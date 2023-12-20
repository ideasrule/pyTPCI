/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtMeanIon print mean ionization fractions or temperatures for all elements */
#include "cddefines.h"
#include "prt.h"
#include "geometry.h"
#include "dense.h"
#include "mean.h"
#include "elementnames.h"

void PrtMeanIon(
	/* this is either 't' or 'i', for mean ionization or temperature */
	char chType,
	/* true include density, false do not */
	bool lgDensity,
	/* this is stream used for io, is stdout when called by final,
	 * is save unit when save output generated */
	FILE *ioMEAN )
{
	long int i, 
	  limit, 
	  n, 
	  nelem;
	/* >>chng 04 dec 31, had been static, not possible to ever set false if 
	 * ever true, rm static */
	bool lgPrtLots=false;

	realnum aa[LIMELM+1];

	DEBUG_ENTRY( "PrtMeanIon()" );

	/* print mean ionization or temperature for the computed model
	 * the ionization means are relative to the total element */

	const char* type[3] = { "radius", "area", "volume" };

	for( int d=2; d >= 0; --d )
	{
		/* only print spherical if not plane parallel */
		/* >>>chng 99 may 01, test included one for sphere being set,
		 * so no printout when not set.  now print even when sphere not
		 * set but geo is not plane parallel */
		if( geometry.lgGeoPP && d > 0 )
			continue;

		/* get means for hydrogen */
		mean.MeanIon(chType,ipHYDROGEN,d,&n,aa,lgDensity);

		/* now print hydrogen, special since part of title goes here */
		fprintf( ioMEAN, "\n Hydrogen  " );
		for( i=0; i < 3; i++ )
		{
			fprintf( ioMEAN, "%7.3f", aa[i] );
		}
		fprintf(ioMEAN," (H2)");
		if( chType=='i' && lgDensity )
		{
			fprintf( ioMEAN, 
				 "         Log10 Mean Ionisation (over %s*electron density)\n", type[d] );
		}
		else if( chType=='i' )
		{
			fprintf( ioMEAN, 
				 "                 Log10 Mean Ionisation (over %s)\n", type[d] );
		}
		else if( chType=='t' && lgDensity )
		{
			fprintf( ioMEAN, 
				 "          Log10 Mean Temperature (over %s*electron density)\n", type[d] );
		}
		else if( chType=='t' )
		{
			fprintf( ioMEAN, 
				 "                  Log10 Mean Temperature (over %s)\n", type[d] );
		}
		else
		{
			fprintf( ioQQQ, " PrtMeanIon called with insane job: %c\n", chType );
			TotalInsanity();
		}

		/* ionization fractions for remaining elements */
		for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				mean.MeanIon(chType,nelem,d,&n,aa,lgDensity);
				limit = MIN2(17,n);
				fprintf( ioMEAN, " %10.10s", elementnames.chElementName[nelem]  );

				for( i=0; i < limit; i++ )
				{
					fprintf( ioMEAN, "%7.3f", aa[i] );
				}
				fprintf( ioMEAN, "\n" );

				if( n > 17 )
				{
					lgPrtLots = true;
					fprintf( ioMEAN, "           " );
					for( i=17; i < n; i++ )
					{
						fprintf( ioMEAN, "%7.3f", aa[i] );
					}
					fprintf( ioMEAN, "\n" );
				}
			}
		}

		fprintf( ioMEAN, "\n         " );
		for( i=1; i <= 17; i++ )
		{
			fprintf( ioMEAN, "%7ld", i );
		}
		fprintf( ioMEAN, "\n" );
		if( lgPrtLots )
		{
			fprintf( ioMEAN, "         " );
			for( i=18; i <= LIMELM; i++ )
			{
				fprintf( ioMEAN, "%7ld", i );
			}
			fprintf( ioMEAN, "\n" );
		}
	}
	return;
}
