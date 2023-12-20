/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonNelem ionization balance for elements without special cases */
#include "cddefines.h"
#include "trace.h"
#include "ionbal.h"
#include "dense.h"

void IonNelem(
  /* this is debug flag */
  bool lgPrintIt,
  /* nelem is the atomic number on the C scale, 0 for H */
  long int nelem)
{
	DEBUG_ENTRY( "IonNelem()" );

	if( dense.lgElmtOn[nelem] )
	{
		ion_zero(nelem);

		ion_photo(nelem,lgPrintIt);

		// collisional ionization rates
		ion_collis(nelem);

		// charge exchange
		ion_CX( nelem );

		// total recombination
		ion_recomb(lgPrintIt,nelem);

		/* solve for ionization balance */
		ion_solver(nelem,lgPrintIt);

		if( trace.lgTrace && trace.lgHeavyBug /*|| nelem==ipNICKEL*/ )
		{
			fprintf( ioQQQ, "     IonNelem nelem\t%li\tfnzone\t%6.2f\tfrac\t",
					nelem, fnzone );
			for( int i=0; i < nelem+2; i++ )
			{
				fprintf( ioQQQ, "\t%10.3e", dense.xIonDense[nelem][i]/
				  dense.gas_phase[nelem] );
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	return;
}
