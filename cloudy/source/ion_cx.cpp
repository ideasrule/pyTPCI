/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*IonNelem ionization balance for elements without special cases */
#include "cddefines.h"
#include "ionbal.h"
#include "atmdat.h"
#include "iso.h"
#include "dense.h"

/** charge exchange ionization / recombination */
void ion_CX( long nelem )
{
	DEBUG_ENTRY( "ion_CX()" );

	ASSERT( nelem < LIMELM);
	ASSERT( nelem > 1 );

	long limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);

	for( long ion=0; ion < dense.IonLow[nelem]; ion++ )
		ionbal.CX_recomb_rate_used[nelem][ion] = 0;

	// array goes up to [nelem+1] since no ion stage above bare nucleus
	for( long ion=limit+1; ion <= nelem; ion++ )
		ionbal.CX_recomb_rate_used[nelem][ion] = 0;

	for( long ion=dense.IonLow[nelem]; ion <= limit; ion++ )
	{
		/* number of bound electrons of the ion after recombination,
		 * for an atom (ion=0) this is equal to nelem+1,
		 * the element on the physical scale, since nelem is
		 * on the C scale, being 5 for carbon */

		ionbal.CX_recomb_rate_used[nelem][ion] = 0.0;
		for (long nelem1=0; nelem1<t_atmdat::NCX; ++nelem1)
		{
			long ipISO=nelem1;
			ionbal.CX_recomb_rate_used[nelem][ion] +=
				/* nelem1^0 + ion charge transfer recombination */
				atmdat.CharExcRecTo[nelem1][nelem][ion]*
				/* following is density [cm-3] of nelem1^0 */
				iso_sp[ipISO][nelem1].st[0].Pop();
		}
	}

	return;
}
