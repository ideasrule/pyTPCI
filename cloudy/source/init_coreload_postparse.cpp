/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitCoreloadPostparse initialization at start, called from cloudy
* after parser one time per core load */
#include "cddefines.h" 
#include "init.h"
#include "dense.h"
#include "iso.h"
#include "species.h"

/*InitCoreloadPostparse initialization at start, called from cloudy
* after parser, one time per core load */
void InitCoreloadPostparse( void )
{

	static int nCalled = 0;

	DEBUG_ENTRY( "InitCoreloadPostparse()" );

	/* only do this once per coreload */
	if( nCalled > 0 )
	{
		return;
	}

	/* this is first call, increment the nCalled counter so we never do this again */
	++nCalled;

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem<LIMELM; ++nelem)
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				iso_update_num_levels( ipISO, nelem );
				ASSERT( iso_sp[ipISO][nelem].numLevels_max > 0 );
				iso_ctrl.nLyman_alloc[ipISO] = iso_ctrl.nLyman[ipISO];
				iso_ctrl.nLyman_max[ipISO] = iso_ctrl.nLyman[ipISO];
				// resolved and collapsed levels
				long numLevels = iso_sp[ipISO][nelem].numLevels_max;
				// "extra" Lyman lines
				numLevels += iso_ctrl.nLyman_alloc[ipISO] - 2;
				// satellite lines (one for doubly-excited continuum)
				if( iso_ctrl.lgDielRecom[ipISO] )
					numLevels += 1;
				iso_sp[ipISO][nelem].st.init( makeChemical( nelem, nelem-ipISO ).c_str(), numLevels );
			}
		}
	}

	return;
}
