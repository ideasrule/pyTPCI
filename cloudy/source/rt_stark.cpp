/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_stark compute stark broadening escape probabilities using Puetter formalism */
#include "cddefines.h"
#include "iso.h"
#include "dense.h"
#include "phycon.h"
#include "rt.h"

STATIC realnum strkar( long nLo, long nHi );

void RT_stark(void)
{
	long int ipLo, 
	  ipHi;

	double aa , ah, 
	  stark, 
	  strkla;

	DEBUG_ENTRY( "RT_stark()" );

	/* only evaluate one time per zone */
	static long int nZoneEval=-1;
	if( nzone==nZoneEval )
		return;
	nZoneEval = nzone;

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( long nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( !dense.lgElmtOn[nelem] )
				continue;

			t_iso_sp* sp = &iso_sp[ipISO][nelem];

			if( !rt.lgStarkON )
			{
				for( ipHi=0; ipHi < sp->numLevels_max; ipHi++ )
				{
					for( ipLo=0; ipLo < sp->numLevels_max; ipLo++ )
					{
						sp->ex[ipHi][ipLo].pestrk = 0.;
						sp->ex[ipHi][ipLo].pestrk_up = 0.;
					}
				}
				continue;
			}

			/* evaluate Stark escape probability from 
			 * >>ref Puetter Ap.J. 251, 446. */

			/* coefficients for Stark broadening escape probability
			 * to be Puetters AH, equation 9b, needs factor of (Z^-4.5 * (nu*nl)^3 * xl) */
			ah = 6.9e-6*1000./1e12/(phycon.sqrte*phycon.te10*phycon.te10*
			  phycon.te03*phycon.te01*phycon.te01)*dense.eden;

			/* include Z factor */
			ah *= powpq( (double)(nelem+1), -9, 2 );

			/* coefficient for all lines except Ly alpha */
			/* equation 10b, except missing tau^-0.6 */
			stark = 0.264*pow(ah,0.4);

			/* coefficient for Ly alpha */
			/* first few factors resemble equation 13c...what about the rest? */
			strkla = 0.538*ah*4.*9.875*(phycon.sqrte/phycon.te10/phycon.te03);

			long ipHi = iso_ctrl.nLyaLevel[ipISO];
			/* Lyman lines always have outer optical depths */
			/* >>chng 02 mar 31, put in max, crashed on some first iteration 
			 * with negative optical depths,
			 * NB did not understand why neg optical depths started */
			aa = (realnum)max(0.,sp->trans(ipHi,0).Emis().TauIn());
			aa = powpq( aa, -3, 4 );
			sp->ex[ipHi][0].pestrk_up = strkla/2.*MAX2(1.,aa);

			/**\todo	2	- Stark is disabled for now since Lya escape causes density dependent
			 * feedback on the radiative transfer.  Would need to redo the escape
			 * probs every time the electron density is updated - see blr89.in for an 
			 * example */
			sp->ex[ipHi][0].pestrk_up =
					MIN2(.01,sp->ex[ipHi][0].pestrk_up);
			sp->ex[ipHi][0].pestrk_up = 0.;
			sp->ex[ipHi][0].pestrk = sp->ex[ipHi][0].pestrk_up *
					sp->trans(ipHi,0).Emis().Aul();

			/* >>chng 06 aug 28, from numLevels_max to _local. */
			for( ipHi=3; ipHi < sp->numLevels_local; ipHi++ )
			{
				if( sp->trans(ipHi,0).ipCont() <= 0 )
					continue;

				sp->ex[ipHi][0].pestrk_up = stark / 2. * 
					strkar( sp->st[0].n(), sp->st[ipHi].n() ) *
					powpq(MAX2(1.,sp->trans(ipHi,0).Emis().TauIn()),-3,4);

				sp->ex[ipHi][0].pestrk_up = MIN2(.01,sp->ex[ipHi][0].pestrk_up);
				sp->ex[ipHi][0].pestrk = sp->trans(ipHi,0).Emis().Aul()*
					sp->ex[ipHi][0].pestrk_up;
			}

			/* zero out rates above iso.numLevels_local */
			for( ipHi=sp->numLevels_local; ipHi < sp->numLevels_max; ipHi++ )
			{
				sp->ex[ipHi][0].pestrk_up = 0.;
				sp->ex[ipHi][0].pestrk = 0.;
			}

			/* all other lines */
			for( ipLo=ipH2s; ipLo < (sp->numLevels_local - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < sp->numLevels_local; ipHi++ )
				{
					if( sp->trans(ipHi,ipLo).ipCont() <= 0 )
						continue;

					aa = stark *
						strkar( sp->st[ipLo].n(), sp->st[ipHi].n() ) *
						powpq(MAX2(1.,sp->trans(ipHi,ipLo).Emis().TauIn()),-3,4);
					sp->ex[ipHi][ipLo].pestrk_up = 
						(realnum)MIN2(.01,aa);

					sp->ex[ipHi][ipLo].pestrk = sp->trans(ipHi,ipLo).Emis().Aul()*
						sp->ex[ipHi][ipLo].pestrk_up;
				}
			}

			/* zero out rates above iso.numLevels_local */
			for( ipLo=(sp->numLevels_local - 1); ipLo<(sp->numLevels_max - 1); ipLo++ )
			{
				for( ipHi=ipLo + 1; ipHi < sp->numLevels_max; ipHi++ )
				{
					sp->ex[ipHi][ipLo].pestrk_up = 0.;
					sp->ex[ipHi][ipLo].pestrk = 0.;
				}
			}
		}
	}

	return;
}

STATIC realnum strkar( long nLo, long nHi )
{
	return (realnum)pow((realnum)( nLo * nHi ),(realnum)1.2f);
}

