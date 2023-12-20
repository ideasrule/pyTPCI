/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_trim raise or lower most extreme stages of ionization considered,
 * called by ConvBase - ion limits were originally set by  */
#include "cddefines.h"
#include "ion_trim.h"

#include "heavy.h"
#include "conv.h"
#include "rfield.h"
#include "phycon.h"
#include "mole.h"
#include "thermal.h"
#include "iso.h"
#include "struc.h"
#include "ionbal.h"
#include "dense.h"
#include "dynamics.h"

void ion_trim_untrim( long nelem )
{
	DEBUG_ENTRY( "ion_trim_untrim()" );
	dense.IonLow[nelem] = 0;
	dense.IonHigh[nelem] = nelem+1;
}
void ion_trim_invalidate( long nelem )
{
	DEBUG_ENTRY( "ion_trim_invalidate()" );
	dense.IonLow[nelem] = -1;
	dense.IonHigh[nelem] = -1;
}

STATIC void ion_trim_from_set( long nelem);

void ion_trim_init()
{
	DEBUG_ENTRY( "ion_trim_init()" );
	for(long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		if( !dense.lgElmtOn[nelem] )
		{
			ion_trim_invalidate(nelem);
			continue;
		}
		
		if( dense.lgSetIoniz[nelem] )
		{
			/* >>chng 04 jan 13, add this test, caught by Orly Gnat */
			/* check on actual zero abundances of lower stages - this will only 
			 * happen when ionization is set with element ionization command */
			ion_trim_from_set(nelem);
		}
		else
		{
			// IonHigh[n] is the highest stage of ionization present
			// the IonHigh array index is on the C scale, so [0] is hydrogen
			// the value is also on the C scale, so element [nelem] can range
			// from 0 to nelem+1 
			dense.IonHigh[nelem] = nelem + 1;
			
			dense.IonLow[nelem] = 0;
			// for very intense radiation fields very heavy elements, N>Fe, will fail
			// in ion_solver due to ill conditioned matrix.  all populations are in
			// fully stripped state.  Start will fully stripped ion distribution in this case.
			if( rfield.uh > 1e15 )
			{
				//trim down highest stage to be within incident radiation field
				while ( ionbal.lgTrimhiOn &&
						  rfield.anu(Heavy.ipHeavy[nelem][dense.IonHigh[nelem]-1]) > 
						  rfield.anu(rfield.nflux) && dense.IonHigh[nelem]>1 ) 
					--dense.IonHigh[nelem];
				
				if (ionbal.lgTrimloOn)
					dense.IonLow[nelem] = max( 0 , dense.IonHigh[nelem]-1 );
			}
		}				
		
		// make low-stage populations zero
		for( long ion=0; ion<dense.IonLow[nelem]; ++ion )
		{
			dense.xIonDense[nelem][dense.IonLow[nelem]] += dense.xIonDense[nelem][ion];
			dense.xIonDense[nelem][ion] = 0.;
		}
		for( long ion=nelem+1; ion>dense.IonHigh[nelem]; --ion )
		{
			dense.xIonDense[nelem][dense.IonHigh[nelem]] += dense.xIonDense[nelem][ion];
			dense.xIonDense[nelem][ion] = 0.;
		}
	}
}

void ion_trim_from_set()
{
	DEBUG_ENTRY( "ion_trim_from_set()" );
	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		if( dense.lgSetIoniz[nelem] )
		{
			ion_trim_from_set(nelem);
		}			  
	}
}

STATIC void ion_trim_from_set( long nelem )
{
	DEBUG_ENTRY( "ion_trim_from_set()" );
	dense.IonLow[nelem] = 0;
	dense.IonHigh[nelem] = nelem + 1;
	while( ionbal.lgTrimloOn && dense.SetIoniz[nelem][dense.IonLow[nelem]] < dense.density_low_limit )
		++dense.IonLow[nelem];
	while( ionbal.lgTrimhiOn && dense.SetIoniz[nelem][dense.IonHigh[nelem]] < dense.density_low_limit )
		--dense.IonHigh[nelem];
}

void ion_trim_small (long nelem, double abund_total)
{
	DEBUG_ENTRY( "ion_trim_small()" );
	while( ionbal.lgTrimhiOn && dense.IonHigh[nelem] > dense.IonLow[nelem] && 
		dense.xIonDense[nelem][dense.IonHigh[nelem]] < 1e-25*abund_total &&
		dense.IonHigh[nelem] > 1)
	{
		ASSERT( dense.xIonDense[nelem][dense.IonHigh[nelem]] >= 0. );
		/* zero out abundance and heating due to stage of ionization we are about to zero out */
		dense.xIonDense[nelem][dense.IonHigh[nelem]] = 0.;
		thermal.setHeating(nelem,dense.IonHigh[nelem]-1,0.);
		/* decrement counter */
		--dense.IonHigh[nelem];
	}
}

void mole_ion_trim(void)
{
	DEBUG_ENTRY( "mole_ion_trim()" );
	for (ChemNuclideList::iterator atom = nuclide_list.begin();
		  atom != nuclide_list.end(); ++atom)
	{
		if( !(*atom)->lgHasLinkedIon())
			continue;
		const long int nelem=(*atom)->el()->Z-1;
		if( !dense.lgElmtOn[nelem] )
			continue;

		for (long int ion=0;ion<nelem+2;ion++) 
		{
			if ((*atom)->ipMl[ion] != -1)
			{
				if ( ionbal.lgTrimloOn && dense.IonLow[nelem] > ion )
				{
					if( dense.xIonDense[nelem][ion] > ( ionbal.trimlo * dense.gas_phase[nelem] ) )
					{
						dense.IonLow[nelem] = ion;
					}
					else
					{
						dense.xIonDense[nelem][dense.IonLow[nelem]] += dense.xIonDense[nelem][ion];
						dense.xIonDense[nelem][ion] = 0.;
					}
				}

				if ( ionbal.lgTrimhiOn && dense.IonHigh[nelem] < ion )
				{
					if( dense.xIonDense[nelem][ion] > ( ionbal.trimhi * dense.gas_phase[nelem] ) )
					{
						dense.IonHigh[nelem] = ion;
					}
					else
					{
						dense.xIonDense[nelem][dense.IonHigh[nelem]] += dense.xIonDense[nelem][ion];
						dense.xIonDense[nelem][ion] = 0.;
					}
				}
			}
		}
	}
}

void ion_trim(
	/* nelem is on the C scale, 0 for H, 5 for C */
	long int nelem )
{

	/* this will remember that higher stages trimed up or down */
	bool lgHi_Down = false;
	bool lgHi_Up = false;
	bool lgHi_Up_enable;
	/* this will remember that lower stages trimmed up or own*/
	bool lgLo_Up = false;
	bool lgLo_Down = false;
	long int ion_lo_save = dense.IonLow[nelem],
		ion_hi_save = dense.IonHigh[nelem];
	long int ion;
	realnum trimhi , trimlo;

	/* this is debugging code that turns on print after certain number of calls */
	/*realnum xlimit = SMALLFLOAT *10.;*/
	/*static int ncall=0;
	if( nelem==5 )
		++ncall;*/

	DEBUG_ENTRY( "ion_trim()" );

	/* this routine should not be called early in the simulation, before
	 * things have settled down */
	ASSERT( conv.nTotalIoniz>2 );

	/*confirm all ionization stages are within their range of validity */
	ASSERT( nelem >= ipHELIUM && nelem < LIMELM );
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonHigh[nelem] <= nelem+1 );
	/* IonHigh can be equal to IonLow */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	/* during search phase of mostly neutral matter the electron density
	 * can be vastly too large, and the ionization suppressed as a result.
	 * during search do not trim down or up as much */
	if( conv.lgSearch )
	{
		trimhi = (realnum)ionbal.trimhi * 1e-4f;
		trimlo = (realnum)ionbal.trimlo * 1e-4f;
	}
	else
	{
		trimhi = (realnum)ionbal.trimhi;
		trimlo = (realnum)ionbal.trimlo;
	}

	/* helium is special case since abundance is so high, and He+ CT with 
	 * CO is the dominant CO destruction process in molecular regions */
	if( nelem == ipHELIUM )
	{
		/* never want to trip up a lower stage of ionization */
		trimlo = SMALLFLOAT;

		/* if He+ is highest stage of ionization, probably want to keep it
		 * since important for CO chemistry in molecular clouds */
		if( dense.IonHigh[ipHELIUM] == 1 )
		{
			trimhi = MIN2( trimhi , 1e-20f );
		}
		else if( dense.IonHigh[ipHELIUM] == 2 )
		{
			if( conv.lgSearch )
			{
				/* during search phase we may be quite far from solution, in the
				 * case of a PDR sim the electron density may be vastly higher 
				 * than it will be with stable solution.  He++ can be very important
				 * for the chemistry and we want to consider it.  Make the
				 * threshold for ignoring He++ much higher than normal to prevent
				 * a premature removal of the ion */
				trimhi = MIN2( trimhi , 1e-17f ); 
			}
			else
			{
				/* similar smaller upper limit for ion*/
				trimhi = MIN2( trimhi , 1e-12f ); 
			}
		}
	}

	/* logic for PDRs, for elements included in chemistry, need stable solutions, 
	 * keep 3 ion stages in most cases - added by NA to do HII/PDR sims */
	if( !mole_global.lgNoMole )
	{
		trimlo = SMALLFLOAT;
		if(dense.IonHigh[nelem] ==2)
		{
			trimhi = MIN2(trimhi, 1e-20f);
		}
	}

	/* raise or lower highest and lowest stages of ionization to
	 * consider only significant stages
	 * IonHigh[nelem]  stage of ionization  */

	/* this is a special block for initialization only - it checks on absurd abundances
	 * and can trim multiple stages of ionization at one time. */
	if( conv.lgSearch )
	{
		/* special - trim down higher stages if they have essentially zero abundance */
		while( ionbal.lgTrimhiOn && 
			(dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem] < dense.density_low_limit || 
			dense.xIonDense[nelem][dense.IonHigh[nelem]] < dense.density_low_limit ) && 
			/* >>chng 02 may 12, rm +1 since this had effect of not allowing fully atomic */
			dense.IonHigh[nelem] > dense.IonLow[nelem] )
		{
			/* dense.xIonDense[nelem][i] is density of i+1 th ionization stage (cm^-3)
			 * the -1 is correct for the heating, -1 since heating is caused by ionization of stage below it */
			dense.xIonDense[nelem][dense.IonHigh[nelem]-1] +=
				dense.xIonDense[nelem][dense.IonHigh[nelem]];
			dense.xIonDense[nelem][dense.IonHigh[nelem]] = 0.;
			thermal.setHeating(nelem,dense.IonHigh[nelem]-1,0.);
			if( dense.IonHigh[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonHigh[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
					iso_sp[ipISO][nelem].st[level].Pop() = 0.;
			}

			/* decrement high stage limit */
			--dense.IonHigh[nelem];
			ASSERT( dense.IonHigh[nelem] >= 0);
			/* remember this happened */
			lgHi_Down = true;
		}

		/* special - trim up lower stages trim if they have essentially zero abundance */
		while( ionbal.lgTrimloOn && 
			(dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] < dense.density_low_limit || 
			dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit ) && 
			dense.IonLow[nelem] < dense.IonHigh[nelem] - 1 )
		{
			/* dense.xIonDense[nelem][i] is density of ith ionization stage (cm^-3)
			 * there is no-1 since we are removing the agent that heats */
			dense.xIonDense[nelem][dense.IonLow[nelem]+1] += 
				dense.xIonDense[nelem][dense.IonLow[nelem]];
			dense.xIonDense[nelem][dense.IonLow[nelem]] = 0.;
			/* >>chng 01 aug 04, remove -1 which clobbers thermal.heating when IonLow == 0 */
			thermal.setHeating(nelem,dense.IonLow[nelem],0.);
			if( dense.IonLow[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonLow[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
					iso_sp[ipISO][nelem].st[level].Pop() = 0.;
			}

			/* increment low stage limit */
			++dense.IonLow[nelem];
			lgLo_Up = true;
		}
	}

	/* sanity check */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	/* this checks on error condition where the gas is stupendously highly ionized,
	 * the low limit is already one less than high, and we need to raise the low
	 * limit further */
	if( dense.IonHigh[nelem] > 1 &&
		(dense.IonLow[nelem]==dense.IonHigh[nelem]-1) &&
		(dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit) )
	{
		dense.xIonDense[nelem][dense.IonLow[nelem]] = (realnum)dense.density_low_limit;
		return;
	}

	/* trim down high stages that have too small an abundance */
	while( ionbal.lgTrimhiOn && 
		dense.IonHigh[nelem] > dense.IonLow[nelem]  && 
		( (dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem] <= 
		trimhi ) ||
		(dense.xIonDense[nelem][dense.IonHigh[nelem]] <= dense.density_low_limit ) )
		) 
	{
		/* >>chng 03 sep 30, the atom and its first ion are a special case
		 * since we want to compute even trivial ions in molecular clouds */
		if( dense.IonHigh[nelem]>1 ||
			(dense.IonHigh[nelem]==1&&dense.xIonDense[nelem][1]<100.*dense.density_low_limit) )
		{
			dense.xIonDense[nelem][dense.IonHigh[nelem]-1] += 
				dense.xIonDense[nelem][dense.IonHigh[nelem]];
			dense.xIonDense[nelem][dense.IonHigh[nelem]] = 0.;
			thermal.setHeating(nelem,dense.IonHigh[nelem]-1,0.);
			if( dense.IonHigh[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonHigh[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
					iso_sp[ipISO][nelem].st[level].Pop() = 0.;
			}
			--dense.IonHigh[nelem];
			lgHi_Down = true;
		}
		else
		{
			break;
		}
	}

	/* trim up highest stages - will this be done? */
	lgHi_Up_enable = true;
	/* trimming can be turned off */
	if( !ionbal.lgTrimhiOn )
		lgHi_Up_enable = false;
	/* high stage is already fully stripped */
	if( dense.IonHigh[nelem]>=nelem+1 )
		lgHi_Up_enable = false;
	/* we have previously trimmed this stage down */
	if( lgHi_Down )
		lgHi_Up_enable = false;
	/* we are in neutral gas */
	if( dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN]<0.9  )
		lgHi_Up_enable = false;

	if( lgHi_Up_enable )
	{
		realnum abundold = 0;
		if( nzone>2 ) 
			abundold = struc.xIonDense[nzone-3][nelem][dense.IonHigh[nelem]]/
			SDIV( struc.gas_phase[nzone-3][nelem]);
		realnum abundnew = dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem];
		/* only raise highest stage if ionization potential of next highest stage is within
		 * continuum array and the abundance of the highest stage is significant */
		if( (Heavy.Valence_IP_Ryd[nelem][dense.IonHigh[nelem]] < rfield.anu(rfield.nflux-1) ||
			Heavy.Valence_IP_Ryd[nelem][dense.IonHigh[nelem]] < phycon.te_ryd*10.) &&
			dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem] > trimhi &&
			/* this checks that abundance of highest stage is increasing */
			abundnew > abundold*1.01 )
		{
			/*fprintf(ioQQQ,"uuppp %li %li \n", nelem, dense.IonHigh[nelem] );*/
			/* raise highest level of ionization */
			++dense.IonHigh[nelem];
			lgHi_Up = true;
			/* set abundance to almost zero so that sanity check elsewhere will be ok */
			dense.xIonDense[nelem][dense.IonHigh[nelem]] = (dense.density_low_limit*10.);
			dense.xIonDense[nelem][dense.IonHigh[nelem]-1] -= 
				dense.xIonDense[nelem][dense.IonHigh[nelem]];
			long int ipISO = nelem - dense.IonHigh[nelem];
			if (ipISO >= 0 && ipISO < NISO)
				iso_sp[ipISO][nelem].st[0].Pop() = dense.xIonDense[nelem][dense.IonHigh[nelem]];
		}
	}

	/* sanity check */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	/* lower lowest stage of ionization if we have significant abundance at current lowest */
	if( dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] > trimlo &&
		dense.IonLow[nelem] > 0 )
	{
		/* lower lowest level of ionization */
		--dense.IonLow[nelem];
		lgLo_Down = true;
		/* >>chng 01 aug 02, set this to zero so that sanity check elsewhere will be ok */
		dense.xIonDense[nelem][dense.IonLow[nelem]] = (realnum)dense.density_low_limit;
		dense.xIonDense[nelem][dense.IonLow[nelem]+1] -= 
			dense.xIonDense[nelem][dense.IonLow[nelem]];
		long int ipISO = nelem - dense.IonLow[nelem];
		if (ipISO >= 0 && ipISO < NISO)
			iso_sp[ipISO][nelem].st[0].Pop() = dense.xIonDense[nelem][dense.IonLow[nelem]];
	}

	/* raise lowest stage of ionization, but only if we are near illuminated face of cloud*/
	else if( ionbal.lgTrimloOn && nzone < 10 &&
		(dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] <= (realnum)trimlo) && 
		(dense.IonLow[nelem] < (dense.IonHigh[nelem] - 2) ) )
	{
		/* raise lowest level of ionization */
		dense.xIonDense[nelem][dense.IonLow[nelem]+1] += 
			dense.xIonDense[nelem][dense.IonLow[nelem]];
		dense.xIonDense[nelem][dense.IonLow[nelem]] = 0.;
		/* no minus one in below since this is low bound, already bounds at atom */
		thermal.setHeating(nelem,dense.IonLow[nelem],0.);
		if( dense.IonLow[nelem] == nelem+1-NISO )
		{
			long int ipISO = nelem - dense.IonLow[nelem];
			ASSERT( ipISO>=0 && ipISO<NISO );
			for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
				iso_sp[ipISO][nelem].st[level].Pop() = 0.;
		}
		++dense.IonLow[nelem];
		lgLo_Up = true;
	}
	/* test on zero */
	else if( ionbal.lgTrimloOn && (dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit) && 
		/*>>chng 06 may 24, from IonLow < IonHigh to <IonHigh-1 because
		 * old test would allow low to be raised too high, which then
		 * leads to insanity */
		(dense.IonLow[nelem] < dense.IonHigh[nelem]-1) )
	{
		while(dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit && 
			/* >>chng 07 feb 27 from < IonHigh to < IonHigh-1 to prevent raising
			 * low to high */
			(dense.IonLow[nelem] < dense.IonHigh[nelem]-1) )
		{
			/* raise lowest level of ionization */
			dense.xIonDense[nelem][dense.IonLow[nelem]+1] += 
				dense.xIonDense[nelem][dense.IonLow[nelem]];
			dense.xIonDense[nelem][dense.IonLow[nelem]] = 0.;
			/* no minus one in below since this is low bound, already bounds at atom */
			thermal.setHeating(nelem,dense.IonLow[nelem],0.);
			if( dense.IonLow[nelem] == nelem+1-NISO )
			{
				long int ipISO = nelem - dense.IonLow[nelem];
				ASSERT( ipISO>=0 && ipISO<NISO );
				for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
					iso_sp[ipISO][nelem].st[level].Pop() = 0.;
			}
			++dense.IonLow[nelem];
			lgLo_Up = true;
		}
	}

	/* sanity check */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	/* these are standard bounds checks that appear throughout this routine
	 * dense.xIonDense[IonLow] is first one with positive abundances
	 * 
	 * in case where lower ionization stage was just lowered the abundance
	 * was set to dense.density_low_limit so test must be < dense.density_low_limit */
	ASSERT( dense.IonLow[nelem] <= 1 ||
		dense.xIonDense[nelem][dense.IonLow[nelem]-1] == 0. );

	ASSERT( (dense.IonLow[nelem]==0 && dense.IonHigh[nelem]==0 ) || lgLo_Up ||
			  ! ionbal.lgTrimloOn ||
			  dense.xIonDense[nelem][dense.IonLow[nelem]] >= dense.density_low_limit ||
			  dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] >= dense.density_low_limit ||
			  /*>>chng 06 may 24, include this to cover case where cap is set by not allowing
				* lower limit to come up to upper limit */
			  (dense.IonLow[nelem]==dense.IonHigh[nelem]-1 ));

	{
		/* option to print out what has happened so far .... */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nelem == ipHELIUM )
		{
			if( lgHi_Down ||lgHi_Up ||lgLo_Up ||lgLo_Down )
			{
				fprintf(ioQQQ,"DEBUG TrimZone\t%li\t",nzone );
				if(  lgHi_Down )
				{
					fprintf(ioQQQ,"high dn %li to %li",
						ion_hi_save , 
						dense.IonHigh[nelem] );
				}
				if( lgHi_Up )
				{
					fprintf(ioQQQ,"high up %li to %li",
						ion_hi_save , 
						dense.IonHigh[nelem] );
				}
				if( lgLo_Up )
				{
					fprintf(ioQQQ,"low up %li to %li",
						ion_lo_save , 
						dense.IonLow[nelem] );
				}
				if( lgLo_Down )
				{
					fprintf(ioQQQ,"low dn %li to %li",
						ion_lo_save , 
						dense.IonLow[nelem] );
				}
				fprintf(ioQQQ,"\n" );
			}
		}
	}

	/* option to print if any trimming occurred */
	if( lgHi_Down || lgHi_Up || lgLo_Up || lgLo_Down )
	{
		conv.lgIonStageTrimed = true;
		{
			/* option to print out what has happened so far .... */
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC && nelem==ipHELIUM )
			{
				fprintf(ioQQQ,"DEBUG ion_trim zone\t%.2f \t%li\t", fnzone, nelem);
				if( lgHi_Down  )
					fprintf(ioQQQ,"\tHi_Down");
				if( lgHi_Up )
					fprintf(ioQQQ,"\tHi_Up");
				if( lgLo_Up )
					fprintf(ioQQQ,"\tLo_Up");
				if( lgLo_Down )
					fprintf(ioQQQ,"\tLo_Down");
				fprintf(ioQQQ,"\n");
			}
		}
	}

	/* asserts that that densities of ion stages that are now turned 
	 * off are zero */
	for( ion=0; ion<dense.IonLow[nelem]; ++ion )
	{
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}
	for( ion=dense.IonHigh[nelem]+1; ion<nelem+1; ++ion )
	{
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}

	for( ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
	{
		/* in case where ionization stage was changed the
		 * density is zero, we set it to SMALLFLOAT since a density
		 * of zero is not possible */
		dense.xIonDense[nelem][ion] = MAX2(dense.xIonDense[nelem][ion], SMALLFLOAT);
	}
	return;
}

// Try winnowing down ion_trim
void ion_trim2(
	/* nelem is on the C scale, 0 for H, 5 for C */
	long int nelem,
	long int oldIonRange[2]
	)
{

	/* this will remember that higher stages trimed up or down */
	bool lgHi_Down = false;
	/* this will remember that lower stages trimmed up or own*/
	bool lgLo_Up = false;
	long int ion_lo_save = dense.IonLow[nelem],
		ion_hi_save = dense.IonHigh[nelem];
	long int ion;
	realnum trimhi , trimlo, trimcharge;
	long int const minstates = 3;
	
	/* this is debugging code that turns on print after certain number of calls */
	/*realnum xlimit = SMALLFLOAT *10.;*/
	/*static int ncall=0;
	if( nelem==5 )
		++ncall;*/

	DEBUG_ENTRY( "ion_trim2()" );

	/*confirm all ionization stages are within their range of validity */
	ASSERT( nelem >= ipHYDROGEN && nelem < LIMELM );
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonHigh[nelem] <= nelem+1 );
	/* IonHigh can be equal to IonLow */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	// Trimming limit on fraction of total charge in system
	trimcharge = 1e-6;

	/* during search phase of mostly neutral matter the electron density
	 * can be vastly too large, and the ionization suppressed as a result.
	 * during search do not trim down or up as much */
	if( conv.lgSearch )
	{
		trimhi = (realnum)ionbal.trimhi * 1e-4f;
		trimlo = (realnum)ionbal.trimlo * 1e-4f;
	}
	else
	{
		trimhi = (realnum)ionbal.trimhi;
		trimlo = (realnum)ionbal.trimlo;
	}

	/* helium is special case since abundance is so high, and He+ CT with 
	 * CO is the dominant CO destruction process in molecular regions */
	if( nelem == ipHELIUM )
	{
		/* never want to trip up a lower stage of ionization */
		trimlo = SMALLFLOAT;

		/* if He+ is highest stage of ionization, probably want to keep it
		 * since important for CO chemistry in molecular clouds */
		if( dense.IonHigh[ipHELIUM] == 1 )
		{
			trimhi = MIN2( trimhi , 1e-20f );
		}
		else if( dense.IonHigh[ipHELIUM] == 2 )
		{
			if( conv.lgSearch )
			{
				/* during search phase we may be quite far from solution, in the
				 * case of a PDR sim the electron density may be vastly higher 
				 * than it will be with stable solution.  He++ can be very important
				 * for the chemistry and we want to consider it.  Make the
				 * threshold for ignoring He++ much higher than normal to prevent
				 * a premature removal of the ion */
				trimhi = MIN2( trimhi , 1e-17f ); 
			}
			else
			{
				/* similar smaller upper limit for ion*/
				trimhi = MIN2( trimhi , 1e-12f ); 
			}
		}
	}

	/* logic for PDRs, for elements included in chemistry, need stable solutions, 
	 * keep 3 ion stages in most cases - added by NA to do HII/PDR sims */
	if( !mole_global.lgNoMole )
	{
		trimlo = SMALLFLOAT;
		if(dense.IonHigh[nelem] ==2)
		{
			trimhi = MIN2(trimhi, 1e-20f);
		}
	}

	/* raise or lower highest and lowest stages of ionization to
	 * consider only significant stages
	 * IonHigh[nelem]  stage of ionization  */

	/* sanity check */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	/* this checks on error condition where the gas is stupendously highly ionized,
	 * the low limit is already one less than high, and we need to raise the low
	 * limit further */
	if( dense.IonHigh[nelem] > 1 &&
		(dense.IonLow[nelem]==dense.IonHigh[nelem]-1) &&
		(dense.xIonDense[nelem][dense.IonLow[nelem]] < dense.density_low_limit) )
	{
		dense.xIonDense[nelem][dense.IonLow[nelem]] = (realnum)dense.density_low_limit;
		return;
	}

	// Some hysteresis in trimming criteria should help prevent oscillations
	double threshextend = 1.0, threshkeep = 0.125;
	double itrim =  ( dense.IonHigh[nelem] > oldIonRange[1] ) ? threshextend : threshkeep;
	/* trim down high stages that have too small an abundance */
	if ( ionbal.lgTrimhiOn && 
	     dense.IonHigh[nelem] > dense.IonLow[nelem] )
	{
		while (
			( ( dense.xIonDense[nelem][dense.IonHigh[nelem]]/dense.gas_phase[nelem] <= trimhi*itrim &&
			    dense.xIonDense[nelem][dense.IonHigh[nelem]] <= trimcharge*dense.EdenTrue*itrim ) ||
			  dense.xIonDense[nelem][dense.IonHigh[nelem]] <= dense.density_low_limit*itrim )				 
			)
		{
			if ( dense.IonHigh[nelem] <= dense.IonLow[nelem] + minstates-1 )
				break;
			itrim = threshkeep;
			dense.xIonDense[nelem][dense.IonHigh[nelem]-1] += 
				dense.xIonDense[nelem][dense.IonHigh[nelem]];
			dense.xIonDense[nelem][dense.IonHigh[nelem]] = 0.;
			
			if( dynamics.lgAdvection || dynamics.lgTimeDependentStatic )
			{
				dynamics.Source[nelem][dense.IonHigh[nelem]-1] += 
					dynamics.Source[nelem][dense.IonHigh[nelem]];
				dynamics.Source[nelem][dense.IonHigh[nelem]] = 0.;
			}

			thermal.setHeating(nelem,dense.IonHigh[nelem]-1,0.);
			long int ipISO = nelem - dense.IonHigh[nelem];
			if ( ipISO>=0 && ipISO<NISO )
			{
				for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
					iso_sp[ipISO][nelem].st[level].Pop() = 0.;
			}
			--dense.IonHigh[nelem];
		}
		lgHi_Down = dense.IonHigh[nelem] < oldIonRange[1];
	}

	/* sanity check */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	itrim = ( dense.IonLow[nelem] < oldIonRange[0] ) ? threshextend : threshkeep ;
	/* raise lowest stage of ionization, but only if we are near illuminated face of cloud*/
	if ( ionbal.lgTrimloOn )
	{
		while( (dense.xIonDense[nelem][dense.IonLow[nelem]]/dense.gas_phase[nelem] <= (realnum)trimlo*itrim) && 
		       (dense.IonLow[nelem] < (dense.IonHigh[nelem] - (minstates - 1)) ) )
		{
			itrim = threshkeep;
			/* raise lowest level of ionization */
			dense.xIonDense[nelem][dense.IonLow[nelem]+1] += 
				dense.xIonDense[nelem][dense.IonLow[nelem]];
			dense.xIonDense[nelem][dense.IonLow[nelem]] = 0.;

			if( dynamics.lgAdvection || dynamics.lgTimeDependentStatic )
			{
				dynamics.Source[nelem][dense.IonLow[nelem]+1] += 
					dynamics.Source[nelem][dense.IonLow[nelem]];
				dynamics.Source[nelem][dense.IonLow[nelem]] = 0.;
			}

			/* no minus one in below since this is low bound, already bounds at atom */
			thermal.setHeating(nelem,dense.IonLow[nelem],0.);
			long int ipISO = nelem - dense.IonLow[nelem];
			if ( ipISO>=0 && ipISO<NISO )
			{
				for( long level = 0; level < iso_sp[ipISO][nelem].numLevels_max; level++ )
					iso_sp[ipISO][nelem].st[level].Pop() = 0.;
			}
			++dense.IonLow[nelem];
		}
		lgLo_Up = dense.IonLow[nelem] > oldIonRange[0];
	}

	/* sanity check */
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );

	/* these are standard bounds checks that appear throughout this routine
	 * dense.xIonDense[IonLow] is first one with positive abundances
	 * 
	 * in case where lower ionization stage was just lowered the abundance
	 * was set to dense.density_low_limit so test must be < dense.density_low_limit */
	ASSERT( dense.IonLow[nelem] <= 1 ||
		dense.xIonDense[nelem][dense.IonLow[nelem]-1] == 0. );

	{
		/* option to print out what has happened so far .... */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nelem == ipHELIUM )
		{
			if( lgHi_Down ||lgLo_Up )
			{
				fprintf(ioQQQ,"DEBUG TrimZone\t%li\t",nzone );
				if(  lgHi_Down )
				{
					fprintf(ioQQQ,"high dn %li to %li",
						ion_hi_save , 
						dense.IonHigh[nelem] );
				}
				if( lgLo_Up )
				{
					fprintf(ioQQQ,"low up %li to %li",
						ion_lo_save , 
						dense.IonLow[nelem] );
				}
				fprintf(ioQQQ,"\n" );
			}
		}
	}

	/* option to print if any trimming occurred */
	if( lgHi_Down || lgLo_Up )
	{
		conv.lgIonStageTrimed = true;
		{
			/* option to print out what has happened so far .... */
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC && nelem==ipHELIUM )
			{
				fprintf(ioQQQ,"DEBUG ion_trim zone\t%.2f \t%li\t", fnzone, nelem);
				if( lgHi_Down  )
					fprintf(ioQQQ,"\tHi_Down");
				if( lgLo_Up )
					fprintf(ioQQQ,"\tLo_Up");
				fprintf(ioQQQ,"\n");
			}
		}
	}

	/* asserts that that densities of ion stages that are now turned 
	 * off are zero */
	for( ion=0; ion<dense.IonLow[nelem]; ++ion )
	{
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}
	for( ion=dense.IonHigh[nelem]+1; ion<nelem+1; ++ion )
	{
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}

	for( ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
	{
		/* in case where ionization stage was changed the
		 * density is zero, we set it to SMALLFLOAT since a density
		 * of zero is not possible */
		dense.xIonDense[nelem][ion] = MAX2(dense.xIonDense[nelem][ion], SMALLFLOAT);
	}
	return;
}

void ion_widen(long nelem)
{
	DEBUG_ENTRY( "ion_widen()" );
	if ( dense.lgSetIoniz[nelem] || 
		 ( dense.IonLow[nelem] == 0 && dense.IonHigh[nelem] == nelem+1 ) )
	{
		return;
	}
	double abund = 0.0;
	for (long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion)
	{
		abund += dense.xIonDense[nelem][ion];
	}
	double abundnew = abund;
	if (dense.IonLow[nelem] > 0 )
	{
		--dense.IonLow[nelem];
		dense.xIonDense[nelem][dense.IonLow[nelem]] = (dense.density_low_limit*10.);
		abundnew += (dense.density_low_limit*10.);
		long ipISO = nelem - dense.IonLow[nelem];
		if (ipISO >= 0 && ipISO < NISO)
			iso_sp[ipISO][nelem].st[0].Pop() = dense.xIonDense[nelem][dense.IonLow[nelem]];
	}
	if ( dense.IonHigh[nelem] < nelem+1)
	{
		++dense.IonHigh[nelem];
		dense.xIonDense[nelem][dense.IonHigh[nelem]] = (dense.density_low_limit*10.);
		abundnew += (dense.density_low_limit*10.);
		long ipISO = nelem - dense.IonHigh[nelem];
		if (ipISO >= 0 && ipISO < NISO)
			iso_sp[ipISO][nelem].st[0].Pop() = dense.xIonDense[nelem][dense.IonHigh[nelem]];
	}
	double frac = abund/abundnew;
	for (long ion=dense.IonLow[nelem];ion<=dense.IonHigh[nelem];++ion)
	{
		dense.xIonDense[nelem][ion] *= frac;
	}
}

#ifdef NDEBUG
void ion_trim_validate(long, bool)
{
	(void)0;
}
#else
void ion_trim_validate(long nelem, bool lgIonizTrimCalled)
{
	DEBUG_ENTRY( "ion_trim_validate()" );
	for( long ion=0; ion<dense.IonLow[nelem]; ++ion )
	{
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}
	/*if( nelem==5 ) fprintf(ioQQQ,"carbbb\t%li\n", dense.IonHigh[nelem]);*/
	for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
	{
		/* >>chng 02 feb 06, had been > o., chng to > SMALLFLOAT to
		 * trip over VERY small floats that failed on alphas, but not 386
		 * 
		 * in case where lower ionization stage was just lowered or
		 * trimmed down the abundance
		 * was set to SMALLFLOAT so test must be < SMALLFLOAT */
		/* >>chng 02 feb 19, add check for search phase.  During this search
		 * models with extreme ionization (all neutral or all ionized) can
		 * have extreme but non-zero abundances far from the ionization peak for
		 * element with lots of electrons.  These will go away once the model
		 * becomes stable */
		/* >>chng 03 dec 01, add check on whether ion trim was called 
		 * conserve.in threw assert when iontrim not called and abund grew small */
		ASSERT( conv.lgSearch || !lgIonizTrimCalled ||
				  /* this can happen if all C is in the form of CO -- 
					  generalize test to being that at least one ion must survive */
				  dense.IonLow[nelem] == dense.IonHigh[nelem] ||
				  dense.xIonDense[nelem][ion] >= SMALLFLOAT || 
				  dense.xIonDense[nelem][ion]/dense.gas_phase[nelem] >= SMALLFLOAT );
	}
	for( long ion=dense.IonHigh[nelem]+1; ion<=nelem+1; ++ion )
	{
		ASSERT( ion >= 0 );
		ASSERT( dense.xIonDense[nelem][ion] == 0. );
	}
}
#endif
