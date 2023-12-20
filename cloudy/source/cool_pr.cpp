/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*coolpr stores coolants before block printed, when printing cooling agents */
#include "cddefines.h"
#include "cooling.h"
#include "thermal.h"

#define	NCOLSAV	100

void coolpr(
	FILE * io,
	/* the line label */
	const char *chLabel, 
	/* the line wavelength */
	realnum lambda, 
	/* the ratio of cooling to total, negative if a heat source  */
	double ratio, 
	/* the job to do, one of "ZERO", "DOIT", or "DONE"  */
	const char *chJOB
	)
{
	static char chLabsv[NCOLSAV][NCOLNT_LAB_LEN+1];

	static char chSig[NCOLSAV];

	long int i, 
	  ipAr[NCOLSAV], 
	  j, 
	  limit;

	static long int nCoolant = 0; 
	static realnum sav[NCOLSAV];

	realnum SavMax, 
	  scratch[NCOLSAV];

	static realnum csav[NCOLSAV];

	DEBUG_ENTRY( "coolpr()" );

	/* routine is called with two flags, "ZERO" and "DONE" to 
	 * initialize and complete the printout.  Any other label is
	 * interpreted as a line label */
	if( strcmp(chJOB,"ZERO") == 0 )
	{
		/* nCoolant is the counter through the array of coolants,
		 * zero it if new job to do */
		nCoolant = 0;
		for( i=0; i<NCOLSAV; ++i )
		{
			scratch[i] = FLT_MAX;
			ipAr[i] = LONG_MAX;
		}
	}

	else if( strcmp(chJOB,"DOIT") == 0 )
	{
		strcpy( chLabsv[nCoolant], chLabel );

		if( lambda < 10000. )
		{
			sav[nCoolant] = lambda;
		}
		else
		{
			sav[nCoolant] = lambda/10000.f;
		}

		csav[nCoolant] = (realnum)ratio;
		/* is this coolant really cooling (+) or a heat source? */
		if( ratio < 0. )
		{
			chSig[nCoolant] = 'n';
		}
		else
		{
			chSig[nCoolant] = ' ';
		}

		/* increment the counter, so this is the number actually in the stack */
		++nCoolant;

		/* this is limit to how much we can save */
		if( nCoolant >= NCOLSAV )
		{
			fprintf( ioQQQ, "  coolpr ran out of room, increase NCOLSAV.\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( strcmp(chJOB,"DONE") == 0 )
	{
		/* want to print sorted list of coolants sorted from strongest to faintest */
		for( i=0; i < nCoolant; i++ )
		{
			/* save abs val so we pick up both heating and cooling */
			scratch[i] = (realnum)fabs(csav[i]);
		}

		for( i=0; i < nCoolant; i++ )
		{
			SavMax = 0.;
			/* following will be reset in following loop */
			ipAr[i] = -LONG_MAX;

			/* find largest of remaining coolants */
			for( j=0; j < nCoolant; j++ )
			{
				if( scratch[j] > SavMax )
				{
					SavMax = scratch[j];
					/* ipAr will point to coolant within saved stack */
					ipAr[i] = j;
				}
			}

			ASSERT( i >= 0 && i < NCOLSAV );
			ASSERT( ipAr[i] >=0 && ipAr[i] < NCOLSAV );
			/* set it to zero so we can look for next strongest */
			scratch[ipAr[i]] = 0.;
		}

		/* now print this stack in order or strength, seven across a line */
		for( j=0; j < nCoolant; j += 7 )
		{
			limit = MIN2(nCoolant,j+7);
			fprintf( io, "     " );
			for( i=j; i < limit; i++ )
			{
				ASSERT( i < NCOLSAV );

				fprintf( io, 
					" %s %.2f%c%6.3f", 
					/* label for the coolant, like "C  4" */
					chLabsv[ipAr[i]], 
					/* wavelength */
					sav[ipAr[i]], 
					/* usually space, but n if negative coolant */
					chSig[ipAr[i]],
					/* fraction of total cooling */
					csav[ipAr[i]] );
			}
			fprintf( io, " \n" );
		}
	}

	else 
	{
		fprintf( ioQQQ, "  coolpr called with insane job =%s=\n",chJOB );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

#undef NCOLSAV
