/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvIoniz called by ConvEdenIonz, it calls ConvBase until converged */
#include "cddefines.h"
#include "thermal.h"
#include "trace.h"
#include "conv.h"

/* this routine is called by ConvEdenIoniz, it calls ConvBase
 * until it converges or overruns the loop limit */
void ConvIoniz()
{
	DEBUG_ENTRY( "ConvIoniz()" );

	/* expand limit to number of calls to ConvBase during search phase */
	int LoopLimit = conv.lgSearch ? 20 : 10;

	/* this is ionization/electron density convergence loop
	 * keep calling ConvBase until lgIonDone is true */
	for( int i=0; i < LoopLimit; ++i )
	{
		/* compute the current ionization, ots rates, secondary ionization rates */
		ConvBase(i);

		if( trace.nTrConvg >= 4 )
		{
			/* cooling has not been evaluated yet */
			fprintf( ioQQQ, "    ConvIoniz4 %d heat: %.2e cool: %.2e ", 
				 i, thermal.htot , thermal.ctot );

			/* this is flag saying whether or not ionization/eden has converged */
			if( conv.lgConvIoniz() )
			{
				fprintf( ioQQQ, " ioniz converged\n" );
			}
			else
			{
				fprintf( ioQQQ, " ioniz no conv: %s old %.4e new %.4e OscilOTS %c\n", 
							conv.chConvIoniz() , 
							conv.convIonizOldVal() ,
							conv.convIonizNewVal() ,
				  TorF(conv.lgOscilOTS));
			}
		}

		if( conv.lgConvIoniz() )
			break;
		
	}

	if( trace.nTrConvg>=4 )
	{
		if (! conv.lgConvIoniz())
		{
			fprintf( ioQQQ, 
						"    ConvIoniz4>>>>>>>>>>exit without converging after %i tries!!!!\n", LoopLimit);
		}
		/* if trace convergence is in operation and we did not converge, give warning */
		//ConvFail("ioni","");
		//return 1;
	}
}
