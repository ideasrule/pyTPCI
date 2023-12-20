/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "conv.h"
#include "trace.h"
#include "grainvar.h"
#include "rfield.h"
#include "mole.h"
#include "dense.h"
#include "iso.h"
#include "deuterium.h"

// eden_sum: sum all contributions to the electron density, sets variable dense.EdenTrue
// called by ConvBase - ConvEdenIoniz actually updates the electron density dense.eden
// here we allow dense.EdenTrue to become negative, ConvEdenIoniz will deal with this
// and will also assure that dense.eden is always positive.
void eden_sum()
{
	DEBUG_ENTRY( "eden_sum()" );

	if( dense.EdenSet > 0.f )
	{
		/* electron density set with set eden command */
		dense.EdenTrue = dense.EdenSet;
		dense.eden_from_metals = 1.;

		if( trace.lgTrace || trace.lgESOURCE )
			fprintf( ioQQQ, "     eden_sum zn: %.2f eden set to: %.4e\n", fnzone, dense.EdenSet );
	}
	else if( dense.EdenFraction > 0.f )
	{
		/* electron fraction set with set eden fraction command */
		dense.EdenTrue = dense.EdenFraction*dense.gas_phase[ipHYDROGEN];
		dense.eden_from_metals = 1.;

		if( trace.lgTrace || trace.lgESOURCE )
		{
			fprintf( ioQQQ, "     eden_sum zn: %.2f eden ratio set to: %.4e, eden is: %.4e\n",
				fnzone, dense.EdenFraction, dense.EdenTrue );
		}
	}
	else
	{
		/* EdenExtra is normally zero, set with EDEN command, to add extra e- */
		dense.EdenTrue = dense.EdenExtra;

		/* sum over all ions */
		double eden_ions[LIMELM];
		double sum_all_ions = 0.;
		double sum_metals = 0.;
		for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
		{
			eden_ions[nelem] = 0.;
			for( long ion=1; ion <= nelem+1; ion++ )
				eden_ions[nelem] += ion*dense.xIonDense[nelem][ion];

			sum_all_ions += eden_ions[nelem];
			if( nelem >= ipLITHIUM )
				sum_metals += eden_ions[nelem];
		}
		sum_all_ions += deut.xIonDense[1];
		dense.EdenTrue += sum_all_ions;

		/* electrons contributed by heavy molecules */
		// include electrons lost to negative atomic ions, e.g. H-, but
		// not +ve which are treated by the ionization ladders
		mole.elec = 0.;
		for( long i=0; i < mole_global.num_calc; i++ )
		{
			if( mole_global.list[i]->isIsotopicTotalSpecies() && (!mole_global.list[i]->isMonatomic() || mole_global.list[i]->charge < 0) )
				mole.elec += mole.species[i].den*mole_global.list[i]->charge;
		}

		dense.EdenTrue += mole.elec;
		
		/* gv.lgGrainElectrons - should grain electron source/sink be included in overall electron sum?
		 * default is true, set false with no grain electrons command */
		dense.EdenTrue += gv.TotalEden*gv.lgGrainElectrons;

		/* fraction of electrons from ions heavier than helium */
		dense.eden_from_metals = safe_div( sum_metals, dense.EdenTrue, 1. );

		if( trace.lgTrace || trace.lgESOURCE )
		{
			fprintf( ioQQQ, 
				 "     eden_sum zn: %.2f current: %.4e new true: %.4e ions: %.4e mole: %.4e"
				 " grain: %.4e extra: %.4e LaOTS: %.4e\n",
				 fnzone ,
				 dense.eden , 
				 dense.EdenTrue , 
				 sum_all_ions ,
				 mole.elec ,
				 gv.TotalEden*gv.lgGrainElectrons,
				 dense.EdenExtra ,
				 rfield.otslin[iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont()-1] );

			if( trace.lgNeBug )
			{
				for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
				{
					if( nelem == 0 )
						fprintf( ioQQQ, "      eden_sum H -Ne:" );
					else if( nelem == 10 )
						fprintf( ioQQQ, "      eden_sum Na-Ca:" );
					else if( nelem == 20 )
						fprintf( ioQQQ, "      eden_sum Sc-Zn:" );
					fprintf( ioQQQ, " %.4e", eden_ions[nelem] );
					if( nelem%10 == 9 )
						fprintf( ioQQQ, "\n" );
				}
			}
		}
	}

	/* case where electron density is set with set eden command, make sure we use it */
	ASSERT( dense.EdenSet <= 0.f || fp_equal((realnum)dense.EdenTrue, dense.EdenSet) );
}
