/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*EdenChange - update electron density and its dependent quantities */
#include "cddefines.h"
#include "dense.h"
#include "rfield.h"
#include "phycon.h"
#include "conv.h"
#include "rt.h"

void EdenChange( double EdenNew )
{
	static double EdenOld=-1;

	DEBUG_ENTRY( "EdenChange()" );

	// don't let electron density be zero - logs are taken
	ASSERT( EdenNew > 0. );

	// init EdenOld on first sweep through this sim
	if( conv.nTotalIoniz==0 )
		EdenOld = dense.eden;

	// this confirms that eden is only changed in this routine
	if( !fp_equal( dense.eden, EdenOld ) )
		TotalInsanity();

	dense.eden = EdenNew;
	EdenOld = EdenNew;

	dense.EdenHCorr = dense.eden +
		/* dense.HCorrFac is unity by default and changed with the set HCOR command */
		dense.xIonDense[ipHYDROGEN][0]*1.7e-4 * dense.HCorrFac;
	dense.EdenHCorr_f = (realnum)dense.EdenHCorr;
	dense.edensqte = dense.EdenHCorr/phycon.sqrte;
	dense.cdsqte = dense.edensqte*COLL_CONST;
	dense.SqrtEden = sqrt(dense.eden);

	/* evaluate the plasma frequency one time per zone to avoid PF moving across
	 * a line during convergence loops */
	if( nzone != rfield.nZonePlsFrqEval || conv.lgSearch )
	{
		rfield.nZonePlsFrqEval = nzone;
		rfield.plsfrq = (realnum)((ELEM_CHARGE_ESU/sqrt(PI*ELECTRON_MASS)/FR1RYD)*sqrt(dense.eden));

		if( rfield.ipPlasma > 0 )
		{
			/* increase index for plasma frequency until within proper cell */
			while( rfield.plsfrq > rfield.anumax(rfield.ipPlasma) )
				++rfield.ipPlasma;

			/* decrease index for plasma frequency until within proper cell */
			while( rfield.ipPlasma > 2 && rfield.plsfrq < rfield.anumin(rfield.ipPlasma) )
				--rfield.ipPlasma;
		}

		/* also remember the largest plasma frequency we encounter */
		rfield.plsfrqmax = MAX2(rfield.plsfrqmax, rfield.plsfrq);

		/* is plasma frequency within energy grid? */
		if( rfield.plsfrq > rfield.anu(0) )
		{
			rfield.lgPlasNu = true;
		}
	}

	// if plasma frequency has changed we need to update transitions - those
	// below plasma frequency do not exist - emission rate set to smallfloat
	// only do this in search phase since plasma frequency moving across transition
	// would create discontinuous changes in cooling that would present noise
	// to the solver.  The electron density should not change by much during the
	// solution for a zone
	static double EdenEval=-1;
	if( conv.lgSearch && !fp_equal(EdenEval,dense.eden) )
	{
		EdenEval = dense.eden;
		RT_line_all_escape ( NULL );
	}

	return;
}
