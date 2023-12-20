/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveSpecial generate output for the save special command */
#include "cddefines.h"
#include "save.h"
#include "wind.h"
#include "opacity.h"
#include "dense.h"
#include "radius.h"

/*SaveSpecial generate output for the save special command */
void SaveSpecial(FILE* ioPUN , 
  const char *chTime)
{
	/*long int i;*/

	DEBUG_ENTRY( "SaveSpecial()" );

	if( strncmp(chTime,"LAST",4) == 0 )
	{
		/* code to execute only after last zone */
#		if 0
		long ipISO , nelem , limit , i;
		double EdenAbund , fach;
#		include "physconst.h"
#		include "hydrogenic.h"
		ipISO = ipHYDROGEN;
		nelem = ipHYDROGEN;

		/* in all following the factor of two is because a single
		 * decay produces two photons */
		EdenAbund = iso_sp[ipH_LIKE][nelem].st[ipH2s].Pop*8.226*powi(1.+nelem,6);
		fprintf(ioPUN," 2s = %.3e\n", EdenAbund);

		/* upper limit to H-like 2-phot is energy of La, which is in ipCont-1 cell */
		limit = iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).ipCont()-1;
		/* remember sum of rates, this will add up to twice the real rate since
		 * each transition makes two photons */
		for( i=0; i < limit; i++ )
		{
			/*>>chng 01 jan 23, previous change had doubled cross section for H two-photon,
			 * so here we divide by 2 to get old answer */
			/** \todo	2	this most likely needs to be changed in light of new 2nu treatment	*/
			fach = iso_sp[ipISO][nelem].TwoNu[0].As2nu[i]/2.f;
			fach *= rfield.anu2(i)/rfield.widflx(i)*EN1RYD;
			fprintf(ioPUN,"%.3e\t%.3e\t%.3e\n", 
				RYDLAM/1e4/rfield.anu(i) , fach , fach*(realnum)EdenAbund );
		}
#		endif

	}
	else
	{
		/* code to do for every zone */
		fprintf(ioPUN,"%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
			radius.Radius ,
			wind.AccelCont ,
			wind.fmul ,
			opac.opacity_sct[1000],
			dense.eden , 
			dense.xMassDensity,
			dense.gas_phase[ipHYDROGEN] );
	}

	return;
}
