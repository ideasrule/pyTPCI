/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtAllTau master routine controlling printout of optical depths at end of calculation */
#include "cddefines.h"
#include "iso.h"
#include "opacity.h"
#include "dense.h"
#include "elementnames.h"
#include "geometry.h"
#include "prt.h"
#include "mole.h"
#include "freebound.h"

/*PrtAllTau master routine controlling printout of optical depths at end of calculation */
void PrtAllTau(void)
{
	long int i, 
		n,
	  nelem;
	realnum fcon, 
	  flin;

	DEBUG_ENTRY( "PrtAllTau()" );

	/* optical depths used by code are total through model,
	 * when sphere is set, this is twice optical depth through
	 * computed structure */
	if( geometry.lgSphere )
	{
		fcon = 2.;
		if( geometry.lgStatic )
		{
			flin = 2.;
		}
		else
		{
			flin = 1.;
		}
	}
	else
	{
		fcon = 1.;
		flin = 1.;
	}

	/* continuum optical depths and column densities */
	fprintf( ioQQQ, "\n Contin Optical Depths: COMP:");

	fprintf( ioQQQ,PrintEfmt("%9.2e", opac.telec));
	fprintf( ioQQQ, "    H-:");
	fprintf( ioQQQ,PrintEfmt("%9.2e",opac.thmin ));

	/* R(1300) is Rayleigh scattering */
	fprintf( ioQQQ, " R(1300):");
	fprintf( ioQQQ,PrintEfmt("%9.2e", findspecieslocal("H")->column*6.71e-24));

	fprintf( ioQQQ, "  H2+:");
	fprintf( ioQQQ,PrintEfmt("%9.2e", findspecieslocal("H2+")->column*7e-18));

	fprintf( ioQQQ, "  Bra:");
	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_local >= 5 )
	{
		long ip5p = iso_sp[ipH_LIKE][ipHYDROGEN].QN2Index(5, 1, 2);
		ASSERT( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip5p,ipH4s).ipCont() > 0 );
		PrintE82( ioQQQ , opac.TauTotalGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip5p,ipH4s).ipCont()-1]/fcon);
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}
	fprintf( ioQQQ, "\n" );

	fprintf( ioQQQ, "                          Pa:");

	if( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_local > ipH4p )
	{
		fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH3s).ipCont()-1]/fcon));
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}

	fprintf( ioQQQ, "    Ba:");

	if( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_local > 3 )
	{
		fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH3p,ipH2s).ipCont()-1]/fcon));
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}

	fprintf( ioQQQ, "      Hb:");

	if( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_local > 4 )
	{
		fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH2s).ipCont()-1]/fcon));
	}
	else
	{
		PrintE82( ioQQQ , 0.);
	}

	fprintf( ioQQQ, "   La:");
	fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauTotalGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont()-1]/fcon));

	fprintf( ioQQQ, "     1r:");
	PrintE93( ioQQQ , opac.TauTotalGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1]/fcon);

	if( dense.lgElmtOn[ipHELIUM] )
	{
		fprintf( ioQQQ, "  1.8:");
		PrintE82( ioQQQ , opac.TauTotalGeo[0][iso_sp[ipHE_LIKE][ipHELIUM].fb[ipHe1s1S].ipIsoLevNIonCon-1]/fcon);

		fprintf( ioQQQ, " 4.:");
		PrintE93( ioQQQ , opac.TauTotalGeo[0][iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon-1]/fcon);
	}
	fprintf( ioQQQ, "\n");

	/* print optical depths of some metal lines */
	prtmet();

	/* print H-like H, He+ optical depths */
	for( nelem=ipHYDROGEN; nelem<=ipHELIUM; ++nelem )
	{
		/* helium may be turned off */
		if( dense.lgElmtOn[nelem] )
		{
			const int NUMB_PER_LINE = 8;
			fprintf( ioQQQ, "\n Old, new %s%2li continuum optical depths:\n",
				elementnames.chElementSym[nelem] ,
				nelem+1);
			/* absorption continuum optical depths are energies of the h-like ion continua
			 * loop over old, then new optical depths */
			for( i=1; i>=0; --i )
			{
				/* print ground, skip t2, then do 2p */
				for( n=ipH1s; n < iso_sp[ipH_LIKE][nelem].numLevels_max - iso_sp[ipH_LIKE][nelem].nCollapsed_max; n++ )
				{
					if( n==ipH2s )
						continue;
					if( n%NUMB_PER_LINE ==1)
						fprintf(ioQQQ,"\n");
					/* this, combined with "continue" above, ensures that we print
					 * 1 (1s), 2(tot 2), then 3 */
					fprintf( ioQQQ , "%6ld",MAX2(1,n));
					fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauAbsGeo[i][iso_sp[ipH_LIKE][nelem].fb[n].ipIsoLevNIonCon-1]/fcon));
				}
				fprintf( ioQQQ, "\n" );
			}

			/* now do h-like line optical depths */
			fprintf( ioQQQ, "\n Old, new %s%2li mean line optical depths:\n",
				elementnames.chElementSym[nelem] ,
				nelem+1);
			/* Lya is a special case due to 2s-2p resolution - explicitly print it first */
			fprintf( ioQQQ, "%3i-%2i",2, 1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot()*SQRTPI/flin ));
			/* total optical depth in 3-2s and 3-2p, is total of both so 2-1 is correct for 3-2*/
			/* 06 aug 28, from numLevels_max to _local. */
			for( n=3; n <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_local; n++ )
			{
				if( n%NUMB_PER_LINE ==1)
					fprintf(ioQQQ,"\n");
				fprintf( ioQQQ, "%3ld-%2ld",n, n-1 );
				fprintf( ioQQQ,PrintEfmt("%9.2e", 
					/* just do nP - n'S,
					 *>>chng 12 jul25, from line center to mean line optical depths */
							 iso_sp[ipH_LIKE][nelem].trans( iso_sp[ipH_LIKE][nelem].QN2Index(n, 1, 2),
								iso_sp[ipH_LIKE][nelem].QN2Index(n-1, 0, 2) ).Emis().TauTot()*SQRTPI/flin ));
			}
			if( prt.lgPrnIsoCollapsed )
			{
				/* above flag set with print line iso collapsed command
				 * not done by default */
				for( n=iso_sp[ipH_LIKE][nelem].numLevels_local - iso_sp[ipH_LIKE][nelem].nCollapsed_local; n < iso_sp[ipH_LIKE][nelem].numLevels_local; n++ )
				{
					if( iso_sp[ipH_LIKE][nelem].st[n].n() % NUMB_PER_LINE ==1)
						fprintf(ioQQQ,"\n");
					fprintf( ioQQQ, "%3ld-%2ld", iso_sp[ipH_LIKE][nelem].st[n].n(), iso_sp[ipH_LIKE][nelem].st[n-1].n() );
					fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipH_LIKE][nelem].trans(n,n-1).Emis().TauTot()*SQRTPI/flin ));
				}
			}

			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "%3i-%2i",2, 1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn()*SQRTPI/flin ));

			for( n=3; n <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_local; n++ )
			{
				if( n%NUMB_PER_LINE ==1)
					fprintf(ioQQQ,"\n");
				fprintf( ioQQQ, "%3ld-%2ld",n, n-1 );
				fprintf( ioQQQ,PrintEfmt("%9.2e", 
					/* just do nP - n'S */
							 iso_sp[ipH_LIKE][nelem].trans( iso_sp[ipH_LIKE][nelem].QN2Index(n, 1, 2),
								iso_sp[ipH_LIKE][nelem].QN2Index(n-1, 0, 2) ).Emis().TauIn()*SQRTPI/flin ));
			}
			if( prt.lgPrnIsoCollapsed )
			{
				for( n=iso_sp[ipH_LIKE][nelem].numLevels_local - iso_sp[ipH_LIKE][nelem].nCollapsed_local; n < iso_sp[ipH_LIKE][nelem].numLevels_local; n++ )
				{
					if( iso_sp[ipH_LIKE][nelem].st[n].n() % NUMB_PER_LINE ==1)
						fprintf(ioQQQ,"\n");
					fprintf( ioQQQ, "%3ld-%2ld", iso_sp[ipH_LIKE][nelem].st[n].n(), iso_sp[ipH_LIKE][nelem].st[n-1].n() );
					fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipH_LIKE][nelem].trans(n,n-1).Emis().TauIn()*SQRTPI/flin ));
				}
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* ================================================================================ */

	/* print He I lines if helium exists */
	if( dense.lgElmtOn[ipHELIUM] )
	{
		fprintf( ioQQQ, "\n Old He Is optical depths:" );
		for( i=0; i < 5; i++ )
		{
			fprintf( ioQQQ, "%5ld", i+1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauAbsGeo[1][iso_sp[ipHE_LIKE][ipHELIUM].fb[i].ipIsoLevNIonCon-1]/fcon) );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " New He Is optical depths:" );
		for( i=0; i < 5; i++ )
		{
			fprintf( ioQQQ, "%5ld", i+1 );
			fprintf( ioQQQ,PrintEfmt("%9.2e", opac.TauAbsGeo[0][iso_sp[ipHE_LIKE][ipHELIUM].fb[i].ipIsoLevNIonCon-1]/fcon ));
		}
		fprintf( ioQQQ, "\n" );

		/* ================================================================================*/

		fprintf( ioQQQ, "          Old He Is Lines:" );
		fprintf( ioQQQ, " %4d",584 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe2p1P,ipHe1s1S).Emis().TauTot()*SQRTPI/flin ));
		fprintf( ioQQQ, " %4d",3889 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe3p3P,ipHe2s3S).Emis().TauTot()*SQRTPI/flin ));
		fprintf( ioQQQ, " %4d",5016 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe3p1P,ipHe2s1S).Emis().TauTot()*SQRTPI/flin ));
		fprintf( ioQQQ, " %4d",5876 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe3d3D,ipHe2p3P2).Emis().TauTot()*SQRTPI/flin ));
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "          New He Is Lines:" );
		fprintf( ioQQQ, " %4d",584 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe2p1P,ipHe1s1S).Emis().TauIn()*SQRTPI/flin ));
		fprintf( ioQQQ, " %4d",3889 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe3p3P,ipHe2s3S).Emis().TauIn()*SQRTPI/flin ));
		fprintf( ioQQQ, " %4d",5016 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe3p1P,ipHe2s1S).Emis().TauIn()*SQRTPI/flin ));
		fprintf( ioQQQ, " %4d",5876 );
		fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHe3d3D,ipHe2p3P2).Emis().TauIn()*SQRTPI/flin ));
		fprintf( ioQQQ, "\n" );

		/* ================================================================================*/
	}
	return;
}
