/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HydroLevel calls iso_level to solve for ionization balance 
 * and level populations of model hydrogen atom */
/*PrtHydroTrace1 print trace info for hydrogen-like species */
#include "cddefines.h"
#include "iso.h"
#include "secondaries.h"
#include "trace.h"
#include "phycon.h"
#include "ionbal.h"
#include "hydrogenic.h"
#include "freebound.h"
#include "two_photon.h"
#include "dense.h"

/*PrtHydroTrace1a print trace info for hydrogen-like species */
STATIC void PrtHydroTrace1a( long ipISO, long nelem )
{
	double colfrc, 
	  phtfrc, 
	  secfrc;

	DEBUG_ENTRY( "PrtHydroTrace1a()" );

	/* identify how atom is ionized for full trace */
	if( iso_sp[ipISO][nelem].xIonSimple > 0. )
	{
		/* fraction of ionization due to photoionization */
		phtfrc = iso_sp[ipISO][nelem].fb[ipH1s].gamnc/((dense.eden*(iso_sp[ipISO][nelem].RadRec_effec + 
			ionbal.CotaRate[nelem]) )*
			iso_sp[ipISO][nelem].xIonSimple);

		/* fraction of ionization due to collisional ionization */
		colfrc = (iso_sp[ipISO][nelem].fb[ipH1s].ColIoniz*dense.EdenHCorr )/
			((dense.eden*(iso_sp[ipISO][nelem].RadRec_effec + 
			ionbal.CotaRate[0]) )*
			iso_sp[ipISO][nelem].xIonSimple);

		/* fraction of ionization due to secondary ionization */
		secfrc = secondaries.csupra[nelem][nelem]/((dense.eden*(iso_sp[ipISO][nelem].RadRec_effec + 
			ionbal.CotaRate[0]) )*
			iso_sp[ipISO][nelem].xIonSimple);
	}
	else
	{
		phtfrc = 0.;
		colfrc = 0.;
		secfrc = 0.;
	}

	fprintf( ioQQQ, "     HydroLevel Z=%2ld called, simple II/I=",nelem);
	PrintE93( ioQQQ, iso_sp[ipISO][nelem].xIonSimple);
	fprintf( ioQQQ," PhotFrc:");
	PrintE93( ioQQQ,phtfrc);
	fprintf(ioQQQ," ColFrc:");
	PrintE93( ioQQQ,colfrc);
	fprintf( ioQQQ," SecFrc");
	PrintE93( ioQQQ, secfrc);
	fprintf( ioQQQ,"  Te:");
	PrintE93( ioQQQ,phycon.te);
	fprintf( ioQQQ," eden:");
	PrintE93( ioQQQ,dense.eden);
	fprintf( ioQQQ,"\n"); 
	return;
}

/*PrtHydroTrace1 print trace info for hydrogen-like species */
STATIC void PrtHydroTrace1( long ipISO, long nelem )
{
	long int ipHi , ipLo , i;

	DEBUG_ENTRY( "PrtHydroTrace1()" );

	fprintf( ioQQQ, 
		"       HydroLevel%3ld finds arrays, with optical depths defined? %li induced 2ph=%12.3e\n", 
		nelem, iteration, iso_sp[ipISO][nelem].TwoNu[0].induc_up );
	/* 06 aug 28, from numLevels_max to _local. */
	for( ipHi=ipH2s; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
	{
		fprintf( ioQQQ, "up:%2ld", ipHi );
		fprintf( ioQQQ, "lo" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ, "%9ld", ipLo );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " A*esc" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul()*
				iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pesc() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " A*ees" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul()*
				iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pelec_esc() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " tauin" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " t tot" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauTot() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " Esc  " );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pesc() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " Eesc " );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pelec_esc() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " Dest " );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pdest()) );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " A*dst" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul()*
				iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pdest() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " StrkE" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipISO][nelem].ex[ipHi][ipLo].pestrk_up ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " B(ul)" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().pump()*
				iso_sp[ipISO][nelem].st[ipLo].g()/iso_sp[ipISO][nelem].st[ipHi].g() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " tcont" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauCon() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "%3ld", ipHi );
		fprintf( ioQQQ, " C(ul)" );
		for( ipLo=ipH1s; ipLo < ipHi; ipLo++ )
		{
			fprintf( ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].trans(ipHi,ipLo).Coll().ColUL( colliders ) ));
		}
		fprintf( ioQQQ, "\n" );

		if( ipHi == 2 && ipISO==ipH_LIKE && nelem==ipHYDROGEN )
		{
			fprintf( ioQQQ, "    FeIIo");
			fprintf( ioQQQ,PrintEfmt("%9.2e", 
				hydro.dstfe2lya* iso_sp[ipISO][nelem].trans(ipH2p,ipH1s).Emis().Aul() ));
			fprintf( ioQQQ, "\n");
		}
	}

	fprintf( ioQQQ, "         " );
	/* 06 aug 28, from numLevels_max to _local. */
	for( i=1; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
	{
		fprintf( ioQQQ, "%9ld", i );
	}
	fprintf( ioQQQ, "\n" );
	return;
}

/*HydroLevel calls iso_level to solve for ionization balance 
 * and level populations of model hydrogen atom */
void HydroLevel(long ipISO, long nelem)
{
	long int i; 

	DEBUG_ENTRY( "HydroLevel()" );

	/* check that we were called with valid charge */
	ASSERT( nelem >= 0);
	ASSERT( nelem < LIMELM );

	/* option to print some rates */
	if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO]) && (nelem == trace.ipIsoTrace[ipISO]) )
		PrtHydroTrace1(ipISO, nelem);

	if( trace.lgHBug && trace.lgTrace )
		PrtHydroTrace1a(ipISO, nelem);

	/* this is main trace h-like printout */
	if( (trace.lgIsoTraceFull[ipISO] && trace.lgTrace) && (nelem == trace.ipIsoTrace[ipISO]) )
	{
		fprintf( ioQQQ, "       HLEV HGAMNC" );
		PrintE93( ioQQQ, iso_sp[ipISO][nelem].fb[ipH1s].gamnc );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH2s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].fb[i].gamnc ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV TOTCAP" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=1; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].fb[i].RateCont2Level/dense.eden ));
		}
		fprintf( ioQQQ," tot");
		fprintf( ioQQQ,PrintEfmt("%10.2e", ionbal.RateRecomTot[nelem][nelem-ipISO]/dense.eden ) );
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV IND Rc" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].fb[i].RecomInducRate*iso_sp[ipISO][nelem].fb[i].PopLTE ));
		}
		fprintf( ioQQQ, "\n" );

		/* incuded recombination rate coefficients */
		fprintf( ioQQQ, "       IND Rc LTE " );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e",
				iso_sp[ipISO][nelem].fb[i].gamnc*iso_sp[ipISO][nelem].fb[i].PopLTE ));
		}
		fprintf( ioQQQ, "\n" );

		/* LTE level populations */
		fprintf( ioQQQ, "       HLEV   HLTE" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].fb[i].PopLTE ));
		}
		fprintf( ioQQQ, "\n" );

		/* fraction of total ionization due to collisions from given level */
		fprintf( ioQQQ, "       HLEVfr cion" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", 
				iso_sp[ipISO][nelem].fb[i].ColIoniz*dense.EdenHCorr/MAX2(1e-30,iso_sp[ipISO][nelem].fb[i].RateLevel2Cont) ) );
		}
		fprintf( ioQQQ, "\n" );

		/* fraction of total ionization due to photoionization from given level */
		if( ionbal.RateRecomTot[nelem][nelem]> 0. )
		{
			fprintf( ioQQQ, "       HLEVfrPhIon" );
			/* 06 aug 28, from numLevels_max to _local. */
			for( i=ipH1s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
			{
				fprintf(ioQQQ,PrintEfmt("%9.2e", 
					iso_sp[ipISO][nelem].fb[i].gamnc/MAX2(1e-30,iso_sp[ipISO][nelem].fb[i].RateLevel2Cont) ) );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, "       HLEV     HN" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].st[i].Pop() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV   b(n)" );
		/* 06 aug 28, from numLevels_max to _local. */
		for( i=ipH1s; i < iso_sp[ipISO][nelem].numLevels_local; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].st[i].DepartCoef() ));
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "       HLEV   X12tot");
		fprintf(ioQQQ,PrintEfmt("%9.2e", secondaries.x12tot ));
		fprintf( ioQQQ," Grn dest:");
		fprintf(ioQQQ,PrintEfmt("%9.2e",
		  ionbal.RateIoniz[nelem][nelem][nelem+1] ));
		fprintf(ioQQQ, "\n"); 
	}

	if( trace.lgTrace )
	{
		/* iso.RecomTotal[nelem] is gross rec coef, computed here while filling in matrix
		 * elements, all physical processes included. 
		 * RadRec_effec is total effective radiative only */
		fprintf( ioQQQ, "       HydroLevel Z:%2ld return %s te=",
			nelem,
			iso_sp[ipISO][nelem].chTypeAtomUsed );
		PrintE93( ioQQQ,phycon.te);
		fprintf( ioQQQ," density=%.4e", dense.xIonDense[nelem][nelem-ipISO] );

		fprintf( ioQQQ," simple=%.4e",iso_sp[ipISO][nelem].xIonSimple);

		fprintf( ioQQQ," b1=%.2e",iso_sp[ipISO][nelem].st[ipH1s].DepartCoef());

		fprintf( ioQQQ," ion rate=%.4e",ionbal.RateIonizTot(nelem,nelem-ipISO) );

		fprintf( ioQQQ," TotRec=%.4e",ionbal.RateRecomTot[nelem][nelem-ipISO]);

		fprintf( ioQQQ," RadRec=%.4e",iso_sp[ipISO][nelem].RadRec_effec);
		fprintf( ioQQQ, "\n");
	}
	return;
}
