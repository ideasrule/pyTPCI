/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_recomb generate recombination coefficients for any species */
/*ion_recombAGN generate recombination coefficients for AGN table */
#include "cddefines.h"
#include "phycon.h"
#include "heavy.h"
#include "grainvar.h"
#include "conv.h"
#include "thermal.h"
#include "iso.h"
#include "abund.h"
#include "save.h"
#include "elementnames.h"
#include "atmdat.h"
#include "ionbal.h"
#include "dense.h"

/*ion_recomb generate recombination coefficients for any species */
void ion_recomb(
  /* this is debug flag */
  bool lgPrintIt,
  /* nelem is the atomic number on the C scale, 0 for H */
  long int nelem/*, 
  double tlow[]*/)
{
	long int i, 
	  ion, 
	  limit;

	DEBUG_ENTRY( "ion_recomb()" );

	ASSERT( nelem < LIMELM);
	ASSERT( nelem > 1 );

	/* check that range of ionization is correct */
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonLow[nelem] <= nelem+1 );

	atmdat.nsbig = MAX2(dense.IonHigh[nelem]+1,atmdat.nsbig);

	/* this routine only does simple two-level species, 
	 * loop over ions will be <= limit, IonHigh is -1 since very
	 * highest stage of ionization is not recombined into.  
	 * for Li, will do only atom, since ions are H and He like,
	 * so limit is zero */

	/* zero-out loop comes before main loop since there are off-diagonal
	 * elements in the main ionization loop, due to multi-electron processes */
	/* >>chng 00 dec 07, limit changed to identical to ion_solver */
	for( ion=0; ion <= dense.IonHigh[nelem]-1; ion++ )
	{
		ionbal.RateRecomTot[nelem][ion] = 0.;
	}

	for( long ipISO=ipH_LIKE; ipISO<MIN2(NISO,nelem+1); ipISO++ )
	{
		if( ((dense.IonHigh[nelem] >= nelem - ipISO) &&
			  (dense.IonLow[nelem] <= nelem - ipISO)) )
		{
			ionbal.RateRecomTot[nelem][nelem-ipISO] = ionbal.RateRecomIso[nelem][ipISO];
		}
	}

	limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);
	ASSERT( limit >= -1 );

	/* these are counted elsewhere and must not be added here */
	Heavy.xLyaHeavy[nelem][nelem] = 0.;
	Heavy.xLyaHeavy[nelem][nelem-1] = 0.;

	/* IonLow is 0 for the atom, limit chosen to NOT include iso sequences  */
	for( ion=dense.IonLow[nelem]; ion <= limit; ion++ )
	{
		/* number of bound electrons of the ion after recombination,
		 * for an atom (ion=0) this is
		 * equal to nelem+1, the element on the physical scale, since nelem is
		 * on the C scale, being 5 for carbon */
		ASSERT(ionbal.DR_Badnell_rate_coef[nelem][ion] >= 0);
		ASSERT(ionbal.RR_rate_coef_used[nelem][ion] >= 0);

		/* sum of recombination rates [units s-1] for radiative, three body, charge transfer */
		ionbal.RateRecomTot[nelem][ion] = 
			dense.eden* (
			ionbal.RR_rate_coef_used[nelem][ion] + 
			ionbal.DR_Badnell_rate_coef[nelem][ion] +
			ionbal.CotaRate[ion] );

		/* >>chng 01 jun 30, FRAC_LINE was 0.1, not 1, did not include anything except
		 * radiative recombination, the radrec term */
		static const double FRAC_LINE = 1.;
		/* was 0.1 */
		/*Heavy.xLyaHeavy[nelem][ion] = (realnum)(dense.eden*radrec*FRAC_LINE );*/
		Heavy.xLyaHeavy[nelem][ion] = (realnum)(dense.eden*
			(ionbal.RR_rate_coef_used[nelem][ion]+ionbal.DR_Badnell_rate_coef[nelem][ion])*FRAC_LINE );
	}

	/* option to save rec coefficients */
	if( save.lgioRecom || lgPrintIt )
	{
		/* >>chng 04 feb 22, make option to print ions for single element */
		FILE *ioOut;
		if( lgPrintIt )
			ioOut = ioQQQ;
		else
			ioOut = save.ioRecom;

		/* print name of element */
		fprintf( ioOut, 
			" %s recombination coefficients fnzone:%.2f \tte\t%.4e\tne\t%.4e\n", 
			elementnames.chElementName[nelem] , fnzone , phycon.te , dense.eden );

		/*limit = MIN2(11,dense.IonHigh[nelem]);*/
		/* >>chng 05 sep 24, just print one long line - need info */
		limit = dense.IonHigh[nelem];
		// give ion stage
		for( i=0; i<limit; ++i )
			fprintf( ioOut, "%10ld",i+1);
		fprintf( ioOut, "\n");
		
		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RR_rate_coef_used[nelem][i] );
		}
		fprintf( ioOut, " radiative used vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RR_Verner_rate_coef[nelem][i] );
		}
		fprintf( ioOut, " old Verner vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RR_Badnell_rate_coef[nelem][i] );
		}
		fprintf( ioOut, " new Badnell vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			/* >>chng 06 jan 19, from div by eden to div by H0 - want units of cm3 s-1 but
			 * no single collider does this so not possible to get rate coefficient easily
			 * H0 is more appropriate than electron density */
			fprintf( ioOut, "%10.2e", ionbal.CX_recomb_rate_used[nelem][i]/SDIV(dense.xIonDense[ipHYDROGEN][0]) );
		}
		fprintf( ioOut, " CT/n(H0)\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.CotaRate[ion] );
		}
		fprintf( ioOut, " 3body vs Z /ne\n" );

		/* note different upper limit - this routine does grain rec for all ions */
		for( i=0; i < dense.IonHigh[nelem]; i++ )
		{
			fprintf( ioOut, "%10.2e", gv.GrainChTrRate[nelem][i+1][i]/dense.eden );
		}
		fprintf( ioOut, " Grain vs Z /ne\n" );
		fprintf( ioOut, " old Nussbaumer Storey DR vs Z\n" );

		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.DR_Badnell_rate_coef[nelem][i] );
		}
		fprintf( ioOut, " new Badnell DR vs Z\n" );

		/* total recombination rate, with density included - this goes into the matrix */
		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RateRecomTot[nelem][i] );
		}
		fprintf( ioOut, 
			" total rec rate (with density) for %s\n", 
			elementnames.chElementSym[nelem] );
		for( i=0; i < limit; i++ )
		{
			fprintf( ioOut, "%10.2e", ionbal.RateRecomTot[nelem][i]/dense.eden );
		}
		fprintf( ioOut, 
			" total rec rate / ne for %s\n\n", 
			elementnames.chElementSym[nelem] );

		/* spill over to next line for many stages of ionization */
		if( dense.IonHigh[nelem] > 11 )
		{
			limit = MIN2(29,dense.IonHigh[nelem]);
			fprintf( ioOut, " R " );
			for( i=11; i < limit; i++ )
			{
				fprintf( ioOut, "%10.2e", dense.eden*ionbal.CotaRate[ion] );
			}
			fprintf( ioOut, "\n" );

			fprintf( ioOut, "   " );
			for( i=11; i < limit; i++ )
			{
				fprintf( ioOut, "%10.2e", ionbal.RateRecomTot[nelem][i] );
			}
			fprintf( ioOut, "\n\n" );
		}
	}

	/* >>chng 02 nov 09, from -2 to -NISO */
	/*limit = MIN2(nelem-2,dense.IonHigh[nelem]-1);*/
	limit = MIN2(nelem-NISO,dense.IonHigh[nelem]-1);
	for( i=dense.IonLow[nelem]; i <= limit; i++ )
	{
		ASSERT( Heavy.xLyaHeavy[nelem][i] > 0. );
		ASSERT( ionbal.RateRecomTot[nelem][i] > 0. );
	}
	return;
}

/*ion_recombAGN generate recombination coefficients for AGN table */
void ion_recombAGN( FILE * io )
{
	static const int N1LIM = 3;
	static const int N2LIM = 4;
	double te1[N1LIM]={ 5000., 10000., 20000.};
	double te2[N2LIM]={ 20000.,50000.,100000.,1e6};
	/* this is boundary between two tables */
	double BreakEnergy = 100./13.0;
	long int nelem, ion , i;
	/* this will hold element symbol + ionization */
	char chString[100],
		chOutput[100];
	/* save temp here	*/
	double TempSave = phycon.te;
	/* save ne here	*/
	double EdenSave = dense.eden;

	DEBUG_ENTRY( "ion_recombAGN()" );

	EdenChange( 1. );
	/*atmdat_readin();*/

	/* first put header on file */
	fprintf(io,"X+i\\Te");
	for( i=0; i<N1LIM; ++i )
	{
		phycon.te = te1[i];
		fprintf(io,"\t%.0f K",phycon.te);
	}
	fprintf(io,"\n");

	/* now do loop over temp, but add elements */
	for( nelem=ipLITHIUM; nelem<LIMELM; ++nelem )
	{
		/* this list of elements included in the AGN tables is defined in zeroabun.c */
		if( dense.lgElmtOn[nelem] && abund.lgAGN[nelem] )
		{
			for( ion=0; ion<=nelem; ++ion )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][ion] > 0.05 );

				if( Heavy.Valence_IP_Ryd[nelem][ion] > BreakEnergy )
					break;

				/* print chemical symbol */
				sprintf(chOutput,"%s", 
					elementnames.chElementSym[nelem]);
				/* some elements have only one letter - this avoids leaving a space */
				if( chOutput[1]==' ' )
					chOutput[1] = chOutput[2];
				/* now ionization stage */
				if( ion==0 )
				{
					sprintf(chString,"0 ");
				}
				else if( ion==1 )
				{
					sprintf(chString,"+ ");
				}
				else
				{
					sprintf(chString,"+%li ",ion);
				}
				strcat( chOutput , chString );
				fprintf(io,"%5s",chOutput );

				for( i=0; i<N1LIM; ++i )
				{
					TempChange(te1[i] , false);
					dense.IonLow[nelem] = 0;
					dense.IonHigh[nelem] = nelem+1;
					ConvBase(0);
					fprintf(io,"\t%.2e",ionbal.RateRecomTot[nelem][ion]);
				}
				fprintf(io,"\n");
			}
			fprintf(io,"\n");
		}
	}

	/* second put header on file */
	fprintf(io,"X+i\\Te");
	for( i=0; i<N2LIM; ++i )
	{
		TempChange(te2[i] , false);
		fprintf(io,"\t%.0f K",phycon.te);
	}
	fprintf(io,"\n");

	/* now do same loop over temp, but add elements */
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		/* this list of elements included in the AGN tables is defined in zeroabun.c */
		if( dense.lgElmtOn[nelem] && abund.lgAGN[nelem] )
		{
			for( ion=0; ion<=nelem; ++ion )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][ion] > 0.05 );

				if( Heavy.Valence_IP_Ryd[nelem][ion] <= BreakEnergy )
					continue;

				/* print chemical symbol */
				fprintf(io,"%s", 
					elementnames.chElementSym[nelem]);
				/* now ionization stage */
				if( ion==0 )
				{
					fprintf(io,"0 ");
				}
				else if( ion==1 )
				{
					fprintf(io,"+ ");
				}
				else
				{
					fprintf(io,"+%li",ion);
				}

				for( i=0; i<N2LIM; ++i )
				{
					TempChange(te2[i] , false);
					dense.IonLow[nelem] = 0;
					dense.IonHigh[nelem] = nelem+1;
					ConvBase(0);
					fprintf(io,"\t%.2e",ionbal.RateRecomTot[nelem][ion]);
				}
				fprintf(io,"\n");
			}
			fprintf(io,"\n");
		}
	}

	TempChange(TempSave , true);
	EdenChange( EdenSave );
	return;
}
