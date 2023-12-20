/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*AbundancesPrt print all abundances, both gas phase and grains */
/*AbundancesSet sets initial abundances after parameters are entered by reading input */
/*AbundancesTable interpolate on table of points to do 'element table' command, */
/*PrtElem print chemical composition at start of calculation */
#include "cddefines.h"
#include "phycon.h"
#include "called.h"
#include "stopcalc.h"
#include "thermal.h"
#include "trace.h"
#include "elementnames.h"
#include "dense.h"
#include "radius.h"
#include "grainvar.h"
#include "abund.h"
#include "deuterium.h"
#include "physconst.h"

/*PrtElem print chemical composition at start of calculation */
STATIC void PrtElem(
  /* the job to do, the options are "init", "fill", "flus" */
  const char *chJob, 
  /* label for the element */
  const char *chLabl, 
  /* its abundance */
  double abund_prt);

/*AbundancesPrt print all abundances, both gas phase and grains */
void AbundancesPrt( void )
{
	long int i;
	double GrainNumRelHydrSilicate ,
		GrainNumRelHydrCarbonaceous ,
		GrainNumRelHydr_PAH,
		GrainMassRelHydrSilicate,
		GrainMassRelHydrCarbonaceous,
		GrainMassRelHydr_PAH;

	DEBUG_ENTRY( "AbundancesPrt()" );

	/* this is main loop to print abundances of each element */
	if( called.lgTalk )
	{
		PrtElem("initG","  ",0.);/* initialize print routine for gas*/
		for( i=0; i < LIMELM; i++ )
		{
			if( dense.lgElmtOn[i] )
			{
				/* fill in print buffer with abundances */
				PrtElem("fill",(char*)elementnames.chElementSym[i],
				  abund.GasPhase[i]);
			}
		}

		/* flush the print buffer */
		PrtElem("flus","  ",0.);
		/* final carriage return */
		fprintf( ioQQQ, " \n" );

		/* now grains if present */
		if( gv.lgDustOn() )
		{
			/* we will first print the total abundances of each element locked up in grains */
			/* initialize print routine for dust*/
			PrtElem("initD","  ",0.);
			double sum_grains = 0.;
			for( i=0; i < LIMELM; i++ )
			{
				if( gv.elmSumAbund[i]>SMALLFLOAT )
				{
					/* fill in print buffer with abundances */
					PrtElem("fill",(char*)elementnames.chElementSym[i],
						gv.elmSumAbund[i]/dense.gas_phase[ipHYDROGEN]);
					sum_grains += gv.elmSumAbund[i]/dense.gas_phase[ipHYDROGEN];
				}
			}
			/* flush the print buffer */
			PrtElem("flus","  ",0.);

			double sum_depl = abund.SumDepletedAtoms();
			double offset = 0.;
			if( sum_depl > 0. )
			{
				offset = (sum_depl / sum_grains - 1.) * 100.;
				sum_depl = log10( sum_depl );
			}
			else
			{
				offset = - 100.;
				sum_depl = -30.;
			}
			fprintf(ioQQQ, "                                                   "
					"Sum Atoms Depleted: %.4f\n", sum_depl );
			fprintf(ioQQQ, "                                                   "
					"Depletion/Grains-1: %4.2f%%\n", offset );
			/* final carriage return */
			fprintf( ioQQQ, " \n" );

			/* this is used to store grain number density per hydrogen */
			GrainNumRelHydrSilicate = 0.;
			GrainNumRelHydrCarbonaceous = 0;
			GrainNumRelHydr_PAH = 0.;
			GrainMassRelHydrSilicate = 0.;
			GrainMassRelHydrCarbonaceous = 0;
			GrainMassRelHydr_PAH = 0.;

			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{

				/* number density of grains per hydrogen, the ratio
				 * gv.bin[nd].IntVol/gv.bin[nd].AvVol is the number of grain particles 
				 * per H at standard grain abundance*/
				realnum DensityNumberPerHydrogen = 
					(gv.bin[nd].IntVol/gv.bin[nd].AvVol)*gv.bin[nd].dstAbund / 
					gv.bin[nd].GrnDpth;
				/* mass of grains per hydrogen */
				realnum DensityMassPerHydrogen = 
					gv.bin[nd].IntVol*gv.bin[nd].dustp[0]*gv.bin[nd].dstAbund/
					(realnum)ATOMIC_MASS_UNIT / gv.bin[nd].GrnDpth;

				/* >>chng 06 mar 05, fix expression for calculating grain number density, PvH */
				if( gv.bin[nd].matType == MAT_CAR || gv.bin[nd].matType == MAT_CAR2 ||
				    gv.bin[nd].matType == MAT_SIC )
				{
					/* carbonaceous grains */
					GrainNumRelHydrCarbonaceous += DensityNumberPerHydrogen;
					GrainMassRelHydrCarbonaceous += DensityMassPerHydrogen;
				}
				else if( gv.bin[nd].matType == MAT_SIL || gv.bin[nd].matType == MAT_SIL2 )
				{
					/* silicate grains */
					GrainNumRelHydrSilicate += DensityNumberPerHydrogen;
					GrainMassRelHydrSilicate += DensityMassPerHydrogen;
				}
				else if( gv.bin[nd].matType == MAT_PAH || gv.bin[nd].matType == MAT_PAH2 )
				{
					/* PAHs - full abundance - remove possible factor accounting for
					 * variation of abundances with physical conditions - this will 
					 * be the PAH abundance with scale factor of unity */
					GrainNumRelHydr_PAH += DensityNumberPerHydrogen;
					GrainMassRelHydr_PAH += DensityMassPerHydrogen;
				}
				else
					TotalInsanity();
			}

			/* now print total number of grains of each type */
			fprintf(ioQQQ,"              Number of grains per hydrogen (scale=1)                         Mass of grains per hydrogen (scale=1)\n");
			fprintf(ioQQQ,"        Carbonaceous: %.3f  Silicate: %.3f  PAH: %.3f         Carbonaceous: %.3f  Silicate: %.3f  PAH: %.3f\n\n" ,
				log10( MAX2( 1e-30, GrainNumRelHydrCarbonaceous ) ) ,
				log10( MAX2( 1e-30, GrainNumRelHydrSilicate ) ) ,
				log10( MAX2( 1e-30, GrainNumRelHydr_PAH ) ) ,
				log10( MAX2( 1e-30, GrainMassRelHydrCarbonaceous ) ) ,
				log10( MAX2( 1e-30, GrainMassRelHydrSilicate ) ) ,
				log10( MAX2( 1e-30, GrainMassRelHydr_PAH ) )  );
		}
	}
	return;
}

/*AbundancesSet print all abundances, both gas phase and grains */
void AbundancesSet(void)
{
	long int i, 
	  nelem;
	double fac;
	static bool lgFirstCall=true;
	static bool lgElOnOff[LIMELM];

	DEBUG_ENTRY( "AbundancesSet()" );

	/* if this is the first call to this routine in this core load, 
	 * save the state of the lgElmOn array, so that it is possible
	 * to turn off elements in later models, but not turn on an
	 * element that was initially turned off.  This is necessary since
	 * the Create... routines that create space for elements will
	 * not be revisited in later models.  You can turn off an initially
	 * enabled element, but not turn a disabled one on.  */

	if( lgFirstCall )
	{
		/* first call - save the initial state of the lgElmtOn vector */
		for( i=0; i<LIMELM; ++i )
		{
			lgElOnOff[i] = dense.lgElmtOn[i];
		}
	}
	lgFirstCall = false;

	/* make sure that initially false elements remain off, while letting 
	 * enabled elements be turned off */
	for( i=ipHYDROGEN; i<LIMELM; ++i )
	{
		dense.lgElmtOn[i] = lgElOnOff[i] && dense.lgElmtOn[i];
	}

	/* rescale so that abundances are H=1 */
	for( i=ipHYDROGEN; i < LIMELM; i++ )
	{
		abund.GasPhase[i] = abund.ReferenceAbun[i] / abund.ReferenceAbun[ipHYDROGEN];
	}

	/* set current abundances to "solar" times metals scale factor
	 * and grain depletion factor */
	abund.GasPhase[ipHELIUM] *= abund.DepletionScaleFactor[ipHELIUM]*abund.ScaleElement[ipHELIUM];

	/* option for density or abundance variations, this flag is true by default,
	 * set in zero, but set false if variations are enabled AND these
	 * are not density variations, but rather abundances */
	if( dense.lgDenFlucOn )
	{
		/* usual case - either density fluctuations or none at all */
		fac = 1.;
	}
	else
	{
		/* abundance fluctuations enabled, set initial value */
		fac = dense.cfirst*cos(dense.flcPhase) + dense.csecnd;
	}

	for( i=ipLITHIUM; i < LIMELM; i++ )
	{
		abund.GasPhase[i] *= (realnum)(abund.ScaleMetals*abund.DepletionScaleFactor[i]*
					abund.ScaleElement[i]*fac);
	}

	/* now fix abundance of any element with element table set */
	if( abund.lgAbTaON )
	{
		for( nelem=ipHELIUM; nelem < LIMELM; ++nelem )
		{
			if( abund.lgAbunTabl[nelem] )
			{
				abund.GasPhase[nelem] = (realnum)(AbundancesTable(radius.Radius,
				  radius.depth,nelem+1));
			}
		}
	}

	/* dense.gas_phase[nelem] contains total abundance of element */
	/* the density of hydrogen itself has already been set at this point -
	 * it is set when commands parsed, most likely by the hden command -
	 * set all heavier elements */
	/* if abund.GasPhase[ipHYDROGEN] == 1, consistency doesn't hurt that much */
	for( nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		/* this implements the element off limit xxx command, where
		 * xxx is the limit to the smallest n(A)/n(H) that will remain on */
		if( abund.GasPhase[nelem] < dense.AbundanceLimit )
			dense.lgElmtOn[nelem] = false;

		if( dense.lgElmtOn[nelem] )
		{
			dense.SetGasPhaseDensity( nelem, abund.GasPhase[nelem]*dense.gas_phase[ipHYDROGEN] );
			if( dense.gas_phase[nelem] <= 0. )
			{
				fprintf( ioQQQ, " Abundances must be greater than zero.  "
					"Check entered abundance for element%3ld  = %2.2s, value=%.2e\n",
					nelem, elementnames.chElementSym[nelem] , dense.gas_phase[nelem]);
				fprintf( ioQQQ, " Use the ELEMENT <NAME> OFF command to disable an element.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			else if( dense.gas_phase[nelem] < SMALLFLOAT )
			{
				fprintf(ioQQQ," Abundance for %s is %.2e, less than lower "
					"limit of %.3e, so turning element off.\n",
					elementnames.chElementSym[nelem],
					dense.gas_phase[nelem],
					SMALLFLOAT );
				dense.lgElmtOn[nelem] = false;
			}
			else if( dense.gas_phase[nelem] > MAX_DENSITY )
			{
				fprintf(ioQQQ," Abundance for %s is %.2e.  This version of Cloudy does not "
					"permit densities greater than %e cm-3.\n",
					elementnames.chElementSym[nelem],
					dense.gas_phase[nelem],
					MAX_DENSITY );
				cdEXIT(EXIT_FAILURE);
			}
		}
		if( !dense.lgElmtOn[nelem] )
		{
			/* >>chng 04 apr 20, set to zero if element is off */
			dense.SetGasPhaseDensity( nelem, 0. );
		}

		/* Set all neutral ions to maintain invariant */
		dense.xIonDense[nelem][0] = dense.gas_phase[nelem];
		for( long int ion=1; ion < LIMELM+1; ion++ )
		{
			dense.xIonDense[nelem][ion] = 0.;
		}

		if( deut.lgElmtOn )
			InitDeuteriumIonization();
	}

	SumDensities();

	/* if stop temp set below default then we are going into cold and possibly 
	 * molecular gas - check some parameters in this case */
	if( called.lgTalk && (StopCalc.TempLoStopZone < phycon.TEMP_STOP_DEFAULT || 
		/* thermal.ConstTemp def is zero, set pos when used */
		(thermal.ConstTemp > 0. && thermal.ConstTemp < phycon.TEMP_STOP_DEFAULT ) ) )
	{

		/* print warning if temperature set below default but C > O */
		if( dense.lgElmtOn[ipOXYGEN] && dense.gas_phase[ipCARBON]/SDIV( dense.gas_phase[ipOXYGEN]) >= 1. )
		{
			fprintf( ioQQQ, "\n >>> \n"
							" >>> The simulation is going into possibly molecular gas but the carbon/oxygen abundance ratio is greater than unity.\n" );
			fprintf( ioQQQ, " >>> Standard interstellar chemistry networks are designed for environments with C/O < 1.\n" );
			fprintf( ioQQQ, " >>> The chemistry network may (or may not) collapse deep in molecular regions where CO is fully formed.\n" );
			fprintf( ioQQQ, " >>> \n\n\n\n\n" );
		}
	}

	if( trace.lgTrace )
	{
		realnum sumx , sumy , sumz = 0.;

		sumx = dense.gas_phase[ipHYDROGEN]*dense.AtomicWeight[ipHYDROGEN];
		sumy = dense.gas_phase[ipHELIUM]*dense.AtomicWeight[ipHELIUM];

		fprintf( ioQQQ, "\n AbundancesSet sets following densities (cm^-3); \n" );
		for( i=0; i<3; i++ )
		{
			for( nelem=i*10; nelem < i*10+10; nelem++ )
			{
				fprintf( ioQQQ, " %2.2s", elementnames.chElementSym[nelem] );
				PrintE82( ioQQQ, dense.gas_phase[nelem] );
				if( nelem>ipHELIUM )
					sumz += dense.gas_phase[nelem]*dense.AtomicWeight[nelem];
			}
			fprintf( ioQQQ, " \n" );
		}
		fprintf( ioQQQ, "\n AbundancesSet sets following abundances rel to H; \n" );
		for( i=0; i<3; i++ )
		{
			for( nelem=i*10; nelem < i*10+10; nelem++ )
			{
				fprintf( ioQQQ, " %2.2s", elementnames.chElementSym[nelem] );
				PrintE82( ioQQQ, dense.gas_phase[nelem]/dense.gas_phase[ipHYDROGEN] );
			}
			fprintf( ioQQQ, " \n" );
		}
		fprintf( ioQQQ, " \n" );
		fprintf(ioQQQ," Gas-phase mass fractions, X:%.3e Y:%.3e Z:%.3e\n\n",
			sumx/SDIV(sumx+sumy+sumz) , 
			sumy/SDIV(sumx+sumy+sumz) , 
			sumz/SDIV(sumx+sumy+sumz) );
	}
	return;
}

/* this is number of elements across one line */
#define	NELEM1LINE	9

/*PrtElem print chemical composition at start of calculation */
STATIC void PrtElem(
  /* the job to do, the options are "init", "fill", "flus" */
  const char *chJob, 
  /* label for the element */
  const char *chLabl, 
  /* its abundance */
  double abund_prt)
{
	static char chAllLabels[NELEM1LINE][14];/* buffer where elements will be stored*/
	long int i, 
	  noffset;
	static long int nelem;  /* counter for number of elements read in*/

	DEBUG_ENTRY( "PrtElem()" );

	if( strcmp(chJob,"initG") == 0 )
	{
		/* gas phase abundances */
		nelem = 0;
		fprintf( ioQQQ, 
			"                                                  Gas Phase Chemical Composition\n" );
	}
	else if( strcmp(chJob,"initD") == 0 )
	{
		/* abundances in grains */
		nelem = 0;
		fprintf( ioQQQ, 
			"                                                    Grain Chemical Composition\n" );
	}

	else if( strcmp(chJob,"fill") == 0 )
	{
		/* print log of abundance to avoid exponential output */
		abund_prt = log10( abund_prt );
		/* stuff in labels and abundances */
		sprintf( chAllLabels[nelem], "  %2.2s:%8.4f", chLabl, abund_prt );
		if( nelem == NELEM1LINE-1 )
		{
			/* we hit as many as it will hold - print it out and reset*/
			fprintf( ioQQQ, "      " );
			for( i=0; i < NELEM1LINE; i++ )
			{
				fprintf( ioQQQ, "%13.13s", chAllLabels[i] );
			}
			fprintf( ioQQQ, "\n" );
			/* reset counter to zero */
			nelem = 0;
		}
		else
		{
			/* just increment */
			++nelem;
		}
	}

#	if 0
	/* Do this if you want to know about PAH number abundance */
	else if( strcmp(chJob,"fillp") == 0 )
	{
		/* print log of abundance to avoid exponential output */
		abund_prt = log10( abund_prt );

		/* stuff in labels and abundances */
		sprintf( chAllLabels[nelem], "  %2.2s:%8.4f", chLabl, abund_prt );
		if( nelem == NELEM1LINE-1 )
		{
			/* we hit as many as it will hold - print it out and reset*/
			fprintf( ioQQQ, "      " );
			for( i=0; i < NELEM1LINE; i++ )
			{
				fprintf( ioQQQ, "%13.13s", chAllLabels[i] );
			}
			fprintf( ioQQQ, "\n" );
			/* reset counter to zero */
			nelem = 0;
		}
		else
		{
			/* just increment */
			++nelem;
		}
	}
#	endif

	else if( strcmp(chJob,"flus") == 0 )
	{
		/* flush the stack */
		i = NELEM1LINE - (nelem - 2);
		noffset = i/2-1;
		/* make format pretty */
		fprintf( ioQQQ, "      " );

		for(i=0; i < noffset; i++)
		{
			/* skip out this many fields */
			fprintf( ioQQQ, "             " );
		}

		/* if nelem is even we need to space out another 8 */
		if( !(nelem%2) && nelem > 0)
			fprintf( ioQQQ,"        ");

		for( i=0; i < nelem; i++ )
		{
			fprintf( ioQQQ, "%13.13s", chAllLabels[i] );
		}

		fprintf( ioQQQ, "\n" );
	}
	else
	{
		fprintf( ioQQQ, " PrtElem does not understand job=%4.4s\n", 
		  chJob );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}


/*AbundancesTable interpolate on table of points to do 'element table' command, */
double AbundancesTable(double r0, 
  double depth, 
  long int iel)
{
	bool lgHit;
	long int j;
	double frac, 
	  tababun_v, 
	  x;

	DEBUG_ENTRY( "AbundancesTable()" );
	/* interpolate on table of points to do 'element table' command, based 
	 * on code by K Volk, each line is log radius and abundance. */

	/* interpolate on radius or depth? */
	if( abund.lgAbTaDepth[iel-1] )
	{
		/* depth key appeared = we want depth */
		x = log10(depth);
	}
	else
	{
		/* use radius */
		x = log10(r0);
	}

	/* this will be reset below, but is here as a safety check */
	tababun_v = -DBL_MAX;

	if( x < abund.AbTabRad[0][iel-1] || x >= abund.AbTabRad[abund.nAbunTabl-1][iel-1] )
	{
		fprintf( ioQQQ, " requested radius outside range of AbundancesTable\n" );
		fprintf( ioQQQ, " radius was%10.2e min, max=%10.2e%10.2e\n", 
		  x, abund.AbTabRad[0][iel-1], abund.AbTabRad[abund.nAbunTabl-1][iel-1] );
		cdEXIT(EXIT_FAILURE);
	}

	else
	{
		lgHit = false;
		j = 1;

		while( !lgHit && j <= abund.nAbunTabl - 1 )
		{
			if( abund.AbTabRad[j-1][iel-1] <= (realnum)x && 
				abund.AbTabRad[j][iel-1] > (realnum)x )
			{
				frac = (x - abund.AbTabRad[j-1][iel-1])/(abund.AbTabRad[j][iel-1] - 
				  abund.AbTabRad[j-1][iel-1]);
				tababun_v = abund.AbTabFac[j-1][iel-1] + frac*
				  (abund.AbTabFac[j][iel-1] - abund.AbTabFac[j-1][iel-1]);
				lgHit = true;
			}
			++j;
		}

		if( !lgHit )
		{
			fprintf( ioQQQ, " radius outran dlaw table scale, requested=%6.2f largest=%6.2f\n", 
			  x, abund.AbTabRad[abund.nAbunTabl-1][iel-1] );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* got it, now return value, not log of density */
	tababun_v = exp10(tababun_v);
	return( tababun_v );
}

