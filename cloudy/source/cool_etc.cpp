/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolSum  total cooling from all entries into cooling stack */
/*CoolZero set cooling and heating stack to zero */
/*CoolAdd add coolants to the cooling stack, called in evaluation of cooling function */
#include "cddefines.h"
#include "taulines.h"
#include "thermal.h"
#include "ionbal.h"
#include "cooling.h"

/*CoolAdd add coolants to the cooling stack, called in evaluation of cooling function */
void CoolAdd(
  const char *chLabel, 
  realnum lambda, 
  double cool)
{

	DEBUG_ENTRY( "CoolAdd()" );

	/* this flag indicates (true) that we are between when cooling was set to
	 * zero with call to CoolZero, and when final sum was used.  Any call
	 * after final summation (false) will be ignored and so is fatal error */
	ASSERT( thermal.lgCoolEvalOK );

	/* this can be done with an assert since these results cannot possibly
	 * depend on user input */
	ASSERT( thermal.ncltot < NCOLNT );

	/* copy coolant label into stack */
	ASSERT( strlen( thermal.chClntLab[thermal.ncltot]) < NCOLNT_LAB_LEN );
	strcpy( thermal.chClntLab[thermal.ncltot], chLabel);

	/* now the wavelength */
	thermal.collam[thermal.ncltot] = lambda;

	/* normal line cooling */
	thermal.cooling[thermal.ncltot] = MAX2(0.,cool);

	/* possible line heating - not supposed to be done this way!
	 * this is intrinsic positive number, to be added to heating */
	thermal.heatnt[thermal.ncltot] = MAX2(0.,-cool);

	/* now increment counter, this is the number of coolants entered */
	thermal.ncltot += 1;
	return;
}

/*CoolZero set cooling and heating stack to zero */
void CoolZero(void)
{

	DEBUG_ENTRY( "CoolZero()" );

	thermal.ncltot = 0;
	thermal.dCooldT = 0.;
	thermal.dHeatdT = 0.;

	/* >>chng 03 nov 29, from explicit loop to memset to save time */
	memset(thermal.cooling , 0 , NCOLNT*sizeof(thermal.cooling[0] ) );
	memset(thermal.heatnt , 0 , NCOLNT*sizeof(thermal.heatnt[0] ) );

	/* initialize coolants' cooling data */
	for( int i = 0; i <= LIMELM ; i++ )
	    thermal.elementcool[i] = 0.;
	thermal.dima = 0.;
	
	/* this flag indicates that it is ok to add coolants to cooling
	 * stack since between first zero, and final sum - CoolAdd checks
	 * that this is true */
	thermal.lgCoolEvalOK = true;


	/* now zero out these arrays */
	for( long nelem=0; nelem< LIMELM; ++nelem )
	{
		for( long ion=0; ion<nelem+1; ++ion )
		{
			ionbal.ExcitationGround[nelem][ion] = 0.;
		}
	}

	return;
}

/*CoolSum  total cooling from all entries into cooling stack */
void CoolSum(double *total)
{
	long int i;

	DEBUG_ENTRY( "CoolSum()" );

	/* routine to add together all line heating and cooling */

	*total = 0.;
	thermal.coolheat = 0.;
	/* this is sum of agents that should be coolants
	 * coolheat will be coolants that came out as heat */
	for( i=0; i < thermal.ncltot; i++ )
	{
		*total += thermal.cooling[i];
		thermal.coolheat += thermal.heatnt[i];
	}
	thermal.setHeating(0,12, thermal.coolheat);

	/* make comment if negative cooling ever significant */
	if( thermal.htot > 0. )
	{
		if( thermal.coolheat/thermal.htot > 0.01 )
		{
			/* CoolHeatMax was set to zero at start of calc, we want very biggest */
			for( i=0; i < thermal.ncltot; i++ )
			{
				if( thermal.heatnt[i]/thermal.htot > thermal.CoolHeatMax )
				{
					thermal.CoolHeatMax = (realnum)(thermal.heatnt[i]/thermal.htot);
					thermal.wlCoolHeatMax = thermal.collam[i];
					strcpy( thermal.chCoolHeatMax, thermal.chClntLab[i] );
				}
			}
		}
	}

	/* this sum of lines that were heat sources - this
	 * part was not counted as heating in call to cooling add routine
	 * since atom_level2 and atom_level3 cooling routines separate this out 
	 * into ->cool() and ->heat() - this does
	 * NOT double count line heating */

	thermal.heatl = 0.;
	for( i=0; i < nWindLine; i++ )
	{
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
			thermal.heatl += TauLine2[i].Coll().heat();
	}

	/* line heating added in following, only here */
	thermal.setHeating(0,22, thermal.heatl);
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && thermal.heatl/thermal.ctot > 0.1 )
		{
			double thresh = 0.1;
			fprintf(ioQQQ," all heating lines > %.4f of heatl printed next \n",
				thresh );
			for( i=0; i < nWindLine; i++ )
			{
				if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
				{
					if( TauLine2[i].Coll().heat()/thermal.heatl > thresh )
						DumpLine( TauLine2[i] );
				}
			}

			/*Atomic & Molecular Lines CHIANTI & Leiden lines*/
			for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
			{
				for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
					  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
				{
					if( (*em).Tran().Coll().heat()/thermal.heatl > thresh )
						DumpLine( (*em).Tran() );
				}
			}
		}
	}

	/*begin sanity check */
	if( *total <= 0. )
	{
		fprintf( ioQQQ, " CoolSum finds cooling <= 0%10.2e\n", 
		  *total );
	}
	if( thermal.heatl/thermal.ctot < -1e-15 )
	{
		fprintf( ioQQQ, " CoolSum finds negative heating %10.2e %10.2e\n", 
		  thermal.heatl, thermal.ctot );
	}
	/*end sanity check */

	thermal.lgCoolEvalOK = false;
	return;
}
