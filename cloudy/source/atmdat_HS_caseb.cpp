/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* atmdat_HS_caseB - interpolate on line emissivities from Storey & Hummer tables for hydrogen */
#include "cddefines.h"  
#include "atmdat.h"  

double atmdat_HS_caseB( 
	  /* upper and lower quantum numbers, abort if iHi>25*/
	  long int iHi, 
	  long int iLo, 
	  /* element number, 1 or 2 at this point, but decremented to C scale later */
	  long int nelem, 
	  /* temperature to interpolate, for H between 500-30,000K*/
	  double TempIn, 
	  /* density to interpolate */
	  double DenIn,
	  /* case - 'a' or 'b' */
	  char chCase )

/* general utility to interpolate line emissivities from the Storey & Hummer tables
 of case B emissivities.  
 iHi<25, iLo, the principal quantum
 numbers, and are upper and lower levels in any order.  
 nelem element number on physicial scale, 1 or 2 have data
 TempIn = temperature, and must lie within the range of the table, which depends on
 the ion charge, and is 500 - 30,000K for hydrogen.  
 DenIn is the density and must not exceed the high density limit to the table.  

 routine returns -1 if conditions outside temperature range, or
 above density range of tabulated case b results
 If desired density is low limit, lower limit is used instead 
*/

{
	long int 
		ipTemp, /*pointer to temperature*/
		ipDens, /*pointer to density*/
		ipDensHi ,
		ipTempHi;
	int ipUp , ipLo , ip;
	double x1 , x2 , yy1 , yy2 , z1 , z2 , za , zb ,z;
	int iCase;

	DEBUG_ENTRY( "atmdat_HS_caseB()" );

	/*make sure nelem is 1 or 2*/
	if( nelem<ipHYDROGEN || nelem> HS_NZ )
	{
		printf("atmdat_HS_caseB called with improper nelem, was%li and must be 1 or 2",nelem);
		cdEXIT(EXIT_FAILURE);

	}
	/* now back nelem back one, to be on c scale */
	--nelem;

	/* case A or B? */
	if( chCase == 'a' || chCase=='A' )
	{
		iCase = 0;
	}
	else if(  chCase == 'b' || chCase=='B' )
	{
		iCase = 1;
	}
	else
	{
		printf("atmdat_HS_caseB called with improper case, was %c and must be A or B",chCase);
		cdEXIT(EXIT_FAILURE);
	}

	/*===========================================================*/
	/* following is option to have principal quantum number given in either order, 
	 * final result is that ipUp and ipLo will be the upper and lower levels */
	if( iHi > iLo )
	{
		ipUp = (int)iHi;  ipLo = (int)iLo;
	}
	else if( iHi < iLo )
	{
		ipUp = (int)iLo; ipLo = (int)iHi;
	}
	else
	{ 
		printf("atmdat_HS_caseB called with indices equal, %ld  %ld  \n",iHi,iLo);
		cdEXIT(EXIT_FAILURE);
	}

	/* now check that they are in range of the predicted levels of their model atom*/
	if( ipLo <1 ) 
	{
		printf(" atmdat_HS_caseB called with lower quantum number less than 1, = %i\n",
			ipLo);
		cdEXIT(EXIT_FAILURE);
	}

	if( ipUp >25 ) 
	{
		printf(" atmdat_HS_caseB called with upper quantum number greater than 25, = %i\n",
			ipUp);
		cdEXIT(EXIT_FAILURE);
	}

	/*===========================================================*/
	/*bail if above high density limit */
	if( DenIn > atmdat.Density[iCase][nelem][atmdat.nDensity[iCase][nelem]-1] ) 
	{
		/* this is flag saying bogus results */
		return -1;
	}

	if( DenIn <= atmdat.Density[iCase][nelem][0] )
	{
		/* this case, desired density is below lower limit to table,
		 * just use the lower limit */
		ipDens = 0;
	}
	else
	{
		/* this case find where within table density lies */
		for( ipDens=0; ipDens < atmdat.nDensity[iCase][nelem]-1; ipDens++ )
		{
			if( DenIn >= atmdat.Density[iCase][nelem][ipDens] && 
				DenIn < atmdat.Density[iCase][nelem][ipDens+1] ) break;
		}
	}


	/*===========================================================*/
	/* confirm within temperature range*/
	if( TempIn < atmdat.ElecTemp[iCase][nelem][0] || 
		TempIn > atmdat.ElecTemp[iCase][nelem][atmdat.ntemp[iCase][nelem]-1] ) 
	{
		/* this is flag saying bogus results */
		return -1;
	}

	/* find where within grid this temperature lies */ 
	for( ipTemp=0; ipTemp < atmdat.ntemp[iCase][nelem]-1; ipTemp++ )
	{
		if( TempIn >= atmdat.ElecTemp[iCase][nelem][ipTemp] && 
			TempIn < atmdat.ElecTemp[iCase][nelem][ipTemp+1] ) break;
	}

	/*===========================================================*/
	/*we now have the array indices within the temperature array*/

	if( ipDens+1 > atmdat.nDensity[iCase][nelem]-1 )
	{ 
		/* special case, when cell is highest density point */
		ipDensHi = atmdat.nDensity[iCase][nelem]-1;
	}
	else if( DenIn < atmdat.Density[iCase][nelem][0])
	{
		 /* or density below lower limit to table, set both bounds to 0 */
		ipDensHi = 0;
	}
	else
	{ 
		ipDensHi = ipDens+1;
	}

	/*special case, if cell is highest temperature point*/
	if( ipTemp+1 > atmdat.ntemp[iCase][nelem]-1 )
	{ 
		ipTempHi = atmdat.ntemp[iCase][nelem]-1;
	}
	else
	{ 
		ipTempHi = ipTemp+1;
	}

	x1 = log10( atmdat.Density[iCase][nelem][ipDens] );
	x2 = log10( atmdat.Density[iCase][nelem][ipDensHi] );

	yy1 = log10( atmdat.ElecTemp[iCase][nelem][ipTemp] );
	yy2 = log10( atmdat.ElecTemp[iCase][nelem][ipTempHi] );

	/*now generate the index to the array, expression from Storey code -1 for c*/
	ip = (int)((((atmdat.ncut[iCase][nelem]-ipUp)*(atmdat.ncut[iCase][nelem]+ipUp-1))/2)+ipLo - 1);

	/*pointer must lie within line array*/
	ASSERT( ip < NLINEHS ); 
	ASSERT( ip >= 0 );

	/* interpolate on emission rate*/
	z1 = log10( atmdat.Emiss[iCase][nelem][ipTemp][ipDens][ip]);
	z2 = log10( atmdat.Emiss[iCase][nelem][ipTemp][ipDensHi][ip]);

	if( fp_equal( x2, x1 ) ) 
	{
		za = z2;
	}
	else 
	{
		za = z1 + log10( DenIn / atmdat.Density[iCase][nelem][ipDens] ) * (z2-z1)/(x2-x1);
	}

	z1 = log10( atmdat.Emiss[iCase][nelem][ipTempHi][ipDens][ip]);
	z2 = log10( atmdat.Emiss[iCase][nelem][ipTempHi][ipDensHi][ip]);

	if( fp_equal( x2, x1 ) )
	{
		zb = z2;
	}
	else 
	{
		zb = z1 + log10( DenIn / atmdat.Density[iCase][nelem][ipDens] ) * (z2-z1)/(x2-x1);
	}

	if( fp_equal( yy2, yy1 ) )
	{
		z = zb;
	}
	else 
	{
		z = za + log10( TempIn / atmdat.ElecTemp[iCase][nelem][ipTemp] ) * (zb-za)/(yy2-yy1);
	}

	return ( exp10(z) );
}
