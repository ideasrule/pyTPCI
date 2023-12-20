/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*AgeCheck check various timescales after calculation complete to confirm time steady OK */
#include "cddefines.h"
#include "physconst.h"
#include "prt.h"
#include "mole.h"
#include "struc.h"
#include "warnings.h"
#include "dense.h"
#include "timesc.h"

void AgeCheck(void)
{
	char chLine[INPUT_LINE_LENGTH];
	long int i, 
	  limit;
	double hold,
	  tlong, 
	  tsound;

	DEBUG_ENTRY( "AgeCheck()" );

	/* cloud age of zero means that age command turned off
	 * negative cloud age means was not set */

	/* remember longest timescale */
	tlong = 0.;

	limit = MAX2(1,nzone-1);
	ASSERT( limit < struc.nzlim );

	/* thermal equilibrium timescales */
	timesc.calc_therm_timesc( limit );


	/* NZLIM is the size of the struc vectors - should be dynamic 
	limit = MIN2( limit , NZLIM-1 );*/

	for( i=0; i < limit; i++ )
	{
		if( struc.coolstr[i]>SMALLFLOAT )
		{
			timesc.time_therm_long =
				MAX2( timesc.time_therm_long ,
				struc.DenParticles[i]*BOLTZMANN*1.5*struc.testr[i]/struc.coolstr[i]);
			timesc.time_therm_short =
				MIN2( timesc.time_therm_short ,
				struc.DenParticles[i]*BOLTZMANN*1.5*struc.testr[i]/struc.coolstr[i]);
			/*>>chng 99 feb 01, had div by heating, changed to cooling so constant
			 * temperature models are more realistic */
		}
	}

	tlong = MAX2(tlong,timesc.time_therm_long);
	if( prt.lgPrnAges )
	{
		sprintf( chLine, "   AGE: longest thermal timescale= %.2es.", 
		  timesc.time_therm_long );
		notein(chLine);
	}

	tlong = MAX2(tlong,timesc.TimeH21cm);
	if( prt.lgPrnAges )
	{
		sprintf( chLine, "   AGE: 21 cm equilibrium timescale= %.2es.", 
		  timesc.TimeH21cm );
		notein(chLine);
	}

	if( timesc.CloudAgeSet > 0. && timesc.time_therm_long > timesc.CloudAgeSet )
	{
		sprintf( chLine, " C-AGE: Thermal equilibrium timescale, %.2es, longer than age", 
		  timesc.time_therm_long );
		caunin(chLine);
	}

	/* check soundt travel time if constant pressure */
	if( strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		tsound = timesc.sound;
		if( prt.lgPrnAges )
		{
			sprintf( chLine, "   AGE: sound travel time= %.2es.", 
			  tsound );
			notein(chLine);
		}

		if( timesc.CloudAgeSet > 0. && tsound > timesc.CloudAgeSet )
		{
			sprintf( chLine, " C-AGE: Sound travel time longer than age in constant pressure model = %.2es", 
			  timesc.time_therm_long );
			caunin(chLine);
		}
	}

	else
	{
		/* do not check if not constant pressure */
		tsound = 0.;
	}
	tlong = MAX2(tlong,tsound);

	/* molecule formation timescale */
	/* >>chng 04 sep 17, - if species are negligible will set to negative number
	 * to retain value but not include in timescales */
	if( findspecieslocal("H2")->xFracLim < 1e-2 )
	{
		timesc.time_H2_Dest_longest *= -1.;
		timesc.time_H2_Form_longest *= -1.;
	}
	tlong = MAX2( tlong , timesc.time_H2_Dest_longest );
	tlong = MAX2( tlong , timesc.time_H2_Form_longest );

	if( findspecieslocal("CO")->xFracLim < 1e-2 )
	{
		timesc.BigCOMoleForm *= -1.;
	}
	tlong = MAX2( tlong , timesc.BigCOMoleForm );

	/* >>chng 97 jan 03, don't print if zero */
	if( prt.lgPrnAges && timesc.time_H2_Dest_longest > 0. )
	{
		sprintf( chLine, "   AGE: longest H2 destruction timescale= %.2es.", 
		  timesc.time_H2_Dest_longest );
		notein(chLine);
	}

	if( prt.lgPrnAges && timesc.time_H2_Form_longest > 0. )
	{
		sprintf( chLine, "   AGE: longest H2 formation timescale= %.2es.", 
		  timesc.time_H2_Form_longest );
		notein(chLine);
	}

	if( timesc.CloudAgeSet > 0. && timesc.time_H2_Dest_longest > timesc.CloudAgeSet )
	{
		sprintf( chLine, " C-AGE: H2 destruction timescale longer than age, = %.2es", 
		  timesc.time_H2_Dest_longest );
		caunin(chLine);
	}

	if( timesc.CloudAgeSet > 0. && timesc.time_H2_Form_longest > timesc.CloudAgeSet )
	{
		sprintf( chLine, " C-AGE: H2 formation timescale longer than age, = %.2es", 
		  timesc.time_H2_Form_longest );
		caunin(chLine);
	}

	if( prt.lgPrnAges && timesc.BigCOMoleForm > 0. )
	{
		sprintf( chLine, "   AGE: longest CO formation timescale= %.2es.", 
		  timesc.BigCOMoleForm );
		notein(chLine);
	}

	if( timesc.CloudAgeSet > 0. && timesc.BigCOMoleForm > timesc.CloudAgeSet )
	{
		sprintf( chLine, " C-AGE: CO formation timescale longer than age, = %.2es", 
		  timesc.BigCOMoleForm );
		caunin(chLine);
	}

	/* hydrogen recombination timescale */
	timesc.time_Hrecom_long = 0.;
	timesc.time_Hrecom_short = 0.;
	for( i=0; i < limit; i++ )
	{
		if( struc.ednstr[i]>SMALLFLOAT )
		{
			hold = struc.ednstr[i]*2.90e-10*pow(struc.testr[i],(realnum)-0.77f);
			timesc.time_Hrecom_long = MAX2(timesc.time_Hrecom_long , 1./hold);
			timesc.time_Hrecom_short = MAX2(timesc.time_Hrecom_short , 1./hold);
		}
	}

	tlong = MAX2(tlong,timesc.time_Hrecom_long);
	if( prt.lgPrnAges )
	{
		sprintf( chLine, "   AGE: longest H recombination timescale= %.2es.", 
		  timesc.time_Hrecom_long );
		notein(chLine);
	}

	if( timesc.CloudAgeSet > 0. && timesc.time_Hrecom_long > timesc.CloudAgeSet )
	{
		sprintf( chLine, " C-AGE: Hydrogen recombination timescale longer than age, = %.2es", 
		  timesc.time_Hrecom_long );
		caunin(chLine);
	}

	/* give age in various units, depending on order of magnitude */
	if( timesc.CloudAgeSet < 0. )
	{
		/* CloudAgeSet initially set to -1, if still the case then age not set */
		if( tlong < 3600. )
		{
			/* less than one day, give only seconds */
			sprintf( chLine, "  !AGE: Cloud age was not set.  Longest timescale was %.2e s.", 
			  tlong );
			bangin(chLine);
		}

		else if( tlong < 8.64e4 )
		{
			/* less than one day, give seconds and hours */
			sprintf( chLine, "  !AGE: Cloud age was not set.  Longest timescale was %.2e s = %.2e hours.", 
			  tlong, tlong/3600. );
			bangin(chLine);
		}

		else if( tlong < 3e7/12. )
		{
			/* less than one month, give seconds and days */
			sprintf( chLine, "  !AGE: Cloud age was not set.  Longest timescale was %.2e s = %.2e days.", 
			  tlong, tlong/86400. );
			bangin(chLine);
		}

		else if( tlong < 3e7 )
		{
			/* less than one year, give seconds and months */
			sprintf( chLine, "  !AGE: Cloud age was not set.  Longest timescale was %.2e s = %.2e months.", 
			  tlong, (tlong/3.15569e7)*12. );
			bangin(chLine);
		}

		else
		{
			/* more than one year, give seconds and years */
			sprintf( chLine, "  !AGE: Cloud age was not set.  Longest timescale was %.2e s = %.2e years.", 
			  tlong, tlong/3.15569e7 );
			bangin(chLine);
		}
	}

	else
	{
		/* age set, and passed tests, still say longest */
		if( tlong < 3e7 )
		{
			/* less than one year, give only seconds */
			sprintf( chLine, "   AGE: Cloud age was %.2es, Longest timescale was %.2e s.", 
			  timesc.CloudAgeSet, tlong );
			notein(chLine);
		}

		else
		{
			/* more than one year, give seconds and years */
			sprintf( chLine, "   AGE: Cloud age was %.2e s.  Longest timescale was %.2e s = %.2e years.", 
			  timesc.CloudAgeSet, tlong, tlong/3.15569e7 );
			notein(chLine);
		}
	}
	return;
}
