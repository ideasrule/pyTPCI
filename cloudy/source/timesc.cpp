/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "timesc.h"
#include "physconst.h"
#include "struc.h"
t_timesc timesc;

void t_timesc::zero()
{
	/* age of the cloud, to check for time-steady */
	CloudAgeSet = -1.f;
	/* some timescale for CO and H2 */
	time_H2_Dest_longest = 0.;
	time_H2_Form_longest = 0.;
	/* remains neg if not evaluated */
	time_H2_Dest_here = -1.;
	time_H2_Form_here = 0.;

	BigCOMoleForm = 0.;

	TimeH21cm = 0.;
	sound_speed_isothermal = 0.;	
}

void t_timesc::calc_therm_timesc( long int izone )
{
	if( izone <= 0 )
	{
		izone = MAX2(1,nzone-1);
		ASSERT( izone < struc.nzlim );

		/* NZLIM is the size of the struc vectors - should be dynamic 
		 * izone = MIN2( izone, NZLIM-1 );*/
	}

	time_therm_long = 0;
	time_therm_short = BIGDOUBLE;

	for( long i = 0; i < izone; i++ )
	{
		/* >>chng 99 feb 01, had div by heating, changed to cooling so constant
		 * temperature models are more realistic */
		if( struc.coolstr[i] > SMALLFLOAT )
		{
			double dt = 1.5 * BOLTZMANN * struc.DenParticles[i] * struc.testr[i] / struc.coolstr[i];
			//double dt = 1.5 * BOLTZMANN * struc.DenParticles[i] * safe_div( (double)struc.testr[i] , struc.coolstr[i] );
			time_therm_long = MAX2( time_therm_long , dt );
			time_therm_short= MIN2( time_therm_short, dt );
			//      printf("dt = %g\t long = %g\t short = %g\n", dt, time_therm_long, time_therm_short);}
		}
	//      printf( "*** long = %g\t short = %g\n", time_therm_long, time_therm_short );
	}

	return;
}
