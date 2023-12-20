/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef TIMESC_H_
#define TIMESC_H_

/* timesc.h */

#include "module.h"

/**AgeCheck check various timescales after calculation complete to confirm time steady OK */
void AgeCheck(void);

struct t_timesc : public module {
	/** timescales tracked by the code, some set in checkage, some in radinc */

	void zero();
	void comment(t_warnings&) {}

	const char *chName() const
	{
		return "timesc";
	}

	/** compton equilibrium timescale */
	double tcmptn;

	/** thermal timescales */
	double time_therm_long , 
		time_therm_short;

	/* compute thermal timescales */
	void calc_therm_timesc( long izone );

	/** hydrogen recombination timescale */
	double time_Hrecom_long ,
		time_Hrecom_short; 

	/** sound is sound travel time (sec) updated in radinc */
	double sound;

	/** age of cloud set with age command */
	realnum CloudAgeSet;

	/** local and longest timescales for CO and H2 molecules to form */
	double 
		/* allocated in zero.c */
	  time_H2_Dest_longest, 
	  time_H2_Form_longest, 
	  time_H2_Dest_here, 
	  time_H2_Form_here, 
	  BigCOMoleForm;

	/** isothermal sound speed */
	double sound_speed_isothermal;

	/** adiabatic sound speed assuming monatomic gas - gamnma is 5/3*/
	double sound_speed_adiabatic;

	/** timescale for photoerosion of iron due to gamma rays */
	realnum TimeErode;

	/** H 21 cm line to come into equilibrium */
	double TimeH21cm;

	};
extern t_timesc timesc;


#endif /* TIMESC_H_ */
