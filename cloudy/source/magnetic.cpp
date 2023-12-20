/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseMagnet parse magnetic field command  */
/*Magnetic_init initialize magnetic field parameters */
/*Magnetic_reinit - reinitialized magnetic field at start of new iteration */
/*Magnetic_evaluate evaluate some parameters to do with magnetic field */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "doppvel.h"
#include "optimize.h"
#include "input.h"
#include "wind.h"
#include "magnetic.h"
#include "parser.h"

t_magnetic magnetic;

/* the initial magnetic field */
static double Btangl_init;

/* this is logical var, set in zero, which says whether the magnetic field has been
	* initialized */
static bool lgBinitialized;

/* the current magnetic field */
static double Btangl_here;

/* the initial parallel and tangential fields for ordered case */
static double Bpar_init, Btan_init;

/* the current parallel and tangential fields for ordered case */
static double Bpar_here, Btan_here;

/* this is the gamma power law index, default is 4. / 3. */
static double gamma_mag;

/*Magnetic_evaluate evaluate some parameters to do with magnetic field */
void Magnetic_evaluate(void)
{

	DEBUG_ENTRY( "Magnetic_evaluate()" );

	/* this flag set true if magnetic field is specified */
	if( magnetic.lgB )
	{
		static double density_initial, 
			/* the square of the Alfven velocity at illuminated face */
			v_A;

		/* does magnetic field need to be initialized for this iteration? 
		 * flag is reset false at init of code, and at start of every iteration */
		if( !lgBinitialized )
		{
			lgBinitialized = true;

			/* set initial tangled field */
			Btangl_here = Btangl_init;

			/* set initial ordered field */
			/* mag field angle_wrt_los set when ordered field specified */
			Bpar_here = Bpar_init;
			Btan_here = Btan_init;

			/* XMassDensity was set above, safe to use this on first call */
			density_initial = dense.xMassDensity;

			/* this is used to update tangential field */
			v_A = POW2(Bpar_init) / (PI4 * density_initial );
		}

		/* now update parameters in tangled field case */
		/* magnetic pressure is a function of the local density, use gamma law */
		Btangl_here = Btangl_init * pow(dense.xMassDensity/density_initial, gamma_mag/2.);

		/* ordered components of field - parallel field is always constant - find tangential component -
		 * but only in wind case */
		if( !wind.lgStatic() )
		{
			/* N B - must preserve sign in this equation - will blow if product of wind speeds is equal to v_A) */
			/* wind.windv*wind.windv0 == v_A should not occur since mag pressure goes to inf */
			Btan_here = Btan_init * (POW2(wind.windv0) - v_A)/ (wind.windv*wind.windv0-v_A);
		}

		/* magnetic pressure - sum of two fields - can have negative pressure (tension)
		 * is ordered field dominates */
		magnetic.pressure = POW2(Btangl_here)/PI8 + (POW2(Btan_here) - POW2(Bpar_here)) / PI8;

		/* energy density - this is positive */
		magnetic.energydensity = POW2(Btangl_here)/PI8 + (POW2(Btan_here) + POW2(Bpar_here)) / PI8;

		/* option for turbulence in equipartition with B field */
		if( DoppVel.lgTurbEquiMag )
		{
			/* >>chng 05 jan 26, as per Robin Williams email,
			 * evaluate energydensity above, which is +ve, and use that for
			 * velocity here - had used pressure but could not evaluate when negative */
			/* turbulent velocity is mag pres over density */
			/* >>chng 06 apr 19, use DoppVel.Heiles_Troland_F */
			/*DoppVel.TurbVel = (realnum)sqrt(2.*magnetic.energydensity/dense.xMassDensity);*/
			DoppVel.TurbVel = (realnum)sqrt(6.*magnetic.energydensity/dense.xMassDensity/
				DoppVel.Heiles_Troland_F);
			/* >>chng 06 apr 19, do not double mag pressure, really count turbulence as pressure */
			/* double magnetic pressure to account for ram pressure due to turbulence,
			 * which is not counted elsewhere 
			magnetic.pressure *= 2.;*/
		}

		/* input parser made sure gamma != 1, default magnetic gamma is 4/3 */
		magnetic.EnthalpyDensity = gamma_mag/(gamma_mag-1.) *
			POW2(Btangl_here)/PI8 + (POW2(Btan_here) + POW2(Bpar_here)) / PI4;
	}
	else
	{
		magnetic.pressure = 0.;
		magnetic.energydensity = 0.;
		magnetic.EnthalpyDensity = 0.;
	}
	return;
}

/*Magnetic_reinit - reinitialized magnetic field at start of new iteration */
void Magnetic_reinit(void)
{
	DEBUG_ENTRY( "Magnetic_reinit()" );

	/* this says whether B has been initialized in this run */
	lgBinitialized = false;
	return;
}

/* t_magnetic::zero initialize magnetic field parameters */
void t_magnetic::zero(void)
{

	DEBUG_ENTRY( "t_magnetic::zero()" );

	gamma_mag = 4./3.;
	lgB = false;
	/* this says whether B has been initialized in this run */
	lgBinitialized = false;
	/* the initial tangled and ordered fields */
	Btangl_init = 0.;
	Btangl_here = DBL_MAX;
	pressure = DBL_MAX;
	energydensity = DBL_MAX;
	Bpar_init = 0.;
	Btan_init = 0.;
	Bpar_here = DBL_MAX;
	Btan_here = DBL_MAX;
	EnthalpyDensity = DBL_MAX;
	return;
}

/*ParseMagnet parse magnetic field command  */
void ParseMagnet(Parser &p )
{
	bool lgTangled;
	double angle_wrt_los=-1. , Border_init=-1.;

	DEBUG_ENTRY( "ParseMagnet()" );

	/* flag saying B turned on */
	magnetic.lgB = true;

	/* check whether ordered is specified - if not this is tangled */
	if( p.nMatch("ORDE") )
	{
		/* ordered field case */
		lgTangled = false;

		/* field is ordered, not isotropic - need field direction - spcified
		 * by angle away from beam of light - 0 => parallel to beam */

		/* this is the log of the ordered field strength */
		Border_init = p.getNumberCheckAlwaysLog("ordered field");

		/* this is the angle (default in degrees) with respect to the line of sight */
		angle_wrt_los = p.getNumberCheck("LOS angle");

		/* if radians is on the line then the number already is in radians - 
		 * else convert to radians */
		if( !p.nMatch("RADI") )
			angle_wrt_los *= PI2 / 360.;

		/* now get initial components that we will work with */
		Bpar_init = Border_init*cos(angle_wrt_los);
		Btan_init = Border_init*sin(angle_wrt_los);

	}
	else
	{
		/* tangled field case */
		lgTangled = true;
		/* this is the log of the tangled field strength */
		Btangl_init = p.getNumberCheckAlwaysLog("tangled field");

		/* optional gamma for dependence on pressure */
		gamma_mag = p.getNumberDefault("field gamma law", 4./3.);

		if( gamma_mag != 0. && gamma_mag <= 1. )
		{
			/* impossible value for gamma */
			fprintf( ioQQQ, 
				" This value of gamma (%.3e) is impossible.  Must have gamma = 0 or > 1.\n Sorry.\n",
				gamma_mag );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		/* number of parameters */
		optimize.nvarxt[optimize.nparm] = 2;
		if( lgTangled )
		{
			/* tangled field case */
			// keyword LOG is not needed, but its presence is checked elsewhere
			strcpy( optimize.chVarFmt[optimize.nparm], "MAGNETIC FIELD TANGLED= %f LOG GAMMA= %f" );
			optimize.vparm[0][optimize.nparm] = (realnum)log10( Btangl_init );
			/* this is not varied */
			optimize.vparm[1][optimize.nparm] = (realnum)gamma_mag;
		}
		else
		{
			/* ordered field case */
			// keyword LOG is not needed, but its presence is checked elsewhere
			strcpy( optimize.chVarFmt[optimize.nparm], "MAGNETIC FIELD ORDERED= %f LOG ANGLE RADIANS= %f" );
			optimize.vparm[0][optimize.nparm] = (realnum)log10( Border_init );
			/* this is not varied */
			optimize.vparm[1][optimize.nparm] = (realnum)angle_wrt_los;
		}

		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vincr[optimize.nparm] = 0.5;
		++optimize.nparm;
	}
	return;
}
