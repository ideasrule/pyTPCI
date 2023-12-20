/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "conv.h"
t_conv conv;

void t_conv::zero()
{
	/* this counts number of times ionize is called by PressureChange, in current zone
	 * these are reset here, so that we count from first zone not search */
	nPres2Ioniz = 0;

	/* clear flag indicating the ionization convergence failures 
	 * have ever occurred in current zone
	lgConvIonizThisZone = false; */
	lgFirstSweepThisZone = false;
	lgLastSweepThisZone = false;

	resetCounters();

	/* cooling tolerance heating tolerance - allowed error in heating -  cooling balance */
	/*HeatCoolRelErrorAllowed = 0.02f;*/
	/* >>chng 04 sep 25, change te tolerance from 0.02 to 4x smaller, 0.005, drove instabilities
	 * in chemistry */
	HeatCoolRelErrorAllowed = 0.005f;

	/* this is the default allowed relative error in the electron density */
	EdenErrorAllowed = 1e-2;

	IonizErrorAllowed = 1e-2;

	dCmHdT = 0.;

	LimFail = 20;
	lgMap = false;

	/* this counts how many times ionize is called in this model after startr,
	 * and is flag used by ionize to understand it is being called the first time*/
	nTotalIoniz = 0;
	/* these are variables to remember the biggest error in the
	 * electron density, and the zone in which it occurred */
	BigEdenError = 0.;
	AverEdenError = 0.;
	BigHeatCoolError = 0.;
	AverHeatCoolError = 0.;
	BigPressError = 0.;
	AverPressError = 0.;
	strcpy( chSolverEden, "SECA" );
	strcpy( chSolverTemp, "vWDB" );
	strcpy( chNotConverged, "none" );
	resetConvIoniz();
	/* iterate to convergence flag */
	lgAutoIt = false;
	/* convergence criteria */
	autocv = 0.20f;
	lgAllTransitions = false;
	lgConvTemp = true;
	lgConvPres = true;
	lgConvEden = true;
	lgUpdateCouplings = false;
	/* >>chng 04 jan 25, only set lgConvIoniz true where used in ConvXXX path */
	/*lgConvIoniz = true;*/

	/* this is the default allowed relative error in the pressure */	
	PressureErrorAllowed = 0.01f;

	MaxFractionalDensityStepPerIteration = 0.2;

	/* default error in total gas-phase density of each element, including molecules */
	GasPhaseAbundErrorAllowed = 1e-5f;

	/* this is abort option set with SET PRESIONIZ command */
	limPres2Ioniz = 3000;

	nTeFail = 0;
	nTotalFailures = 0;
	nPreFail = 0;
	failmx = 0.;
	nIonFail = 0;
	nPopFail = 0;
	nNeFail = 0;
	nGrainFail = 0;
	dCmHdT = 0.;
}

ConvergenceCounter t_conv::register_(const string name)
{
	DEBUG_ENTRY( "t_conv::register_()" );
	for( size_t i=0; i<m_counters.size(); ++i )
	{
		if (name == m_labels[i])
		{
			fprintf(ioQQQ,"Internal error: convergence counter '%s' already"
					  " registered\n", name.c_str());
			cdEXIT(EXIT_FAILURE);
		}
	}
	size_t newtype = m_counters.size();
	m_counters.push_back(0);
	m_counters_zone.push_back(0);
	m_labels.push_back(name);
	return ConvergenceCounter(this, newtype);
}

