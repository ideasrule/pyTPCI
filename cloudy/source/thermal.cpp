/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "thermal.h"
#include "cooling.h"
t_thermal thermal;

void t_thermal::zero()
{
	DEBUG_ENTRY( "t_thermal::zero()" );
	for( long int i=0; i<LIMELM; i++)
		heavycollcool[i] = 0.;

	lgCNegChk = true;
	CoolHeatMax = 0.;
	wlCoolHeatMax = 0;
	totcol = 0.;
	heatl = 0.;
	coolheat = 0.;
	lgCExtraOn = false;
	CoolExtra = 0.;
	ctot = 1.;

	htot = 1.;
	power = 0.;
	FreeFreeTotHeat = 0.;
	char_tran_cool = 0.;
	char_tran_heat = 0.;
	
	HeatZero();
	CoolZero();
}
