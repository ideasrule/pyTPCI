/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "hextra.h"
t_hextra hextra;

void t_hextra::zero()
{
	DEBUG_ENTRY( "t_hextra::zero()" );
	/* effects of fast neutrons */
	hextra.frcneu = 0.;
	hextra.effneu = 1.;
	hextra.totneu = 0.;
	hextra.lgNeutrnHeatOn = false;
	hextra.CrsSecNeutron = 4e-26;	
}

