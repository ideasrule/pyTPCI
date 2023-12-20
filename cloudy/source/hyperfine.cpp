/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "hyperfine.h"
t_hyperfine hyperfine;

void t_hyperfine::zero()
{
	DEBUG_ENTRY( "t_hyperfine::zero()" );
	lgLya_pump_21cm = true;	
	LyaSourceFunctionShape_assumed = EXCITATION;
}
