/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*wcnint initialize stack or warnings, cautions, notes */
/*warnin enter warnings at the end of the calculations into large stack */
/*notein enter a note about calculation into comment array */
/*bangin called by routine comment to enter surprise into comment stack */
/*caunin called by comment to enter caution into comment stack */
#include "cddefines.h"
#include "warnings.h"

t_warnings warnings;

void t_warnings::zero(void)
{
	DEBUG_ENTRY( "t_warnings::zero()" );

	/* this sub is called first, to initialize the variables */
	chRgcln.clear();
	chWarnln.clear();
	chCaunln.clear();
	chBangln.clear();
	chNoteln.clear();
	lgWarngs = false;
	lgCautns = false;
}
