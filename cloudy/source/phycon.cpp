/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "phycon.h"

t_phycon phycon;

void t_phycon::zero()
{
	DEBUG_ENTRY( "t_phycon::zero()" );
	lgPhysOK = true;
	/* largest relative changes in Te, ne, H+, H2, and CO in structure
	 * this is computed as part of prtcomment so does not exist when code not talking,
	 * set to zero in zero and still zero if prtcomment not called */
	BigJumpTe = 0.;
	BigJumpne = 0.;
	BigJumpH2 = 0.;
	BigJumpCO = 0.;
	
}
