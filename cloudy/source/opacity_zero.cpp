/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*OpacityZero zero out opacity save arrays, save old opacity in OldOpacSave array */
#include "cddefines.h"
#include "rfield.h"
#include "opacity.h"

void OpacityZero(void)
{
	long int i;

	DEBUG_ENTRY( "OpacityZero()" );


	for( i=0; i < rfield.nflux_with_check; i++ )
	{
		opac.opacity_sct[i] = 0.;
		/* save the current opacities */
		opac.OldOpacSave[i] = opac.opacity_abs[i];
		opac.opacity_abs[i] = 0.;
	}

	/* only zero out the static array if we are going to
	 * totally redo the static part */
	if( opac.lgRedoStatic )
	{
		/*fprintf(ioQQQ," OpacityZero is zeroing out the static opacities\n");*/
		for( i=0; i < rfield.nflux_with_check; i++ )
		{
			opac.OpacStatic[i] = 0.;
		}
	}
	return;
}

/* set old opac array to current versin during search phase */
void OpacityZeroOld(void)
{
	long int i;

	DEBUG_ENTRY( "OpacityZeroOld()" );


	for( i=0; i < rfield.nflux_with_check; i++ )
	{
		/* save the current opacities */
		opac.OldOpacSave[i] = opac.opacity_abs[i];
	}
	return;
}
