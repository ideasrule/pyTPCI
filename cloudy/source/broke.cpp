/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"

#include "broke.h"

#include "version.h"
#include "warnings.h"

t_broke broke;

void t_broke::comment(t_warnings& w)
{
	DEBUG_ENTRY( "t_broke::comment()" );

	/* this flag set by call to fixit routine,
	 * to show that parts of the code need repair. 
	 * lgRelease is true if this is release version */
	if( lgFixit && !t_version::Inst().lgRelease )
	{
		w.bangin( "  !The code needs to be fixed - search for fixit()." );
	}

	/* this flag set by call to CodeReview routine,
	* to show that parts of the code need to be reviewed. 
	* lgRelease is true if this is release version */
	if( lgCheckit  && !t_version::Inst().lgRelease )
	{
		w.bangin( "  !New code needs to be reviewed - search for CodeReview()." );
	}
	/* this flag set by call to routine broken ( ); 
	 * and show that the code is broken. */
	if( broke.lgBroke )
	{
		w.warnin( " W-The code is broken - search for broken()." );
	}
}
	
