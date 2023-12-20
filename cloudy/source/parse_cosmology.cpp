/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCosmology parse cosmological parameters and options */
#include "cddefines.h"
#include "cosmology.h"
#include "parser.h"
#include "rfield.h"

void ParseCosmology( Parser &p )
{
	realnum param;

	DEBUG_ENTRY( "ParseCosmology()" );

	cosmology.lgDo = true;
	
	// Use Sobolev optical depths	
	strcpy( rfield.chDffTrns, "SOB" );

	/* set an omega parameter */
	if( p.nMatch("OMEG") )
	{
		param = (realnum)p.FFmtRead();

		if( p.nMatch( "BARY") )
		{
			cosmology.omega_baryon = param;
		}
		else if( p.nMatch( "RADI") )
		{
			cosmology.omega_rad = param;
		}
		else if( p.nMatch( "MATT") )
		{
			cosmology.omega_matter = param;
		}
		else if( p.nMatch( "LAMB") )
		{
			cosmology.omega_lambda = param;
		}
		else if( p.nMatch( " K ") || p.nMatch( "CURV") )
		{
			cosmology.omega_k = param;
		}
		else
		{
			fixit("must terminate clause safely.  what should it do?");
			TotalInsanity();
		}
	}
	else if( p.nMatch("HUBB") )
	{
		param = (realnum)p.FFmtRead();

		cosmology.h = param;

		if( param <= 0 || param > 1.1 )
		{
			fprintf( ioQQQ," This sets the variable h, which has units 100 km/s/Mpc, and is typically 0.71\n" );
		}
	}

	return;
}
