/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseRadius parse the radius command */
#include "cddefines.h"
/*#define	PARSCL	18.489396*/
#include "physconst.h"
#include "optimize.h"
#include "radius.h"
#include "iterations.h"
#include "input.h"
#include "parser.h"

void ParseRadius(Parser &p)
{
	bool 
	  lgR2set, 
	  lgRLog;
	double a, 
	  convl,
	  r;

	DEBUG_ENTRY( "ParseRadius()" );

	/* log of inner and outer radii, default second=infinity,
	 * if R2<R1 then R2=R1+R2
	 * there is an optional keyword, "PARSEC" on the line, to use PC as units */
	convl = p.nMatch("PARS") ? log10( PARSEC ) : 0.;

	/* if linear appears on line, then radius is linear, otherwise, log */
	lgRLog = !p.nMatch("LINE");

	r = p.FFmtRead();
	if( p.lgEOL() )
		p.NoNumb("radius");

	/* option for linear or log radius */
	if( lgRLog )
	{
		r += convl;
	}
	else
	{
		if( r > 0. )
		{
			r = log10(r) + convl;
		}
		else
		{
			fprintf(ioQQQ,"The first radius is negative and linear is set - this is impossible.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( r > 37. || r < -37. )
	{
		fprintf(ioQQQ,"WARNING - the log of the radius is %e - this is too big.\n", r );
		fprintf(ioQQQ," Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	radius.Radius = exp10(r);
	radius.lgRadiusKnown = true;

	/* check for second number, which indicates thickness or outer radius of model */
	a = p.FFmtRead();
	if( p.lgEOL() )
	{
		/* not set */
		lgR2set = false;
	}
	else
	{
		/* outer radius is set, */
		lgR2set = true;

		/* log or linear option is still in place */
		if( lgRLog )
		{
			a += convl;
		}
		else
		{
			/* linear radius - convert to log but first make sure that a is > 0 */
			if( a > 0. )
			{
				a = log10(a) + convl;
			}
			else
			{
				fprintf(ioQQQ,"The second radius is negative and linear is set - this is impossible.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		if( a > 37. || a < -37. )
		{
			fprintf(ioQQQ,"WARNING - the log of the second radius is %e - this is too big.\n", a );
			/* flush buffers since we shall soon throw an fpe */
			fflush( ioQQQ );
		}
		a = exp10(a);
		/* check whether it was thickness or outer radius,
		 * we want thickness to be total thickness of modeled region,
		 * NOT outer radius */
		if( a > radius.Radius )
			iterations.StopThickness[0] = a - radius.Radius;
		else
			iterations.StopThickness[0] = a;

		for( long int i=1; i < iterations.iter_alloc; i++ )
		{
			iterations.StopThickness[i] = iterations.StopThickness[0];
		}
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;

		/* flag saying second outer radius or thickness was set */
		if( lgR2set )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "RADIUS %f depth or outer R %f LOG" );
			optimize.nvarxt[optimize.nparm] = 2;
			/* second number is thickness or outer radius */
			optimize.vparm[1][optimize.nparm] = (realnum)log10(a);
			fprintf(ioQQQ,
				" WARNING - outer radius or thickness was set with a variable radius.\n");
			fprintf(ioQQQ,
				" The interpretation of the second number can change from radius to depth as radius changes.\n");
			fprintf(ioQQQ,
				" Do not use the second parameter unless you are certain that you know what you are doing.\n");
			fprintf(ioQQQ,
				" Consider using the STOP THICKNESS or STOP RADIUS command instead.\n");
		}
		else
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "RADIUS= %f LOG" );
			optimize.nvarxt[optimize.nparm] = 1;
		}

		/* log of radius is first number */
		optimize.vparm[0][optimize.nparm] = (realnum)log10(radius.Radius);
		optimize.vincr[optimize.nparm] = 0.5;
		++optimize.nparm;
	}
	return;
}
