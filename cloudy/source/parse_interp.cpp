/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseInterp parse parameters on interpolate command */
#include "cddefines.h"
#include "called.h"
#include "rfield.h"
#include "trace.h"
#include "input.h"
#include "parser.h"

void ParseInterp(Parser &p)
{
	DEBUG_ENTRY( "ParseInterp()" );

	/* 
	 * this sub reads in the "interpolate" command, and has
	 * logic for the "continue" lines as well.
	 * OUTPUT:
	 * TNU is vector of energies where the grid is defined
	 * TSLOP initially is vector of log fnu at each freq
	 * converted into slopes here too
	 */

	if( rfield.nShape >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many spectra entered.  Increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	// we reserve space for 1000 table points, no worries if that is too small though...
	ASSERT( rfield.tNu[rfield.nShape].empty() );
	rfield.tNu[rfield.nShape].reserve( 1000 );
	ASSERT( rfield.tslop[rfield.nShape].empty() );
	rfield.tslop[rfield.nShape].reserve( 1000 );
	ASSERT( rfield.tFluxLog[rfield.nShape].empty() );
	rfield.tFluxLog[rfield.nShape].reserve( 1000 );

	strcpy( rfield.chSpType[rfield.nShape], "INTER" );

	/* read all of the numbers on a line */
	long npairs = 0;

	/* this is flag saying that all numbers are in */
	bool lgDONE = false;

	/* this flag says we hit end of command stream */
	p.m_lgEOF = false;
	while( !lgDONE && !p.m_lgEOF )
	{
		/* keep scanning numbers until we hit eol for current line image */
		while( true )
		{
			Energy E( p.FFmtRead() );
			realnum FluxLog = (realnum)p.FFmtRead();
			if( p.lgEOL() )
				break;
			rfield.tNu[rfield.nShape].push_back( E );
			rfield.tFluxLog[rfield.nShape].push_back( FluxLog );
			++npairs;
		}

		do
		{
			/* read a new line, checking for EOF */
			p.getline();
			
			/* option to ignore all comment lines. 
			 * >>chng 06 sep 04 use routine to check for comments */
		}
		while( !p.m_lgEOF && p.isComment());

		/* print the line, but only if it is a continue line */
		if( called.lgTalk && p.hasCommand("CONTINUE") )
		{
			p.echo();
		}

		/* is this a continue line? */
		if( ! p.hasCommand("CONTINUE") )
		{
			/* we have a line but next command, not continue */
			lgDONE = true;
		}

		/* this is another way to hit end of input stream - blank lines */
		if( p.last() )
			p.m_lgEOF = true;
	}

	/* if valid next line, backup one line */
	if( lgDONE )
	{
		--input.nRead;
	}

	/* done reading all of the possible lines */
	if( npairs < 2 )
	{
		fprintf( ioQQQ, "There must be at least 2 pairs to interpolate,\nSorry\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( rfield.tNu[rfield.nShape][0].Ryd() == 0. )
	{
		/* special case - if first energy is zero then it is low energy */
		if( rfield.tNu[rfield.nShape][1].Ryd() > 0. )
		{
			/* second energy positive, assume linear Ryd */
			rfield.tNu[rfield.nShape][0].set(rfield.emm());
		}
		else if( rfield.tNu[rfield.nShape][1].Ryd() < 0. )
		{
			/* second energy negative, assume log of Ryd */
			rfield.tNu[rfield.nShape][0].set(log10(rfield.emm()));
		}
		else
		{
			/* second energy zero, not allowed */
			fprintf( ioQQQ, 
				"An energy of zero was entered for element%3ld in INTERPOLATE and is not allowed.\nSorry\n", 
			  rfield.nShape );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* convert from log(Hz) to Ryd if first nu>5 */
	if( rfield.tNu[rfield.nShape][0].Ryd() >= 5. )
	{
		for( long i=0; i < npairs; i++ )
		{
			if( rfield.tNu[rfield.nShape][i].Ryd()>300.)
			{
				fprintf(ioQQQ,"DISASTER in INTERPPLATE command syntax\n");
				fprintf(ioQQQ,"log of photon energy in Hz was entered, but the value, %g is too large.\nSorry.\n",
						rfield.tNu[rfield.nShape][i].Ryd());
				cdEXIT(EXIT_FAILURE);
			}
			rfield.tNu[rfield.nShape][i].set( 
				exp10(rfield.tNu[rfield.nShape][i].Ryd())/FR1RYD);
		}
	}
	else if( rfield.tNu[rfield.nShape][0].Ryd() < 0. )
	{
		/* energies entered as logs */
		for( long i=0; i < npairs; i++ )
		{
			rfield.tNu[rfield.nShape][i].set(
				exp10((double)rfield.tNu[rfield.nShape][i].Ryd()));
		}
	}
	else
	{
		/* numbers are linear Rydbergs */
		for( long i=0; i < npairs; i++ )
		{
			if( rfield.tNu[rfield.nShape][i].Ryd() == 0. )
			{
				fprintf( ioQQQ, "An energy of zero was entered for pair %3ld in INTERPOLATE and is not allowed.\nSorry\n",
				  i );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* option to print debugging file then exit */
	{
		/* following should be set true to print information */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			for( long i=0; i < npairs; i++ )
			{
				fprintf(ioQQQ,"%.4e\t%.3e\n",
					rfield.tNu[rfield.nShape][i].Ryd(),
					exp10((double)rfield.tFluxLog[rfield.nShape][i]) * rfield.tNu[rfield.nShape][i].Ryd());
			}
			cdEXIT(EXIT_SUCCESS);
		}
	}

	for( long i=0; i < npairs-1; i++ )
	{
		/* check that frequencies are monotonically increasing */
		if( rfield.tNu[rfield.nShape][i+1].Ryd() <= rfield.tNu[rfield.nShape][i].Ryd() )
		{
			fprintf( ioQQQ, "The energies MUST be in increasing order.  Energy #%3ld=%10.2e Ryd was greater than or equal to the next one.\nSorry.\n", 
				i, rfield.tNu[rfield.nShape][i].Ryd() );
			cdEXIT(EXIT_FAILURE);
		}

		/* TFAC is energy, and TSLOP is slope in f_nu; not photons */
		rfield.tslop[rfield.nShape].push_back(
			(realnum)((rfield.tFluxLog[rfield.nShape][i+1]-rfield.tFluxLog[rfield.nShape][i])/
				  log10(rfield.tNu[rfield.nShape][i+1].Ryd()/rfield.tNu[rfield.nShape][i].Ryd()))
			);
	}
	rfield.tslop[rfield.nShape].push_back( 0. );

	/* now check that array is defined over all energies */
	if( rfield.tNu[rfield.nShape][0].Ryd() > rfield.emm() )
	{
		/* not defined over low energy part of array */
		fprintf( ioQQQ, 
			"\n NOTE The incident continuum was not defined over the entire energy range. Some energies are set to zero.\n" );
		fprintf( ioQQQ, 
			" NOTE You may be making a BIG mistake.\n\n" );
	}

	/* check on IR */
	if( rfield.tNu[rfield.nShape][0].Ryd() > rfield.emm() )
		rfield.lgMMok = false;

	if( rfield.tNu[rfield.nShape][0].Ryd() > 1/36. )
		rfield.lgHPhtOK = false;

	/* gamma ray, EGAMRY is roughly 100MeV */
	if( rfield.tNu[rfield.nShape][npairs-1].Ryd() < rfield.egamry() )
		rfield.lgGamrOK = false;

	/* EnerGammaRay is roughly 100keV; high is gamma ray */
	if( rfield.tNu[rfield.nShape][npairs-1].Ryd() < rfield.EnerGammaRay )
		rfield.lgXRayOK = false;

	/* renormalize continuum so that flux we will interpolate upon passes through unity
	 * at near 1 Ryd.  but first we must find 1 Ryd in the array.
	 * find 1 Ryd, npairs is one less than number of continuum pairs */
	/* If no flux is defined at 1 Ryd, use the nearest endpoint of the supplied spectrum */
	long n;
	for( n=0; n < npairs; n++ )
	{
		if( rfield.tNu[rfield.nShape][n].Ryd() > 1. )
			break;
	}
	/* if present, n is now the table point where rfield.tNuRyd[n-1] is <= 1 ryd,
	 * and rfield.tNuRyd[n] > 1 ryd; max(n-1,0) is the nearest endpoint otherwise */
	realnum fac = rfield.tFluxLog[rfield.nShape][max(n-1,0)];

	for( long i=0; i < npairs; i++ )
	{
		rfield.tFluxLog[rfield.nShape][i] -= fac;
		/*fprintf(ioQQQ,"DEBUG parse %li %e \n", i , rfield.tFluxLog[rfield.nShape][i] );*/
	}

	/* option to print out results at this stage - "trace continuum" */
	if( trace.lgConBug && trace.lgTrace )
	{
		fprintf( ioQQQ, " Table for this continuum;\ni\tTNU\tTFAC\tTSLOP, npairs=%li\n",
			 npairs );
		for( long i=0; i < npairs-1; i++ )
		{
			fprintf( ioQQQ, "%li\t%.4e\t%.4e\t%.4e\n", 
				 i , rfield.tNu[rfield.nShape][i].Ryd(), 
				 rfield.tFluxLog[rfield.nShape][i], rfield.tslop[rfield.nShape][i] );
		}
		fprintf( ioQQQ, "%li\t%.4e\t%.4e\n", 
			 npairs , rfield.tNu[rfield.nShape][npairs].Ryd(), 
			 rfield.tFluxLog[rfield.nShape][npairs]);
	}

	/* finally check that we are within dynamic range of this cpu */
	double cmin = log10( FLT_MIN );
	double cmax = log10( FLT_MAX );
	bool lgHit = false;
	for( long i=0; i < npairs; i++ )
	{
		if( rfield.tFluxLog[rfield.nShape][i] <= cmin || 
		    rfield.tFluxLog[rfield.nShape][i] >= cmax )
		{
			fprintf(ioQQQ,
				" The log of the flux specified in interpolate pair %li is not within dynamic range of this CPU - please rescale.\n",i);
			fprintf(ioQQQ,
				" The photon energy is %f Rydberg and the log of the flux is %f.\n\n",
				rfield.tNu[rfield.nShape][i].Ryd() , 
				rfield.tFluxLog[rfield.nShape][i]);
			lgHit = true;
		}
	}
	if( lgHit )
	{
		fprintf(ioQQQ,"\n NOTE The log of the flux given in an interpolate command is outside the range of this cpu.\n");
		fprintf(ioQQQ," NOTE I will try to renormalize it to be within the range of this cpu, but if I crash, this is a likely reason.\n");
		fprintf(ioQQQ," NOTE Note that the interpolate command only is used for the shape of the continuum.\n");
		fprintf(ioQQQ," NOTE The order of magnitude of the flux is not used in any way.\n");
		fprintf(ioQQQ," NOTE For safety this could be of order unity.\n\n");
	}

	rfield.ncont[rfield.nShape] = npairs;

	ASSERT( rfield.tNu[rfield.nShape].size() == (size_t)npairs );
	ASSERT( rfield.tFluxLog[rfield.nShape].size() == (size_t)npairs );
	ASSERT( rfield.tslop[rfield.nShape].size() == (size_t)npairs );

	/* now increment number of continua */
	++rfield.nShape;

	return;
}
