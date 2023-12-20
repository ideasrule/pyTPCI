/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ffun evaluate total flux for sum of all continuum sources */
/*ffun1 derive flux at a specific energy, for one continuum */
/*ReadTable called by TABLE READ to read in continuum from PUNCH TRANSMITTED CONTINUUM */
#include "cddefines.h"
#include "rfield.h"
#include "continuum.h"

double PlanckFunction( double Temp, double E_Ryd );

/* evaluate sum of all individual continua at one energy, return is
* continuum intensity */
double ffun(
			/* the energy in Rydbergs where the continuum will be evaluated */
			double anu )
{
	double frac_beam_time;
	/* fraction of beamed continuum that is constant */
	double frac_beam_const;
	/* fraction of continuum that is isotropic */
	double frac_isotropic;
	double a;

	DEBUG_ENTRY( "ffun()" );

	/* call real function */
	a = ffun( anu , &frac_beam_time , &frac_beam_const , &frac_isotropic );
	return a;
}

/* evaluate sum of all individual continua at one energy, return is
 * continuum intensity */
double ffun(
	/* the energy in Rydbergs where the continuum will be evaluated */
	double anu , 
	/* fraction of beamed continuum that is varies with time */
	double *frac_beam_time,
	/* fraction of beamed continuum that is constant */
	double *frac_beam_const,
	/* fraction of continuum that is isotropic */
	double *frac_isotropic )
{
	double ffun_v;
	static bool lgWarn = false;
	double flx_beam_time , flx_beam_const , flx_isotropic;

	DEBUG_ENTRY( "ffun()" );

	/* This routine, ffun, returns the sum of photons per unit time, area, energy,
	 * for all continua in the calculation.  ffun1 is called and returns 
	 * a single continuum.  We loop over all nShape continuum sources - 
	 * ipspec points to each */
	ffun_v = 0.;
	flx_beam_time = 0.;
	flx_beam_const = 0.;
	flx_isotropic = 0.;
	for( rfield.ipSpec=0; rfield.ipSpec < rfield.nShape; rfield.ipSpec++ )
	{
		double one = ffun1(anu)*rfield.spfac[rfield.ipSpec];
		ffun_v += one;

		/* find fraction of total that is constant vs variable and 
		 * isotropic vs beamed */
		if( rfield.lgBeamed[rfield.ipSpec] )
		{
			if( rfield.lgTimeVary[rfield.ipSpec] )
				flx_beam_time += one;
			else
				flx_beam_const += one;
		}
		else
			flx_isotropic += one;
	}

	/* at this point rfield.flux is the sum of the following three continua
	 * now keep track of the three different types */
	if( ffun_v < SMALLDOUBLE )
	{
		*frac_beam_time = 0.;
		*frac_beam_const = 1.;
		*frac_isotropic = 0.;
		ffun_v = 0.;
	}
	else
	{
		/* fraction of beamed continuum that varies with time */
		*frac_beam_time = flx_beam_time / ffun_v;
		/* part of beamed continuum that is constant */
		*frac_beam_const = flx_beam_const / ffun_v;
		/* the constant isotropic continuum */
		*frac_isotropic = flx_isotropic / ffun_v;
	}
	ASSERT( *frac_beam_time >=0. && *frac_beam_time<=1.+3.*DBL_EPSILON );
	ASSERT( *frac_beam_const >=0.&& *frac_beam_const<=1.+3.*DBL_EPSILON );
	ASSERT( *frac_isotropic >=0. && *frac_isotropic<=1.+3.*DBL_EPSILON );
	ASSERT( fabs( 1.-*frac_beam_time-*frac_beam_const-*frac_isotropic)<
		10.*DBL_EPSILON);

	if( ffun_v > BIGFLOAT && !lgWarn )
	{
		lgWarn = true;
		fprintf( ioQQQ, " FFUN:  The net continuum is very intense.\n" );
		fprintf( ioQQQ, " I will try to press on, but may have problems.\n" );
	}
	return ffun_v;
}

/*ffun1 derive flux at a specific energy, for one continuum */
double ffun1(double xnu)
{
	char chKey[6];
	long int i;
	double fac, 
	  ffun1_v, 
	  y;
	static bool lgWarn = false;

	DEBUG_ENTRY( "ffun1()" );


	/* confirm that pointer is within range */
	ASSERT( rfield.ipSpec >= 0 );
	ASSERT( rfield.ipSpec < rfield.nShape );

	/* FFUN1 returns photons per unit time, area, energy, for one continuum
	 * ipspec is the pointer to the continuum source, in the order
	 * entered on the command lines */

	/*begin sanity check */
	ASSERT( xnu > rfield.emm() );
	ASSERT( xnu < rfield.egamry() );
	/*end sanity check */

	strcpy( chKey, rfield.chSpType[rfield.ipSpec] );

	if( strcmp(chKey,"AGN  ") == 0 )
	{
		/* power law with cutoff at both ends
		 * nomenclature all screwed up - slope is cutoff energy in Ryd,
		 * cutoff[1][i] is ratio of two continua from alpha ox
		 * cutoff[2][i] is slope of bb temp */
		ffun1_v = pow(xnu,-1. + rfield.cutoff[rfield.ipSpec][1])*
		  sexp(xnu/rfield.slope[rfield.ipSpec])*sexp(0.01/
		  xnu);
		/* only add on x-ray for energies above 0.1 Ryd */
		if( xnu > 0.1 )
		{
			if( xnu < 7350. )
			{
				/* cutoff is actually normalization constant
				 * below 100keV continuum is nu-1 */
				ffun1_v += pow(xnu,rfield.cutoff[rfield.ipSpec][2] - 
				  1.)*rfield.cutoff[rfield.ipSpec][0]*sexp(1./
				  xnu);
			}
			else
			{
				ffun1_v += pow(7350.,rfield.cutoff[rfield.ipSpec][2] - 
				  1.)*rfield.cutoff[rfield.ipSpec][0]/
				  POW3(xnu/7350.);
			}
		}

	}
	else if( strcmp(chKey,"POWER") == 0 )
	{
		/* power law with cutoff at both ends */
		ffun1_v = pow(xnu,-1. + rfield.slope[rfield.ipSpec])*
		  sexp(xnu/rfield.cutoff[rfield.ipSpec][0]+rfield.cutoff[rfield.ipSpec][1]/
		  xnu);

	}
	else if( strcmp(chKey,"DISKB") == 0 )
	{
		long numSteps = 100;
		double TempHi, TempLo;
		// Let temps take any values.  If equal, just do blackbody...
		if( fp_equal( rfield.slope[rfield.ipSpec], rfield.cutoff[rfield.ipSpec][0] ) )
		{
			ffun1_v = PlanckFunction( rfield.slope[rfield.ipSpec], xnu );
		}
		else
		{
			if( rfield.slope[rfield.ipSpec] > rfield.cutoff[rfield.ipSpec][0] )
			{
				TempHi = rfield.slope[rfield.ipSpec];
				TempLo = rfield.cutoff[rfield.ipSpec][0];
			}
			else
			{
				TempLo = rfield.slope[rfield.ipSpec];
				TempHi = rfield.cutoff[rfield.ipSpec][0];
			}
			ASSERT( TempLo < TempHi );
			double LogDeltaT = (log10(TempHi) - log10(TempLo))/(numSteps-1.);
			ffun1_v = 0.;
			for( long i=0; i<numSteps; i++ )
			{
				double Temp = exp10( log10(TempHi) - double(i) * LogDeltaT );
				double relativeWeight = powpq( TempHi/Temp, 8, 3 ) * exp10( LogDeltaT );
				ffun1_v += PlanckFunction( Temp, xnu ) * relativeWeight;
			}
		}
	}
	else if( strcmp(chKey,"BLACK") == 0 )
	{
		ffun1_v = PlanckFunction( rfield.slope[rfield.ipSpec], xnu );
	}
	else if( strcmp(chKey,"INTER") == 0 )
	{
		/* interpolate on tabulated input spectrum, factor of 1.0001 to 
		 * make sure that requested energy is within bounds of array */
		if( xnu >= rfield.tNu[rfield.ipSpec][0].Ryd()*1.000001 )
		{
			/* loop starts at second array element, [1], since want to 
			 * find next continuum energy greater than desired point */
			for( i=1; i < rfield.ncont[rfield.ipSpec]; ++i )
			{
				if( xnu < rfield.tNu[rfield.ipSpec][i].Ryd() )
				{
					/* the energy xnu is between points rfield.tNuRyd[rfield.ipSpec][i-1]
					 * and rfield.tNuRyd[rfield.ipSpec][i] - do linear 
					 * interpolation in log log space */
					y = rfield.tFluxLog[rfield.ipSpec][i-1] + 
					  rfield.tslop[rfield.ipSpec][i-1]*
					  log10(xnu/rfield.tNu[rfield.ipSpec][i-1].Ryd());

					/* return value is photon density, div by energy */
					ffun1_v = exp10(y);

					/* this checks that overshoots did not occur - interpolated
					 * value must be between lowest and highest point */
#					ifndef NDEBUG
					double ys1 = MIN2( rfield.tFluxLog[rfield.ipSpec][i-1],rfield.tFluxLog[rfield.ipSpec][i]);
					double ys2 = MAX2( rfield.tFluxLog[rfield.ipSpec][i-1],rfield.tFluxLog[rfield.ipSpec][i]);
					ys1 = exp10( ys1 );
					ys2 = exp10( ys2 );
					ASSERT( ffun1_v >= ys1/(1.+100.*FLT_EPSILON) );
					ASSERT( ffun1_v <= ys2*(1.+100.*FLT_EPSILON) );
#					endif
					/* return value is photon density, div by energy */
					return( ffun1_v/xnu );
				}
			}
			/* energy above highest in table */
			ffun1_v = 0.;
		}
		else
		{
			/* energy below lowest on table */
			ffun1_v = 0.;
		}
	}
	else if( strcmp(chKey,"BREMS") == 0 )
	{
		/* brems continuum, rough gaunt factor */
		fac = TE1RYD*xnu/rfield.slope[rfield.ipSpec];
		ffun1_v = sexp(fac)/pow(xnu,1.2);

	}
	else if( strcmp(chKey,"LASER") == 0 )
	{
		const double BIG = 1.e10;
		/* a laser, mostly one frequency */
		/* >>chng 01 jul 01, was hard-wired 0.05 rel frac, change to optional
		 * second parameter, with default of 0.05 */
		/*if( xnu > 0.95*rfield.slope[rfield.ipSpec] && xnu < 
		  1.05*rfield.slope[rfield.ipSpec] )*/
		if( xnu > (1.-rfield.cutoff[rfield.ipSpec][0])*rfield.slope[rfield.ipSpec] && 
			xnu < (1.+rfield.cutoff[rfield.ipSpec][0])*rfield.slope[rfield.ipSpec] )
		{
			ffun1_v = BIG;
		}
		else
		{
			ffun1_v = 0.;
		}

	}
	else if( strcmp(chKey,"READ ") == 0 )
	{
		/* use array of values read in on TABLE READ command */
		// this table is already on the Cloudy mesh, so no need to interpolate...
		// this table is not defined for rfield.anu(rfield.nflux) (i.e. the cell for the unit test)
		// we return 0 in that case since the flux value in that cell is irrelevant anyway
		if( xnu >= rfield.anumax(rfield.nflux-1) )
		{
			ffun1_v = 0.;
		}
		else
		{
			i = rfield.ipointC(xnu);
			ASSERT( rfield.ncont[rfield.ipSpec] == rfield.nflux );
			ffun1_v = exp10((double)rfield.tFluxLog[rfield.ipSpec][i])/rfield.anu(i);
		}
	}
	else if( strcmp(chKey,"VOLK ") == 0 )
	{
		/* use array of values from the TABLE STAR or TABLE HM12 command */
		// this table is already on the Cloudy mesh, so no need to interpolate...
		i = rfield.ipointC(xnu);
		ASSERT( rfield.ncont[rfield.ipSpec] == rfield.nflux_with_check );
		/* bug fixed Jul 9 93: FFUN1 = TSLOP(IPSPEC,I) / ANU(I) / ANU(I) */
		ffun1_v = rfield.tslop[rfield.ipSpec][i]/rfield.anu(i);
	}
	else
	{
		fprintf( ioQQQ, " ffun1: I do not understand continuum label \"%s\" for continuum %li.\n", 
		  chKey , rfield.ipSpec);
		cdEXIT(EXIT_FAILURE);
	}

	if( ffun1_v > 1e35 && !lgWarn )
	{
		lgWarn = true;
		fprintf( ioQQQ, " FFUN1:  Continuum %ld is very intense.\n", 
		  rfield.ipSpec );
		fprintf( ioQQQ, " I will try to press on, but may have problems.\n" );
	}
	return ffun1_v;
}

double PlanckFunction( double Temp, double E_Ryd )
{
	const double db_log = log(DBL_MAX);
	double ffun1_v;
	double fac;

	/* black body */
	fac = TE1RYD*E_Ryd/Temp;
	/* >>>chng 00 apr 13 from 80 to log(dbl_max) */
	if( fac > db_log )
	{
		ffun1_v = 0.;
	}
	else
	{
		ffun1_v = E_Ryd*E_Ryd/expm1(fac);
	}

	return ffun1_v;
}
/*outsum sum outward continuum beams */
void outsum(double *outtot, double *outin, double *outout)
{
	long int i;

	DEBUG_ENTRY( "outsum()" );

	*outin = 0.;
	*outout = 0.;
	for( i=0; i < rfield.nflux; i++ )
	{
		/* N.B. in following en1ryd prevents overflow */
		*outin += rfield.anu(i)*(rfield.flux[0][i]*EN1RYD);
		*outout += rfield.anu(i)*(rfield.outlin[0][i] + rfield.outlin_noplot[i] +rfield.ConInterOut[i])*
		  EN1RYD;
	}

	*outtot = *outin + *outout;
	return;
}
