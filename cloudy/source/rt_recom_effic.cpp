/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_recom_effic generate escape probability function for continua, */
#include "cddefines.h"
#include "rt.h"
#include "rfield.h"
#include "phycon.h"
#include "opacity.h"
#include "rt_escprob.h"
#include "cosmology.h"

double RT_recom_effic(long int ip)
{
	long int i;
	double dEner, 
	  denom, 
	  escin, 
	  escout, 
	  hnukt, 
	  receff_v, 
	  sum, 
	  tin, 
	  tout;

	DEBUG_ENTRY( "RT_recom_effic()" );

	/* escape probability function for continua,
	 * formally correct for photoelectric absorption only */

	ASSERT( ip > 0 && ip <= rfield.nflux_with_check );

	if( ip > rfield.nflux )
	{
		/* >>chng 01 dec 18, return had been zero, but this did not
		 * work for case where gas much hotter than continuum, as in a
		 * coronal plasma.  change to return of unity */
		receff_v = 1.;
		return( receff_v );
	}

	/* bug in following statement unocvered June 93 S. Schaefer */
	hnukt = TE1RYD*rfield.anu(ip-1)/phycon.te;

	/* rfield.chDffTrns = "OU2" by default */
	/* inward optical depth and escape prob */
	if( strcmp(rfield.chDffTrns,"OSS") == 0 )
	{
		/* which side of Lyman limit is this? */
		if( rfield.anu(ip) > 0.99 )
		{
			/* this is a simple OTS, turned on with DIFFUSE OTS SIMPLE */
			receff_v = SMALLFLOAT;
		}
		else
		{
			receff_v = 1.;
		}
	}
	else if( strcmp(rfield.chDffTrns,"OTS") == 0 )
	{
		tin = opac.TauAbsGeo[0][ip-1];
		if( tin < 5. )
		{
			escin = esccon(tin,hnukt);
		}
		else
		{
			escin = 1e-4;
		}

		/* outward optical depth */
		tout = opac.TauAbsGeo[1][ip-1] - tin;

		if( iteration > 1 )
		{
			/* check whether we have overrun the optical depth scale */
			if( tout > 0. )
			{
				/* good optical depths in both directions, take mean */
				if( tout < 5. )
				{
					escout = esccon(tout,hnukt);
				}
				else
				{
					escout = 1e-4;
				}
				receff_v = 0.5*(escin + escout);
			}
			else
			{
				/* >>chng 91 apr add logic to prevent big change in
				 * esc prob, resulting in terminal oscillations, when optical
				 * depth scale overrun
				 * tau was negative, use 5% of inward optical depth */
				escout = esccon(tin*0.05,hnukt);
				receff_v = 0.5*(escin + escout);
			}
		}
		else
		{
			receff_v = escin;
		}
	}
	else if( strcmp(rfield.chDffTrns,"OU1") == 0 )
	{
		receff_v = opac.ExpZone[ip+1];
	}
	else if( strcmp(rfield.chDffTrns,"OU2") == 0 )
	{
		/* this is the default rt method, as set in zero
		 * E2TauAbsFace is optical depth to illuminated face */
		receff_v = opac.E2TauAbsFace[ip+1];
	}
	else if( strcmp(rfield.chDffTrns,"OU3") == 0 )
	{
		receff_v = 1.;
	}
	else if( strcmp(rfield.chDffTrns,"OU4") == 0 )
	{
		/* this cannot happen, was the former outward treat
		 * optical depth for this zone */
		if( rfield.ContBoltz[ip-1] > 0. )
		{
			i = ip;
			dEner = phycon.te/TE1RYD*0.5;
			sum = 0.;
			denom = 0.;
			while( rfield.ContBoltz[i-1] > 0. &&
			   rfield.anu(i-1)-rfield.anu(ip-1) < (realnum)dEner &&
			   i <= rfield.nflux )
			{
				sum += rfield.ContBoltz[i-1]*opac.tmn[i-1];
				denom += rfield.ContBoltz[i-1];
				i += 1;
			}
			receff_v = sum/denom;
		}
		else
		{
			receff_v = opac.tmn[ip-1];
		}
	}
	else if( strcmp(rfield.chDffTrns,"SOB") == 0 )
	{
		long int ipRecombEdgeFine = rfield.ipnt_coarse_2_fine[ip];
		double OpacityEffective, EffectiveThickness;
		realnum tau;

		/* find line center opacity - use fine opacity if array indices are OK */
		if( ipRecombEdgeFine>=0 && ipRecombEdgeFine<rfield.nfine && rfield.lgOpacityFine )
		{
			/* use fine opacities fine grid fine mesh to get optical depth 
			 * to continuum source */
			/* total line center optical depth, all overlapping lines included */
			OpacityEffective = rfield.fine_opac_zone[ipRecombEdgeFine];
		}
		else
		{
			OpacityEffective =  opac.opacity_abs[ip];
		}

		if( cosmology.lgDo )
		{
			/* dv/dr (s-1), equal to dv/dt / v */
			/* in this case, dv/dr is just the Hubble factor */
			realnum dvdr = GetHubbleFactor(cosmology.redshift_current);
			EffectiveThickness = SPEEDLIGHT / dvdr;
			tau = (realnum)(OpacityEffective * EffectiveThickness);
		}
		else
			tau = opac.taumin;

		tau = MAX2((double)opac.taumin,tau);

		ASSERT( tau >= 0. );

		if( tau < 1e-5 )
			receff_v = (1. - tau/2.);
		else
			receff_v = (1. - dsexp( tau ) )/ tau;
		ASSERT( receff_v >= 0.f );
		ASSERT( receff_v <= 1.f );
	}
	else
	{
		fprintf( ioQQQ, " RECEFF does not understand the transfer method=%3.3s\n", 
		  rfield.chDffTrns );
		cdEXIT(EXIT_FAILURE);
	}

	receff_v = MAX2((double)opac.otsmin,receff_v);
	/* can get epsilon above unity on cray */
	receff_v = MIN2(1.,receff_v);
	return( receff_v );
}
