/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ZoneEnd last routine called after all zone calculations, before iter_end_check,
 * upon exit radiation field is for outer edge of current zone */
/*ZoneStart set variables that change with each zone, like radius, depth,
 * upon exit flux will be flux at center of zone about to be computed */
#include "cddefines.h"
#include "zones.h"
#include "phycon.h"
#include "opacity.h"
#include "rfield.h"
#include "struc.h"
#include "thermal.h"
#include "dense.h"
#include "h2.h"
#include "geometry.h"
#include "iterations.h"
#include "conv.h"
#include "save.h"
#include "dynamics.h"
#include "radius.h"
#include "doppvel.h"
#include "mole.h"
#include "cosmology.h"
/* this is number of zones to include in guess for next temperatures */
#define	IOFF	3

void ZoneStart(const char *chMode)
{
	bool lgNeOscil,
		lgTeOscil;
	long int i;

	double 
	  /* this is correction factor for dilution of radiation
	   * across current zone, set in ZoneStart to correct
		* flux for center of current zone, and undone in ZoneDone */
	  DilutionCorrec ,
		drFac ,
		dTauThisZone,
	  outwrd, 
	  ratio, 
	  rm, 
	  rad_middle_zone,
	  vin, 
	  vout;

	static double DTaver , DEaver, 
	  dt2 , de2;

	DEBUG_ENTRY( "ZoneStart()" );

	/** \todo 2 is this the best place for this? */
	/* this is a turbulent velocity power law.  */
	if( DoppVel.lgTurbLawOn )
	{
		DoppVel.TurbVel = DoppVel.TurbVelZero * 
			pow( (realnum)(radius.Radius/radius.rinner), DoppVel.TurbVelLaw );
	}

	/* this sub can be called in two ways, 'incr' means increment
	 * radius variables, while 'init' means only initialize rest */
	/* called at start of a zone to update all variables
	 * having to do with starting this zone */

	/* first establish current filling factor (normally 1) since used within
	 * following branches */
	geometry.FillFac = (realnum)(geometry.fiscal*
		pow( radius.Radius/radius.rinner ,(double)geometry.filpow));
	geometry.FillFac = (realnum)MIN2(1.,geometry.FillFac);

	if( strcmp(chMode,"init") == 0 )
	{
		/* initialize the variables at start of calculation, after this exits
		 * radius is equal to the initial radius, not outer edge of zone */

		/* depth to illuminated face */
		radius.depth = 1e-30;

		/* integral of depth times filling factor */
		radius.depth_x_fillfac = radius.depth*geometry.FillFac;
		/* effective radius */
		radius.drad_x_fillfac = radius.drad*geometry.FillFac;

		/* reset radius to inner radius, at this point equal to illuminated face radius */
		radius.Radius = radius.rinner;
		radius.Radius_mid_zone = radius.rinner + radius.dRadSign*radius.drad/2.;

		/* thickness of next zone */
		radius.drNext = radius.drad;

		/* depth to illuminated face */
		radius.depth_mid_zone = radius.drad/2.;

		/* depth variable for globule */
		radius.glbdst = radius.glbrad;

		/* this is case where outer radius is set */
		if( iterations.StopThickness[iteration-1] < 5e29 )
		{
			radius.Depth2Go = iterations.StopThickness[iteration-1];
		}
		else if( iteration > 1 )
		{
			/* this case second or later iteration, use previous thickness */
			radius.Depth2Go = iterations.StopThickness[iteration-2];
		}
		else
		{
			/* first iteration and depth not set, so we cannot estimate it */
			radius.Depth2Go = 0.;
		}
	}
	else if( strcmp(chMode,"incr") == 0 )
	{
		/* update radius variables - called by cloudy at start of this zone's calcs */
		radius.drad_mid_zone = (radius.drad+radius.drNext)/2.;
		radius.drad = radius.drNext;
		/* >>chng 06 mar 21, remember largest and smallest dr over this iteration
		 * for use in recombination time dep logic */
		radius.dr_min_last_iter = MIN2( radius.dr_min_last_iter , radius.drNext );
		radius.dr_max_last_iter = MAX2( radius.dr_max_last_iter , radius.drNext );

		ASSERT( radius.drad > 0. );
		radius.depth += radius.drad;
		radius.depth_mid_zone = radius.depth - radius.drad/2.;

		/* effective radius */
		radius.drad_x_fillfac = radius.drad*geometry.FillFac;

		/* integral of depth times filling factor */
		radius.depth_x_fillfac += radius.drad_x_fillfac;

		/* decrement Depth2Go but do not let become negative */
		radius.Depth2Go = MAX2( 0.,radius.Depth2Go - radius.drad );

		/* increment the radius, during the calculation Radius is the
		 * outer radius of the current zone*/
		radius.Radius += radius.drad*radius.dRadSign;

		/* Radius is now outer edge of this zone since incremented above, 
		 * so need to add drad to it */
		radius.Radius_mid_zone = radius.Radius - radius.dRadSign*radius.drad/2.;

		/***********************************************************************
		 *
		 * attenuate rfield to center of this zone
		 *
		 ***********************************************************************/

		 /* radius was just incremented above, so this resets continuum to
		  * flux at center of zone we are about to compute.  radius will be
		  * incremented immediately following this.  this correction is undone
		  * when ZoneDone called */

		/* this will be the optical thickness of the next zone
		 * AngleIllumRadian is 1/COS(theta), is usually 1, reset with illuminate command,
		 * option for illumination of slab at an angle */

		/* radius.dRNeff is next dr with filling factor, this will only be
		 * used to get approximate correction for attenuation
		 * of radiation in this zone,
		 * equations derived in hazy, Netzer&Ferland 83, for factor accounting
		 * any changes here should be checked with both sphonly.in and pphonly*/
		/* honlyotspp seems most sensitive to the 1.35 */
		drFac = radius.drad*geometry.FillFac*geometry.DirectionalCosin*1.35;

		if( cosmology.lgDo )
			drFac = 0.;

		/* dilute radiation so that rfield is in center of zone, this
		 * goes into tmn at start of zone here, but is later divided out
		 * when ZoneEnd called */
		DilutionCorrec = 1./POW2(
			(radius.Radius-radius.dRadSign*radius.drad/2.)/(radius.Radius-radius.dRadSign*radius.drad) );

		/* note that this for loop is through <= nflux, to include the [nflux]
		 * unit integration verification token.  The token was set in
		 * in MadeDiffuse and propagated out in metdif.  the opacity at that energy is
		 * zero, so only the geometrical part of tmn will affect things.  The final
		 * sum is verified in PrtComment */
		for( i=0; i <= rfield.nflux; i++ )
		{
			dTauThisZone = opac.opacity_abs[i]*drFac;
			/* TMN is array of scale factors which account for attenuation
			 * of continuum across the zone (necessary to conserve energy
			 * at the 1% - 5% level.) sphere effects in
			 * drNext was set by NEXTDR and will be next dr */

			if( dTauThisZone < 1e-4 )
			{
				/* this small tau limit is the most likely branch, about 60% in parispn.in */
				opac.tmn[i] = 1.f;
			}
			else if( dTauThisZone < 5. )
			{
				/* this middle tau limit happens in the remaining 40% */
				opac.tmn[i] = (realnum)((1. - exp(-dTauThisZone))/(dTauThisZone));
			}
			else
			{
				/* this happens almost not at all */
				opac.tmn[i] = (realnum)(1./(dTauThisZone));
			}

			/* now add on correction for dilution across this zone */
			opac.tmn[i] *= (realnum)DilutionCorrec;

			rfield.flux_beam_const[i] *= opac.tmn[i];
			rfield.flux_beam_time[i] *= opac.tmn[i];
			rfield.flux_isotropic[i] *= opac.tmn[i];
			rfield.flux[0][i] = rfield.flux_beam_const[i] + rfield.flux_beam_time[i] +
				rfield.flux_isotropic[i];

			/* >>chng 03 nov 08, update SummedCon here since flux changed */
			rfield.SummedCon[i] = rfield.flux[0][i] + rfield.SummedDif[i];
			rfield.SummedOcc[i] = rfield.SummedCon[i]*rfield.convoc[i];
		}

		/* following is distance to globule, usually does not matter */
		radius.glbdst -= (realnum)radius.drad;

		/* if gt 5th zone, and not constant pressure, and not oscillating te
		 * then guess next temperature : option to predict next temperature
		 * NZLIM is dim of structure variables saving temp, do data if nzone>NZLIM
		 * IOFF is number of zones to look over, set set to 3 in the define above */
		/* >>chng 01 mar 12, add option to not do this, set with no tepred command */
		if( (strcmp(dense.chDenseLaw,"CPRE") != 0) && 
			thermal.lgPredNextTe && 
			(nzone > IOFF + 1)  )
		{
			phycon.TeInit = phycon.te;
			phycon.EdenInit = dense.eden;
			lgTeOscil = false;
			lgNeOscil = false;
			dt2 = 0.;
			de2 = 0.;
			DTaver = 0.;
			DEaver = 0.;
			for( i=nzone - IOFF-1; i < (nzone - 1); i++ )
			{
				double dt1 = dt2;
				double de1 = de2;
				/* this will get the average change in temperature for the 
				 * past IOFF zones */
				dt2 = struc.testr[i-1] - struc.testr[i];
				de2 = struc.ednstr[i-1] - struc.ednstr[i];
				DTaver += dt2;
				DEaver += de2;
				if( dt1*dt2 < 0. )
				{
					lgTeOscil = true;
				}
				if( de1*de2 < 0. )
				{
					lgNeOscil = true;
				}
			}

			/* option to guess next electron temperature if no oscillations */
			if( !lgTeOscil )
			{
				DTaver /= (double)(IOFF);
				/* don't want to over correct, use smaller of two */
				dt2 = fabs(dt2);
				rm = fabs(DTaver);
				DTaver = sign(MIN2(rm,dt2),DTaver);
				/* do not let it change by more than 5% of current temp */
				/* >>chng 03 mar 18, from 5% to 1% - convergence is much better
				 * now, and throwing the next Te off by 5% could easily disturb
				 * the solvers - primal.in was a good case */
				DTaver = sign(MIN2(rm,0.01*phycon.te),DTaver);
				/* this actually changes the temperature */
				TempChange(phycon.te - DTaver , true);
			}
			else
			{
				/* temp was oscillating - too dangerous to do anything */
				DTaver = 0.;
			}
			/* this is used to remember the proposed temperature, for output
			 * in the save predictor command */
			phycon.TeProp = phycon.te;

			/* option to guess next electron density if no oscillations */
			if( !lgNeOscil )
			{
				DEaver /= IOFF;
				de2 = fabs(de2);
				rm = fabs(DEaver);
				/* don't want to over correct, use smaller of two */
				DEaver = sign(MIN2(rm,de2),DEaver);
				/* do not let it change by more than 5% of current temp */
				DEaver = sign(MIN2(rm,0.05*dense.eden),DEaver);
				/* this actually changes the temperature */
				EdenChange( dense.eden - DEaver );
			}
			else
			{
				/* temp was oscillating - too dangerous to do anything */
				DEaver = 0.;
			}
			/* this is used to remember the proposed temperature, for output
			 * in the save predictor command */
			phycon.EdenProp = dense.eden;
			/* must call TempChange since ionization has changed, there are some
			* terms that affect collision rates (H0 term in electron collision) */
			TempChange(phycon.te , false);
		}
	}

	else
	{
		fprintf( ioQQQ, " PROBLEM ZoneStart called with insane argument, %4.4s\n", 
		  chMode );
		/* TotalInsanity exits after announcing a problem */
		TotalInsanity();
	}

	/* do advection if enabled */
	if( dynamics.lgAdvection || dynamics.lgTimeDependentStatic )
		DynaStartZone();

	/* clear flag indicating the ionization convergence failures 
	 * have ever occurred in current zone
	conv.lgConvIonizThisZone = false; */
	conv.resetCountersZone();

	/* this will say how many times the large H2 molecule has been called in this zone -
	 * if not called (due to low H2 abundance) then not need to update its line arrays */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom)
		(*diatom)->nCall_this_zone = 0;

	/* check whether zone thickness is too small relative to radius */
	if( strcmp(dense.chDenseLaw,"GLOB") == 0 )
	{
		ratio = radius.drad/(radius.glbdst + radius.drad);
		/* this only matters for globule command */
		if( ratio < 1e-14 )
		{
			radius.lgdR2Small = true;
		}
		else
		{
			radius.lgdR2Small = false;
		}
	}

	/* area factor, ratio of radius to outer edge of this zone
	 * relative to radius of illuminated face of cloud
	 * is unity in intensity case, per unit area */
	radius.r1r0sq = POW2( radius.Radius/radius.rinner );
	if( radius.lgPredLumin )
	{
		if( save.lgLuminosityOld )
		{
			/* use old style, save continuum luminosity per unit cloud area
			 * at inner radius, when SET SAVE LUMINOSITY OLD used */
			radius.PI4_Radius_sq = radius.r1r0sq;
			radius.PI4_rinner_sq = 1.;
		}
		else
		{
			radius.PI4_Radius_sq = PI4 * POW2( radius.Radius );
			radius.PI4_rinner_sq = PI4 * POW2( radius.rinner );
		}
	}
	else
	{
		/* intensity case with possible sphericity.  Report
		 * intensity per unit inner cloud area
		 */
		radius.PI4_Radius_sq = radius.r1r0sq;
		radius.PI4_rinner_sq = 1.;
	}

	/* following only approximate, used for center of next zone */
	radius.dRNeff = radius.drNext*geometry.FillFac;

	/* this is the middle of the zone */
	rad_middle_zone = radius.Radius - radius.dRadSign*radius.drad/2.;

	/* this is used for long slit integrals */
	if( radius.drad/radius.Radius > 0.01 )
	{
		double Ropp = radius.Radius - radius.dRadSign*radius.drad;
		radius.darea_x_fillfac = PI*abs(pow2(radius.Radius) - pow2(Ropp))*geometry.FillFac;
	}					     
	else
		radius.darea_x_fillfac = PI2*rad_middle_zone*radius.drad*geometry.FillFac;

	/* Radius is outer radius of this zone, so radius - drad is inner radius
	 * rinner is inner radius of nebula; dVeffVol is the effective vol element
	 * dVeffVol is vol rel to inner radius, so has units cm
	 * for plane parallel dVeffVol is dReff */
	if( radius.drad/radius.Radius > 0.01 )
	{
		double r1 = radius.Radius - radius.dRadSign*radius.drad;
		double rin_zone = min(r1,radius.Radius);
		double rout_zone = max(r1,radius.Radius);
		vin = pow2(rin_zone/radius.rinner)*rin_zone/3.;
		if( rin_zone > radius.CylindHigh )
		{
			// the volume of the cap of a sphere is given here:
			// http://en.wikipedia.org/wiki/Spherical_cap
			// we need two of these...
			double h = rin_zone-radius.CylindHigh;
			double v2cap = pow2(h/radius.rinner)*(rin_zone - h/3.)/2.;
			vin -= v2cap;
		}
		vout = pow2(rout_zone/radius.rinner)*rout_zone/3.;
		if( rout_zone > radius.CylindHigh )
		{
			double h = rout_zone-radius.CylindHigh;
			double v2cap = pow2(h/radius.rinner)*(rout_zone - h/3.)/2.;
			vout -= v2cap;
		}
		/* this is the usual case the full volume, just the difference in the two vols */
		/* this will become the true vol when multiplied by 4pi*rinner^2 */
		radius.dVeffVol = (vout - vin)*geometry.FillFac;
		/* this is correction for slit projected onto resolved object -
		 * this only happens when aperture command is entered.
		 * when slit and cylinder are both used it is assumed that the slit
		 * is oriented along the rotational axis of the cylinder
		 * the case where the slit is perpendicular to this axis is identical
		 * to the normal spherical case, so can be done by removing the cylinder command
		 * slits oriented at any other angle can be done by dividing the cylinder height
		 * by the cosine of the angle */
		/* default of iEmissPower is 2, set to 0 is just aperture beam, 
		 * and to 1 if aperture long slit set */
		if( geometry.iEmissPower == 2 )
		{
			radius.dVeffAper = radius.dVeffVol;
		}
		else if( geometry.iEmissPower == 1 )
		{
			double ain = (rin_zone/radius.rinner)*rin_zone/2.;
			if( rin_zone > radius.CylindHigh )
			{
				// the area of a circular segment is given here:
				// http://en.wikipedia.org/wiki/Circular_segment
				// we need two of those...
				double Theta = 2.*acos(min(radius.CylindHigh/rin_zone,1.));
				ain *= 1. - max(Theta - sin(Theta),0.)/PI;
			}
			double aout = (rout_zone/radius.rinner)*rout_zone/2.;
			if( rout_zone > radius.CylindHigh )
			{
				double Theta = 2.*acos(min(radius.CylindHigh/rout_zone,1.));
				aout *= 1. - max(Theta - sin(Theta),0.)/PI;
			}
			radius.dVeffAper = (aout - ain)*geometry.FillFac;
		}
		else if( geometry.iEmissPower == 0 )
		{
			radius.dVeffAper = radius.drad*geometry.FillFac;
		}
	}
	else
	{
		/* thin cell limit */
		/* rad_middle_zone is the middle of the zone */
		radius.dVeffVol = (rad_middle_zone/radius.rinner)*radius.drad*geometry.FillFac*
			(min(rad_middle_zone,radius.CylindHigh)/radius.rinner);
		if( geometry.iEmissPower == 2 )
		{
			radius.dVeffAper = radius.dVeffVol;
		}
		else if( geometry.iEmissPower == 1 )
		{
			radius.dVeffAper = (rad_middle_zone/radius.rinner)*radius.drad*geometry.FillFac;
			if( rad_middle_zone > radius.CylindHigh )
			{
				double Theta = 2.*acos(min(radius.CylindHigh/rad_middle_zone,1.));
				double q = sqrt(max(1.-pow2(radius.CylindHigh/rad_middle_zone),0.))*rad_middle_zone;
				radius.dVeffAper *= 1. - max(Theta - sin(Theta),0.)/PI - max(1. - cos(Theta),0.)*radius.CylindHigh/(PI*q);
			}
		}
		else if( geometry.iEmissPower == 0 )
		{
			radius.dVeffAper = radius.drad*geometry.FillFac;
		}
	}

	/* covering factor, default is unity */
	radius.dVeffVol *= geometry.covgeo;
	radius.dVeffAper *= geometry.covaper;

	/* these are needed for line intensities in outward beam
	 * lgSphere set, geometry.covrt usually 1, 0 when not lgSphere
	 * so outward is 1 when lgSphere set 1/2 when not set, */
	outwrd = (1. + geometry.covrt)/2.;

	/*>>>chng 99 apr 23, from above to below */
	/*radius.dVolOutwrd = outwrd*POW2( (radius.Radius-radius.drad_x_fillfac/2.)/radius.Radius) * 
	  radius.drad;*/
	/* this includes covering fact, the effective vol,, and 1/r^2 dilution,
	 * dReff includes filling factor.  the covering factor part is 1 for sphere,
	 * 1/2 for open */
	/*radius.dVolOutwrd = outwrd*radius.Radius*radius.drad_x_fillfac/(radius.Radius + 
	  2.*radius.drad);*/
	/* dReff from above, includes filling factor */
	radius.dVolOutwrd = outwrd*radius.drad_x_fillfac;
	ASSERT( radius.dVolOutwrd > 0. );

	/* following will be used to "uncorrect" for this in lines when predicting continua
	radius.GeoDil = radius.Radius/(radius.Radius + 2.*radius.drad);*/

	/* this should multiply the line to become the net inward.  geo part is 1/2 for
	 * open geometry, 0 for closed.  for case of isotropic emitter only this and dVolOutwrd
	 * above are used */
	radius.dVolReflec = (1. - outwrd)*radius.drad_x_fillfac*radius.r1r0sq;

	if( geometry.lgSphere )
	{
		/* inward beams do not go in when lgSphere set since geometry symmetric */
		radius.BeamInIn = 0.;
		radius.BeamInOut = radius.Radius*radius.drad_x_fillfac/(radius.Radius + 
		  2.*radius.drad);
	}
	else
	{
		radius.BeamInIn = radius.drad_x_fillfac*radius.r1r0sq;

		/* inward beams do not go out */
		radius.BeamInOut = 0.;
	}
	/* factor for outwardly directed part of line */
	radius.BeamOutOut = radius.Radius*radius.drad_x_fillfac/(radius.Radius + 
	  2.*radius.drad);
	return;
}

void ZoneEnd(void)
{
	DEBUG_ENTRY( "ZoneEnd()" );

	/***********************************************************************
	 *
	 *  correct rfield for attenuation from center of zone to inner edge
	 *
	 ***********************************************************************/

	// verify sanity of rfield.nPositive, ==0 is case with no light
	//fprintf(ioQQQ,"DEBUGGG, %li %li %.3e\n",rfield.nPositive, rfield.nflux,
	//		rfield.SummedCon[rfield.nPositive-1]) ;
	if( rfield.nPositive>0 )
	{
		ASSERT( rfield.SummedCon[rfield.nPositive-1]>0. );
		if (rfield.FluxFaint > 0.0)
		{
			for( long i=rfield.nPositive; i<rfield.nflux; ++i )
				ASSERT( rfield.SummedCon[i]==0.);
		}
	}
	//fprintf(ioQQQ,"DEBUGGG, %li %li \n",rfield.nflux,rfield.nflux_with_check) ;
	ASSERT( rfield.nflux+1 == rfield.nflux_with_check) ;
	 /* radius is outer radius of this zone, this resets continuum to
	  * flux at illuminated face of zone we have already computed. */

	/* opac.tmn defined in ZoneStart */
	/* undilute radiation so that rfield is at face of zone */
	/* NB, upper limit of sum includes [nflux] to test unit verification cell */
	for( long i=0; i <= rfield.nflux; i++ )
	{
		rfield.flux_beam_const[i] /= opac.tmn[i];
		rfield.flux_beam_time[i] /= opac.tmn[i];
		rfield.flux_isotropic[i] /= opac.tmn[i];
		rfield.flux[0][i] = rfield.flux_beam_const[i] + rfield.flux_beam_time[i] +
			rfield.flux_isotropic[i];
		/* >>chng 03 nov 08, update SummedCon here since flux changed */
		rfield.SummedCon[i] = rfield.flux[0][i] + rfield.SummedDif[i];
		rfield.SummedOcc[i] = rfield.SummedCon[i]*rfield.convoc[i];
	}

	rfield.setTrimming();

	/* do advection if enabled */
	if( dynamics.lgAdvection || dynamics.lgTimeDependentStatic )
		DynaEndZone();

	if (0)
		mole_rk_bigchange();

	return;
}
