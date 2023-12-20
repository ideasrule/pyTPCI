/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*OpacityAdd1Subshell add opacity due to single shell to main opacity array*/
/*OpacityAdd1SubshellInduc add opacity of individual species, including stimulated emission */
#include "cddefines.h"
#include "rfield.h"
#include "hydrogenic.h"
#include "opacity.h"

void OpacityAdd1Subshell(
	/*ipOpac is opacity index within opac opacity offset for this species */
	long int ipOpac, 
	/* lower freq limit to opacity range on energy mesh */
	long int ipLowLim, 
	/* upper limit to opacity range on energy mesh */
	long int ipUpLim, 
	/* abundance, we bail if zero */
	realnum abundance,
	/* either static 's' or volitile 'v' */
	char chStat )
{
	long int i ,
	  ipOffset, 
	  limit;

	DEBUG_ENTRY( "OpacityAdd1Subshell()" );

	/* code spends roughly 20% of its time in this loop*/

	ASSERT( chStat == 's' || chStat == 'v' );

	ipOffset = ipOpac - ipLowLim;
	ASSERT( ipLowLim > 0 );
	/* >>chng 02 aug 13, negative offset is ok, it is only this plus ipLowLim that matters */
	/*ASSERT( ipOffset >= 0 );*/

	limit = MIN2(ipUpLim,rfield.nflux);

	/* do nothing if abundance is zero, or if static opacities do not
	 * need to be redone */
	if( abundance <= 0. || (chStat=='s' && !opac.lgRedoStatic) )
	{ 
		return;
	}

	/* volative (outer shell, constantly reevaluated) or static opacity? */
	if( chStat=='v' )
	{
		for( i=ipLowLim-1; i < limit; i++ )
		{
			opac.opacity_abs[i] += opac.OpacStack[i+ipOffset]*abundance;
		}
	}
	else
	{
		for( i=ipLowLim-1; i < limit; i++ )
		{
			opac.OpacStatic[i] += opac.OpacStack[i+ipOffset]*abundance;
		}
	}
	return;
}

/*OpacityAdd1SubshellInduc add opacity of individual species, including stimulated emission */
void OpacityAdd1SubshellInduc(
	/* pointer to opacity within opacity stack */
	long int ipOpac, 
	/* pointer to low energy in continuum array for this opacity band */
	long int ipLowEnergy, 
	/* pointer to high energy in continuum array for this opacity */
	long int ipHiEnergy, 
	/* this abundance of this species, may be zero */
	double abundance, 
	/* the departure coef, may be infinite or zero */
	double DepartCoef ,
	/* either 'v' for volitile or 's' for static opacities */
	char chStat )
{
	long int i, 
	  iup, 
	  k;

	DEBUG_ENTRY( "OpacityAdd1SubshellInduc()" );

	/* add opacity of individual species, including stimulated emission
	 * abundance is the density of the lower level (cm^-3)
	 * DepartCoef is its departure coefficient, can be zero */

	/* this is opacity offset, must be positive */
	ASSERT( ipOpac > 0 );

	/* check that chStat is either 'v' or 's' */
	ASSERT( chStat == 'v' || chStat == 's' );

	/* do nothing if abundance is zero, or if static opacities do not
	 * need to be redone */
	if( abundance <= 0. || (chStat=='s' && !opac.lgRedoStatic) )
	{ 
		return;
	}

	k = ipOpac - ipLowEnergy;

	/* DepartCoef is dep coef, rfield.lgInducProcess is turned off with 'no induced' command */
	if( (DepartCoef > 1e-35 && rfield.lgInducProcess) && hydro.lgHInducImp )
	{
		iup = MIN2(ipHiEnergy,rfield.nflux);
		/* >>>chng 99 apr 29, following was present, caused pdr to make opac at impossible energy*/
		/*iup = MAX2(ipLowEnergy,iup);*/
		double DepartCoefInv = 1./DepartCoef;
		if( chStat == 'v' )
		{
			/* volitile opacities, always reevaluate */
			for( i=ipLowEnergy-1; i < iup; i++ )
			{
				opac.opacity_abs[i] += opac.OpacStack[i+k]*abundance*
					MAX2(0. , 1.-  rfield.ContBoltz[i]*DepartCoefInv);
			}
		}
		else
		{
			/* static opacities, save in special array */
			for( i=ipLowEnergy-1; i < iup; i++ )
			{
				opac.OpacStatic[i] += opac.OpacStack[i+k]*abundance*
					MAX2(0. , 1.-  rfield.ContBoltz[i]*DepartCoefInv);
			}
		}
	}

	else
	{
		/* DepartCoef is the departure coef, which can go to zero, 
		 * neglect stimulated emission in this case */
		iup = MIN2(ipHiEnergy,rfield.nflux);
		/* >>>chng 99 apr 29, following was present, caused pdr to make opac at impossible energy*/
		/*iup = MAX2(ipLowEnergy,iup);*/
		if( chStat == 'v' )
		{
			for( i=ipLowEnergy-1; i < iup; i++ )
			{
				opac.opacity_abs[i] += opac.OpacStack[i+k]*abundance;
			}
		}
		else
		{
			for( i=ipLowEnergy-1; i < iup; i++ )
			{
				opac.OpacStatic[i] += opac.OpacStack[i+k]*abundance;
			}
		}
	}

	return;
}
