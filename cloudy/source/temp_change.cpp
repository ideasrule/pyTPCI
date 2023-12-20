/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*tfidle update some temperature dependent variables */
/*tauff compute optical depth where cloud is thin to free-free and plasma freq */
#include "cddefines.h"
#include "conv.h"
#include "opacity.h"
#include "dense.h"
#include "phycon.h"
#include "stopcalc.h"
#include "trace.h"
#include "rfield.h"
#include "doppvel.h"
#include "radius.h"
#include "wind.h"
#include "thermal.h"
#include "vectorize.h"
#include "taulines.h"

/*tauff compute optical depth where cloud is thin to free-free and plasma freq */
STATIC void tauff();
/* On first run, fill GauntFF with gaunt factors	*/

/**tfidle update some temperature dependent variables 
\param lgForceUpdate option to force update of all variables 
*/
STATIC void tfidle(bool lgForceUpdate);

/**TempChange change kinetic temperature, calls tfidle
*/
void TempChange(double TempNew ,
		/* option to force update of all variables */
		bool lgForceUpdate)
{
	DEBUG_ENTRY( "TempChange()" );

	const char* message = "\n This is often due to mixing luminosity/intensity commands, having"
		" an unphysical SED, or too high a flux of cosmic rays.\n"
		" To find the problem, try setting a constant temperature (and STOP ZONE 1),"
		" rerun the model, and check the output to find what is happening.\n\n";

	/* set new temperature */
	if( TempNew > phycon.TEMP_LIMIT_HIGH )
	{
		/* temp is too high */
		fprintf(ioQQQ,"\n\n PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is above the upper limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_HIGH );
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n%s", message);
		throw cloudy_abort("the kinetic temperature is above the upper limit");
	}
	else if( TempNew < phycon.TEMP_LIMIT_LOW )
	{
		/* temp is too low */
		fprintf(ioQQQ,"\n\n PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is below the lower limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_LOW );
		fprintf(ioQQQ," Consider setting a lowest temperature with the SET TEMPERATURE FLOOR command.\n");
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n%s", message);	
		throw cloudy_abort("the kinetic temperature is below the lower limit");
	}
	else if( TempNew < StopCalc.TeFloor )
	{
		if( trace.lgTrace || trace.nTrConvg>=2 )
			fprintf(ioQQQ,"temp_change: temp change floor hit, TempNew=%.3e TeFloor=%.3e, "
					"setting constant temperature, nTotalIoniz=%li\n",
					 TempNew , StopCalc.TeFloor , conv.nTotalIoniz);
		/* temperature floor option  -
		 * go to constant temperature calculation if temperature
		 * falls below floor */
		thermal.lgTemperatureConstant = true;
		thermal.ConstTemp = (realnum)StopCalc.TeFloor;
		phycon.te = thermal.ConstTemp;
		/*fprintf(ioQQQ,"DEBUG TempChange hit temp floor, setting const temp to %.3e\n",
			phycon.te );*/
	}
	else
	{
		/* temp is within range */
		phycon.te = TempNew;
	}

	/* now update related variables */
	if ( conv.lgSearch )
	{
		/* in search phase, force species2.cpp level trimming to be
		 * re-tested */
		for( long ipSpecies=0; ipSpecies<nSpecies; ++ipSpecies )
		{
			dBaseSpecies[ipSpecies].numLevels_local = dBaseSpecies[ipSpecies].numLevels_max;
		}		
	}
	tfidle(lgForceUpdate);
	return;
}
/**TempChange change kinetic temperature, calls tfidle
 * but does not update extensive variables or check for temperature floor,
 * intended for use by routines that are sanity checks rather than real calculation */
void TempChange(double TempNew)
{

	DEBUG_ENTRY( "TempChange()" );

	/* set new temperature */
	if( TempNew > phycon.TEMP_LIMIT_HIGH )
	{
		/* temp is too high */
		fprintf(ioQQQ," PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is above the upper limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_HIGH );
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n");
		throw cloudy_abort("the kinetic temperature is above the upper limit");
	}
	else if( TempNew < phycon.TEMP_LIMIT_LOW )
	{
		/* temp is too low */
		fprintf(ioQQQ," PROBLEM DISASTER - the kinetic temperature, %.3eK,"
			" is below the lower limit of the code, %.3eK.\n",
			TempNew , phycon.TEMP_LIMIT_LOW );
		fprintf(ioQQQ," Consider setting a lowest temperature with the SET TEMPERATURE FLOOR command.\n");
		fprintf(ioQQQ," This calculation is aborting.\n Sorry.\n");
		throw cloudy_abort("the kinetic temperature is below the lower limit");
	}
	else
	{
		/* temp is within range */
		phycon.te = TempNew;
	}

	/* now update related variables */
	tfidle( false );
	return;
}

void tfidle(
	/* option to force update of all variables */
	bool lgForceUpdate)
{
	static double tgffused=-1.;
	static double ttused = 0.;
	static bool lgZLogSet = false;

	DEBUG_ENTRY( "tfidle()" );

	/* called with lgForceUpdate true in zero.c, when we must update everything */
	if( lgForceUpdate )
	{
		ttused = -1.;
		tgffused = -1.;
	}

	/* check that eden not negative */
	if( dense.eden <= 0. )
	{
		fprintf( ioQQQ, "tfidle called with a zero or negative electron density,%10.2e\n", 
		  dense.eden );
		TotalInsanity();
	}

	/* check that temperature not negative */
	if( phycon.te <= 0. )
	{
		fprintf( ioQQQ, "tfidle called with a negative electron temperature,%10.2e\n", 
		  phycon.te );
		TotalInsanity();
	}

	/* one time only, set up array of logs of charge squared */
	if( !lgZLogSet )
	{
		for( long nelem=0; nelem<LIMELM; ++nelem )
		{
			/* this array is used to modify the log temperature array
			 * defined below, for hydrogenic species of charge nelem+1 */
			phycon.sqlogz[nelem] = log10( POW2(nelem+1.) );
		}
		lgZLogSet = true;
	}

	if( ! fp_equal( phycon.te, ttused ) )
	{
		ttused = phycon.te;
		thermal.te_update = phycon.te;
		/* current temperature in various units */
		phycon.te_eV = phycon.te/EVDEGK;
		phycon.te_ryd = phycon.te/TE1RYD;
		phycon.te_wn = phycon.te / T1CM;

		phycon.tesqrd = POW2(phycon.te);
		phycon.sqrte = sqrt(phycon.te);
		thermal.halfte = 0.5/phycon.te;
		thermal.tsq1 = 1./phycon.tesqrd;
		phycon.te32 = phycon.te*phycon.sqrte;
		phycon.teinv = 1./phycon.te;

		phycon.alogte = log10(phycon.te);
		phycon.alnte = log(phycon.te);

		phycon.telogn[0] = phycon.alogte;
		for( int i=1; i < 7; i++ )
		{
			phycon.telogn[i] = phycon.telogn[i-1]*phycon.telogn[0];
		}

		phycon.te10 = pow(phycon.te,0.10);
		phycon.te20 = phycon.te10 * phycon.te10;
		phycon.te30 = phycon.te20 * phycon.te10;
		phycon.te40 = phycon.te30 * phycon.te10;
		phycon.te70 = phycon.sqrte * phycon.te20;
		phycon.te90 = phycon.te70 * phycon.te20;

		phycon.te01 = pow(phycon.te,0.01);
		phycon.te02 = phycon.te01 * phycon.te01;
		phycon.te03 = phycon.te02 * phycon.te01;
		phycon.te04 = phycon.te02 * phycon.te02;
		phycon.te05 = phycon.te03 * phycon.te02;
		phycon.te07 = phycon.te05 * phycon.te02;

		phycon.te001 = pow(phycon.te,0.001);
		phycon.te002 = phycon.te001 * phycon.te001;
		phycon.te003 = phycon.te002 * phycon.te001;
		phycon.te004 = phycon.te002 * phycon.te002;
		phycon.te005 = phycon.te003 * phycon.te002;
		phycon.te007 = phycon.te005 * phycon.te002;
		/*>>>chng 06 June 30 -Humeshkar Nemala*/
		phycon.te0001 = pow(phycon.te ,0.0001);
		phycon.te0002 = phycon.te0001 * phycon.te0001;
		phycon.te0003 = phycon.te0002 * phycon.te0001;
		phycon.te0004 = phycon.te0002 * phycon.te0002;
		phycon.te0005 = phycon.te0003 * phycon.te0002;
		phycon.te0007 = phycon.te0005 * phycon.te0002;

	}

	/* >>>chng 99 nov 23, removed this line, so back to old method of h coll */
	/* used for hydrogenic collisions */
	/* 
	 * following electron density has approximate correction for neutrals
	 * corr of hi*1.7e-4 accounts for col ion by HI; 
	 * >>refer	H0	correction for collisional contribution		Drawin, H.W. 1969, Zs Phys 225, 483.
	 * also quoted in Dalgarno & McCray 1972
	 * extensive discussion of this in 
	 *>>refer	H0	collisions	Lambert, D.L. 
	 * used EdenHCorr instead
	 * edhi = eden + hi * 1.7e-4
	 */
	dense.EdenHCorr = dense.eden + 
		/* dense.HCorrFac is unity by default and changed with the set HCOR command */
		dense.xIonDense[ipHYDROGEN][0]*1.7e-4 * dense.HCorrFac;
	dense.EdenHCorr_f = (realnum)dense.EdenHCorr;
	
	/*>>chng 93 jun 04,
	 * term with hi added June 4, 93, to account for warm pdr */
	/* >>chng 05 jan 05, Will Henney noticed that 1.e-4 used here is not same as
	 * 1.7e-4 used for EdenHCorr, which had rewritten the expression.
	 * change so that edensqte uses EdenHCorr rather than reevaluating */
	/*dense.edensqte = ((dense.eden + dense.xIonDense[ipHYDROGEN][0]/1e4)/phycon.sqrte);*/
	dense.edensqte = dense.EdenHCorr/phycon.sqrte;
	dense.cdsqte = dense.edensqte*COLL_CONST;
	dense.SqrtEden = sqrt(dense.eden);

	/* rest have to do with radiation field and frequency mesh which may not be defined yet */
	if( !lgRfieldAllocated || !rfield.lgMeshSetUp() )
		return;

	/* correction factors for induced recombination, 
	 * also used as Boltzmann factors
	 * check for zero is because ContBoltz is zeroed out in initialization
	 * of code, its possible this is a constant density grid of models
	 * in which the code is called as a subroutine */
	/* >>chng 01 aug 21, must also test on size of continuum nflux because 
	 * conintitemp can increase nflux then call this routine, although 
	 * temp may not have changed */
	if( ! fp_equal(tgffused, phycon.te) || rfield.ContBoltz[0] <= 0. )
	{
		tgffused = phycon.te;
		for( long i=0; i < rfield.nflux_with_check; ++i )
			rfield.vexp_arg[i] = -rfield.anu(i)/phycon.te_ryd;
		/* atom_level2 uses ContBoltz to see whether the rates will be significant.
		 * If the numbers did not agree then this test would be flawed, resulting in
		 * div by zero */
		vexp( rfield.vexp_arg.data(), rfield.ContBoltz.data(), 0, rfield.nflux_with_check );
		for( long i=0; i < rfield.nflux_with_check; ++i )
			rfield.vexp_arg[i] = -rfield.anumin(i)/phycon.te_ryd;
		/* this is Boltzmann factor averaged over the width of the cell */
		vexp( rfield.vexp_arg.data(), rfield.ContBoltzAvg.data(), 0, rfield.nflux_with_check );
		for( long i=0; i < rfield.nflux_with_check; ++i )
			rfield.vexp_arg[i] = -rfield.widflx(i)/phycon.te_ryd;
		vexpm1( rfield.vexp_arg.data(), rfield.ContBoltzHelp1.data(), 0, rfield.nflux_with_check );
		for( long i=0; i < rfield.nflux_with_check; ++i )
		{
			rfield.ContBoltzHelp1[i] /= rfield.vexp_arg[i];
			rfield.ContBoltzAvg[i] *= rfield.ContBoltzHelp1[i];
		}
		vexp( rfield.vexp_arg.data(), rfield.ContBoltzHelp2.data(), 0, rfield.nflux_with_check );
		/* ipMaxBolt is number of non-zero cells, so non-zero up through ipMaxBolt-1 */
		rfield.ipMaxBolt = rfield.nflux_with_check;
		for( long i=rfield.nflux_with_check-1; i >= 0; --i )
		{
			if( rfield.ContBoltz[i] > 0. )
				break;
			rfield.ipMaxBolt = i;
		}
	}

	/* find frequency where thin to bremsstrahlung or plasma frequency */
	tauff();
}

/*tauff compute optical depth where cloud is thin to free-free and plasma freq */
STATIC void tauff()
{
	realnum fac;

	/* simply return if space not yet allocated */
	if( !lgOpacAllocated )
		return;

	DEBUG_ENTRY( "tauff()" );

	if( !conv.nTotalIoniz )
		rfield.ipEnergyBremsThin = 0;

	/* routine sets variable ipEnergyBremsThin, index for lowest cont cell that is optically thin */
	/* find frequency where continuum thin to free-free */
	while( rfield.ipEnergyBremsThin < rfield.nflux && 
		opac.TauAbsGeo[0][rfield.ipEnergyBremsThin] >= 1. )
	{
		++rfield.ipEnergyBremsThin;
	}

	/* TFF will be frequency where cloud becomes optically thin to bremsstrahlung
	 * >>chng 96 may 7, had been 2, change as per Kevin Volk bug report */
	if( rfield.ipEnergyBremsThin > 1 && opac.TauAbsGeo[0][rfield.ipEnergyBremsThin] > 0.001 )
	{
		/* tau can be zero when plasma frequency is within energy grid, */
		fac = (1.f - opac.TauAbsGeo[0][rfield.ipEnergyBremsThin-1])/(opac.TauAbsGeo[0][rfield.ipEnergyBremsThin] - 
		  opac.TauAbsGeo[0][rfield.ipEnergyBremsThin-1]);
		fac = MAX2(fac,0.f);
		rfield.EnergyBremsThin = rfield.anu(rfield.ipEnergyBremsThin-1) + rfield.widflx(rfield.ipEnergyBremsThin-1)*fac;
	}
	else
	{
		rfield.EnergyBremsThin = 0.f;
	}

	/* did not include plasma freq before
	 * function returns larger of these two frequencies */
	rfield.EnergyBremsThin = MAX2(rfield.plsfrq,rfield.EnergyBremsThin);

	/* now increment ipEnergyBremsThin still further, until above plasma frequency */
	while( rfield.ipEnergyBremsThin < rfield.nflux && 
		rfield.anu(rfield.ipEnergyBremsThin) <= rfield.EnergyBremsThin )
	{
		++rfield.ipEnergyBremsThin;
	}
	return;
}

realnum GetDopplerWidth( realnum massAMU )
{
	ASSERT( massAMU > 0. );
	// force a fairly conservative upper limit
	ASSERT( massAMU < 10000. );

	/* usually TurbVel =0, reset with turbulence command
	 * cm/s here, but was entered in km/s with command */
	double turb2 = POW2(DoppVel.TurbVel);

	/* this is option to dissipate the turbulence.  DispScale is entered with
	 * dissipate keyword on turbulence command.  The velocity is reduced here,
	 * by an assumed exponential scale, and also adds heat */
	if( DoppVel.DispScale > 0. )
	{
		/* square of exp depth dependence */
		turb2 *= sexp( 2.*radius.depth / DoppVel.DispScale );
	}

	/* in case of D-Critical flow include initial velocity as
	 * a component of turbulence */
	if( ! ( wind.lgBallistic() || wind.lgStatic() ) )
	{
		turb2 += POW2(wind.windv0);
	}

	realnum width = (realnum)sqrt(2.*BOLTZMANN/ATOMIC_MASS_UNIT*phycon.te/massAMU+turb2);
	ASSERT( width > 0.f );
	return width;
}

realnum GetAveVelocity( realnum massAMU )
{
#if 0
	/* usually TurbVel =0, reset with turbulence command
	 * cm/s here, but was entered in km/s with command */
	double turb2 = POW2(DoppVel.TurbVel);

	/* this is option to dissipate the turbulence.  DispScale is entered with
	 * dissipate keyword on turbulence command.  The velocity is reduced here,
	 * by an assumed exponential scale, and also adds heat */
	if( DoppVel.DispScale > 0. )
	{
		/* square of exp depth dependence */
		turb2 *= sexp( 2.*radius.depth / DoppVel.DispScale );
	}

	/* in case of D-Critical flow include initial velocity as
	 * a component of turbulence */
	if( ! ( wind.lgBallistic() || wind.lgStatic() ) )
	{
		turb2 += POW2(wind.windv0);
	}
#endif
	DEBUG_ENTRY( "GetAveVelocity()" );

	/* this is average (NOT rms) particle speed for Maxwell distribution, Mihalas 70, 9-70 */
	fixit("turbulence was included here for molecules but not ions.  Now neither. Resolve.");
	return (realnum)sqrt(8.*BOLTZMANN/PI/ATOMIC_MASS_UNIT*phycon.te/massAMU);
}
