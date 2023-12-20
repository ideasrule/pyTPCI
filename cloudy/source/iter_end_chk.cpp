/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iter_end_check after each zone by Cloudy, determines whether model is complete */
#include "cddefines.h"
/*  */
#ifdef EPS
#	undef EPS
#endif
#define	EPS	1.00001
#include "lines.h"
#include "mole.h"
#include "conv.h"
#include "iterations.h"
#include "trace.h"
#include "dense.h"
#include "colden.h"
#include "taulines.h"
#include "hmi.h"
#include "prt.h"
#include "phycon.h"
#include "geometry.h"
#include "stopcalc.h"
#include "opacity.h"
#include "thermal.h"
#include "cooling.h"
#include "predcont.h"
#include "pressure.h"
#include "radius.h"
#include "called.h"
#include "wind.h"
#include "hcmap.h"
#include "rfield.h"
#include "flux.h"

/*dmpary print all coolants for some zone, as from print cooling command */
STATIC void dmpary(void);

int iter_end_check(void)
{
	bool lgDone, 
	  lgEndFun_v, 
	  lgPrinted;
	long int i;
	double oxy_in_grains;

	DEBUG_ENTRY( "iter_end_check()" );

	/* >>chng 05 nov 22 - NPA.  Stop calculation when fraction of oxygen frozen
	 * out on grains gets too high - 
	 * NB this test is not used since StopCalc.StopDepleteFrac is set to > 1 */
	oxy_in_grains = 0.0f;
	for(i=0;i<mole_global.num_calc;++i)
	{
		/* define the abundance of oxygen frozen out on grain surfaces */
		if( ! mole_global.list[i]->lgGas_Phase && mole_global.list[i]->isIsotopicTotalSpecies() )
			oxy_in_grains += mole.species[i].den*mole_global.list[i]->nElement(ipOXYGEN);
	}
	/*fprintf(ioQQQ, "DEBUG oxy in grains %.2e %e %e\n", 
		oxy_in_grains ,
		oxy_in_grains/MAX2(SMALLFLOAT,dense.gas_phase[ipOXYGEN]) , StopCalc.StopDepleteFrac );*/

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " iter_end_check called, zone %li.\n" , nzone);
	}
	ASSERT( hcmap.MapZone >= 00 || !conv.lgSearch );

	/* >>chng 97 jun 09, now called before first zone with nzone 0 */
	if( nzone == 0 )
	{
		lgEndFun_v = false;

		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " iter_end_check returns, doing nothing since zone 0.\n" );
		}
		return( lgEndFun_v );
	}

	/* check that species is valid, if we may stop for this reason */
	static bool lgSpeciesKnownGood=false;
	static molezone *SpeciesCurrent;
	if( StopCalc.lgStopSpeciesColumn && !lgSpeciesKnownGood )
	{
		SpeciesCurrent =
			findspecieslocal_validate(StopCalc.chSpeciesColumn.c_str());
		lgSpeciesKnownGood = true;
	}

	/* check whether trace is needed for this zone and iteration */
	if( (nzone >= trace.nznbug && iteration >= trace.npsbug) && trace.lgTrOvrd )
	{
		if( trace.nTrConvg==0 )
		{
			geometry.nprint = 1;
			trace.lgTrace = true;
		}
		else
			/* trace convergence entered = but with negative flag = make positive,
			 * abs and not mult by -1 since may trigger more than one time */
			trace.nTrConvg = abs( trace.nTrConvg );
	}

	/* option to turn printout on after certain zone; only the master rank talks */
	if( prt.lgPrtStart && prt.nstart == nzone )
	{
		called.lgTalk = cpu.i().lgMPI_talk();
	}

	/* check whether model is done, various criteria used. */
	lgDone = false;

	/* this is flag to check whether stopping reason was bad */
	conv.lgBadStop = false;

	/* set temperature floor option  -
	 * go to constant temperature calculation if temperature
	 * falls below floor */
	if( phycon.te < StopCalc.TeFloor )
	{
		thermal.lgTemperatureConstant = true;
		thermal.ConstTemp = (realnum)StopCalc.TeFloor;
		phycon.te = thermal.ConstTemp;
		TempChange(thermal.ConstTemp , false);
	}

	/* check on radiation pressure - constant pressure unstable if dominated
	 * by radiation pressure */
	if( (pressure.lgPres_radiation_ON && pressure.pbeta > 1.0) && 
		(strcmp(dense.chDenseLaw ,"CPRE") == 0) && 
		/* >>chng 03 aug 20, check on extreme values of pbeta, and abort if true */
		(iterations.lgLastIt||(pressure.pbeta>1000.)) &&
		/* >>chng 03 aug 19, add check on pbeta, and abort even if "no abort"
		 * was specified, since rad pres dominated limit can lead to VERY
		 * small H densities and crash due to underflow */
		(pressure.lgRadPresAbortOK||(pressure.pbeta>1000.)) )
	{
		/* const total pres model; if RadPres>PGAS, then unstable, stop */
		if( called.lgTalk )
		{
			fprintf( ioQQQ, "\n STOP since P(rad)/P(gas)=%7.3f >1.0\n", 
			  pressure.pbeta );

			fprintf( ioQQQ, " Line contributors to radiation pressure follows:\n" );
			PrtLinePres(ioQQQ);
		}
		lgDone = true;
		conv.lgBadStop = true;
		strncpy( StopCalc.chReasonStop, "of radiation pressure.", sizeof(StopCalc.chReasonStop) );
	}

	/* radius and resulting volume too large for this cpu */
	if( radius.drad_x_fillfac*radius.r1r0sq > BIGFLOAT/10.)
	{
		/* too big */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "volume too large for this cpu.", sizeof(StopCalc.chReasonStop) );
	}
	/* supersonic outflowing wind, initial velocity, windv0, was > 0,
	 * but it has coasted to a stop - lgVelPos is false */
	else if( !wind.lgVelPos && wind.lgBallistic() )
	{
		/* wind solution with negative velocity */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "wind veloc too small.", sizeof(StopCalc.chReasonStop) );
	}

	else if( !wind.lgStatic() && fabs(wind.windv) < StopCalc.StopVelocity )
	{
		/* stop if absolute value of velocity falls below value */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "wind V too small.", sizeof(StopCalc.chReasonStop) );
	}

	/* this flag says that 21cm line optical depth is the stop quantity */
	else if( StopCalc.lgStop21cm && (HFLines[0].Emis().TauCon() >=  StopCalc.tauend) )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "21 cm optical depth.", sizeof(StopCalc.chReasonStop) );
	}

	else if( rfield.extin_mag_V_extended >= StopCalc.AV_extended )
	{
		/* stop at specified AV for (1-g) in scattering opacity */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "A_V reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( rfield.extin_mag_V_point >= StopCalc.AV_point )
	{
		/* stop at specified AV without (1-g) in scattering opacity */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "A_V reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( StopCalc.xMass!=0. &&
		log10(SDIV(dense.xMassTotal))+1.0992099+ 2.*log10(radius.rinner) >= StopCalc.xMass )
	{
		/* stop at specified AV without (1-g) in scattering opacity */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "mass reached.", sizeof(StopCalc.chReasonStop) );
	}

	/* >>chg 02 may 31, added pressure.lgSonicPoint logic */
	/* WJH 19 May 2004: for some models, we want to get through the
	 * sonic point and out the other side */ 
	else if( pressure.lgSonicPoint && pressure.lgSonicPointAbortOK )
	{
		/* D-critical solution reached sonic point */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "sonic point reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( dense.EdenTrue==0 )
	{
		/* calculation failed */
		conv.lgBadStop = true;
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "zero electron density.", sizeof(StopCalc.chReasonStop) );
	}

	else if( radius.lgdR2Small )
	{
		lgDone = true;
		conv.lgBadStop = true;
		strncpy( StopCalc.chReasonStop, "DR small rel to thick.", sizeof(StopCalc.chReasonStop) );
	}

	else if( dense.eden < StopCalc.StopElecDensity )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "lowest EDEN reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( dense.eden/dense.gas_phase[ipHYDROGEN] < StopCalc.StopElecFrac )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "low electron fraction.", sizeof(StopCalc.chReasonStop) );
	}

	/* >>chng 05 nov 22, NA add this stop condition - stop when too many molecules
	 * are ices or solids on grains - the limit StopCalc.StopDepleteFrac is changed
	 * with the stop molecular depletion command */
	else if( dense.lgElmtOn[ipOXYGEN] &&
		(oxy_in_grains/MAX2(SMALLFLOAT,dense.gas_phase[ipOXYGEN])) > StopCalc.StopDepleteFrac )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "freeze out fraction.", sizeof(StopCalc.chReasonStop) );
	}

	/*else if( 2.*findspecieslocal("H2")->den/dense.gas_phase[ipHYDROGEN] < StopCalc.StopH2MoleFrac )*/
	else if( 2.*hmi.H2_total/dense.gas_phase[ipHYDROGEN] > StopCalc.StopH2MoleFrac )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "large H_2/H fraction.", sizeof(StopCalc.chReasonStop) );
	}

	else if( dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN] < 
		StopCalc.StopHPlusFrac )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "low H_+/H fraction.", sizeof(StopCalc.chReasonStop) );
	}

	else if( radius.lgDrMinUsed )
	{
		/* dr became too small */
		conv.lgBadStop = true;
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "DRAD small- set DRMIN.", sizeof(StopCalc.chReasonStop) );
	}

	else if( radius.depth >= iterations.StopThickness[iteration-1]/EPS )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "outer radius reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( StopCalc.iptnu >= 0 && 
		opac.TauAbsGeo[0][StopCalc.iptnu-1] >= StopCalc.tauend/EPS )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "optical depth reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( StopCalc.lgStopSpeciesColumn && 
		SpeciesCurrent->column >= StopCalc.col_species/EPS )
	{
		/* StopCalc.col_species default set to COLUMN_INIT == 1e30 */
		lgDone = true;
		sprintf( StopCalc.chReasonStop, "%s column dens reached.", StopCalc.chSpeciesColumn.c_str() );
		//strncpy( StopCalc.chReasonStop, "H column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( colden.colden[ipCOL_HTOT] >= StopCalc.HColStop/EPS )
	{
		/* StopCalc.HColStop default set to COLUMN_INIT == 1e30 */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "H column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( findspecieslocal("H+")->column >= StopCalc.colpls/EPS )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "H+ column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( (findspecieslocal("H2")->column+findspecieslocal("H2*")->column) >= StopCalc.col_h2/EPS )
	{
		/* >>chng 03 apr 15, add molecular hydrogen */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "H2 column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( (2.*(findspecieslocal("H2")->column+findspecieslocal("H2*")->column) + findspecieslocal("H")->column) >= StopCalc.col_h2_nut/EPS )
	{
		/* >>chng 04 feb 10, stopping command for H2 + H I */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "2 H2+H0 column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( colden.H0_ov_Tspin >= StopCalc.col_H0_ov_Tspin/EPS )
	{
		/* >>chng 05 jan 09, stopping command for N(H0) / Tspin */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "N(H0)/Tspin column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( findspecieslocal("CO")->column >= StopCalc.col_monoxco/EPS )
	{
		/* >>chng 03 oct 27--Nick Abel, add carbon monoxide */
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "CO column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( findspecieslocal("H")->column >= StopCalc.colnut/EPS )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "H0 column dens reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( phycon.te > StopCalc.TempHiStopZone )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "highest Te reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( phycon.te < StopCalc.TempLoStopZone )
	{
		lgDone = true;
		strncpy( StopCalc.chReasonStop, "lowest Te reached.", sizeof(StopCalc.chReasonStop) );
	}

	else if( nzone >= iterations.nend[iteration-1] )
	{
		lgDone = true;
		geometry.lgZoneTrp = true;
		strncpy( StopCalc.chReasonStop, "NZONE reached.", sizeof(StopCalc.chReasonStop) );
	}

	/* option to stop calculation when line intensity ratio reaches certain value,
	 * nstpl is number of stop line commands entered */
	else if( StopCalc.nstpl > 0 || StopCalc.ContIndex.size() > 0 )
	{
		/* line ratio exceeded maximum permitted value
		 * do not consider case where norm line has zero intensity */
		for( i=0; i < StopCalc.nstpl; i++ )
		{
			/* the second line is always set to something, default is H beta */
			if( LineSave.lines[StopCalc.ipStopLin2[i]].SumLine(StopCalc.nEmergent[i]) > 0. )
			{
				if( LineSave.lines[StopCalc.ipStopLin1[i]].SumLine(StopCalc.nEmergent[i])/
					 LineSave.lines[StopCalc.ipStopLin2[i]].SumLine(StopCalc.nEmergent[i]) > 
					StopCalc.stpint[i] )
				{
					lgDone = true;
					sprintf( StopCalc.chReasonStop, "line %s reached", 
								LineSave.lines[StopCalc.ipStopLin1[i]].label().c_str() );
				}
			}
		}
		/* continuum flux exceeded maximum permitted value */
		for( size_t k=0; k < StopCalc.ContIndex.size(); ++k )
		{
			// there are 4 entries for each wavelength: nFnu, nInu, InwT, InwC
			long ind = t_PredCont::Inst().offset() + 4*StopCalc.ContIndex[k];
			double nFnu_model = LineSave.lines[ind].SumLine(0) * radius.Conv2PrtInten;
			if( nFnu_model >= StopCalc.ContNFnu[k].get("erg/s/cm2") )
			{
				lgDone = true;
				sprintf( StopCalc.chReasonStop, "flux %s reached", 
							LineSave.lines[ind].label().c_str() );
			}
		}
	}

	else if( radius.drNext <= 0. )
	{
		/* this cant happen */
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " drNext=%10.2e STOP\n", radius.drNext );
		}
		conv.lgBadStop = true;
		strncpy( StopCalc.chReasonStop, "internal error - DRAD.", sizeof(StopCalc.chReasonStop) );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	lgPrinted = false;
	if( lgDone )
	{
		/* flag to call it quits */
		lgEndFun_v = true;
		PrtZone();
		lgPrinted = true;
	}

	else
	{
		/* passed all the tests, keep going */
		/* check whether this zone should be printed */
		if( (nzone % geometry.nprint == 0 || 
		  nzone == 1) || trace.nTrConvg )
		{
			PrtZone();
			lgPrinted = true;
		}
		/* flag to keep going */
		lgEndFun_v = false;
	}

	/* dump cooling arrays for this zone? */
	if( prt.nzdump == nzone || prt.nzdump == 0 )
		dmpary();

	/* do map of cooling function if desired, and not yet done */
	/* >>chng 02 may 29, MapZone < = to <= 0 - map 0 did not work */
	/* >>chng 04 jun 16, change to MapZone = 0 for map of first zone then quit,
	 * -1 is not set, positive, do map of that zone */
	if( !hcmap.lgMapDone && (hcmap.MapZone == 0 || nzone == hcmap.MapZone) )
	{
		/* print last zone if not already done */
		if( !lgPrinted )
		{
			PrtZone();
		}

		/* say that we are doing a map */
		hcmap.lgMapBeingDone = true;

		/* save old output file then redirect to map file */
		/* >>chng 01 mar 28, ioMAP may not be initialized, PvH */
		if( ioMAP != NULL )
			map_do(ioMAP, " map");
		else
			map_do(ioQQQ, " map");

		/* stop after doing map */
		lgEndFun_v = true;
		strncpy( StopCalc.chReasonStop, "MAP command used-stop.", sizeof(StopCalc.chReasonStop) );

		/* >>chng 03 jun 06, reset iterations since we want to stop even if
		 * iterate xx is specified, bug caught by Joop Schaye */
		iterations.itermx = 0;

		/* make really sure that the string contained in StopCalc.chReasonStop is properly terminated */
		StopCalc.chReasonStop[sizeof(StopCalc.chReasonStop)-1] = '\0';

		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " iter_end_check returns after map.\n" );
		}
		return( lgEndFun_v );
	}

	if( lgEndFun_v && prt.lgOnlyZone )
	{
		cdEXIT(EXIT_FAILURE);
	}

	/* the string contained in StopCalc.chReasonStop must be properly 
	 * terminated -this can't fail - strlen returns the number of characters
	 * in str, excluding the terminal NULL.*/
	if( strlen( StopCalc.chReasonStop ) >= nCHREASONSTOP-1 )
		TotalInsanity();

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " iter_end_check bottom return.\n" );
	}
	return( lgEndFun_v );
}

#ifdef EPS
#	undef EPS
#endif
#define	EPS	0.005
/*dmpary print all coolants for some zone, as from print cooling command */
STATIC void dmpary(void)
{
	long int i;
	realnum ratio;

	DEBUG_ENTRY( "dmpary()" );

	fprintf( ioQQQ, 
		" This is a print out of the cooling array for zone number %3ld\n", 
	  nzone );

	fprintf( ioQQQ, 
		" The total heating was HTOT=%10.2e erg/s/cm3, the total cooling was CTOT=%10.2e, and the temperature was%10.3eK.\n", 
	  thermal.htot, thermal.ctot, phycon.te );

	fprintf( ioQQQ, 
		" All coolants greater than%6.2f%% of the total will be printed.\n", 
	  EPS*100. );

	/* flag all significant coolants */
	coolpr(ioQQQ,"ZERO",1,0.,"ZERO");
	for( i=0; i < thermal.ncltot; i++ )
	{
		ratio = (realnum)(thermal.cooling[i]/thermal.ctot);
		if( fabs(ratio) > EPS )
		{
			coolpr(ioQQQ,thermal.chClntLab[i],thermal.collam[i],
			  ratio,"DOIT");
		}

		ratio = (realnum)(thermal.heatnt[i]/thermal.ctot);
		if( fabs(ratio) > EPS )
		{
			coolpr(ioQQQ,thermal.chClntLab[i],thermal.collam[i],
			  ratio,"DOIT");
		}
	}
	coolpr(ioQQQ,"DONE",1,0.,"DONE");
	return;
}
#undef EPS
