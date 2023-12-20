/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitSimPostparse initialization at start of simulation, 
 * called from cloudy after parser, 
 * called one time per sim in grid 
 * not called for every iteration */
#include "cddefines.h" 
#include "init.h" 
#include "dense.h"
#include "struc.h"
#include "elementnames.h"
#include "atoms.h"
#include "iterations.h"
#include "pressure.h"
#include "trace.h"
#include "radius.h"
#include "thermal.h"
#include "heavy.h"
#include "wind.h"
#include "iso.h"
#include "h2.h"
#include "monitor_results.h"
#include "taulines.h"
#include "mole.h"
#include "rfield.h"
#include "continuum.h"
#include "mean.h"

/*InitSimPostparse initialization after parser, 
 * called one time for each simulation in a grid, 
 * only before first iteration
 * sets initial or zeros values before start of a simulation */
void InitSimPostparse( void )
{
	DEBUG_ENTRY( "InitSimPostparse()" );
	static double one=1.0;

	struc.nzonePreviousIteration = -1;

	thermal.thist = 0.;
	thermal.tlowst = 1e20f;
	thermal.nUnstable = 0;
	thermal.lgUnstable = false;

	colliders.init();

	mole_global.init();

	mole.set_ion_locations();

	findspecieslocal("e-")->location = &(dense.eden);
	findspecieslocal("grn")->location = &(mole.grain_area);
	findspecieslocal("PHOTON")->location = &one;
	findspecieslocal("CRPHOT")->location = &one;
	findspecieslocal("CRP")->location = &one;

	/* >> chng 06 Nov 24 rjrw: Move reaction definitions here so they can be controlled by parsed commands */
	mole_create_react();

	//mole_cmp_num_in_out_reactions();

	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_Reset();

	/* read LTE radiative cooling, if small model in use */
	if( ! h2.lgEnabled )
		h2.H2_Read_LTE_cooling_per_H2();

	/* zero out diffuse Lya emission */
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		for( long ion=0; ion<=nelem; ++ion )
			Heavy.xLyaHeavy[nelem][ion] = 0;

	/* convert STOP RADIUS command to STOP THICKNESS command now that inner radius is known */
	if( iterations.StopRadius[0] > 0. )
	{
		for( long j=0; j < iterations.iter_alloc; j++ )
		{
			if( iterations.StopRadius[j] > radius.Radius )
				iterations.StopThickness[j] = iterations.StopRadius[j] - radius.Radius;
			else
			{
				fprintf( ioQQQ, " PROBLEM stopping radius is <= inner radius. Bailing out.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* this term applies centrifugal acceleration if DISK option on wind command */
	wind.DiskRadius = 0;
	if( wind.lgDisk )
		wind.DiskRadius = radius.Radius;

	iterations.lgLastIt = false;

	rfield_opac_zero( 0 , rfield.nflux_with_check );

	rfield.lgUSphON = false;

	/* where cloud is optically thick to brems */
	rfield.ipEnergyBremsThin = 0;
	rfield.EnergyBremsThin = 0.;
	rfield.comtot = 0.;
	rfield.cmcool = 0.;
	rfield.cinrat = 0.;
	rfield.extin_mag_B_point = 0.;
	rfield.extin_mag_V_point = 0.;
	rfield.extin_mag_B_extended = 0.;
	rfield.extin_mag_V_extended = 0.;
	rfield.EnerGammaRay = 7676.;

	for( vector<TransitionList>::iterator it = AllTransitions.begin(); it != AllTransitions.end(); ++it )
	{
		for( TransitionList::iterator tr = it->begin(); tr != it->end(); ++tr )
		{
			tr->Emis().TauTrack().clear();
		}
	}

	/* usually computed in pressure_change, but not called before first zone */
	if (wind.comass == 0.0)
	{
		wind.AccelGravity = 0.f;
	}
	else
	{
		if (radius.Radius < 0.0)
		{
			fprintf(ioQQQ,"PROBLEM: must have positive initial radius for central gravity model\n");
			cdEXIT(EXIT_FAILURE);
		}
		wind.AccelGravity = (realnum)(GRAV_CONST*wind.comass*SOLAR_MASS/POW2(radius.Radius)*
			(1.-wind.DiskRadius/radius.Radius) );
	}
	if( trace.lgTrace && trace.lgWind )
	{
		fprintf(ioQQQ," InitSimPostparse sets AccelGravity %.3e lgDisk?%c\n",
			wind.AccelGravity , 
			TorF(wind.lgDisk) );
	}

	/* pressure related variables */
	pressure.PresInteg = 0.;
	pressure.PresIntegElecThin = 0.;
	pressure.pinzon = 0.;

	/* update iso sequence level information */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			iso_sp[ipISO][nelem].Reset();
			/* these elements that are turned off */
			if( !dense.lgElmtOn[nelem] )
			{
				iso_sp[ipISO][nelem].numLevels_max = 0;
				iso_sp[ipISO][nelem].nCollapsed_max = 0;
				iso_sp[ipISO][nelem].n_HighestResolved_max = 0;

				iso_sp[ipISO][nelem].numLevels_local = 0;
				iso_sp[ipISO][nelem].nCollapsed_local = 0;
				iso_sp[ipISO][nelem].n_HighestResolved_local = 0;
			}
			else
			{
				iso_update_num_levels( ipISO, nelem );
				/* must always have at least one collapsed level, unless compiling recomb data file. */
				ASSERT( iso_sp[ipISO][nelem].nCollapsed_max > 0 || iso_ctrl.lgCompileRecomb[ipISO] );
			}
		}
	}

	if( iso_ctrl.lgPrintNumberOfLevels )
	{
		fprintf(ioQQQ,"\n\n Number of levels in ions treated by iso sequences.\n");
		fprintf(ioQQQ," ISO   Element  hi-n(l-resolved) #(l-resolved)  n(collapsed)\n");
		/* option to print number of levels for each element */
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( long nelem=ipISO; nelem<LIMELM; ++nelem )
			{
				/* print number of levels */
				fprintf(ioQQQ," %s  %s    %4li            %4li           %4li \n",
					iso_ctrl.chISO[ipISO] ,
					elementnames.chElementSym[nelem],
					iso_sp[ipISO][nelem].n_HighestResolved_max,
					iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max,
					iso_sp[ipISO][nelem].nCollapsed_max);
			}
		}
	}
	atoms.p2nit = 0.;
	atoms.d5200r = 0.;
	atoms.rateMg2 = 0.;

	MonitorResults.SumErrorCaseMonitor = 0.;
	MonitorResults.nSumErrorCaseMonitor = 0;

	mean.setup_molecules();

	return;
}
